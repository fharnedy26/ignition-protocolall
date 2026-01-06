"""
Genetic algorithm for evolving fuel molecules.

Provides functions for mutation, crossover, and population evolution.
"""

import random
from typing import List, Tuple, Optional, Dict, Union, Iterator
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Silence RDKit warnings (handle different RDKit versions)
try:
    from rdkit.Chem import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        # RDLogger not available in this RDKit version, warnings will show
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)

from molecular_generator.fitness import (
    calculate_molecular_properties,
    compute_fuel_fitness,
    compute_home_synthesis_fitness,
    estimate_synthesis_feasibility
)


def is_valid_fuel(smiles: str, allowed_atoms: set = {'C', 'H', 'O'}) -> bool:
    """
    Check if SMILES is a valid molecule and contains only allowed atoms.
    
    Args:
        smiles: SMILES string to validate
        allowed_atoms: Set of allowed atomic symbols
    
    Returns:
        True if valid fuel molecule, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
        return atoms.issubset(allowed_atoms)
    except:
        return False


def mutate_smiles(smiles: str) -> str:
    """
    Randomly modify a SMILES string to produce a new variant.
    
    Mutation types:
    - insert: Add a random token
    - delete: Remove a token
    - swap: Swap two tokens
    
    Args:
        smiles: Original SMILES string
    
    Returns:
        Mutated SMILES string
    """
    tokens = list(smiles)
    if len(tokens) < 3:
        return smiles  # Too short to mutate
    
    mutation_type = random.choice(['insert', 'delete', 'swap'])
    idx = random.randint(0, len(tokens) - 1)
    
    if mutation_type == 'insert':
        insert_token = random.choice(['C', 'O', '=', '(', ')'])
        tokens.insert(idx, insert_token)
    
    elif mutation_type == 'delete' and len(tokens) > 3:
        tokens.pop(idx)
    
    elif mutation_type == 'swap' and len(tokens) > 3:
        idx2 = random.randint(0, len(tokens) - 1)
        tokens[idx], tokens[idx2] = tokens[idx2], tokens[idx]
    
    return ''.join(tokens)


def evolve_population(
    initial_population: pd.DataFrame,
    generations: int = 10,
    top_k: int = 10,
    children_per_parent: int = 5,
    verbose: bool = True
) -> Tuple[pd.DataFrame, List[float]]:
    """
    Run evolution over several generations.
    
    Args:
        initial_population: DataFrame with columns ['SMILES', 'Fitness', ...]
        generations: Number of generations to evolve
        top_k: Number of top molecules to select as parents
        children_per_parent: Number of children to generate per parent
        verbose: If True, print progress
    
    Returns:
        Tuple of (final_population_dataframe, fitness_history)
    """
    best_fitness = []
    population = initial_population.copy()
    
    for gen in range(generations):
        # Select top performers
        top = population.sort_values(by='Fitness', ascending=False).head(top_k)
        new_smiles = []
        
        # Generate children through mutation
        for smiles in top['SMILES']:
            for _ in range(children_per_parent):
                child = mutate_smiles(smiles)
                if is_valid_fuel(child):
                    new_smiles.append(child)
        
        # Calculate properties for new generation
        new_data = []
        for s in new_smiles:
            props = calculate_molecular_properties(s)
            if props:
                props['Fitness'] = compute_fuel_fitness(props)
                new_data.append(props)
        
        if new_data:
            population = pd.DataFrame(new_data)
        
        # Track best fitness
        if len(population) > 0:
            best = population.sort_values(by='Fitness', ascending=False).iloc[0]
            best_fitness.append(best['Fitness'])
            
            if verbose:
                print(f"Gen {gen+1}: Best fitness = {best['Fitness']:.2f} | "
                      f"SMILES = {best['SMILES']}")
        else:
            if verbose:
                print(f"Gen {gen+1}: No valid molecules generated")
            best_fitness.append(0.0)
    
    return population, best_fitness


def evolve_with_fragments(
    fragments: Union[List[str], Iterator[str], str],
    generations: int = 100,
    population_size: int = 50,
    survivors: int = 10,
    n_frags_per_molecule: int = 3,
    verbose: bool = True,
    fragment_cache_size: int = 100000
) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Evolve molecules using fragment-based assembly.
    
    Args:
        fragments: Fragment source - can be:
            - List[str]: List of fragment SMILES strings (backward compatible)
            - Iterator[str]: Iterator yielding fragment SMILES strings
            - str: File path to fragment file (supports .tgz, .gz, .smi, .txt)
        generations: Number of generations
        population_size: Size of population each generation
        survivors: Number of top molecules to keep
        n_frags_per_molecule: Number of fragments to combine
        verbose: If True, print progress
        fragment_cache_size: Cache size when using file paths or iterators (default: 100000)
    
    Returns:
        Tuple of (final_population_dataframe, evolution_history)
    """
    from molecular_generator.fragments import (
        merge_fragments,
        build_fragment_cache,
        is_fuel_relevant_fragment
    )
    from molecular_generator.fitness import compute_fuel_fitness_v3
    
    # Handle different fragment input types
    if isinstance(fragments, str):
        # File path - build cache with fuel-relevance filter
        if verbose:
            print(f"Loading fragments from file: {fragments}")
        fragments = build_fragment_cache(
            fragments,
            cache_size=fragment_cache_size,
            filter_func=is_fuel_relevant_fragment
        )
        if verbose:
            print(f"Loaded {len(fragments):,} fuel-relevant fragments from cache")
    elif hasattr(fragments, '__iter__') and not isinstance(fragments, (str, list)):
        # Iterator (but not a list or string) - build cache from iterator using reservoir sampling
        if verbose:
            print("Building fragment cache from iterator...")
        cache = []
        count = 0
        for smiles in fragments:
            if is_fuel_relevant_fragment(smiles):
                count += 1
                if len(cache) < fragment_cache_size:
                    cache.append(smiles)
                else:
                    # Reservoir sampling: replace with probability fragment_cache_size/count
                    j = random.randint(0, count - 1)
                    if j < fragment_cache_size:
                        cache[j] = smiles
        fragments = cache
        if verbose:
            print(f"Loaded {len(fragments):,} fuel-relevant fragments from iterator")
    # If fragments is already a list, use it directly (backward compatible)
    
    if not fragments or len(fragments) == 0:
        raise ValueError("No fragments available for evolution")
    
    history = []
    population = []
    seen_smiles = set()
    
    # Initialize population
    attempts = 0
    while len(population) < population_size and attempts < population_size * 10:
        selected = random.sample(fragments, k=min(n_frags_per_molecule, len(fragments)))
        if len(selected) >= 2:
            merged = merge_fragments(selected[0], selected[1])
            if merged and merged not in seen_smiles and is_valid_fuel(merged):
                score = compute_fuel_fitness_v3(merged)
                population.append((merged, score))
                seen_smiles.add(merged)
        attempts += 1
    
    # Evolution loop with early stopping and diversity injection
    best_fitness_history = []
    plateau_count = 0
    plateau_threshold = 50  # Stop if no improvement for 50 generations
    
    for gen in range(generations):
        # Sort by fitness
        population.sort(key=lambda x: x[1], reverse=True)
        
        # Track best
        if population:
            best_smiles, best_score = population[0]
            best_fitness_history.append(best_score)
            
            # Early stopping: check for plateau
            if len(best_fitness_history) > plateau_threshold:
                recent_best = max(best_fitness_history[-plateau_threshold:])
                if abs(recent_best - best_score) < 0.01:  # No improvement within threshold
                    plateau_count += 1
                    if plateau_count >= 5 and verbose:
                        print(f"Gen {gen+1}: Plateau detected (no improvement for {plateau_threshold} gens). Stopping early.")
                        break
                else:
                    plateau_count = 0
            
            mol = Chem.MolFromSmiles(best_smiles)
            if mol:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                h_donors = Descriptors.NumHDonors(mol)
                rings = rdMolDescriptors.CalcNumRings(mol)
                
                atoms = {}
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    atoms[symbol] = atoms.get(symbol, 0) + 1
                o_c = atoms.get('O', 0) / max(atoms.get('C', 1), 1)
                
                history.append({
                    'Generation': gen + 1,
                    'Best_SMILES': best_smiles,
                    'Fitness': best_score,
                    'MolWt': mw,
                    'logP': logp,
                    'H_Donors': h_donors,
                    'Rings': rings,
                    'O/C': o_c,
                })
        
        # Select survivors
        survivors_list = population[:survivors]
        
        # Generate new population
        new_population = survivors_list.copy()
        attempts = 0
        
        # Diversity injection: occasionally add random new molecules
        # Increase diversity rate if plateau detected
        base_diversity_rate = 0.25  # 25% base diversity injection (increased from 15% for better exploration)
        if plateau_count > 0:
            diversity_rate = min(0.4, base_diversity_rate + plateau_count * 0.05)  # Increase up to 40%
        else:
            diversity_rate = base_diversity_rate
        
        # Also keep some lower-fitness but unique molecules for diversity
        diversity_keep = max(2, int(population_size * 0.1))  # Keep 10% for diversity
        
        while len(new_population) < population_size and attempts < population_size * 10:
            # Diversity injection: random exploration
            if random.random() < diversity_rate and len(fragments) >= 2:
                # Create completely new molecule from random fragments
                selected = random.sample(fragments, k=min(n_frags_per_molecule, len(fragments)))
                if len(selected) >= 2:
                    merged = merge_fragments(selected[0], selected[1])
                    if merged and merged not in seen_smiles and is_valid_fuel(merged):
                        score = compute_fuel_fitness_v3(merged)
                        new_population.append((merged, score))
                        seen_smiles.add(merged)
            else:
                # Crossover: merge two random survivors (or random from population for diversity)
                if len(survivors_list) >= 2:
                    # Occasionally pick from broader population for diversity
                    if random.random() < 0.2 and len(population) > survivors:
                        # Pick one from survivors, one from broader population
                        p1 = random.choice(survivors_list)
                        p2 = random.choice(population[survivors:min(survivors+diversity_keep, len(population))])
                    else:
                        # Standard: both from survivors
                        p1, p2 = random.sample(survivors_list, 2)
                    
                    merged = merge_fragments(p1[0], p2[0])
                    if merged and merged not in seen_smiles and is_valid_fuel(merged):
                        score = compute_fuel_fitness_v3(merged)
                        new_population.append((merged, score))
                        seen_smiles.add(merged)
            attempts += 1
        
        # Ensure we keep some diversity by adding unique lower-fitness molecules
        if len(new_population) < population_size and len(population) > survivors:
            diversity_candidates = population[survivors:min(survivors+diversity_keep*2, len(population))]
            for candidate in diversity_candidates:
                if len(new_population) >= population_size:
                    break
                if candidate[0] not in seen_smiles:
                    new_population.append(candidate)
                    seen_smiles.add(candidate[0])
        
        population = new_population
        
        if verbose and (gen + 1) % 10 == 0:
            if population:
                print(f"Gen {gen+1}: Best fitness = {population[0][1]:.2f}")
    
    # Convert to DataFrame
    df = pd.DataFrame(history)
    return df, history


def evolve_home_synthesizable_molecules(
    fragments: List[str],
    generations: int = 100,
    population_size: int = 50,
    survivors: int = 10,
    n_frags_per_molecule: int = 3,
    min_feasibility: float = 60.0,
    synthesis_weight: float = 0.5,
    fuel_weight: float = 0.3,
    verbose: bool = True
) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Evolve molecules with synthesis feasibility constraint.
    
    This function generates molecules optimized for both fuel performance
    and home synthesis feasibility.
    
    Args:
        fragments: List of fragment SMILES strings (preferably from
                   home_synthesis_fragments)
        generations: Number of generations
        population_size: Size of population each generation
        survivors: Number of top molecules to keep
        n_frags_per_molecule: Number of fragments to combine
        min_feasibility: Minimum synthesis feasibility score (0-100)
        synthesis_weight: Weight of synthesis in fitness (0-1)
        fuel_weight: Weight of fuel properties in fitness (0-1)
        verbose: If True, print progress
        
    Returns:
        Tuple of (final_population_dataframe, evolution_history)
        DataFrame includes synthesis metadata columns
    """
    from molecular_generator.fragments import merge_fragments
    from molecular_generator.synthesis import suggest_synthesis_route
    
    history = []
    population = []
    seen_smiles = set()
    
    # Initialize population
    attempts = 0
    while len(population) < population_size and attempts < population_size * 10:
        selected = random.sample(fragments, k=min(n_frags_per_molecule, len(fragments)))
        if len(selected) >= 2:
            merged = merge_fragments(selected[0], selected[1])
            if merged and merged not in seen_smiles and is_valid_fuel(merged):
                # Use home synthesis fitness
                score = compute_home_synthesis_fitness(
                    merged,
                    synthesis_weight=synthesis_weight,
                    fuel_weight=fuel_weight
                )
                population.append((merged, score))
                seen_smiles.add(merged)
        attempts += 1
    
    # Evolution loop
    for gen in range(generations):
        # Sort by fitness
        population.sort(key=lambda x: x[1], reverse=True)
        
        # Track best
        if population:
            best_smiles, best_score = population[0]
            mol = Chem.MolFromSmiles(best_smiles)
            if mol:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                h_donors = Descriptors.NumHDonors(mol)
                rings = rdMolDescriptors.CalcNumRings(mol)
                
                atoms = {}
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    atoms[symbol] = atoms.get(symbol, 0) + 1
                o_c = atoms.get('O', 0) / max(atoms.get('C', 1), 1)
                
                # Get synthesis metadata
                synth = estimate_synthesis_feasibility(best_smiles)
                
                history.append({
                    'Generation': gen + 1,
                    'Best_SMILES': best_smiles,
                    'Fitness': best_score,
                    'MolWt': mw,
                    'logP': logp,
                    'H_Donors': h_donors,
                    'Rings': rings,
                    'O/C': o_c,
                    'Synthesis_Feasibility': synth['feasibility'],
                    'Synthesis_Complexity': synth['complexity'],
                })
        
        # Select survivors
        survivors_list = population[:survivors]
        
        # Generate new population
        new_population = survivors_list.copy()
        attempts = 0
        while len(new_population) < population_size and attempts < population_size * 10:
            # Crossover: merge two random survivors
            if len(survivors_list) >= 2:
                p1, p2 = random.sample(survivors_list, 2)
                merged = merge_fragments(p1[0], p2[0])
                if merged and merged not in seen_smiles and is_valid_fuel(merged):
                    # Use home synthesis fitness
                    score = compute_home_synthesis_fitness(
                        merged,
                        synthesis_weight=synthesis_weight,
                        fuel_weight=fuel_weight
                    )
                    new_population.append((merged, score))
                    seen_smiles.add(merged)
            attempts += 1
        
        population = new_population
        
        if verbose and (gen + 1) % 10 == 0:
            if population:
                best_smiles, best_score = population[0]
                synth = estimate_synthesis_feasibility(best_smiles)
                print(f"Gen {gen+1}: Best fitness = {best_score:.2f} | "
                      f"Feasibility = {synth['feasibility']:.1f} | "
                      f"Complexity = {synth['complexity']:.1f}")
    
    # Convert to DataFrame and filter by feasibility
    df = pd.DataFrame(history)
    
    # Filter by minimum feasibility
    if len(df) > 0:
        df = df[df['Synthesis_Feasibility'] >= min_feasibility]
        
        # Add synthesis route information
        if len(df) > 0:
            routes = []
            for smiles in df['Best_SMILES']:
                route = suggest_synthesis_route(smiles)
                routes.append(route)
            df['Synthesis_Route'] = routes
    
    return df, history

