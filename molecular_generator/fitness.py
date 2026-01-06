"""
Fitness functions for evaluating fuel molecules.

Provides functions to calculate molecular properties and assign fitness scores
based on fuel-relevant criteria.
"""

from typing import Dict, Optional, List
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


def calculate_molecular_properties(smiles: str) -> Optional[Dict]:
    """
    Calculate basic molecular descriptors for a given SMILES string.
    
    Args:
        smiles: SMILES string of the molecule
    
    Returns:
        Dictionary with molecular properties, or None if invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        
        # Count atoms
        atoms = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atoms[symbol] = atoms.get(symbol, 0) + 1
        
        num_C = atoms.get('C', 0)
        num_O = atoms.get('O', 0)
        num_H = atoms.get('H', 0)
        
        # Calculate O/C ratio
        o_c_ratio = num_O / num_C if num_C > 0 else 0.0
        
        return {
            'SMILES': smiles,
            'MW': round(mw, 2),
            'C': num_C,
            'O': num_O,
            'H': num_H,
            'H_donors': h_donors,
            'logP': round(logp, 2),
            'rings': rings,
            'O/C': round(o_c_ratio, 3),
        }
    except:
        return None


def compute_fuel_fitness(properties: Dict) -> float:
    """
    Compute fitness score for a fuel molecule based on its properties.
    
    Fitness considers:
    - Energy potential (molecular weight, logP)
    - Oxygen balance (O/C ratio)
    - Volatility (logP)
    - H-bonding (reactivity)
    
    Args:
        properties: Dictionary with molecular properties (from calculate_molecular_properties)
    
    Returns:
        Fitness score (higher is better)
    """
    mw = properties.get('MW', 0)
    C = properties.get('C', 0)
    O = properties.get('O', 0)
    logp = properties.get('logP', 0)
    h_donors = properties.get('H_donors', 0)
    
    # Energy score (simplified: MW weighted by logP)
    energy_score = (mw / 10) + (logp * 2)
    
    # Oxygen balance - want ~10-20% oxygen content (O/C ~0.1-0.2)
    o_c_ratio = properties.get('O/C', 0)
    if o_c_ratio < 0.05:
        oxygen_score = -2
    elif o_c_ratio > 0.5:
        oxygen_score = -2
    else:
        oxygen_score = 2 - abs(o_c_ratio - 0.2) * 10  # peak fitness near 0.2
    
    # Volatility/handling: penalize extreme logP
    volatility_penalty = -abs(logp - 1.5) * 2  # ideal ~1.5
    
    # H-bonding: minor penalty for reactivity
    h_penalty = -h_donors * 0.5
    
    # Total fitness
    fitness = energy_score + oxygen_score + volatility_penalty + h_penalty
    return round(fitness, 2)


def compute_fuel_fitness_v3(smiles: str) -> float:
    """
    Compute fitness score directly from SMILES string (v3 algorithm).
    
    This is a more refined version that considers:
    - Molecular weight (target ~150)
    - logP (target ~1.5)
    - H-bond donors (penalized)
    - O/C ratio (target ~0.2)
    - Ring count (penalized)
    
    Args:
        smiles: SMILES string of the molecule
    
    Returns:
        Fitness score (higher is better, minimum 0)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        
        # Count atoms
        atoms = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atoms[symbol] = atoms.get(symbol, 0) + 1
        
        num_C = atoms.get('C', 1)  # Avoid division by zero
        num_O = atoms.get('O', 0)
        o_c = num_O / num_C
        
        # Fitness calculation
        score = 50.0
        score -= abs(mw - 150) * 0.05  # Penalize deviation from ideal MW
        score -= abs(logp - 1.5) * 2  # Penalize deviation from ideal logP
        score -= h_donors * 1.5  # Penalize H-bond donors
        score -= abs(o_c - 0.2) * 15  # Penalize deviation from ideal O/C
        score -= rings * 0.3  # Slight penalty for rings
        
        return max(score, 0.0)
    except:
        return 0.0


def count_functional_groups(mol) -> Dict[str, int]:
    """
    Count different functional groups in a molecule.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dictionary with counts of functional groups
    """
    try:
        ester_pattern = Chem.MolFromSmiles('C(=O)O')
        ether_pattern = Chem.MolFromSmiles('COC')
        alcohol_pattern = Chem.MolFromSmiles('CO')
        acid_pattern = Chem.MolFromSmiles('C(=O)O')
        aldehyde_pattern = Chem.MolFromSmiles('C=O')
        ketone_pattern = Chem.MolFromSmiles('CC(=O)C')
        
        return {
            'esters': len(mol.GetSubstructMatches(ester_pattern)),
            'ethers': len(mol.GetSubstructMatches(ether_pattern)),
            'alcohols': len(mol.GetSubstructMatches(alcohol_pattern)),
            'acids': len(mol.GetSubstructMatches(acid_pattern)),
            'aldehydes': len(mol.GetSubstructMatches(aldehyde_pattern)),
            'ketones': len(mol.GetSubstructMatches(ketone_pattern)),
        }
    except:
        return {}


def has_ester_group(mol) -> bool:
    """Check if molecule has an ester group."""
    try:
        ester_pattern = Chem.MolFromSmiles('C(=O)O')
        return mol.HasSubstructMatch(ester_pattern)
    except:
        return False


def has_ether_group(mol) -> bool:
    """Check if molecule has an ether group."""
    try:
        ether_pattern = Chem.MolFromSmiles('COC')
        return mol.HasSubstructMatch(ether_pattern)
    except:
        return False


def has_carboxylic_acid(mol) -> bool:
    """Check if molecule has a carboxylic acid group."""
    try:
        # Carboxylic acid: C(=O)O with H attached
        acid_pattern = Chem.MolFromSmiles('C(=O)O')
        return mol.HasSubstructMatch(acid_pattern)
    except:
        return False


def has_alcohol_group(mol) -> bool:
    """Check if molecule has an alcohol group."""
    try:
        alcohol_pattern = Chem.MolFromSmiles('CO')
        return mol.HasSubstructMatch(alcohol_pattern)
    except:
        return False


def estimate_synthesis_feasibility(smiles: str) -> Dict[str, any]:
    """
    Estimate if a molecule can be synthesized at home.
    
    Returns:
        Dictionary with:
        - feasibility: Score 0-100 (higher = easier to make)
        - complexity: Rating 1-10 (lower = simpler)
        - reactions: List of reaction types needed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {'feasibility': 0.0, 'complexity': 10.0, 'reactions': []}
        
        feasibility = 50.0  # Base score
        complexity = 5.0
        reactions = []
        
        # Get molecular properties
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        
        # Reward smaller molecules
        if mw < 200:
            feasibility += 20.0
            complexity -= 1.0
        elif mw > 300:
            feasibility -= 20.0
            complexity += 2.0
        
        # Check for ester groups (easy: acid + alcohol)
        if has_ester_group(mol):
            feasibility += 30.0
            complexity -= 2.0
            reactions.append('esterification')
        
        # Check for ether groups (moderate: Williamson)
        if has_ether_group(mol):
            feasibility += 20.0
            complexity -= 1.0
            reactions.append('etherification')
        
        # Check for alcohol groups (can be starting material or product)
        if has_alcohol_group(mol):
            feasibility += 10.0
        
        # Penalize complex structures
        if rings > 2:
            feasibility -= 30.0
            complexity += 2.0
        elif rings > 0:
            feasibility -= 10.0
            complexity += 1.0
        
        # Count functional groups
        func_groups = count_functional_groups(mol)
        total_groups = sum(func_groups.values())
        if total_groups > 3:
            feasibility -= 20.0
            complexity += 1.0
        
        # Penalize unsaturation (harder to control)
        double_bonds = sum(1 for bond in mol.GetBonds() 
                          if bond.GetBondType() == Chem.BondType.DOUBLE)
        if double_bonds > 2:
            feasibility -= 15.0
            complexity += 1.0
        
        # Penalize triple bonds (very difficult)
        triple_bonds = sum(1 for bond in mol.GetBonds() 
                          if bond.GetBondType() == Chem.BondType.TRIPLE)
        if triple_bonds > 0:
            feasibility -= 25.0
            complexity += 2.0
        
        return {
            'feasibility': max(0.0, min(100.0, feasibility)),
            'complexity': max(1.0, min(10.0, complexity)),
            'reactions': list(set(reactions))  # Remove duplicates
        }
    except:
        return {'feasibility': 0.0, 'complexity': 10.0, 'reactions': []}


def estimate_oxidation_susceptibility(smiles: str) -> float:
    """
    Estimate oxidation susceptibility based on molecular structure.
    Higher values = more prone to oxidation.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Oxidation susceptibility score (0-10, higher = more prone)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
        
        oxidation_score = 0.0
        
        # Count double bonds (unsaturation prone to oxidation)
        double_bonds = sum(1 for bond in mol.GetBonds() 
                          if bond.GetBondType() == Chem.BondType.DOUBLE)
        oxidation_score += double_bonds * 1.5
        
        # Count aldehyde groups (very prone to oxidation)
        aldehyde_pattern = Chem.MolFromSmiles('C=O')
        aldehydes = len(mol.GetSubstructMatches(aldehyde_pattern))
        oxidation_score += aldehydes * 3.0
        
        # Count allylic positions (C=C-C, prone to autoxidation)
        # Simplified: count double bonds with adjacent carbons
        oxidation_score += double_bonds * 0.5
        
        # Count benzylic positions (aromatic ring with CH2, prone to oxidation)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        oxidation_score += aromatic_rings * 1.0
        
        # Alcohols can oxidize to aldehydes/ketones
        if has_alcohol_group(mol):
            oxidation_score += 1.0
        
        return min(10.0, oxidation_score)
    except:
        return 0.0


def compute_home_synthesis_fitness(
    smiles: str,
    synthesis_weight: float = 0.5,
    fuel_weight: float = 0.3
) -> float:
    """
    Fitness function that prioritizes:
    1. Synthesis feasibility (can be made at home)
    2. Fuel properties (good performance)
    3. Simplicity (inverse of complexity)
    
    Args:
        smiles: SMILES string of the molecule
        synthesis_weight: Weight for synthesis feasibility (default 0.5)
        fuel_weight: Weight for fuel properties (default 0.3)
    
    Returns:
        Combined fitness score (higher is better)
    """
    # Get standard fuel fitness
    fuel_fitness = compute_fuel_fitness_v3(smiles)
    
    # Get synthesis feasibility
    synth = estimate_synthesis_feasibility(smiles)
    synth_score = synth['feasibility']
    
    # Simplicity score (inverse of complexity, normalized to 0-100)
    simplicity_score = (10.0 - synth['complexity']) * 10.0
    
    # Combined score
    total_fitness = (
        synth_score * synthesis_weight +
        fuel_fitness * fuel_weight +
        simplicity_score * (1.0 - synthesis_weight - fuel_weight)
    )
    
    return max(0.0, total_fitness)


def compute_stability_test_fitness(
    smiles: str,
    target_oxidation: Optional[float] = None
) -> float:
    """
    Generate molecules specifically for stability testing.
    
    Can target:
    - High oxidation risk (for testing)
    - Low oxidation risk (for controls)
    - Intermediate polarity (for phase separation tests)
    
    Args:
        smiles: SMILES string of the molecule
        target_oxidation: Target oxidation susceptibility (None = moderate)
    
    Returns:
        Fitness score for stability testing
    """
    oxidation = estimate_oxidation_susceptibility(smiles)
    
    if target_oxidation is not None:
        # Reward molecules near target oxidation level
        score = 50.0 - abs(oxidation - target_oxidation) * 5.0
    else:
        # Default: reward moderate oxidation risk (good for testing)
        score = 50.0 - abs(oxidation - 3.0) * 5.0
    
    # Also consider synthesis feasibility
    synth = estimate_synthesis_feasibility(smiles)
    score += synth['feasibility'] * 0.3
    
    return max(0.0, score)


def estimate_phase_separation_risk(
    components: List[Dict],
    vol_fractions: List[float]
) -> float:
    """
    Estimate risk of phase separation based on polarity differences.
    
    Args:
        components: List of component dictionaries with 'logP' or similar
        vol_fractions: Volume fractions for each component
    
    Returns:
        Phase separation risk score (0-10, higher = more risk)
    """
    try:
        # Extract logP values (or estimate from molecular properties)
        logp_values = []
        for comp in components:
            if 'logP' in comp:
                logp_values.append(comp['logP'])
            elif 'SMILES' in comp:
                # Estimate logP from SMILES
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(comp['SMILES'])
                if mol:
                    logp_values.append(Descriptors.MolLogP(mol))
        
        if len(logp_values) < 2:
            return 0.0
        
        # Calculate weighted average logP
        weighted_logp = sum(logp * frac for logp, frac in zip(logp_values, vol_fractions))
        
        # Calculate variance in logP (higher variance = more separation risk)
        variance = sum(frac * (logp - weighted_logp)**2 
                      for logp, frac in zip(logp_values, vol_fractions))
        
        # Large logP differences indicate phase separation risk
        logp_range = max(logp_values) - min(logp_values)
        
        # Risk increases with logP range
        risk = min(10.0, logp_range * 0.5 + variance * 2.0)
        
        return risk
    except:
        return 0.0


def estimate_ph_change(smiles: str) -> Dict[str, float]:
    """
    Estimate potential pH change from molecular structure.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Dictionary with:
        - initial_ph_estimate: Estimated initial pH
        - acid_groups: Number of acid groups
        - base_groups: Number of base groups
        - ph_change_potential: Potential for pH change (0-10)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {
                'initial_ph_estimate': 7.0,
                'acid_groups': 0,
                'base_groups': 0,
                'ph_change_potential': 0.0
            }
        
        # Count carboxylic acids
        acid_pattern = Chem.MolFromSmiles('C(=O)O')
        acid_groups = len(mol.GetSubstructMatches(acid_pattern))
        
        # Count amines (basic groups)
        amine_pattern = Chem.MolFromSmiles('CN')
        base_groups = len(mol.GetSubstructMatches(amine_pattern))
        
        # Estimate initial pH
        if acid_groups > 0:
            initial_ph = 3.0 + acid_groups * 0.5  # Acids lower pH
        elif base_groups > 0:
            initial_ph = 11.0 - base_groups * 0.5  # Bases raise pH
        else:
            initial_ph = 7.0  # Neutral
        
        # pH change potential increases with number of acid/base groups
        ph_change_potential = min(10.0, (acid_groups + base_groups) * 2.0)
        
        return {
            'initial_ph_estimate': max(0.0, min(14.0, initial_ph)),
            'acid_groups': acid_groups,
            'base_groups': base_groups,
            'ph_change_potential': ph_change_potential
        }
    except:
        return {
            'initial_ph_estimate': 7.0,
            'acid_groups': 0,
            'base_groups': 0,
            'ph_change_potential': 0.0
        }