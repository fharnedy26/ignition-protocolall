"""
Bridge module: Convert generated molecules to component format for blend optimizer.

This module converts SMILES molecules into the component CSV format required
by the fuel blend optimization engine.
"""

import hashlib
from typing import Dict, List, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from molecular_generator.fitness import (
    calculate_molecular_properties,
    estimate_synthesis_feasibility,
    estimate_oxidation_susceptibility
)
from molecular_generator.synthesis import suggest_synthesis_route


def estimate_fuel_properties(smiles: str) -> Dict[str, float]:
    """
    Estimate fuel properties from molecular structure.
    
    Uses empirical correlations and molecular descriptors to estimate:
    - RON, MON (octane numbers)
    - Cetane number
    - LHV (lower heating value)
    - PMI, TSI (sooting indices)
    - Density
    - Distillation temperatures (T10, T50, T90)
    - Vapor pressure (RVP)
    
    Args:
        smiles: SMILES string of the molecule
    
    Returns:
        Dictionary with estimated fuel properties
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        # Basic descriptors
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
        o_c_ratio = num_O / num_C if num_C > 0 else 0.0
        
        # Estimate RON (Research Octane Number)
        # Simplified correlation: higher MW and lower O/C generally higher RON
        # Branched molecules have higher RON
        branching = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ron = 90.0 + (mw - 100) * 0.1 - o_c_ratio * 20 + branching * 2
        ron = max(50.0, min(120.0, ron))  # Clamp to reasonable range
        
        # MON (Motor Octane Number) is typically 5-10 points lower than RON
        mon = ron - 7.0
        mon = max(50.0, min(115.0, mon))
        
        # Cetane number (for diesel)
        # Higher MW, lower branching, higher O/C -> higher cetane
        cetane = 50.0 + (mw - 100) * 0.15 - branching * 1.5 + o_c_ratio * 10
        cetane = max(20.0, min(80.0, cetane))
        
        # LHV (Lower Heating Value) in MJ/kg
        # Based on elemental composition (Dulong formula approximation)
        # LHV ≈ 33.9*C + 121.4*H - 15.3*O (in MJ/kg, simplified)
        # More accurate: use molecular formula
        total_atoms = num_C + num_H + num_O
        if total_atoms > 0:
            # Simplified: assume ~42 MJ/kg for hydrocarbons, reduce with O
            lhv_mj_kg = 42.0 - o_c_ratio * 15.0
        else:
            lhv_mj_kg = 42.0
        lhv_mj_kg = max(20.0, min(50.0, lhv_mj_kg))
        
        # PMI (Particulate Matter Index) - sooting tendency
        # Higher aromatic content and MW -> higher PMI
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        pmi = 10.0 + (mw - 100) * 0.1 + aromatic_rings * 5.0 + rings * 2.0
        pmi = max(5.0, min(40.0, pmi))
        
        # TSI (Threshold Sooting Index) - similar to PMI
        tsi = pmi * 0.9  # Slightly lower than PMI
        
        # Density (g/mL) - estimate from MW and structure
        # Simplified: assume ~0.75 g/mL base, increase with MW and rings
        density = 0.70 + (mw - 100) * 0.0005 + rings * 0.05
        density = max(0.65, min(0.95, density))
        
        # Distillation temperatures (boiling point estimates)
        # Use logP and MW as proxies
        bp_base = 50.0 + mw * 0.5  # Base boiling point
        bp_correction = -logp * 10  # Lower logP -> higher BP
        
        t10 = bp_base + bp_correction - 20
        t50 = bp_base + bp_correction
        t90 = bp_base + bp_correction + 30
        
        # Clamp to reasonable ranges
        t10 = max(30.0, min(100.0, t10))
        t50 = max(60.0, min(200.0, t50))
        t90 = max(100.0, min(350.0, t90))
        
        # Vapor pressure (RVP in kPa)
        # Higher MW and lower logP -> lower vapor pressure
        rvp = 100.0 - mw * 0.3 + logp * 5.0
        rvp = max(1.0, min(150.0, rvp))
        
        return {
            'RON': round(ron, 1),
            'MON': round(mon, 1),
            'cetane': round(cetane, 1),
            'LHV_MJ_kg': round(lhv_mj_kg, 2),
            'PMI': round(pmi, 1),
            'TSI': round(tsi, 1),
            'density_g_ml': round(density, 3),
            'T10_C': round(t10, 1),
            'T50_C': round(t50, 1),
            'T90_C': round(t90, 1),
            'vapor_pressure_kPa': round(rvp, 1),
            'OC_ratio': round(o_c_ratio, 3),
            'ring_count': rings,
        }
    except Exception as e:
        # Return defaults on error
        return {
            'RON': 90.0,
            'MON': 85.0,
            'cetane': 50.0,
            'LHV_MJ_kg': 42.0,
            'PMI': 20.0,
            'TSI': 20.0,
            'density_g_ml': 0.75,
            'T10_C': 50.0,
            'T50_C': 100.0,
            'T90_C': 150.0,
            'vapor_pressure_kPa': 10.0,
            'OC_ratio': 0.0,
            'ring_count': 0,
        }


def molecule_to_component(
    smiles: str,
    name: Optional[str] = None,
    cost_eur_L: float = 2.0,
    max_vol_frac: float = 0.1,
    is_novel: int = 1,
    waste_credit_eur_L: float = 0.0
) -> Dict[str, any]:
    """
    Convert a SMILES molecule to component dictionary format.
    
    Args:
        smiles: SMILES string of the molecule
        name: Optional name for the component (defaults to SMILES)
        cost_eur_L: Cost per liter in EUR
        max_vol_frac: Maximum volume fraction allowed in blends
        is_novel: Whether this is a novel component (0 or 1)
        waste_credit_eur_L: Waste credit per liter (for novel components)
    
    Returns:
        Dictionary with component data in blend optimizer format
    """
    # Generate unique ID from SMILES
    mol_id = hashlib.md5(smiles.encode()).hexdigest()[:8]
    
    # Estimate properties
    props = estimate_fuel_properties(smiles)
    
    # Get basic properties
    basic_props = calculate_molecular_properties(smiles)
    if not basic_props:
        basic_props = {}
    
    # Get synthesis metadata
    synth = estimate_synthesis_feasibility(smiles)
    synthesis_route = suggest_synthesis_route(smiles)
    oxidation_susceptibility = estimate_oxidation_susceptibility(smiles)
    
    # Build component dictionary
    component = {
        'id': f'MOL_{mol_id}',
        'name': name or smiles[:50],  # Truncate long SMILES
        'density_g_ml': props.get('density_g_ml', 0.75),
        'LHV_MJ_kg': props.get('LHV_MJ_kg', 42.0),
        'RON': props.get('RON', 90.0),
        'MON': props.get('MON', 85.0),
        'cetane': props.get('cetane', 50.0),
        'PMI': props.get('PMI', 20.0),
        'TSI': props.get('TSI', 20.0),
        'OC_ratio': props.get('OC_ratio', 0.0),
        'ring_count': props.get('ring_count', 0),
        'vapor_pressure_kPa': props.get('vapor_pressure_kPa', 10.0),
        'T10_C': props.get('T10_C', 50.0),
        'T50_C': props.get('T50_C', 100.0),
        'T90_C': props.get('T90_C', 150.0),
        'cost_eur_L': cost_eur_L,
        'max_vol_frac': max_vol_frac,
        'flags': '',
        'is_novel': is_novel,
        'waste_credit_eur_L': waste_credit_eur_L,
        'SMILES': smiles,  # Store original SMILES for reference
        # Synthesis metadata
        'synthesis_feasibility': round(synth['feasibility'], 1),
        'synthesis_complexity': round(synth['complexity'], 1),
        'synthesis_reaction_type': synthesis_route.get('reaction_type', 'unknown'),
        'synthesis_difficulty': synthesis_route.get('difficulty', 'hard'),
        'synthesis_yield_estimate': synthesis_route.get('yield_estimate', 'unknown'),
        # Stability metadata
        'oxidation_susceptibility': round(oxidation_susceptibility, 2),
    }
    
    # Add stability notes
    stability_notes = []
    if oxidation_susceptibility > 5.0:
        stability_notes.append('High oxidation risk - monitor for color change and degradation')
    elif oxidation_susceptibility > 2.0:
        stability_notes.append('Moderate oxidation risk')
    else:
        stability_notes.append('Low oxidation risk - good control candidate')
    
    if synth['complexity'] > 6.0:
        stability_notes.append('Complex structure - may have stability issues')
    
    component['stability_notes'] = '; '.join(stability_notes)
    
    return component


def molecules_to_component_csv(
    molecules: List[str],
    output_path: str,
    cost_eur_L: float = 2.0,
    max_vol_frac: float = 0.1,
    is_novel: int = 1,
    waste_credit_eur_L: float = 0.0
) -> pd.DataFrame:
    """
    Convert a list of SMILES molecules to a component CSV file.
    
    Automatically deduplicates molecules to ensure unique component IDs.
    
    Args:
        molecules: List of SMILES strings (duplicates will be removed)
        output_path: Path to save the CSV file
        cost_eur_L: Default cost per liter
        max_vol_frac: Default maximum volume fraction
        is_novel: Whether components are novel (0 or 1)
        waste_credit_eur_L: Default waste credit per liter
    
    Returns:
        DataFrame with component data (unique molecules only)
    """
    # Deduplicate molecules while preserving order
    seen = set()
    unique_molecules = []
    for smiles in molecules:
        if smiles not in seen:
            seen.add(smiles)
            unique_molecules.append(smiles)
    
    if len(unique_molecules) < len(molecules):
        print(f"Note: Removed {len(molecules) - len(unique_molecules)} duplicate molecules")
    
    components = []
    
    for smiles in unique_molecules:
        try:
            component = molecule_to_component(
                smiles,
                cost_eur_L=cost_eur_L,
                max_vol_frac=max_vol_frac,
                is_novel=is_novel,
                waste_credit_eur_L=waste_credit_eur_L
            )
            components.append(component)
        except Exception as e:
            print(f"Warning: Failed to convert {smiles}: {e}")
            continue
    
    if not components:
        raise ValueError("No valid components generated")
    
    df = pd.DataFrame(components)
    
    # Verify no duplicate IDs
    if df['id'].duplicated().any():
        # If somehow duplicates exist, add suffix
        duplicates = df[df['id'].duplicated(keep=False)]
        for idx, row in duplicates.iterrows():
            original_id = row['id']
            counter = 1
            while f"{original_id}_{counter}" in df['id'].values:
                counter += 1
            df.at[idx, 'id'] = f"{original_id}_{counter}"
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} components to {output_path}")
    
    return df

