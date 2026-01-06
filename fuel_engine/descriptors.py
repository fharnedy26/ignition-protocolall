"""
Molecular descriptor calculation for fuel components using RDKit and Mordred.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    Chem = None
    Descriptors = None

try:
    from mordred import Calculator, descriptors
    MORDRED_AVAILABLE = True
except ImportError:
    MORDRED_AVAILABLE = False
    Calculator = None
    descriptors = None


def calculate_basic_descriptors(smiles: str) -> Dict[str, float]:
    """
    Calculate basic molecular descriptors from SMILES string using RDKit.
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        Dictionary of molecular descriptors
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is not available. Install with: pip install rdkit-pypi")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    desc = {}
    
    # Basic descriptors
    desc['molecular_weight'] = Descriptors.MolWt(mol)
    desc['logp'] = Descriptors.MolLogP(mol)
    desc['tpsa'] = Descriptors.TPSA(mol)
    desc['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
    desc['num_aromatic_rings'] = Descriptors.NumAromaticRings(mol)
    desc['num_saturated_rings'] = Descriptors.NumSaturatedRings(mol)
    desc['num_heteroatoms'] = Descriptors.NumHeteroatoms(mol)
    desc['num_heavy_atoms'] = Descriptors.HeavyAtomCount(mol)
    desc['num_carbon_atoms'] = Descriptors.NumCarbonAtoms(mol)
    desc['num_oxygen_atoms'] = Descriptors.NumOxygenAtoms(mol)
    desc['num_nitrogen_atoms'] = Descriptors.NumNitrogenAtoms(mol)
    desc['fraction_csp3'] = Descriptors.FpDensityMorgan1(mol)
    
    # Complexity descriptors
    desc['balaban_j'] = Descriptors.BalabanJ(mol)
    desc['bertz_ct'] = Descriptors.BertzCT(mol)
    desc['chi0'] = Descriptors.Chi0(mol)
    desc['chi1'] = Descriptors.Chi1(mol)
    desc['kappa1'] = Descriptors.Kappa1(mol)
    desc['kappa2'] = Descriptors.Kappa2(mol)
    desc['kappa3'] = Descriptors.Kappa3(mol)
    
    return desc


def calculate_mordred_descriptors(smiles: str) -> Dict[str, float]:
    """
    Calculate comprehensive molecular descriptors using Mordred.
    
    Mordred provides 1800+ molecular descriptors including:
    - Topological descriptors
    - Geometric descriptors
    - Electronic descriptors
    - Constitutional descriptors
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        Dictionary of molecular descriptors (only non-NaN values)
    """
    if not MORDRED_AVAILABLE:
        raise ImportError("Mordred is not available. Install with: pip install mordred")
    
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for Mordred. Install with: pip install rdkit-pypi")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # Create calculator with all descriptors
    calc = Calculator(descriptors, ignore_3D=True)
    
    # Calculate descriptors
    desc_dict = calc(mol)
    
    # Convert to dictionary and filter out NaN/inf values
    result = {}
    for key, value in desc_dict.items():
        if isinstance(value, (int, float)):
            if not (np.isnan(value) or np.isinf(value)):
                result[str(key)] = float(value)
    
    return result


def add_descriptors_to_components(comps: pd.DataFrame, 
                                  smiles_col: str = 'smiles',
                                  use_mordred: bool = True) -> pd.DataFrame:
    """
    Add molecular descriptors to component DataFrame.
    
    This function calculates descriptors for all components that have
    SMILES strings and adds them as new columns.
    
    Args:
        comps: Component DataFrame
        smiles_col: Column name containing SMILES strings
        use_mordred: If True, use Mordred (comprehensive), else use RDKit basic
        
    Returns:
        DataFrame with added descriptor columns
    """
    if smiles_col not in comps.columns:
        raise ValueError(f"Column '{smiles_col}' not found in DataFrame")
    
    comps = comps.copy()
    
    # Get all descriptors for first valid molecule to determine columns
    sample_smiles = None
    for idx, row in comps.iterrows():
        if pd.notna(row[smiles_col]) and str(row[smiles_col]).strip():
            sample_smiles = str(row[smiles_col])
            break
    
    if sample_smiles is None:
        return comps  # No SMILES strings available
    
    # Calculate descriptors for sample
    if use_mordred and MORDRED_AVAILABLE:
        sample_desc = calculate_mordred_descriptors(sample_smiles)
    elif RDKIT_AVAILABLE:
        sample_desc = calculate_basic_descriptors(sample_smiles)
    else:
        raise ImportError("Neither Mordred nor RDKit is available")
    
    # Initialize descriptor columns
    for desc_name in sample_desc.keys():
        if desc_name not in comps.columns:
            comps[desc_name] = np.nan
    
    # Calculate descriptors for all components
    for idx, row in comps.iterrows():
        if pd.notna(row[smiles_col]) and str(row[smiles_col]).strip():
            try:
                smiles = str(row[smiles_col])
                if use_mordred and MORDRED_AVAILABLE:
                    desc = calculate_mordred_descriptors(smiles)
                elif RDKIT_AVAILABLE:
                    desc = calculate_basic_descriptors(smiles)
                else:
                    continue
                
                # Add descriptors to row
                for desc_name, desc_value in desc.items():
                    if desc_name in comps.columns:
                        comps.at[idx, desc_name] = desc_value
            except Exception as e:
                # Skip invalid SMILES, keep NaN values
                continue
    
    return comps


def predict_property_from_descriptors(descriptors: pd.DataFrame,
                                      target_property: str,
                                      method: str = 'linear') -> Optional[np.ndarray]:
    """
    Predict component property from molecular descriptors using simple models.
    
    This is a placeholder for future ML-based property prediction.
    For now, it returns None. In the future, this could use:
    - Linear regression
    - Random forest
    - Neural networks
    - etc.
    
    Args:
        descriptors: DataFrame of molecular descriptors
        target_property: Name of property to predict
        method: Prediction method ('linear', 'rf', 'nn', etc.)
        
    Returns:
        Array of predicted property values, or None if not implemented
    """
    # TODO: Implement ML-based property prediction
    return None


def get_descriptor_availability() -> Dict[str, bool]:
    """Check which descriptor libraries are available."""
    return {
        'rdkit': RDKIT_AVAILABLE,
        'mordred': MORDRED_AVAILABLE,
        'fully_available': RDKIT_AVAILABLE and MORDRED_AVAILABLE
    }

