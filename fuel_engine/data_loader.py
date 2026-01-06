"""
Data loading functions for component libraries and baseline blends.
"""

import csv
from pathlib import Path
from typing import Dict

import pandas as pd

from fuel_engine.constants import REQUIRED_DEFAULTS


def load_components(path: str) -> pd.DataFrame:
    """
    Load component library from CSV.
    
    Args:
        path: Path to components CSV file
        
    Returns:
        DataFrame with component data, ensuring all required columns exist
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file cannot be read or is invalid
    """
    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Components file not found: {path}")
    
    if not path_obj.is_file():
        raise ValueError(f"Path is not a file: {path}")
    
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        raise ValueError(f"Components CSV file is empty: {path}")
    except Exception as e:
        raise ValueError(f"Failed to read components CSV: {e}") from e
    
    if df.empty:
        raise ValueError(f"Components CSV file contains no data: {path}")
    
    # Force IDs to string for stable joins and signatures
    if 'id' in df.columns:
        df['id'] = df['id'].astype(str)
    else:
        df['id'] = df.index.astype(str)
    
    # Validate IDs are unique
    if df['id'].duplicated().any():
        duplicates = df[df['id'].duplicated()]['id'].unique()
        raise ValueError(f"Duplicate component IDs found: {', '.join(duplicates[:5])}")
    
    # Ensure all expected columns exist even if missing in the CSV
    for col, default in REQUIRED_DEFAULTS.items():
        if col not in df.columns:
            df[col] = default
        else:
            df[col] = df[col].fillna(default)
    
    # Validate critical numeric columns
    numeric_cols = ['density_g_ml', 'LHV_MJ_kg', 'cost_eur_L', 'max_vol_frac']
    for col in numeric_cols:
        if col in df.columns:
            if (df[col] < 0).any():
                raise ValueError(f"Negative values found in {col} (must be >= 0)")
            if col == 'max_vol_frac' and (df[col] > 1.0).any():
                raise ValueError(f"max_vol_frac values > 1.0 found (must be <= 1.0)")
    
    # Add is_novel flag if not present
    if 'is_novel' not in df.columns:
        df['is_novel'] = 0
    
    # Validate is_novel is binary
    if 'is_novel' in df.columns:
        invalid = df[~df['is_novel'].isin([0, 1])]
        if not invalid.empty:
            raise ValueError(f"is_novel must be 0 or 1, found invalid values")
    
    return df


def load_baseline(path: str) -> Dict[str, float]:
    """
    Load baseline blend from CSV.
    
    Args:
        path: Path to baseline blend CSV file
        
    Returns:
        Dictionary mapping component IDs to volume fractions
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file cannot be read or is invalid
    """
    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Baseline file not found: {path}")
    
    if not path_obj.is_file():
        raise ValueError(f"Path is not a file: {path}")
    
    baseline = {}
    try:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or 'id' not in reader.fieldnames or 'vol_frac' not in reader.fieldnames:
                raise ValueError("Baseline CSV must have 'id' and 'vol_frac' columns")
            
            for row_num, row in enumerate(reader, start=2):  # start=2 because header is row 1
                if not row.get('id') or not row.get('vol_frac'):
                    raise ValueError(f"Row {row_num}: Missing 'id' or 'vol_frac' value")
                
                comp_id = str(row['id']).strip()
                if not comp_id:
                    raise ValueError(f"Row {row_num}: Empty component ID")
                
                try:
                    vol_frac = float(row['vol_frac'])
                except (ValueError, TypeError) as e:
                    raise ValueError(f"Row {row_num}: Invalid vol_frac value '{row['vol_frac']}': {e}") from e
                
                if vol_frac < 0 or vol_frac > 1:
                    raise ValueError(f"Row {row_num}: vol_frac must be between 0 and 1, got {vol_frac}")
                
                if comp_id in baseline:
                    raise ValueError(f"Duplicate component ID in baseline: {comp_id}")
                
                baseline[comp_id] = vol_frac
    except FileNotFoundError:
        raise
    except Exception as e:
        raise ValueError(f"Failed to read baseline CSV: {e}") from e
    
    # Validate that volume fractions sum to approximately 1.0
    total = sum(baseline.values())
    if abs(total - 1.0) > 0.01:  # Allow small tolerance for floating point
        raise ValueError(f"Baseline volume fractions sum to {total:.4f}, expected ~1.0")
    
    return baseline

