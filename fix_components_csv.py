#!/usr/bin/env python3
"""
Fix duplicate IDs in components CSV file.
"""

import pandas as pd
from pathlib import Path

def fix_duplicates(input_file: str, output_file: str = None):
    """Remove duplicate component IDs from CSV."""
    if output_file is None:
        output_file = input_file
    
    df = pd.read_csv(input_file)
    
    print(f"Original file: {len(df)} rows")
    print(f"Unique IDs: {df['id'].nunique()}")
    
    # Remove duplicates, keeping first occurrence
    df_unique = df.drop_duplicates(subset=['id'], keep='first')
    
    print(f"After deduplication: {len(df_unique)} rows")
    
    df_unique.to_csv(output_file, index=False)
    print(f"Saved fixed file to: {output_file}")
    
    return df_unique

if __name__ == "__main__":
    import sys
    
    input_file = "exports/gdb11_maximum_components.csv"
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    
    fix_duplicates(input_file)

