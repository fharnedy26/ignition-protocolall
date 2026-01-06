#!/usr/bin/env python3
"""Check blend results and see if novel molecules are being used."""

import pandas as pd
import sys
from pathlib import Path

def check_blend(run_dir):
    """Check a blend run directory."""
    run_path = Path(run_dir)
    # Try different possible file names
    blends_file = None
    for name in ["blends.csv", "Top-2.csv", "Portfolio.csv"]:
        candidate = run_path / name
        if candidate.exists():
            blends_file = candidate
            break
    
    if blends_file is None:
        print(f"Error: No blend file found in {run_path}")
        print(f"Available files: {list(run_path.glob('*'))}")
        return
    
    df = pd.read_csv(blends_file)
    if len(df) == 0:
        print("No blends found")
        return
    
    blend = df.iloc[0]  # Best blend
    
    print("=" * 70)
    print("Blend Analysis")
    print("=" * 70)
    print(f"\nBest Score: {blend['score']:.4f}")
    print(f"\nComponents Used:")
    
    novel_used = False
    novel_fraction = blend.get('novel_fraction', 0.0)
    
    # Check if components column exists (new format)
    if 'components' in blend:
        components_str = str(blend['components'])
        print(f"  {components_str}")
        novel_used = novel_fraction > 0.001
    else:
        # Old format with frac_ columns
        for col in df.columns:
            if col.startswith('frac_') and blend[col] > 0.001:
                comp_name = col.replace('frac_', '')
                frac = blend[col]
                is_novel = 'novel' if 'MOL_' in comp_name else 'standard'
                if 'MOL_' in comp_name:
                    novel_used = True
                print(f"  {comp_name:<30} {frac:.3f} ({is_novel})")
    
    print(f"\nNovel Fraction: {novel_fraction:.3f}")
    print(f"Novel Penalty: {blend.get('novel_penalty', 0.0):.4f}")
    print(f"Novel Molecules Used: {'YES' if novel_used or novel_fraction > 0.001 else 'NO'}")
    
    # Check if novel components were available
    components_file = run_path / "components.csv"
    if components_file.exists():
        comps_df = pd.read_csv(components_file)
        novel_comps = comps_df[comps_df.get('is_novel', 0) == 1]
        print(f"\nNovel Components Available: {len(novel_comps)}")
        if len(novel_comps) > 0:
            print("\nAvailable Novel Molecules:")
            for idx, row in novel_comps.head(5).iterrows():
                name = str(row.get('name', row.get('id', '')))[:50]
                print(f"  {row.get('id', 'N/A'):<20} RON={row.get('RON', 0):.1f} LHV={row.get('LHV_MJ_kg', 0):.1f} PMI={row.get('PMI', 0):.1f}")
    
    print("=" * 70)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_dir = sys.argv[1]
    else:
        # Find most recent run
        runs_dir = Path("runs")
        if runs_dir.exists():
            runs = sorted(runs_dir.glob("20*"), reverse=True)
            if runs:
                run_dir = runs[0]
            else:
                print("No runs found")
                sys.exit(1)
        else:
            print("Runs directory not found")
            sys.exit(1)
    
    check_blend(run_dir)

