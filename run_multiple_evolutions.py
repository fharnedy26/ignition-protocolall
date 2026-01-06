#!/usr/bin/env python3
"""
Run multiple independent evolution runs with different seeds.

This explores more of the chemical space by running several independent
evolutionary searches and combining the best results.
"""

from molecular_generator import load_fragments, evolve_with_fragments
from molecular_generator import molecules_to_component_csv
import pandas as pd
from pathlib import Path
import time
import random

def main():
    print("=" * 70)
    print("Multiple Independent Evolution Runs")
    print("=" * 70)
    
    # Load fragments once
    print("\n[Loading] Loading fragments from preprocessed GDB-11 file...")
    try:
        fragments = load_fragments(
            "data/fragments/gdb11_fuel_filtered.smi",
            max_fragments=300000  # Increased from 200K for better diversity
        )
        print(f"[OK] Loaded {len(fragments):,} fragments")
    except Exception as e:
        print(f"[ERROR] {e}")
        return 1
    
    # Run multiple independent evolutions
    num_runs = 10  # Increased from 5 to 10 for better diversity
    all_results = []
    all_histories = []
    
    print(f"\n[Evolution] Running {num_runs} independent evolution runs...")
    print("-" * 70)
    
    for run_id in range(num_runs):
        seed = 42 + run_id * 1000  # Different seed for each run
        random.seed(seed)
        
        print(f"\n[Run {run_id + 1}/{num_runs}] Starting evolution (seed={seed})...")
        start_time = time.time()
        
        try:
            evolution_df, history = evolve_with_fragments(
                fragments=fragments,
                generations=150,  # Increased from 100 for better results
                population_size=75,  # Increased from 50 for more diversity
                survivors=15,  # Increased from 10
                n_frags_per_molecule=3,
                verbose=False  # Less verbose for multiple runs
            )
            
            elapsed = time.time() - start_time
            best_fitness = evolution_df.iloc[0]['Fitness'] if len(evolution_df) > 0 else 0
            unique_count = evolution_df['Best_SMILES'].nunique()
            
            print(f"  [OK] Run {run_id + 1} complete: {len(evolution_df)} molecules, "
                  f"best fitness={best_fitness:.2f}, {unique_count} unique, "
                  f"time={elapsed:.1f}s")
            
            # Add run ID to results
            evolution_df['Run_ID'] = run_id + 1
            all_results.append(evolution_df)
            all_histories.extend(history)
            
        except Exception as e:
            print(f"  [ERROR] Run {run_id + 1} failed: {e}")
            continue
    
    if not all_results:
        print("[ERROR] No successful runs!")
        return 1
    
    # Combine all results
    print("\n[Analysis] Combining results from all runs...")
    print("-" * 70)
    
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Get top molecules across all runs
    top_molecules = combined_df.nlargest(30, 'Fitness')
    unique_top = top_molecules.drop_duplicates(subset=['Best_SMILES'])
    
    print(f"\n   Combined Statistics:")
    print(f"   - Total molecules: {len(combined_df)}")
    print(f"   - Unique molecules: {combined_df['Best_SMILES'].nunique()}")
    print(f"   - Top 30 unique: {len(unique_top)}")
    print(f"   - Best fitness: {top_molecules.iloc[0]['Fitness']:.2f}")
    
    print(f"\n   Top 10 molecules across all runs:")
    print(f"   {'Rank':<6} {'Run':<6} {'SMILES':<45} {'Fitness':<10} {'MolWt':<8} {'logP':<8}")
    print(f"   {'-'*6} {'-'*6} {'-'*45} {'-'*10} {'-'*8} {'-'*8}")
    for i, (idx, row) in enumerate(unique_top.head(10).iterrows(), 1):
        smiles = row['Best_SMILES']
        if len(smiles) > 42:
            smiles = smiles[:39] + "..."
        print(f"   {i:<6} {int(row.get('Run_ID', 0)):<6} {smiles:<45} {row['Fitness']:<10.2f} {row.get('MolWt', 0):<8.1f} {row.get('logP', 0):<8.2f}")
    
    # Save results
    print("\n[Saving] Saving combined results...")
    print("-" * 70)
    
    output_dir = Path("exports")
    output_dir.mkdir(exist_ok=True)
    
    # Save combined results
    combined_df.to_csv(output_dir / "gdb11_multiple_runs_combined.csv", index=False)
    print(f"   [OK] Saved combined results: {output_dir / 'gdb11_multiple_runs_combined.csv'}")
    
    # Save top unique molecules
    unique_top.to_csv(output_dir / "gdb11_multiple_runs_top_unique.csv", index=False)
    print(f"   [OK] Saved top unique molecules: {output_dir / 'gdb11_multiple_runs_top_unique.csv'}")
    
    # Convert top molecules to components
    top_smiles = unique_top.head(20)['Best_SMILES'].tolist()
    try:
        components_df = molecules_to_component_csv(
            molecules=top_smiles,
            output_path=str(output_dir / "gdb11_multiple_runs_components.csv"),
            cost_eur_L=2.0,
            max_vol_frac=0.15,
            is_novel=1,
            waste_credit_eur_L=0.1
        )
        print(f"   [OK] Converted {len(components_df)} molecules to components")
        print(f"   [OK] Saved components: {output_dir / 'gdb11_multiple_runs_components.csv'}")
        
        # Show component summary
        print(f"\n   Top 10 component summary:")
        print(f"   {'ID':<20} {'Name':<35} {'RON':<8} {'LHV':<8} {'PMI':<8}")
        print(f"   {'-'*20} {'-'*35} {'-'*8} {'-'*8} {'-'*8}")
        for idx, row in components_df.head(10).iterrows():
            name = str(row.get('name', row.get('id', '')))[:33]
            print(f"   {str(row.get('id', ''))[:18]:<20} {name:<35} {row.get('RON', 0):<8.1f} {row.get('LHV_MJ_kg', 0):<8.1f} {row.get('PMI', 0):<8.1f}")
            
    except Exception as e:
        print(f"   [WARN] Could not convert to components: {e}")
    
    # Summary
    print("\n" + "=" * 70)
    print("[SUCCESS] Multiple runs complete!")
    print("=" * 70)
    print(f"\nResults:")
    print(f"  - Runs completed: {len(all_results)}")
    print(f"  - Total molecules: {len(combined_df)}")
    print(f"  - Unique molecules: {combined_df['Best_SMILES'].nunique()}")
    print(f"  - Best fitness: {top_molecules.iloc[0]['Fitness']:.2f}")
    print(f"\nOutput files:")
    print(f"  - {output_dir / 'gdb11_multiple_runs_combined.csv'}")
    print(f"  - {output_dir / 'gdb11_multiple_runs_top_unique.csv'}")
    print(f"  - {output_dir / 'gdb11_multiple_runs_components.csv'}")
    print(f"\nNext steps:")
    print(f"  1. Review unique molecules in exports/gdb11_multiple_runs_top_unique.csv")
    print(f"  2. Use components in blend optimizer:")
    print(f"     python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

