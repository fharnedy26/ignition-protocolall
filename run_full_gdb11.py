#!/usr/bin/env python3
"""
Full production run using GDB-11 fragments.

Runs a complete molecular evolution with substantial parameters.
"""

from molecular_generator import load_fragments, evolve_with_fragments
from molecular_generator import molecules_to_component_csv
import pandas as pd
from pathlib import Path

def main():
    print("=" * 70)
    print("GDB-11 Full Production Run")
    print("=" * 70)
    
    # Load fragments
    print("\n[1/4] Loading fragments from preprocessed GDB-11 file...")
    print("-" * 70)
    try:
        fragments = load_fragments(
            "data/fragments/gdb11_fuel_filtered.smi",
            max_fragments=50000  # Use 50K fragments for good diversity
        )
        print(f"[OK] Loaded {len(fragments):,} fragments")
        print(f"   Sample: {fragments[:3]}")
        
    except FileNotFoundError:
        print("[ERROR] Preprocessed file not found!")
        print("   Run: python scripts/preprocess_gdb11.py --input gdb11.tgz --output data/fragments/gdb11_fuel_filtered.smi")
        return 1
    except Exception as e:
        print(f"[ERROR] Error loading fragments: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Run evolution
    print("\n[2/4] Running molecular evolution...")
    print("-" * 70)
    print("   Parameters:")
    print("   - Generations: 100")
    print("   - Population size: 50")
    print("   - Survivors: 10")
    print("   - Fragments per molecule: 3")
    print("\n   Starting evolution...\n")
    
    try:
        evolution_df, history = evolve_with_fragments(
            fragments=fragments,
            generations=100,
            population_size=50,
            survivors=10,
            n_frags_per_molecule=3,
            verbose=True
        )
        
        print(f"\n[OK] Evolution complete!")
        print(f"   Generated {len(evolution_df)} molecules")
        
    except Exception as e:
        print(f"[ERROR] Error during evolution: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Analyze results
    print("\n[3/4] Analyzing results...")
    print("-" * 70)
    
    if len(evolution_df) > 0:
        # Get top molecules
        top_molecules = evolution_df.nlargest(20, 'Fitness')
        
        print(f"\n   Top 10 molecules:")
        print(f"   {'Rank':<6} {'SMILES':<40} {'Fitness':<10} {'MolWt':<8} {'logP':<8}")
        print(f"   {'-'*6} {'-'*40} {'-'*10} {'-'*8} {'-'*8}")
        for i, (idx, row) in enumerate(top_molecules.head(10).iterrows(), 1):
            smiles = row['Best_SMILES']
            if len(smiles) > 37:
                smiles = smiles[:34] + "..."
            print(f"   {i:<6} {smiles:<40} {row['Fitness']:<10.2f} {row.get('MolWt', 0):<8.1f} {row.get('logP', 0):<8.2f}")
        
        best = evolution_df.iloc[0]
        print(f"\n   [BEST] Best molecule: {best['Best_SMILES']}")
        print(f"   Fitness: {best['Fitness']:.2f}")
        print(f"   Molecular Weight: {best.get('MolWt', 'N/A'):.1f}")
        print(f"   logP: {best.get('logP', 'N/A'):.2f}")
        print(f"   H-Donors: {best.get('H_Donors', 'N/A')}")
        print(f"   Rings: {best.get('Rings', 'N/A')}")
        print(f"   O/C Ratio: {best.get('O/C', 'N/A'):.3f}")
        
        # Fitness progression
        if history:
            fitness_history = [h['Fitness'] for h in history]
            print(f"\n   Fitness progression:")
            print(f"   Start: {fitness_history[0]:.2f}")
            print(f"   End:   {fitness_history[-1]:.2f}")
            print(f"   Improvement: {fitness_history[-1] - fitness_history[0]:.2f}")
            print(f"   Best generation: {fitness_history.index(max(fitness_history)) + 1}")
        
        # Save results
        print("\n[4/4] Saving results...")
        print("-" * 70)
        
        # Save evolution history
        output_dir = Path("exports")
        output_dir.mkdir(exist_ok=True)
        
        evolution_df.to_csv(output_dir / "gdb11_evolution_results.csv", index=False)
        print(f"   [OK] Saved evolution results to: {output_dir / 'gdb11_evolution_results.csv'}")
        
        # Convert top molecules to components
        top_smiles = top_molecules.head(10)['Best_SMILES'].tolist()
        try:
            components_df = molecules_to_component_csv(
                molecules=top_smiles,
                output_path=str(output_dir / "gdb11_components.csv"),
                cost_eur_L=2.0,
                max_vol_frac=0.1,
                is_novel=1,
                waste_credit_eur_L=0.1
            )
            print(f"   [OK] Converted {len(components_df)} molecules to components")
            print(f"   [OK] Saved components to: {output_dir / 'gdb11_components.csv'}")
            
            # Show component summary
            print(f"\n   Component summary:")
            print(f"   {'ID':<20} {'Name':<30} {'RON':<8} {'LHV':<8} {'PMI':<8}")
            print(f"   {'-'*20} {'-'*30} {'-'*8} {'-'*8} {'-'*8}")
            for idx, row in components_df.head(5).iterrows():
                name = row.get('name', row.get('id', ''))[:28]
                print(f"   {row.get('id', '')[:18]:<20} {name:<30} {row.get('RON', 0):<8.1f} {row.get('LHV_MJ_kg', 0):<8.1f} {row.get('PMI', 0):<8.1f}")
            
        except Exception as e:
            print(f"   [WARN] Could not convert to components: {e}")
        
        # Summary
        print("\n" + "=" * 70)
        print("[SUCCESS] Full run complete!")
        print("=" * 70)
        print(f"\nResults:")
        print(f"  - Fragments used: {len(fragments):,}")
        print(f"  - Generations: 100")
        print(f"  - Molecules generated: {len(evolution_df)}")
        print(f"  - Best fitness: {best['Fitness']:.2f}")
        print(f"\nOutput files:")
        print(f"  - {output_dir / 'gdb11_evolution_results.csv'}")
        print(f"  - {output_dir / 'gdb11_components.csv'}")
        print(f"\nNext steps:")
        print(f"  1. Review results in exports/gdb11_evolution_results.csv")
        print(f"  2. Use components in blend optimizer:")
        print(f"     python engine.py --components components.csv --novel exports/gdb11_components.csv --allow-novel 1")
        print(f"  3. Run longer evolution (200+ generations) for better results")
        print("=" * 70)
        
    else:
        print("[WARN] No molecules generated")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

