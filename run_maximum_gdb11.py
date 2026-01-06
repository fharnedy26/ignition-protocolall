#!/usr/bin/env python3
"""
Maximum-scale production run using GDB-11 fragments.

Runs a comprehensive molecular evolution with maximum parameters for best results.
"""

from molecular_generator import load_fragments, evolve_with_fragments
from molecular_generator import molecules_to_component_csv
import pandas as pd
from pathlib import Path
import time

def main():
    start_time = time.time()
    
    print("=" * 70)
    print("GDB-11 MAXIMUM-SCALE Production Run")
    print("=" * 70)
    
    # Load fragments - use large sample for maximum diversity
    print("\n[1/4] Loading fragments from preprocessed GDB-11 file...")
    print("-" * 70)
    try:
        fragments = load_fragments(
            "data/fragments/gdb11_fuel_filtered.smi",
            max_fragments=500000  # Use 500K fragments for maximum diversity (increased from 200K)
        )
        print(f"[OK] Loaded {len(fragments):,} fragments")
        print(f"   Sample: {fragments[:5]}")
        
    except FileNotFoundError:
        print("[ERROR] Preprocessed file not found!")
        print("   Run: python scripts/preprocess_gdb11.py --input gdb11.tgz --output data/fragments/gdb11_fuel_filtered.smi")
        return 1
    except Exception as e:
        print(f"[ERROR] Error loading fragments: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Run evolution with maximum parameters
    print("\n[2/4] Running molecular evolution (MAXIMUM SCALE)...")
    print("-" * 70)
    print("   Parameters:")
    print("   - Generations: 500")
    print("   - Population size: 100")
    print("   - Survivors: 20")
    print("   - Fragments per molecule: 3")
    print("   - Fragment pool: 200,000")
    print("\n   This will take several minutes...")
    print("   Starting evolution...\n")
    
    evolution_start = time.time()
    try:
        evolution_df, history = evolve_with_fragments(
            fragments=fragments,
            generations=500,  # Maximum generations
            population_size=100,  # Larger population
            survivors=20,  # More survivors for diversity
            n_frags_per_molecule=3,
            verbose=True
        )
        
        evolution_time = time.time() - evolution_start
        print(f"\n[OK] Evolution complete!")
        print(f"   Generated {len(evolution_df)} molecules")
        print(f"   Evolution time: {evolution_time:.1f} seconds ({evolution_time/60:.1f} minutes)")
        
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
        top_molecules = evolution_df.nlargest(30, 'Fitness')
        
        print(f"\n   Top 15 molecules:")
        print(f"   {'Rank':<6} {'SMILES':<45} {'Fitness':<10} {'MolWt':<8} {'logP':<8}")
        print(f"   {'-'*6} {'-'*45} {'-'*10} {'-'*8} {'-'*8}")
        for i, (idx, row) in enumerate(top_molecules.head(15).iterrows(), 1):
            smiles = row['Best_SMILES']
            if len(smiles) > 42:
                smiles = smiles[:39] + "..."
            print(f"   {i:<6} {smiles:<45} {row['Fitness']:<10.2f} {row.get('MolWt', 0):<8.1f} {row.get('logP', 0):<8.2f}")
        
        best = evolution_df.iloc[0]
        print(f"\n   [BEST] Best molecule: {best['Best_SMILES']}")
        print(f"   Fitness: {best['Fitness']:.2f}")
        print(f"   Molecular Weight: {best.get('MolWt', 'N/A'):.1f}")
        print(f"   logP: {best.get('logP', 'N/A'):.2f}")
        print(f"   H-Donors: {best.get('H_Donors', 'N/A')}")
        print(f"   Rings: {best.get('Rings', 'N/A')}")
        print(f"   O/C Ratio: {best.get('O/C', 'N/A'):.3f}")
        
        # Fitness progression analysis
        if history:
            fitness_history = [h['Fitness'] for h in history]
            print(f"\n   Fitness progression:")
            print(f"   Start: {fitness_history[0]:.2f}")
            print(f"   End:   {fitness_history[-1]:.2f}")
            print(f"   Improvement: {fitness_history[-1] - fitness_history[0]:.2f}")
            print(f"   Best generation: {fitness_history.index(max(fitness_history)) + 1}")
            print(f"   Max fitness: {max(fitness_history):.2f}")
            
            # Show improvement over time
            improvements = []
            for i in range(0, len(fitness_history), 50):
                improvements.append((i+1, fitness_history[i]))
            if len(fitness_history) % 50 != 0:
                improvements.append((len(fitness_history), fitness_history[-1]))
            
            print(f"\n   Fitness milestones:")
            for gen, fit in improvements:
                print(f"   Gen {gen:3d}: {fit:.2f}")
        
        # Diversity analysis
        unique_smiles = evolution_df['Best_SMILES'].nunique()
        print(f"\n   Diversity:")
        print(f"   Unique molecules: {unique_smiles} out of {len(evolution_df)} total")
        print(f"   Diversity ratio: {unique_smiles/len(evolution_df)*100:.1f}%")
        
        # Save results
        print("\n[4/4] Saving results...")
        print("-" * 70)
        
        # Save evolution history
        output_dir = Path("exports")
        output_dir.mkdir(exist_ok=True)
        
        evolution_df.to_csv(output_dir / "gdb11_maximum_evolution_results.csv", index=False)
        print(f"   [OK] Saved evolution results to: {output_dir / 'gdb11_maximum_evolution_results.csv'}")
        
        # Save top 30 molecules
        top_30 = top_molecules.head(30)
        top_30.to_csv(output_dir / "gdb11_top30_molecules.csv", index=False)
        print(f"   [OK] Saved top 30 molecules to: {output_dir / 'gdb11_top30_molecules.csv'}")
        
        # Convert top molecules to components (deduplicate first)
        top_smiles = top_molecules.head(20)['Best_SMILES'].tolist()
        # Remove duplicates while preserving order
        seen = set()
        unique_smiles = []
        for smiles in top_smiles:
            if smiles not in seen:
                seen.add(smiles)
                unique_smiles.append(smiles)
        
        try:
            components_df = molecules_to_component_csv(
                molecules=unique_smiles,  # Use deduplicated list
                output_path=str(output_dir / "gdb11_maximum_components.csv"),
                cost_eur_L=2.0,
                max_vol_frac=0.15,  # Allow up to 15% in blend
                is_novel=1,
                waste_credit_eur_L=0.1
            )
            print(f"   [OK] Converted {len(components_df)} molecules to components")
            print(f"   [OK] Saved components to: {output_dir / 'gdb11_maximum_components.csv'}")
            
            # Show component summary
            print(f"\n   Top 10 component summary:")
            print(f"   {'ID':<20} {'Name':<35} {'RON':<8} {'LHV':<8} {'PMI':<8} {'Density':<8}")
            print(f"   {'-'*20} {'-'*35} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
            for idx, row in components_df.head(10).iterrows():
                name = str(row.get('name', row.get('id', '')))[:33]
                print(f"   {str(row.get('id', ''))[:18]:<20} {name:<35} {row.get('RON', 0):<8.1f} {row.get('LHV_MJ_kg', 0):<8.1f} {row.get('PMI', 0):<8.1f} {row.get('density_g_ml', 0):<8.3f}")
            
        except Exception as e:
            print(f"   [WARN] Could not convert to components: {e}")
            import traceback
            traceback.print_exc()
        
        # Total time
        total_time = time.time() - start_time
        
        # Summary
        print("\n" + "=" * 70)
        print("[SUCCESS] Maximum-scale run complete!")
        print("=" * 70)
        print(f"\nRun Statistics:")
        print(f"  - Fragments used: {len(fragments):,}")
        print(f"  - Generations: 500")
        print(f"  - Population size: 100")
        print(f"  - Molecules generated: {len(evolution_df)}")
        print(f"  - Unique molecules: {unique_smiles}")
        print(f"  - Best fitness: {best['Fitness']:.2f}")
        print(f"  - Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print(f"  - Evolution time: {evolution_time:.1f} seconds ({evolution_time/60:.1f} minutes)")
        print(f"\nOutput files:")
        print(f"  - {output_dir / 'gdb11_maximum_evolution_results.csv'}")
        print(f"  - {output_dir / 'gdb11_top30_molecules.csv'}")
        print(f"  - {output_dir / 'gdb11_maximum_components.csv'}")
        print(f"\nNext steps:")
        print(f"  1. Review results in exports/gdb11_maximum_evolution_results.csv")
        print(f"  2. Use components in blend optimizer:")
        print(f"     python engine.py --components components.csv --novel exports/gdb11_maximum_components.csv --allow-novel 1 --max-add-novel 0.15")
        print(f"  3. Compare with previous runs to see improvement")
        print("=" * 70)
        
    else:
        print("[WARN] No molecules generated")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

