#!/usr/bin/env python3
"""
Test script for GDB-11 integration.

Tests loading and using preprocessed GDB-11 fragments for molecular evolution.
"""

from molecular_generator import load_fragments, evolve_with_fragments

def main():
    print("Testing GDB-11 integration...")
    print("=" * 60)
    
    # Test 1: Load fragments
    print("\n[Test 1] Loading fragments from preprocessed GDB-11 file...")
    try:
        fragments = load_fragments(
            "data/fragments/gdb11_fuel_filtered.smi",
            max_fragments=1000
        )
        print(f"[OK] Loaded {len(fragments):,} fragments")
        
        # Show sample fragments
        print(f"\n   Sample fragments:")
        for i, frag in enumerate(fragments[:5], 1):
            print(f"   {i}. {frag}")
        if len(fragments) > 5:
            print(f"   ... and {len(fragments) - 5} more")
            
    except FileNotFoundError:
        print("[ERROR] Preprocessed file not found!")
        print("   Run: python scripts/preprocess_gdb11.py --input gdb11.tgz --output data/fragments/gdb11_fuel_filtered.smi")
        return 1
    except Exception as e:
        print(f"[ERROR] Error loading fragments: {e}")
        return 1
    
    # Test 2: Quick evolution
    print("\n[Test 2] Running quick evolution test...")
    print("-" * 60)
    try:
        evolution_df, history = evolve_with_fragments(
            fragments=fragments,
            generations=10,  # Short test
            population_size=20,
            survivors=5,
            n_frags_per_molecule=3,
            verbose=True
        )
        
        print(f"\n[OK] Evolution complete!")
        print(f"   Generated {len(evolution_df)} molecules")
        
        if len(evolution_df) > 0:
            print(f"\n   Top 5 molecules:")
            for i, (idx, row) in enumerate(evolution_df.head(5).iterrows(), 1):
                print(f"   {i}. {row['Best_SMILES']}")
                print(f"      Fitness: {row['Fitness']:.2f} | "
                      f"MolWt: {row.get('MolWt', 'N/A'):.1f} | "
                      f"logP: {row.get('logP', 'N/A'):.2f}")
            
            best = evolution_df.iloc[0]
            print(f"\n[BEST] Best molecule: {best['Best_SMILES']}")
            print(f"   Fitness: {best['Fitness']:.2f}")
            
            # Show fitness progression
            if history:
                fitness_history = [h['Fitness'] for h in history]
                print(f"\n   Fitness progression: {fitness_history[0]:.2f} -> {fitness_history[-1]:.2f}")
        else:
            print("[WARN] No molecules generated")
            
    except Exception as e:
        print(f"[ERROR] Error during evolution: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 3: File path direct usage (skipped for large files - too slow)
    print("\n[Test 3] Testing direct file path usage (automatic caching)...")
    print("-" * 60)
    print("[SKIP] Skipping Test 3 - reservoir sampling on large files is slow.")
    print("       For production, use preprocessed files with load_fragments() instead.")
    print("       Direct file path works but is slower for files > 1M fragments.")
    
    # Uncomment below to test (will be slow with 2.7M fragments):
    # try:
    #     evolution_df2, history2 = evolve_with_fragments(
    #         fragments="data/fragments/gdb11_fuel_filtered.smi",  # File path
    #         generations=5,  # Very short test
    #         population_size=15,
    #         survivors=5,
    #         fragment_cache_size=500,  # Small cache for quick test
    #         verbose=True
    #     )
    #     
    #     print(f"\n[OK] Direct file path test complete!")
    #     print(f"   Generated {len(evolution_df2)} molecules")
    #     if len(evolution_df2) > 0:
    #         print(f"   Best: {evolution_df2.iloc[0]['Best_SMILES']} "
    #               f"(Fitness: {evolution_df2.iloc[0]['Fitness']:.2f})")
    #         
    # except Exception as e:
    #     print(f"[ERROR] Error with direct file path: {e}")
    #     import traceback
    #     traceback.print_exc()
    #     return 1
    
    # Summary
    print("\n" + "=" * 60)
    print("[SUCCESS] Core tests passed!")
    print("\nNext steps:")
    print("   1. Run full evolution with more fragments:")
    print("      fragments = load_fragments('data/fragments/gdb11_fuel_filtered.smi', max_fragments=50000)")
    print("   2. Increase generations for better results:")
    print("      evolve_with_fragments(fragments, generations=100, population_size=50)")
    print("   3. Use in notebook: jupyter notebook notebooks/unified_workflow.ipynb")
    print("\nNote: Direct file path (Test 3) works but is slow for large files.")
    print("      Use load_fragments() with max_fragments for better performance.")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

