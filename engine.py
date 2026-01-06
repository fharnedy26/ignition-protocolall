#!/usr/bin/env python3
"""
Minimal Fuel Blend Optimisation Engine
Deterministic, rapid, laboratory-ready blend optimisation from fixed component library
"""

import argparse
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    Parallel = None
    delayed = None

from fuel_engine import (
    greedy_then_refine,
    load_components,
    load_baseline,
    blend_props,
    project_simplex_with_caps,
    write_top_n_csv,
    write_portfolio_csv,
    write_run_json,
)
from fuel_engine.constants import SPECS


def selftest():
    """Run self-test with toy data."""
    print("Running self-test...")
    
    # Create toy components
    toy_comps = pd.DataFrame({
        'id': ['A', 'B', 'C', 'D', 'E'],
        'name': ['Gasoline', 'Ethanol', 'MTBE', 'Toluene', 'Isooctane'],
        'density_g_ml': [0.75, 0.79, 0.74, 0.87, 0.69],
        'LHV_MJ_kg': [44.0, 26.8, 35.2, 40.6, 44.3],
        'RON': [87, 108, 118, 120, 100],
        'MON': [82, 89, 101, 105, 100],
        'cetane': [0, 0, 0, 0, 0],
        'PMI': [15, 5, 10, 25, 12],
        'TSI': [15, 5, 10, 25, 12],
        'OC_ratio': [0, 0.5, 0.2, 0, 0],
        'ring_count': [0, 0, 0, 1, 0],
        'vapor_pressure_kPa': [50, 60, 40, 30, 45],
        'T10_C': [45, 50, 40, 60, 42],
        'T50_C': [95, 100, 90, 110, 92],
        'T90_C': [170, 180, 160, 200, 165],
        'cost_eur_L': [1.0, 1.2, 1.5, 1.8, 1.1],
        'max_vol_frac': [1.0, 0.15, 0.20, 0.30, 1.0],
        'flags': ['', '', '', '', ''],
        'is_novel': [0, 0, 0, 0, 0]
    })
    
    # Test 1: Basic optimisation
    results = greedy_then_refine(toy_comps, K=3, fuel_type='gasoline')
    assert len(results) > 0, "No results generated"
    assert results[0]['score'] > 0, "Negative score"
    print("[OK] Basic optimisation test passed")
    
    # Test 2: Determinism
    results1 = greedy_then_refine(toy_comps, K=3, fuel_type='gasoline')
    results2 = greedy_then_refine(toy_comps, K=3, fuel_type='gasoline')
    assert results1[0]['score'] == results2[0]['score'], "Non-deterministic results"
    print("[OK] Determinism test passed")
    
    # Test 3: Component caps (relaxed for now)
    blend = results[0]
    total_vol = sum(comp['vol_frac'] for comp in blend['components'])
    assert abs(total_vol - 1.0) < 1e-6, f"Volume fractions don't sum to 1: {total_vol}"
    print("[OK] Volume fraction sum test passed")
    
    # Test 4: Spec flags
    flags = blend['flags']
    assert any(flags.values()), "No spec flags generated"
    print("[OK] Spec flags test passed")
    
    # Spec rejection: craft a too-high RVP component and assert strict rejection
    too_rvp = toy_comps.copy()
    too_rvp.loc[0,'vapor_pressure_kPa'] = 200.0
    res = greedy_then_refine(too_rvp, K=2, fuel_type='gasoline', strict=True)
    assert len(res) >= 1
    assert res[0]['properties']['RVP'] <= SPECS['gasoline']['RVP_max'] + 1e-9, "Strict spec not enforced"

    # Novel cap test
    nov = toy_comps.copy()
    nov.loc[1,'is_novel'] = 1
    res = greedy_then_refine(nov, K=3, fuel_type='gasoline', strict=True,
                             allow_novel=True, max_novel=0.03)
    assert res[0]['novel_fraction'] <= 0.0300001, "Novel cap exceeded"

    # Delta budget test
    baseline = {'A':0.6,'E':0.4}
    res = greedy_then_refine(toy_comps, baseline=baseline, K=3, fuel_type='gasoline',
                             strict=True, delta_budget=0.05)
    x0 = np.array([baseline.get(i,0.0) for i in toy_comps['id']])
    xv = np.array([next((c['vol_frac'] for c in res[0]['components'] if c['id']==i),0.0)
                   for i in toy_comps['id']])
    assert np.sum(np.abs(xv - x0)) <= 0.0500001, "Delta budget violated"
    print("[OK] Extended tests passed")
    
    # Jet heuristic smoke test (touch jet path and starts)
    jet_comps = toy_comps.copy()
    # Nudge PMI/density to be in jet range
    jet_comps['PMI'] = [18, 10, 12, 22, 14]
    jet_comps['density_g_ml'] = [0.80, 0.78, 0.79, 0.83, 0.81]
    jet_results = greedy_then_refine(jet_comps, K=3, fuel_type='jet', strict=True)
    assert len(jet_results) >= 1, "No jet results generated"
    assert len(jet_results[0]['components']) <= 3 + 1e-9, "Jet result exceeds K components"
    print("[OK] Jet heuristic test passed")

    # Enforce ≤K components in generic path too
    assert len(results[0]['components']) <= 3 + 1e-9, "Result exceeds K components"
    print("[OK] K-sparsify assertion passed")
    
    # Relative score smoke check (baseline reference)
    baseline = {'A':0.6,'E':0.4}
    res = greedy_then_refine(toy_comps, baseline=baseline, K=3, fuel_type='gasoline',
                             strict=True, delta_budget=0.05, reference_mode='baseline')
    assert 'relative_score' in res[0], "Relative score missing"
    assert np.isfinite(res[0]['relative_score']), "Relative score not finite"
    print("[OK] Relative score test passed")
    
    print("[OK] All self-tests passed")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(description="Minimal Fuel Blend Optimisation Engine")
    
    # Required (unless selftest)
    parser.add_argument('--components', help='Components CSV file')
    
    # Optional
    parser.add_argument('--baseline', help='Baseline blend CSV file')
    parser.add_argument('--K', type=int, default=5, help='Max components in blend')
    parser.add_argument('--top', type=int, default=15, help='Number of top results')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--out', default='runs', help='Output directory')
    
    # Novel components
    parser.add_argument('--allow-novel', type=int, default=0, choices=[0, 1], help='Allow novel components')
    parser.add_argument('--novel', help='Novel components CSV file')
    parser.add_argument('--max-add-novel', type=float, default=0.03, help='Max novel fraction')
    
    # Delta mode
    parser.add_argument('--delta-budget', type=float, default=0.10, help='L1 deviation budget')
    
    # Fuel type and specs
    parser.add_argument('--fuel-type', default='gasoline', choices=['gasoline', 'diesel', 'jet'], help='Fuel type')
    parser.add_argument('--strict-spec', type=int, default=1, choices=[0, 1], help='Strict specifications')
    
    # Reference normalization
    parser.add_argument('--reference-mode', default='auto',
                        choices=['auto','baseline','cheapest','topkey'],
                        help='How to pick the reference "current option" for relative scoring')
    parser.add_argument('--reference-M', type=int, default=3,
                        help='Number of components to use in cheapest/topkey equal-split reference')
    
    # Presets and manufacturability
    parser.add_argument('--preset', default='novel3', choices=['baseline','delta','novel3'],
                        help='Preset configuration for convenience')
    parser.add_argument('--min-component-frac', type=float, default=0.005,
                        help='Minimum volume fraction per component to avoid trace amounts')

    # Exploration controls
    parser.add_argument('--restarts', type=int, default=0)
    parser.add_argument('--min-l1-distance', type=float, default=0.02)
    parser.add_argument('--anneal-prob', type=float, default=0.05)
    parser.add_argument('--anneal-temp', type=float, default=0.01)
    parser.add_argument('--theme', default='none', choices=['none','high_ron','low_pmi','high_lhv','low_cost','all'])
    parser.add_argument('--theme-min', type=int, default=2, help='Min feasible results to keep per theme when sweeping')
    parser.add_argument('--time-budget', type=float, default=0.0, help='Overall time budget in seconds (0=off)')
    parser.add_argument('--max-evals', type=int, default=0, help='Max objective evaluations per start (0=unlimited)')
    parser.add_argument('--pmi-cap', type=float)
    parser.add_argument('--cost-cap', type=float)
    parser.add_argument('--rvp-cap', type=float)
    parser.add_argument('--t90-cap', type=float)
    parser.add_argument('--k-explore', type=int, default=0, choices=[0,1])
    parser.add_argument('--k-max', type=int, default=7)
    parser.add_argument('--novel-screen', type=int, default=0, choices=[0,1], help='Enable 1% novel screening pass')
    parser.add_argument('--waste-credit-default', type=float, default=0.0, help='Default waste credit for novel comps if missing')
    parser.add_argument('--clean-runs', type=int, help='Delete runs older than N days after finishing')
    
    # Other
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1], help='Verbose output')
    parser.add_argument('--selftest', action='store_true', help='Run self-test')
    
    args = parser.parse_args()
    
    if args.selftest:
        selftest()
        return 0
    
    # Check required arguments
    if not args.components:
        parser.error("--components is required")
    
    # Validate components file exists before proceeding
    comps_path = Path(args.components)
    if not comps_path.exists():
        # Try alternative location
        alt_path = Path("data") / comps_path.name
        if alt_path.exists():
            print(f"[INFO] Found components file at alternative location: {alt_path}")
            args.components = str(alt_path)
        else:
            print(f"[ERROR] Components file not found: {args.components}")
            print(f"[HINT] Checked locations:")
            print(f"  - {comps_path.absolute()}")
            print(f"  - {alt_path.absolute()}")
            print(f"[HINT] Ensure the file exists or use --components with the correct path")
            return 1
    
    # Validate argument combinations
    if args.allow_novel and not args.novel:
        print("[WARN] --allow-novel is set but --novel file not specified. Novel components will be empty.")
    
    if args.baseline:
        baseline_path = Path(args.baseline)
        if not baseline_path.exists():
            alt_baseline = Path("data") / baseline_path.name
            if alt_baseline.exists():
                print(f"[INFO] Found baseline file at alternative location: {alt_baseline}")
                args.baseline = str(alt_baseline)
            else:
                print(f"[ERROR] Baseline file not found: {args.baseline}")
                print(f"[HINT] Checked locations:")
                print(f"  - {baseline_path.absolute()}")
                print(f"  - {alt_baseline.absolute()}")
                return 1
        if args.delta_budget is None:
            print("[WARN] --baseline specified but --delta-budget not set. Using default 0.10")
    
    if args.novel:
        novel_path = Path(args.novel)
        if not novel_path.exists():
            alt_novel = Path("data") / novel_path.name
            if alt_novel.exists():
                print(f"[INFO] Found novel components file at alternative location: {alt_novel}")
                args.novel = str(alt_novel)
            else:
                print(f"[WARN] Novel components file not found: {args.novel}")
                print(f"[HINT] Checked locations:")
                print(f"  - {novel_path.absolute()}")
                print(f"  - {alt_novel.absolute()}")
                print(f"[HINT] Continuing without novel components...")
                args.novel = None
                args.allow_novel = 0
    
    if args.min_component_frac < 0 or args.min_component_frac > 1:
        parser.error(f"--min-component-frac must be between 0 and 1, got {args.min_component_frac}")
    
    if args.max_add_novel < 0 or args.max_add_novel > 1:
        parser.error(f"--max-add-novel must be between 0 and 1, got {args.max_add_novel}")
    
    if args.delta_budget is not None and (args.delta_budget < 0 or args.delta_budget > 2):
        parser.error(f"--delta-budget must be between 0 and 2, got {args.delta_budget}")
    
    # Apply preset defaults unless explicitly overridden on CLI
    argv = sys.argv[1:]
    if args.preset == 'baseline':
        if '--strict-spec' not in argv:
            args.strict_spec = 1
        if '--allow-novel' not in argv:
            args.allow_novel = 0
        if '--delta-budget' not in argv:
            args.delta_budget = 0.0
    elif args.preset == 'delta':
        if '--strict-spec' not in argv:
            args.strict_spec = 1
        if '--allow-novel' not in argv:
            args.allow_novel = 0
        if '--delta-budget' not in argv:
            args.delta_budget = 0.05
    elif args.preset == 'novel3':
        if '--strict-spec' not in argv:
            args.strict_spec = 1
        if '--allow-novel' not in argv:
            args.allow_novel = 1
        if '--max-add-novel' not in argv:
            args.max_add_novel = 0.03

    # Set random seed for reproducibility
    np.random.seed(args.seed)
    
    start_time = time.time()
    
    # Load components with error handling
    try:
        comps = load_components(args.components)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        print(f"[HINT] Verify the file path is correct and the file exists")
        print(f"[HINT] Common locations: './components.csv' or './data/components.csv'")
        return 1
    except ValueError as e:
        print(f"[ERROR] Invalid components file: {e}")
        print(f"[HINT] Check that the CSV file:")
        print(f"  - Is not empty")
        print(f"  - Has required columns: id, density_g_ml, LHV_MJ_kg, RON, MON, etc.")
        print(f"  - Is properly formatted (no syntax errors)")
        print(f"[HINT] See README.md for required column format")
        return 1
    except Exception as e:
        print(f"[ERROR] Unexpected error loading components: {e}")
        print(f"[HINT] This may indicate a corrupted file or missing dependencies")
        print(f"[HINT] Try: pip install -r requirements.txt")
        import traceback
        if args.verbose:
            traceback.print_exc()
        return 1
    
    # Apply default waste credit for novel comps if requested
    if 'waste_credit_eur_L' not in comps.columns:
        comps['waste_credit_eur_L'] = 0.0
    if args.waste_credit_default and 'is_novel' in comps.columns:
        mask_novel = comps['is_novel'] == 1
        comps.loc[mask_novel, 'waste_credit_eur_L'] = comps.loc[mask_novel, 'waste_credit_eur_L'].fillna(0.0) + float(args.waste_credit_default)

    # Load novel components if specified
    if args.allow_novel and args.novel:
        try:
            novel_comps = load_components(args.novel)
        except FileNotFoundError as e:
            print(f"[ERROR] Novel components file not found: {e}")
            print(f"[HINT] Verify the file path or disable --allow-novel")
            return 1
        except ValueError as e:
            print(f"[ERROR] Invalid novel components file: {e}")
            print(f"[HINT] Novel components must have the same format as regular components")
            print(f"[HINT] Ensure required columns are present and file is properly formatted")
            return 1
        except Exception as e:
            print(f"[ERROR] Unexpected error loading novel components: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            return 1
        novel_comps['is_novel'] = 1
        if 'waste_credit_eur_L' not in novel_comps.columns:
            novel_comps['waste_credit_eur_L'] = 0.0
        if args.waste_credit_default:
            novel_comps['waste_credit_eur_L'] = novel_comps['waste_credit_eur_L'].fillna(0.0) + float(args.waste_credit_default)
        comps = pd.concat([comps, novel_comps], ignore_index=True)
    
    # Fuel-type-aware filtering (prefer explicit fuel_class; otherwise heuristic)
    if 'fuel_class' in comps.columns:
        fc = args.fuel_type.lower()
        mask = comps['fuel_class'].astype(str).str.contains(fc, case=False, na=False)
        comps = comps[mask].reset_index(drop=True)
    else:
        if args.fuel_type == 'gasoline':
            bad = comps['id'].astype(str).str.contains('diesel|kerosene|biodiesel', case=False, na=False) | \
                  comps['name'].astype(str).str.contains('Diesel|Kerosene|Biodiesel', case=False, na=False)
            comps = comps[~bad].reset_index(drop=True)
    
    # Load baseline if specified
    baseline = None
    if args.baseline:
        try:
            baseline = load_baseline(args.baseline)
        except FileNotFoundError as e:
            print(f"[ERROR] {e}")
            return 1
        except ValueError as e:
            print(f"[ERROR] Invalid baseline file: {e}")
            return 1
        except Exception as e:
            print(f"[ERROR] Unexpected error loading baseline: {e}")
            return 1
    
    # Warn if baseline contains IDs not in the components library
    if baseline is not None:
        lib_ids = set(pd.Series(comps['id']).astype(str).tolist())
        base_ids = set(baseline.keys())
        missing = base_ids - lib_ids
        if missing:
            print(f"[WARN] {len(missing)} baseline component(s) not found in components CSV: "
                  f"{', '.join(sorted(list(missing))[:5])}"
                  f"{' ...' if len(missing) > 5 else ''}")
    
    # Compute baseline properties for delta reporting (if provided)
    baseline_props = None
    if baseline is not None:
        ids = list(comps['id'].astype(str))
        x0 = np.array([baseline.get(i, 0.0) for i in ids])
        x0 = project_simplex_with_caps(x0, comps['max_vol_frac'].values)
        baseline_props = blend_props(x0, comps)

    # Novel 1% screening (optional)
    screened_novels: List[str] = []
    if bool(args.novel_screen) and bool(args.allow_novel):
        scr = greedy_then_refine(
            comps, baseline=baseline, K=min(args.K, 3), fuel_type=args.fuel_type,
            strict=bool(args.strict_spec), allow_novel=True, max_novel=0.01,
            delta_budget=args.delta_budget if baseline else None,
            min_component_frac=float(args.min_component_frac), theme='high_ron',
            anneal_prob=0.0, anneal_temp=0.0, pmi_cap=args.pmi_cap, cost_cap=args.cost_cap,
            rvp_cap=args.rvp_cap, t90_cap=args.t90_cap, k_explore=False, k_max=int(args.k_max),
            num_random_starts=1, tag_theme='novel_screen', tag_restart=-1,
            time_deadline=(time.time() + min(2.0, args.time_budget)) if args.time_budget else None,
            max_evals=min(int(args.max_evals or 200), 200),
            reference_mode=args.reference_mode,
            reference_M=args.reference_M
        )
        if scr:
            for c in scr[0]['components']:
                if any((comps['id'].astype(str) == str(c['id'])) & (comps['is_novel'] == 1)):
                    screened_novels.append(str(c['id']))

    # Orchestrate restarts and theme sweeps
    themes = []
    if args.theme == 'all':
        themes = ['high_ron','low_pmi','high_lhv','low_cost']
    else:
        themes = [args.theme]

    aggregated: List[Dict] = []
    restart_count = max(1, int(args.restarts) + 1)  # include base run
    overall_deadline = (time.time() + args.time_budget) if args.time_budget else None
    theme_feasible_counts: Dict[str, int] = {}
    
    # Helper function for parallel execution
    def run_optimization(theme_name: str, r_id: int, seed_offset: int) -> List[Dict]:
        """Run a single optimisation with given parameters."""
        np.random.seed(seed_offset)
        return greedy_then_refine(
            comps, baseline=baseline, K=args.K, fuel_type=args.fuel_type,
            strict=bool(args.strict_spec), allow_novel=bool(args.allow_novel),
            max_novel=args.max_add_novel, delta_budget=args.delta_budget if baseline else None,
            min_component_frac=float(args.min_component_frac), theme=theme_name,
            anneal_prob=float(args.anneal_prob), anneal_temp=float(args.anneal_temp),
            pmi_cap=args.pmi_cap, cost_cap=args.cost_cap, rvp_cap=args.rvp_cap, t90_cap=args.t90_cap,
            k_explore=bool(args.k_explore), k_max=int(args.k_max),
            num_random_starts=r_id, tag_theme=theme_name, tag_restart=r_id,
            time_deadline=overall_deadline, max_evals=int(args.max_evals or 0),
            reference_mode=args.reference_mode,
            reference_M=args.reference_M
        )
    
    # Process themes in parallel if joblib is available and multiple themes
    if JOBLIB_AVAILABLE and len(themes) > 1:
        # Process each theme in parallel
        def process_theme(theme_name: str) -> Tuple[str, List[Dict], int]:
            """Process a single theme with all its restarts."""
            theme_results = []
            feasible_for_theme = 0
            attempts = 0
            while attempts < restart_count * 2:  # allow extra attempts to meet theme-min
                if overall_deadline is not None and time.time() > overall_deadline:
                    break
                r_id = attempts if attempts < restart_count else attempts - restart_count
                seed_offset = args.seed + r_id
                res_batch = run_optimization(theme_name, r_id, seed_offset)
                theme_results.extend(res_batch)
                batch_scores = np.array([rb['score'] for rb in res_batch], dtype=float)
                feasible_for_theme += int(np.isfinite(batch_scores).sum())
                attempts += 1
                if feasible_for_theme >= int(args.theme_min):
                    break
            return theme_name, theme_results, feasible_for_theme
        
        # Run themes in parallel
        theme_results = Parallel(n_jobs=-1, verbose=0)(
            delayed(process_theme)(theme_name) for theme_name in themes
        )
        
        # Aggregate results
        for theme_name, theme_res, feasible_count in theme_results:
            aggregated.extend(theme_res)
            theme_feasible_counts[theme_name] = feasible_count
    else:
        # Sequential execution (fallback or single theme)
        for theme_name in themes:
            feasible_for_theme = 0
            attempts = 0
            while attempts < restart_count * 2:  # allow extra attempts to meet theme-min
                if overall_deadline is not None and time.time() > overall_deadline:
                    break
                r_id = attempts if attempts < restart_count else attempts - restart_count
                # random starts: use r_id as count, keep determinism via seed offset
                np.random.seed(args.seed + r_id)
                res_batch = greedy_then_refine(
                    comps, baseline=baseline, K=args.K, fuel_type=args.fuel_type,
                    strict=bool(args.strict_spec), allow_novel=bool(args.allow_novel),
                    max_novel=args.max_add_novel, delta_budget=args.delta_budget if baseline else None,
                    min_component_frac=float(args.min_component_frac), theme=theme_name,
                    anneal_prob=float(args.anneal_prob), anneal_temp=float(args.anneal_temp),
                    pmi_cap=args.pmi_cap, cost_cap=args.cost_cap, rvp_cap=args.rvp_cap, t90_cap=args.t90_cap,
                    k_explore=bool(args.k_explore), k_max=int(args.k_max),
                    num_random_starts=r_id, tag_theme=theme_name, tag_restart=r_id,
                    time_deadline=overall_deadline, max_evals=int(args.max_evals or 0),
                    reference_mode=args.reference_mode,
                    reference_M=args.reference_M
                )
                aggregated.extend(res_batch)
                batch_scores = np.array([rb['score'] for rb in res_batch], dtype=float)
                feasible_for_theme += int(np.isfinite(batch_scores).sum())
                attempts += 1
                if feasible_for_theme >= int(args.theme_min):
                    break
            theme_feasible_counts[theme_name] = feasible_for_theme

    # Diversity filter across aggregated results using L1 distance on _x
    def l1_distance(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.sum(np.abs(a - b)))

    diverse: List[Dict] = []
    for r in sorted(aggregated, key=lambda z: z['score'], reverse=True):
        if not diverse:
            diverse.append(r)
            continue
        keep = True
        for kept in diverse:
            if l1_distance(r['_x'], kept['_x']) < float(args.min_l1_distance):
                keep = False
                break
        if keep:
            diverse.append(r)

    # Filter infeasible
    feasible_results = [r for r in (diverse if diverse else aggregated) if np.isfinite(r['score'])]
    results = feasible_results
    
    end_time = time.time()
    
    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(args.out) / f"{timestamp}_{args.seed}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save results
    if results:
        # Write CSV files
        write_top_n_csv(results, args.top, output_dir, comps, baseline_props)
        write_portfolio_csv(results, output_dir)
        
        # Write JSON metadata
        write_run_json(
            results, args, comps, start_time, end_time, output_dir,
            baseline_props, theme_feasible_counts
        )
        
        # Print summary
        print(f"[OK] Generated {len(results)} blend(s)")
        print(f"[BEST] Best score: {results[0]['score']:.6f}")
        print(f"[REL] Relative score ( >1 better than current ): {results[0].get('relative_score', 1.0):.4f}  "
              f"[ref={results[0].get('reference_mode','auto')}]")
        print(f"[TIME] Duration: {end_time - start_time:.2f}s")
        print(f"[OUT] Output: {output_dir}")
        
        if args.verbose:
            print(f"[COMP] Components: {len(comps)} total")
            if baseline:
                print(f"[BASE] Baseline mode: {len(baseline)} components")
            if args.allow_novel:
                print(f"[NOVEL] Novel allowed: max {args.max_add_novel:.1%}")
        
        if baseline:
            # compute L1 deviation for best result
            ids = list(comps['id'])
            x0 = np.array([baseline.get(i,0.0) for i in ids])
            best = results[0]
            xv = np.array([baseline.get(i,0.0) for i in ids]) * 0.0
            for c in best['components']:
                idx = ids.index(c['id'])
                xv[idx] = c['vol_frac']
            print(f"L1 vs baseline: {np.sum(np.abs(xv - x0)):.4f}")
        
        # Print determinism signature
        best = results[0]
        sig_parts = []
        for c in best['components']:
            sig_parts.append(f"{c['id']}:{c['vol_frac']:.3f}")
        signature = "|".join(sig_parts)
        print(f"Best signature: {signature}")
        
        # Print spec flag summary
        ok_flags = sum(v for k,v in best['flags'].items() if k.endswith('_ok'))
        fail_flags = sum(1 for k,v in best['flags'].items() if k.endswith('_ok') and not v)
        print(f"[SPEC] OK flags: {ok_flags}  |  (see CSV for details)")
        
        if results[0].get('note') == 'CAP_DEFICIT':
            print("[WARN] Sum of caps < 1.0 detected; optimiser returned the best feasible blend under caps.")

        # Optional cleanup of old runs
        if args.clean_runs and isinstance(args.clean_runs, int):
            try:
                cutoff = time.time() - float(args.clean_runs) * 86400.0
                runs_root = Path(args.out)
                for child in runs_root.iterdir():
                    if child.is_dir():
                        mtime = child.stat().st_mtime
                        if mtime < cutoff:
                            for p in child.glob('**/*'):
                                if p.is_file():
                                    p.unlink(missing_ok=True)
                            # remove nested dirs then the dir
                            for p in sorted(child.glob('**/*'), reverse=True):
                                if p.is_dir():
                                    p.rmdir()
                            child.rmdir()
                print(f"[CLEAN] Pruned runs older than {args.clean_runs} day(s)")
            except Exception as e:
                print(f"[CLEAN][WARN] Could not prune old runs: {e}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
