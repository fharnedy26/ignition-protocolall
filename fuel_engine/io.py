"""
I/O operations for saving optimization results.
"""

import csv
import json
import platform
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from fuel_engine.constants import SPECS


def write_top_n_csv(results: List[Dict], top_n: int, output_dir: Path,
                   comps: pd.DataFrame, baseline_props: Optional[Dict] = None) -> None:
    """
    Write Top-N CSV file with optimization results.
    
    Args:
        results: List of optimization results
        top_n: Number of top results to write
        output_dir: Output directory path
        comps: Component library DataFrame
        baseline_props: Optional baseline properties for delta computation
    """
    top_n = min(top_n, len(results))
    with open(output_dir / f"Top-{top_n}.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['relative_score','score','raw_score','spec_penalty','novel_penalty',
                  'RON','MON','cetane','LHV_mass','LHV_vol','PMI','TSI','OC_ratio','ring_count',
                  'T10','T50','T90','RVP','cost_L','net_cost_L','waste_credit_total',
                  'novel_fraction','components','recipe_mL_1000','recipe_g_1000','theme','restart_id','feasible']
        # Add baseline deltas if available
        if baseline_props is not None:
            header += ['d_RON','d_LHV_vol','d_PMI','d_RVP','d_T10','d_T90','d_cost_L','d_net_cost_L']
        writer.writerow(header)
        
        # Build density mapping for recipes
        id_to_density = {str(row_id): dens for row_id, dens in zip(comps['id'].astype(str), comps['density_g_ml'])}
        
        for i, result in enumerate(results[:top_n]):
            props = result['properties']
            components_str = '; '.join([f"{c['id']}:{c['vol_frac']:.3f}" for c in result['components']])
            
            # Build per-1000 mL recipes
            recipe_ml_parts = []
            recipe_g_parts = []
            for c in result['components']:
                ml = 1000.0 * c['vol_frac']
                dens = id_to_density.get(str(c['id']), 1.0)
                grams = ml * dens
                recipe_ml_parts.append(f"{c['id']}:{ml:.1f}mL")
                recipe_g_parts.append(f"{c['id']}:{grams:.1f}g")
            recipe_ml = ' | '.join(recipe_ml_parts)
            recipe_g = ' | '.join(recipe_g_parts)

            row = [
                f"{result.get('relative_score', 1.0):.6f}",
                f"{result['score']:.6f}",
                f"{result.get('raw_score', result['score']):.6f}",
                f"{result.get('spec_penalty', 0.0):.6f}",
                f"{result.get('novel_penalty', 0.0):.6f}",
                f"{props['RON']:.1f}",
                f"{props['MON']:.1f}",
                f"{props['cetane']:.1f}",
                f"{props['LHV_mass']:.1f}",
                f"{props['LHV_vol']:.1f}",
                f"{props['PMI']:.1f}",
                f"{props['TSI']:.1f}",
                f"{props['OC_ratio']:.3f}",
                f"{props['ring_count']:.1f}",
                f"{props['T10']:.1f}",
                f"{props['T50']:.1f}",
                f"{props['T90']:.1f}",
                f"{props['RVP']:.1f}",
                f"{props['cost_L']:.3f}",
                f"{props.get('net_cost_L', props['cost_L']):.3f}",
                f"{props.get('waste_credit_total', 0.0):.3f}",
                f"{result['novel_fraction']:.3f}",
                components_str,
                recipe_ml,
                recipe_g,
                str(result.get('theme','none')),
                str(result.get('restart_id',0)),
                '1'
            ]
            if baseline_props is not None:
                row += [
                    f"{props['RON'] - baseline_props['RON']:.2f}",
                    f"{props['LHV_vol'] - baseline_props['LHV_vol']:.2f}",
                    f"{props['PMI'] - baseline_props['PMI']:.2f}",
                    f"{props['RVP'] - baseline_props['RVP']:.2f}",
                    f"{props['T10'] - baseline_props['T10']:.2f}",
                    f"{props['T90'] - baseline_props['T90']:.2f}",
                    f"{props['cost_L'] - baseline_props['cost_L']:.3f}",
                    f"{props.get('net_cost_L', props['cost_L']) - baseline_props.get('net_cost_L', baseline_props['cost_L']):.3f}"
                ]
            writer.writerow(row)


def write_portfolio_csv(results: List[Dict], output_dir: Path) -> None:
    """
    Write Portfolio CSV with diverse set of results.
    
    Args:
        results: List of optimization results
        output_dir: Output directory path
    """
    with open(output_dir / "Portfolio.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank','score','relative_score','theme','restart_id','L1_to_best','components'])
        if results:
            x_best = results[0]['_x']
            for idx, r in enumerate(results):
                l1_to_best = float(np.sum(np.abs(r['_x'] - x_best)))
                comps_str = '; '.join([f"{c['id']}:{c['vol_frac']:.3f}" for c in r['components']])
                writer.writerow([
                    idx + 1,
                    f"{r['score']:.6f}",
                    f"{r.get('relative_score', 1.0):.6f}",
                    str(r.get('theme','none')),
                    str(r.get('restart_id',0)),
                    f"{l1_to_best:.3f}",
                    comps_str
                ])


def write_run_json(results: List[Dict], args, comps: pd.DataFrame, start_time: float,
                   end_time: float, output_dir: Path, baseline_props: Optional[Dict] = None,
                   theme_feasible_counts: Optional[Dict[str, int]] = None) -> None:
    """
    Write RUN.json metadata file.
    
    Args:
        results: List of optimization results
        args: Parsed command-line arguments
        comps: Component library DataFrame
        start_time: Start timestamp
        end_time: End timestamp
        output_dir: Output directory path
        baseline_props: Optional baseline properties
        theme_feasible_counts: Optional theme feasible counts
    """
    theme_feasible_counts = theme_feasible_counts or {}
    
    # Per-theme best summary
    theme_bests = {}
    for r in results:
        t = str(r.get('theme','none'))
        if t not in theme_bests or r['score'] > theme_bests[t]['score']:
            theme_bests[t] = {'score': r['score'], 'components': [(c['id'], c['vol_frac']) for c in r['components']]}
    
    run_data = {
        'seed': args.seed,
        'started': datetime.fromtimestamp(start_time).isoformat(),
        'finished': datetime.fromtimestamp(end_time).isoformat(),
        'duration_seconds': end_time - start_time,
        'args': vars(args),
        'counts': {
            'total_components': len(comps),
            'novel_components': len(comps[comps['is_novel'] == 1]) if 'is_novel' in comps.columns else 0,
            'results_generated': len(results)
        },
        'timing': end_time - start_time,
        'K': args.K,
        'top': args.top,
        'fuel_type': args.fuel_type,
        'strict_spec': bool(args.strict_spec),
        'novel_allowed': bool(args.allow_novel),
        'spec_constants': SPECS.get(args.fuel_type, {}),
        'machine_info': {
            'platform': platform.platform(),
            'python': platform.python_version()
        },
        'normalization': {
            'reference_mode': args.reference_mode,
            'reference_M': args.reference_M
        },
        'preset': args.preset,
        'min_component_frac': args.min_component_frac,
        'baseline_deltas': bool(baseline_props is not None),
        'theme_bests': theme_bests,
        'theme_feasible_counts': {str(k): int(v) for k, v in theme_feasible_counts.items()},
        'themes_with_no_feasible': [t for t, cnt in theme_feasible_counts.items() if cnt <= 0]
    }
    
    with open(output_dir / "RUN.json", 'w') as f:
        json.dump(run_data, f, indent=2)


