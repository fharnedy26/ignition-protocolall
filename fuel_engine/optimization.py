"""
Core optimization algorithms for fuel blend optimization.
"""

import time
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from fuel_engine.constants import (
    WEIGHTS, OPT_PARAMS, PENALTY_WEIGHT, NOVEL_PEN_PER_VOL,
    EPSILON, SIMPLEX_TOL, MAX_SIMPLEX_ITER, EPSILON_SMALL,
    DEFAULT_GREEDY_INIT_FRAC, MIN_COMPONENT_THRESHOLD,
    THEME_WEIGHT_MULTIPLIER
)
from fuel_engine.properties import blend_props, spec_penalties, score


def project_simplex_with_caps(x: np.ndarray, caps: np.ndarray,
                              l1_budget: Optional[float] = None,
                              x0: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Project vector onto simplex with component caps.
    
    Ensures x sums to 1.0 and respects max_vol_frac caps.
    Optionally enforces L1 distance budget from baseline.
    
    Args:
        x: Volume fractions to project
        caps: Maximum volume fractions per component
        l1_budget: Maximum L1 distance from x0 (optional)
        x0: Baseline vector for delta mode (optional)
        
    Returns:
        Projected vector on simplex with caps
    """
    x = np.clip(x, 0.0, caps)
    for _ in range(MAX_SIMPLEX_ITER):
        s = x.sum()
        if abs(s - 1.0) < SIMPLEX_TOL:
            break
        if s <= 0:
            mask = (caps > 0).astype(float)
            denom = max(mask.sum(), 1.0)
            x = mask / denom
        else:
            x *= (1.0 / s)
            x = np.minimum(x, caps)
            x = np.maximum(x, 0.0)
    if l1_budget is not None and x0 is not None:
        d = x - x0
        l1 = np.sum(np.abs(d))
        if l1 > l1_budget:
            x = x0 + (l1_budget / max(l1, EPSILON_SMALL)) * d
            x = project_simplex_with_caps(x, caps, None, None)
    return x


def _create_fuel_type_weights(fuel_type: str, theme: str) -> Dict[str, float]:
    """Create fuel-type aware weights with theme emphasis."""
    weights = WEIGHTS.copy()
    if fuel_type == 'gasoline':
        weights['CET'] = 0.0
    elif fuel_type == 'diesel':
        weights['RON'] = 0.0
    elif fuel_type == 'jet':
        weights['RON'] = 0.0
        weights['CET'] = 0.0
    
    # Theme-based emphasis
    if theme in THEME_WEIGHT_MULTIPLIER:
        multiplier = THEME_WEIGHT_MULTIPLIER[theme]
        if theme == 'high_ron':
            weights['RON'] *= multiplier
        elif theme == 'low_pmi':
            weights['PMI'] *= multiplier
        elif theme == 'high_lhv':
            weights['LHV'] *= multiplier
        elif theme == 'low_cost':
            weights['COST'] *= multiplier
    
    return weights


def _filter_components(comps: pd.DataFrame, allow_novel: bool) -> pd.DataFrame:
    """Filter components based on safety and novel flags."""
    active_comps = comps.copy()
    
    # Safety exclusions
    if 'flags' in active_comps.columns:
        unsafe_mask = active_comps['flags'].str.contains('HALOGEN|METAL', case=False, na=False)
        active_comps = active_comps[~unsafe_mask]
    
    # Novel component filtering
    if not allow_novel:
        active_comps = active_comps[active_comps['is_novel'] == 0]
    
    return active_comps


def greedy_then_refine(comps: pd.DataFrame, baseline: Optional[Dict[str, float]] = None,
                      K: int = 5, fuel_type: str = 'gasoline', strict: bool = True,
                      allow_novel: bool = False, max_novel: float = 0.03,
                      delta_budget: Optional[float] = None,
                      min_component_frac: float = 0.0,
                      theme: str = 'none',
                      anneal_prob: float = 0.0,
                      anneal_temp: float = 0.01,
                      pmi_cap: Optional[float] = None,
                      cost_cap: Optional[float] = None,
                      rvp_cap: Optional[float] = None,
                      t90_cap: Optional[float] = None,
                      k_explore: bool = False,
                      k_max: int = 7,
                      num_random_starts: int = 0,
                      tag_theme: Optional[str] = None,
                      tag_restart: Optional[int] = None,
                      time_deadline: Optional[float] = None,
                      max_evals: int = 0,
                      reference_mode: str = 'auto',
                      reference_M: int = 3) -> List[Dict]:
    """
    Greedy forward selection followed by coordinate descent refinement.
    
    This is the main optimization function that performs fuel blend optimization
    using a combination of greedy selection and coordinate descent.
    
    Args:
        comps: Component library DataFrame
        baseline: Optional baseline blend dictionary
        K: Maximum components in final blend
        fuel_type: Type of fuel ('gasoline', 'diesel', 'jet')
        strict: If True, enforce strict specifications
        allow_novel: Allow novel components
        max_novel: Maximum novel fraction
        delta_budget: L1 distance budget from baseline
        min_component_frac: Minimum manufacturable fraction
        theme: Optimization theme ('none', 'high_ron', 'low_pmi', etc.)
        anneal_prob: Probability of simulated annealing perturbation
        anneal_temp: Temperature for simulated annealing
        pmi_cap: Soft cap on PMI
        cost_cap: Soft cap on cost
        rvp_cap: Soft cap on RVP
        t90_cap: Soft cap on T90
        k_explore: Allow exploration up to k_max
        k_max: Maximum components during exploration
        num_random_starts: Number of random starting points
        tag_theme: Tag for theme tracking
        tag_restart: Tag for restart tracking
        time_deadline: Time deadline for optimization
        max_evals: Maximum objective evaluations
        reference_mode: Mode for reference score ('auto', 'baseline', 'cheapest', 'topkey')
        reference_M: Number of components for reference computation
        
    Returns:
        List of optimization results sorted by score
    """
    # Create fuel-type aware weights
    weights = _create_fuel_type_weights(fuel_type, theme)
    
    # Filter components
    active_comps = _filter_components(comps, allow_novel)
    
    n_comps = len(active_comps)
    if n_comps == 0:
        return []
    
    # Initialize baseline vector
    if baseline:
        x0 = np.zeros(n_comps)
        for i, comp_id in enumerate(active_comps['id']):
            if comp_id in baseline:
                x0[i] = baseline[comp_id]
        x0 = project_simplex_with_caps(x0, active_comps['max_vol_frac'].values)
    else:
        x0 = np.zeros(n_comps)
    
    # Pre-compute arrays for efficiency
    caps_arr = active_comps['max_vol_frac'].values
    cost_arr = active_comps['cost_eur_L'].values if 'cost_eur_L' in active_comps.columns else None
    RON_arr = active_comps['RON'].values if 'RON' in active_comps.columns else None
    cetane_arr = active_comps['cetane'].values if 'cetane' in active_comps.columns else None
    PMI_arr = active_comps['PMI'].values if 'PMI' in active_comps.columns else None
    LHV_arr = active_comps['LHV_MJ_kg'].values if 'LHV_MJ_kg' in active_comps.columns else None
    is_novel_arr = active_comps['is_novel'].values
    
    # Helper function for penalized scoring
    eval_count = 0
    def eval_penalized_score(xv: np.ndarray) -> Tuple[float, Dict[str, float], float, float]:
        nonlocal eval_count
        if time_deadline is not None and time.time() > time_deadline:
            return -np.inf, {}, 0.0, 0.0
        if max_evals and eval_count >= max_evals:
            return -np.inf, {}, 0.0, 0.0
        props_v = blend_props(xv, active_comps)
        eval_count += 1
        pen_v, _flags = spec_penalties(props_v, fuel_type, strict)
        if strict and pen_v > 0:
            return -np.inf, props_v, pen_v, 0.0
        # Epsilon caps (soft penalties)
        if pmi_cap is not None and 'PMI' in props_v:
            pen_v += max(0.0, props_v['PMI'] - pmi_cap) * 0.1
        if cost_cap is not None and 'net_cost_L' in props_v:
            pen_v += max(0.0, props_v['net_cost_L'] - cost_cap) * 1.0
        if rvp_cap is not None and 'RVP' in props_v:
            pen_v += max(0.0, props_v['RVP'] - rvp_cap) * 0.1
        if t90_cap is not None and 'T90' in props_v:
            pen_v += max(0.0, props_v['T90'] - t90_cap) * 0.1
        novel_v = NOVEL_PEN_PER_VOL * np.sum(xv * is_novel_arr)
        s_v = score(xv, props_v, weights, x0 if baseline else None, novel_v)
        s_v -= PENALTY_WEIGHT * pen_v
        return s_v, props_v, pen_v, novel_v
    
    # Helper function to enforce hard novel cap
    def enforce_novel_cap(xv: np.ndarray) -> np.ndarray:
        if not allow_novel:
            return xv
        novel_mask = is_novel_arr.astype(bool)
        total_novel = float(np.sum(xv[novel_mask]))
        if total_novel > max_novel + EPSILON and total_novel > 0:
            scale = max_novel / total_novel
            # scale down novel slice proportionally
            xv_novel = xv[novel_mask] * scale
            # put difference back into non-novel proportionally to available headroom
            freed = np.sum(xv[novel_mask]) - np.sum(xv_novel)
            xv[novel_mask] = xv_novel
            non_mask = ~novel_mask
            # distribute freed volume to non-novel within their caps
            headroom = np.maximum(0.0, caps_arr[non_mask] - xv[non_mask])
            if headroom.sum() > 0:
                xv[non_mask] += freed * (headroom / headroom.sum())
            # final safety projection
            xv = project_simplex_with_caps(xv, caps_arr, None, None)
        return xv
    
    # Helper function for scoring arbitrary vectors (used by reference computation)
    def score_vector(xv: np.ndarray) -> float:
        props_v = blend_props(xv, active_comps)
        pen_v, _ = spec_penalties(props_v, fuel_type, strict)
        if strict and pen_v > 0:
            return -np.inf
        novel_v = NOVEL_PEN_PER_VOL * np.sum(xv * is_novel_arr)
        s_v = score(xv, props_v, weights, x0 if baseline else None, novel_v)
        s_v -= PENALTY_WEIGHT * pen_v
        return float(s_v)
    
    # Nested function to run optimization from a given start
    def run_from_start(x_start: np.ndarray) -> Dict:
        # project start
        xv = project_simplex_with_caps(x_start.copy(), caps_arr,
                                       delta_budget, x0 if baseline else None)
        xv = enforce_novel_cap(xv)
        
        x = xv.copy()
        selected = set(np.nonzero(x > EPSILON)[0].tolist())
        
        # Determine exploration K
        K_use = min(k_max, K) if k_explore else K
        
        # Greedy add up to K_use
        for _ in range(min(K_use, n_comps - len(selected))):
            if time_deadline is not None and time.time() > time_deadline:
                break
            best_idx = -1
            s_cur, props_cur, pen_cur, nov_cur = eval_penalized_score(x)
            best_gain = 0.0
            for i in range(n_comps):
                if i in selected: 
                    continue
                xt = x.copy()
                xt[i] = min(DEFAULT_GREEDY_INIT_FRAC, active_comps.iloc[i]['max_vol_frac'])
                xt = project_simplex_with_caps(xt, caps_arr,
                                               delta_budget, x0 if baseline else None)
                if allow_novel:
                    if np.sum(xt * is_novel_arr) > max_novel:
                        continue
                s_new, _, pen_new, _ = eval_penalized_score(xt)
                gain = s_new - s_cur
                if gain > best_gain:
                    best_gain, best_idx = gain, i
            if best_idx >= 0:
                selected.add(best_idx)
                x[best_idx] = min(DEFAULT_GREEDY_INIT_FRAC, active_comps.iloc[best_idx]['max_vol_frac'])
                x = project_simplex_with_caps(x, caps_arr,
                                              delta_budget, x0 if baseline else None)
                x = enforce_novel_cap(x)
            else:
                break

        # Coordinate descent refine
        for delta in OPT_PARAMS['delta_steps']:
            for _ in range(OPT_PARAMS['max_iterations']):
                if time_deadline is not None and time.time() > time_deadline:
                    break
                improved = False
                s_cur, props_cur, pen_cur, nov_cur = eval_penalized_score(x)
                for i in list(selected):
                    # increase
                    xt = x.copy()
                    inc = min(delta, active_comps.iloc[i]['max_vol_frac'] - x[i])
                    if inc > 0:
                        xt[i] += inc
                        xt = project_simplex_with_caps(xt, caps_arr,
                                                       delta_budget, x0 if baseline else None)
                        xt = enforce_novel_cap(xt)
                        s_new, _, pen_new, _ = eval_penalized_score(xt)
                        if s_new > s_cur:
                            x = xt; improved = True; break
                    # decrease
                    xt = x.copy()
                    xt[i] = max(0.0, xt[i] - delta)
                    xt = project_simplex_with_caps(xt, caps_arr,
                                                   delta_budget, x0 if baseline else None)
                    xt = enforce_novel_cap(xt)
                    s_new, _, pen_new, _ = eval_penalized_score(xt)
                    if s_new > s_cur:
                        x = xt; improved = True; break
                # Stochastic 2D perturbation: swap/shift between two components
                if not improved and anneal_prob > 0.0 and np.random.rand() < anneal_prob and len(selected) >= 2:
                    idxs = list(selected)
                    i1, i2 = np.random.choice(idxs, size=2, replace=False)
                    xt = x.copy()
                    shift = min(delta, xt[i2])
                    xt[i1] = min(active_comps.iloc[i1]['max_vol_frac'], xt[i1] + shift)
                    xt[i2] = max(0.0, xt[i2] - shift)
                    xt = project_simplex_with_caps(xt, caps_arr,
                                                   delta_budget, x0 if baseline else None)
                    xt = enforce_novel_cap(xt)
                    s_new, _, _, _ = eval_penalized_score(xt)
                    if s_new > s_cur:
                        x = xt; improved = True
                    else:
                        # Annealed acceptance
                        if anneal_temp > 0:
                            prob = np.exp((s_new - s_cur) / max(anneal_temp, EPSILON_SMALL))
                            if np.random.rand() < prob:
                                x = xt; improved = True
                if not improved:
                    break

        # Post-refine sparsify: keep top-K by volume, zero the rest, then reproject and re-enforce novel cap
        if K < n_comps:
            keep_idx = np.argsort(-x)[:K]
            mask = np.zeros_like(x, dtype=bool)
            mask[keep_idx] = True
            x = x * mask.astype(float)
            x = project_simplex_with_caps(x, caps_arr,
                                          delta_budget, x0 if baseline else None)
            x = enforce_novel_cap(x)
        
        # Enforce minimum manufacturable fraction: zero traces, then reproject
        if min_component_frac > 0.0:
            small_mask = x < min_component_frac
            if np.any(small_mask):
                x[small_mask] = 0.0
                x = project_simplex_with_caps(x, caps_arr,
                                              delta_budget, x0 if baseline else None)
                x = enforce_novel_cap(x)
        
        # finalize
        x = enforce_novel_cap(x)  # final safety check
        s_fin, props_fin, pen_fin, nov_fin = eval_penalized_score(x)
        # compute true raw score without any penalties
        raw_no_pen = score(x, props_fin, weights, x0 if baseline else None, 0.0)
        components = []
        for i in range(n_comps):
            if x[i] > MIN_COMPONENT_THRESHOLD:
                components.append({'id': active_comps.iloc[i]['id'],
                                   'name': active_comps.iloc[i]['name'],
                                   'vol_frac': float(x[i])})
        components.sort(key=lambda c: (-c['vol_frac'], str(c['id'])))
        sum_x = float(np.sum(x))
        note = "CAP_DEFICIT" if sum_x < 1.0 - EPSILON else ""
        ret = {
            'score': float(s_fin),
            'raw_score': float(raw_no_pen),
            'spec_penalty': float(pen_fin),
            'novel_penalty': float(NOVEL_PEN_PER_VOL * np.sum(x * is_novel_arr)),
            'novel_fraction': float(np.sum(x * is_novel_arr)),
            'components': components,
            'properties': props_fin,
            'flags': spec_penalties(props_fin, fuel_type, strict)[1],
            '_x': x,  # for dedupe
            'theme': tag_theme if tag_theme is not None else theme,
            'restart_id': int(tag_restart) if tag_restart is not None else 0
        }
        if note:
            ret['note'] = note
        return ret

    # Reference score computation for relative normalization
    def compute_reference_score() -> Tuple[float, np.ndarray, str]:
        mode = reference_mode
        M = reference_M
        
        # 1) baseline
        if mode in ('auto','baseline') and baseline:
            xr = project_simplex_with_caps(x0.copy(), caps_arr, delta_budget, x0 if baseline else None)
            xr = enforce_novel_cap(xr)
            sr = score_vector(xr)
            if sr != -np.inf:
                return sr, xr, 'baseline'
            if mode == 'baseline':
                return sr, xr, 'baseline'  # may be -inf if infeasible

        # 2) cheapest equal-split
        if mode in ('auto','cheapest'):
            if cost_arr is not None:
                idx = np.argsort(cost_arr)[:max(1, min(M, n_comps))]
                xr = np.zeros(n_comps); xr[idx] = 1.0 / max(len(idx),1)
                xr = project_simplex_with_caps(xr, caps_arr, None, None)
                xr = enforce_novel_cap(xr)
                sr = score_vector(xr)
                if sr != -np.inf:
                    return sr, xr, 'cheapest'

        # 3) top-key equal-split (RON/cetane or low PMI for jet)
        if mode in ('auto','topkey'):
            if fuel_type == 'gasoline' and RON_arr is not None:
                idx = np.argsort(-RON_arr)[:max(1, min(M, n_comps))]
            elif fuel_type == 'diesel' and cetane_arr is not None:
                idx = np.argsort(-cetane_arr)[:max(1, min(M, n_comps))]
            elif fuel_type == 'jet':
                if PMI_arr is not None:
                    idx = np.argsort(PMI_arr)[:max(1, min(M, n_comps))]  # lower is better
                else:
                    idx = np.argsort(-LHV_arr)[:max(1, min(M, n_comps))]
            else:
                idx = np.arange(min(M, n_comps))
            xr = np.zeros(n_comps); xr[idx] = 1.0 / max(len(idx),1)
            xr = project_simplex_with_caps(xr, caps_arr, None, None)
            xr = enforce_novel_cap(xr)
            sr = score_vector(xr)
            if sr != -np.inf:
                return sr, xr, 'topkey'

        # 4) fallback: use best result later to avoid divide-by-zero; placeholder -inf here
        return -np.inf, np.zeros(n_comps), 'fallback'

    # Build deterministic starts
    starts = []

    # baseline start (if any)
    if baseline:
        starts.append(x0.copy())

    # equal split of top-M by key property
    if fuel_type == 'gasoline':
        key = 'RON'
    elif fuel_type == 'diesel':
        key = 'cetane'
    else:  # jet
        # prefer low PMI (soot) first; if PMI missing, fall back to high LHV
        key = 'PMI' if 'PMI' in active_comps.columns else 'LHV_MJ_kg'
    
    if key in active_comps.columns:
        M = min(K, 5)
        if fuel_type == 'jet' and key == 'PMI':
            idx = np.argsort(active_comps['PMI'].values)[:M]   # lowest PMI
        else:
            idx = np.argsort(-active_comps[key].values)[:M]    # highest RON/cetane/LHV
        xs = np.zeros(n_comps); xs[idx] = 1.0 / max(len(idx), 1)
        starts.append(xs)

    # equal split of M cheapest
    if 'cost_eur_L' in active_comps.columns:
        M = min(K, 5)
        idx = np.argsort(active_comps['cost_eur_L'].values)[:M]
        xs = np.zeros(n_comps); xs[idx] = 1.0 / max(len(idx),1)
        starts.append(xs)

    # singleton starts for first P components (stable by id)
    P = min(5, n_comps)
    order = np.argsort(active_comps['id'].astype(str).values)[:P]
    for i in order:
        xs = np.zeros(n_comps)
        xs[i] = min(DEFAULT_GREEDY_INIT_FRAC, active_comps.iloc[i]['max_vol_frac'])
        starts.append(xs)

    # Random Dirichlet starts respecting caps
    for _ in range(max(0, int(num_random_starts))):
        raw = np.random.dirichlet(np.ones(n_comps))
        xs = raw * 1.0
        xs = project_simplex_with_caps(xs, caps_arr,
                                       delta_budget, x0 if baseline else None)
        starts.append(xs)

    # run and dedupe
    seen = {}
    for xs in starts:
        res = run_from_start(xs)
        sig = tuple((res['_x'] * 1000).round().astype(int))
        if sig not in seen or res['score'] > seen[sig]['score']:
            seen[sig] = res

    results = sorted(seen.values(), key=lambda r: r['score'], reverse=True)
    
    # Reference score for relative normalization
    s_ref, x_ref, ref_mode_used = compute_reference_score()
    if s_ref == -np.inf:
        # fallback to best result as reference to avoid div-by-zero; relative ≤ 1 in this case
        s_ref = float(results[0]['score']) if results else 1.0
        ref_mode_used = 'best_fallback'

    for r in results:
        # Guard against tiny or negative reference
        denom = s_ref if abs(s_ref) > EPSILON else 1.0
        r['relative_score'] = float(r['score'] / denom)
        r['reference_mode'] = ref_mode_used
    
    return results


