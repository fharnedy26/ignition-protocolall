"""
Property calculations for fuel blends.
"""

from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import numba
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    jit = None
    prange = None

from fuel_engine.constants import (
    SPECS, WEIGHTS, PENALTY_WEIGHT, NOVEL_PEN_PER_VOL,
    PENALTY_MULTIPLIER_T10, PENALTY_MULTIPLIER_T90, PENALTY_MULTIPLIER_RVP,
    PENALTY_MULTIPLIER_RON, PENALTY_MULTIPLIER_CETANE, PENALTY_MULTIPLIER_DENSITY,
    PENALTY_MULTIPLIER_PMI
)


# Numba-optimized core computation (if numba is available)
if NUMBA_AVAILABLE:
    @jit(nopython=True, parallel=False)
    def _blend_props_core_numba(x, ron, mon, cetane, pmi, tsi, oc_ratio, ring_count,
                                vapor_pressure, cost, lhv_mj_kg, density, t10, t50, t90,
                                waste_credit, is_novel, has_waste_credit):
        """Numba JIT-compiled core computation for maximum performance."""
        n = len(x)
        
        # Compute all properties
        ron_val = 0.0
        mon_val = 0.0
        cetane_val = 0.0
        pmi_val = 0.0
        tsi_val = 0.0
        oc_ratio_val = 0.0
        ring_count_val = 0.0
        rvp_val = 0.0
        cost_val = 0.0
        t10_val = 0.0
        t50_val = 0.0
        t90_val = 0.0
        lhv_vol_val = 0.0
        total_mass_val = 0.0
        
        for i in range(n):
            ron_val += x[i] * ron[i]
            mon_val += x[i] * mon[i]
            cetane_val += x[i] * cetane[i]
            pmi_val += x[i] * pmi[i]
            tsi_val += x[i] * tsi[i]
            oc_ratio_val += x[i] * oc_ratio[i]
            ring_count_val += x[i] * ring_count[i]
            rvp_val += x[i] * vapor_pressure[i]
            cost_val += x[i] * cost[i]
            t10_val += x[i] * t10[i]
            t50_val += x[i] * t50[i]
            t90_val += x[i] * t90[i]
            lhv_vol_val += x[i] * lhv_mj_kg[i] * density[i]
            total_mass_val += x[i] * density[i]
        
        # Net cost calculation
        net_cost_val = cost_val
        waste_credit_total = 0.0
        if has_waste_credit:
            for i in range(n):
                net_cost_val -= x[i] * waste_credit[i] * is_novel[i]
                waste_credit_total += x[i] * waste_credit[i] * is_novel[i]
        
        # LHV mass
        lhv_mass_val = lhv_vol_val / total_mass_val if total_mass_val > 0.0 else 0.0
        
        # Return as tuple (dicts not supported in numba)
        return (ron_val, mon_val, cetane_val, pmi_val, tsi_val, oc_ratio_val,
                ring_count_val, rvp_val, cost_val, net_cost_val, waste_credit_total,
                t10_val, t50_val, t90_val, lhv_vol_val, lhv_mass_val, total_mass_val)


def _blend_props_core(x: np.ndarray, ron: np.ndarray, mon: np.ndarray, cetane: np.ndarray,
                     pmi: np.ndarray, tsi: np.ndarray, oc_ratio: np.ndarray,
                     ring_count: np.ndarray, vapor_pressure: np.ndarray,
                     cost: np.ndarray, lhv_mj_kg: np.ndarray, density: np.ndarray,
                     t10: np.ndarray, t50: np.ndarray, t90: np.ndarray,
                     waste_credit: Optional[np.ndarray] = None,
                     is_novel: Optional[np.ndarray] = None) -> Dict[str, float]:
    """
    Core blend property calculation using numpy arrays (cacheable and JIT-compilable).
    
    This is the inner function that does the actual computation. It's designed to be
    cacheable (using numpy arrays) and JIT-compilable with numba.
    """
    # Use numba-optimized version if available
    if NUMBA_AVAILABLE and waste_credit is not None and is_novel is not None:
        has_waste_credit = True
        result = _blend_props_core_numba(
            x, ron, mon, cetane, pmi, tsi, oc_ratio, ring_count,
            vapor_pressure, cost, lhv_mj_kg, density, t10, t50, t90,
            waste_credit, is_novel, has_waste_credit
        )
        props = {
            'RON': float(result[0]),
            'MON': float(result[1]),
            'cetane': float(result[2]),
            'PMI': float(result[3]),
            'TSI': float(result[4]),
            'OC_ratio': float(result[5]),
            'ring_count': float(result[6]),
            'RVP': float(result[7]),
            'cost_L': float(result[8]),
            'net_cost_L': float(result[9]),
            'waste_credit_total': float(result[10]),
            'T10': float(result[11]),
            'T50': float(result[12]),
            'T90': float(result[13]),
            'LHV_vol': float(result[14]),
            'LHV_mass': float(result[15]),
            'density': float(result[16])
        }
        return props
    elif NUMBA_AVAILABLE:
        # Create dummy arrays for numba
        dummy_waste = np.zeros_like(x)
        dummy_novel = np.zeros_like(x)
        has_waste_credit = False
        result = _blend_props_core_numba(
            x, ron, mon, cetane, pmi, tsi, oc_ratio, ring_count,
            vapor_pressure, cost, lhv_mj_kg, density, t10, t50, t90,
            dummy_waste, dummy_novel, has_waste_credit
        )
        props = {
            'RON': float(result[0]),
            'MON': float(result[1]),
            'cetane': float(result[2]),
            'PMI': float(result[3]),
            'TSI': float(result[4]),
            'OC_ratio': float(result[5]),
            'ring_count': float(result[6]),
            'RVP': float(result[7]),
            'cost_L': float(result[8]),
            'net_cost_L': float(result[8]),  # Same as cost_L
            'waste_credit_total': 0.0,
            'T10': float(result[11]),
            'T50': float(result[12]),
            'T90': float(result[13]),
            'LHV_vol': float(result[14]),
            'LHV_mass': float(result[15]),
            'density': float(result[16])
        }
        return props
    
    # Fallback to pure numpy version
    props = {}
    
    # Volume-weighted properties
    props['RON'] = float(np.sum(x * ron))
    props['MON'] = float(np.sum(x * mon))
    props['cetane'] = float(np.sum(x * cetane))
    props['PMI'] = float(np.sum(x * pmi))
    props['TSI'] = float(np.sum(x * tsi))
    props['OC_ratio'] = float(np.sum(x * oc_ratio))
    props['ring_count'] = float(np.sum(x * ring_count))
    props['RVP'] = float(np.sum(x * vapor_pressure))
    props['cost_L'] = float(np.sum(x * cost))
    
    # Net cost accounts for waste credit applied only to novel components
    if waste_credit is not None and is_novel is not None:
        net_cost_terms = cost - (waste_credit * is_novel)
        props['net_cost_L'] = float(np.sum(x * net_cost_terms))
        props['waste_credit_total'] = float(np.sum(x * waste_credit * is_novel))
    else:
        props['net_cost_L'] = props['cost_L']
        props['waste_credit_total'] = 0.0
    
    # Volatility (piecewise-linear CDF mix)
    props['T10'] = float(np.sum(x * t10))
    props['T50'] = float(np.sum(x * t50))
    props['T90'] = float(np.sum(x * t90))
    
    # Energy density
    lhv_vol = np.sum(x * lhv_mj_kg * density)
    props['LHV_vol'] = float(lhv_vol)
    
    # Mass-weighted LHV
    total_mass = np.sum(x * density)
    props['LHV_mass'] = float(lhv_vol / total_mass) if total_mass > 0 else 0.0
    
    # Blend density
    props['density'] = float(np.sum(x * density))
    
    return props


@lru_cache(maxsize=100000)
def _blend_props_cached(x_tuple: tuple, ron_tuple: tuple, mon_tuple: tuple, cetane_tuple: tuple,
                        pmi_tuple: tuple, tsi_tuple: tuple, oc_ratio_tuple: tuple,
                        ring_count_tuple: tuple, vapor_pressure_tuple: tuple,
                        cost_tuple: tuple, lhv_mj_kg_tuple: tuple, density_tuple: tuple,
                        t10_tuple: tuple, t50_tuple: tuple, t90_tuple: tuple,
                        waste_credit_tuple: Optional[tuple] = None,
                        is_novel_tuple: Optional[tuple] = None) -> Dict[str, float]:
    """
    Cached version of blend_props_core. Converts tuples back to arrays.
    """
    x = np.array(x_tuple)
    ron = np.array(ron_tuple)
    mon = np.array(mon_tuple)
    cetane = np.array(cetane_tuple)
    pmi = np.array(pmi_tuple)
    tsi = np.array(tsi_tuple)
    oc_ratio = np.array(oc_ratio_tuple)
    ring_count = np.array(ring_count_tuple)
    vapor_pressure = np.array(vapor_pressure_tuple)
    cost = np.array(cost_tuple)
    lhv_mj_kg = np.array(lhv_mj_kg_tuple)
    density = np.array(density_tuple)
    t10 = np.array(t10_tuple)
    t50 = np.array(t50_tuple)
    t90 = np.array(t90_tuple)
    
    waste_credit_arr = None
    is_novel_arr = None
    if waste_credit_tuple is not None:
        waste_credit_arr = np.array(waste_credit_tuple)
    if is_novel_tuple is not None:
        is_novel_arr = np.array(is_novel_tuple)
    
    return _blend_props_core(x, ron, mon, cetane, pmi, tsi, oc_ratio, ring_count,
                            vapor_pressure, cost, lhv_mj_kg, density, t10, t50, t90,
                            waste_credit_arr, is_novel_arr)


def blend_props(x: np.ndarray, comps: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate blend properties from volume fractions.
    
    This function extracts component arrays from the DataFrame and calls
    the cached core function for performance.
    
    Args:
        x: Volume fractions array (must sum to ~1.0)
        comps: DataFrame with component properties
        
    Returns:
        Dictionary of blend properties
    """
    # Extract component property arrays
    ron = comps['RON'].values
    mon = comps['MON'].values
    cetane = comps['cetane'].values
    pmi = comps['PMI'].values
    tsi = comps['TSI'].values
    oc_ratio = comps['OC_ratio'].values
    ring_count = comps['ring_count'].values
    vapor_pressure = comps['vapor_pressure_kPa'].values
    cost = comps['cost_eur_L'].values
    lhv_mj_kg = comps['LHV_MJ_kg'].values
    density = comps['density_g_ml'].values
    t10 = comps['T10_C'].values
    t50 = comps['T50_C'].values
    t90 = comps['T90_C'].values
    
    # Extract optional arrays
    waste_credit = None
    is_novel = None
    if 'waste_credit_eur_L' in comps.columns and 'is_novel' in comps.columns:
        waste_credit = comps['waste_credit_eur_L'].values
        is_novel = comps['is_novel'].values
    
    # Convert to tuples for caching (numpy arrays aren't hashable)
    x_tuple = tuple(x)
    ron_tuple = tuple(ron)
    mon_tuple = tuple(mon)
    cetane_tuple = tuple(cetane)
    pmi_tuple = tuple(pmi)
    tsi_tuple = tuple(tsi)
    oc_ratio_tuple = tuple(oc_ratio)
    ring_count_tuple = tuple(ring_count)
    vapor_pressure_tuple = tuple(vapor_pressure)
    cost_tuple = tuple(cost)
    lhv_mj_kg_tuple = tuple(lhv_mj_kg)
    density_tuple = tuple(density)
    t10_tuple = tuple(t10)
    t50_tuple = tuple(t50)
    t90_tuple = tuple(t90)
    
    waste_credit_tuple = None
    is_novel_tuple = None
    if waste_credit is not None:
        waste_credit_tuple = tuple(waste_credit)
    if is_novel is not None:
        is_novel_tuple = tuple(is_novel)
    
    # Call cached version
    return _blend_props_cached(
        x_tuple, ron_tuple, mon_tuple, cetane_tuple, pmi_tuple, tsi_tuple,
        oc_ratio_tuple, ring_count_tuple, vapor_pressure_tuple, cost_tuple,
        lhv_mj_kg_tuple, density_tuple, t10_tuple, t50_tuple, t90_tuple,
        waste_credit_tuple, is_novel_tuple
    )


def _check_range(prop_value: float, min_val: Optional[float], max_val: Optional[float],
                 penalty_multiplier: float) -> Tuple[float, bool, bool]:
    """
    Check if a property is within range and calculate penalty.
    
    Returns:
        Tuple of (penalty, ok_flag, low_flag, high_flag)
    """
    penalty = 0.0
    ok = True
    low = False
    high = False
    
    if min_val is not None and prop_value < min_val:
        penalty = (min_val - prop_value) * penalty_multiplier
        ok = False
        low = True
    elif max_val is not None and prop_value > max_val:
        penalty = (prop_value - max_val) * penalty_multiplier
        ok = False
        high = True
    
    return penalty, ok, low, high


def spec_penalties(props: Dict[str, float], fuel_type: str, strict: bool) -> Tuple[float, Dict[str, bool]]:
    """
    Calculate specification penalties and pass/fail flags.
    
    Args:
        props: Blend properties dictionary
        fuel_type: Type of fuel ('gasoline', 'diesel', or 'jet')
        strict: If True, return -inf for any violation
        
    Returns:
        Tuple of (penalty, flags_dict)
    """
    penalty = 0.0
    flags = {}
    
    if fuel_type not in SPECS:
        return penalty, flags
    
    spec = SPECS[fuel_type]
    
    if fuel_type == 'gasoline':
        # T10 range
        pen, ok, low, high = _check_range(
            props.get('T10', 0.0), spec['T10_min'], spec['T10_max'], PENALTY_MULTIPLIER_T10
        )
        penalty += pen
        if low:
            flags['T10_low'] = False
        elif high:
            flags['T10_high'] = False
        else:
            flags['T10_ok'] = True
        
        # T90 max
        pen, ok, _, high = _check_range(
            props.get('T90', 0.0), None, spec['T90_max'], PENALTY_MULTIPLIER_T90
        )
        penalty += pen
        if high:
            flags['T90_high'] = False
        else:
            flags['T90_ok'] = True
        
        # RVP max
        pen, ok, _, high = _check_range(
            props.get('RVP', 0.0), None, spec['RVP_max'], PENALTY_MULTIPLIER_RVP
        )
        penalty += pen
        if high:
            flags['RVP_high'] = False
        else:
            flags['RVP_ok'] = True
        
        # RON min
        pen, ok, low, _ = _check_range(
            props.get('RON', 0.0), spec['RON_min'], None, PENALTY_MULTIPLIER_RON
        )
        penalty += pen
        if low:
            flags['RON_low'] = False
        else:
            flags['RON_ok'] = True
        
        # Ensure all gasoline flags exist
        for key in ['T10_ok','T10_low','T10_high','T90_ok','T90_high','RVP_ok','RVP_high','RON_ok','RON_low']:
            flags.setdefault(key, False)
    
    elif fuel_type == 'diesel':
        # Cetane min
        pen, ok, low, _ = _check_range(
            props.get('cetane', 0.0), spec['cetane_min'], None, PENALTY_MULTIPLIER_CETANE
        )
        penalty += pen
        if low:
            flags['cetane_low'] = False
        else:
            flags['cetane_ok'] = True
        
        # Density range
        pen, ok, low, high = _check_range(
            props.get('density', 0.0), spec['density_min'], spec['density_max'], PENALTY_MULTIPLIER_DENSITY
        )
        penalty += pen
        if low:
            flags['density_low'] = False
        elif high:
            flags['density_high'] = False
        else:
            flags['density_ok'] = True
        
        # T90 max
        pen, ok, _, high = _check_range(
            props.get('T90', 0.0), None, spec['T90_max'], PENALTY_MULTIPLIER_T90
        )
        penalty += pen
        if high:
            flags['T90_high'] = False
        else:
            flags['T90_ok'] = True
        
        # Ensure all diesel flags exist
        for key in ['cetane_ok','cetane_low','density_ok','density_low','density_high','T90_ok','T90_high']:
            flags.setdefault(key, False)
    
    elif fuel_type == 'jet':
        # Density range
        pen, ok, low, high = _check_range(
            props.get('density', 0.0), spec['density_min'], spec['density_max'], PENALTY_MULTIPLIER_DENSITY
        )
        penalty += pen
        if low:
            flags['density_low'] = False
        elif high:
            flags['density_high'] = False
        else:
            flags['density_ok'] = True
        
        # T90 max
        pen, ok, _, high = _check_range(
            props.get('T90', 0.0), None, spec['T90_max'], PENALTY_MULTIPLIER_T90
        )
        penalty += pen
        if high:
            flags['T90_high'] = False
        else:
            flags['T90_ok'] = True
        
        # PMI max (smoke point proxy)
        pen, ok, _, high = _check_range(
            props.get('PMI', 0.0), None, spec['PMI_max'], PENALTY_MULTIPLIER_PMI
        )
        penalty += pen
        if high:
            flags['PMI_high'] = False
        else:
            flags['PMI_ok'] = True
        
        # Ensure all jet flags exist
        for key in ['density_ok','density_low','density_high','T90_ok','T90_high','PMI_ok','PMI_high']:
            flags.setdefault(key, False)
    
    return penalty, flags


def score(x: np.ndarray, props: Dict[str, float], weights: Dict[str, float], 
          x0: Optional[np.ndarray] = None, novel_penalty: float = 0.0) -> float:
    """
    Calculate single scalar objective score.
    
    Args:
        x: Current volume fractions
        props: Blend properties dictionary
        weights: Weight dictionary for scoring
        x0: Baseline volume fractions (for delta mode)
        novel_penalty: Novel component penalty
        
    Returns:
        Scalar score (higher is better)
    """
    s = 0.0
    
    # Positive terms (with safe access using .get() with defaults)
    s += weights['RON'] * props.get('RON', 0.0)
    s += weights['CET'] * props.get('cetane', 0.0)
    s += weights['LHV'] * props.get('LHV_vol', 0.0)
    
    # Negative terms (with safe access)
    s -= weights['PMI'] * props.get('PMI', 0.0)
    # Use net cost when available
    cost_val = props.get('net_cost_L', props.get('cost_L', 0.0))
    s -= weights['COST'] * cost_val
    
    # Delta mode penalty
    if x0 is not None:
        l1_dev = np.sum(np.abs(x - x0))
        s -= weights['DELTA'] * l1_dev
    
    # Novel component penalty
    s -= novel_penalty
    
    return s

