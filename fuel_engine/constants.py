"""
Constants and configuration for fuel blend optimization.
"""

# Numerical tolerances
EPSILON = 1e-12
EPSILON_SMALL = 1e-6
SIMPLEX_TOL = 1e-12
MAX_SIMPLEX_ITER = 30

# Default weights for scoring
WEIGHTS = {
    'RON': 0.3,      # Gasoline octane
    'CET': 0.3,      # Diesel cetane  
    'LHV': 0.2,      # Energy density
    'PMI': 0.1,      # Soot proxy (negative)
    'COST': 0.05,    # Cost penalty (negative)
    'DELTA': 0.05    # Delta mode penalty (negative)
}

# Penalty weights
PENALTY_WEIGHT = 1.0          # scales soft spec penalty when strict==0
NOVEL_PEN_PER_VOL = 0.02      # penalty per absolute novel vol-fraction (reduced from 0.1 to make novel molecules more competitive)

# Fuel type specifications
SPECS = {
    'gasoline': {
        'T10_min': 40, 'T10_max': 70,
        'T90_max': 180,
        'RVP_max': 60,  # kPa, season-dependent
        'RON_min': 87
    },
    'diesel': {
        'cetane_min': 40,
        'density_min': 0.82, 'density_max': 0.86,
        'T90_max': 370
    },
    'jet': {
        'density_min': 0.775, 'density_max': 0.84,
        'T90_max': 300,
        'PMI_max': 25  # Smoke point proxy
    }
}

# Optimization parameters
OPT_PARAMS = {
    'delta_steps': [0.10, 0.02, 0.005],
    'convergence_threshold': 1e-6,
    'max_iterations': 1000
}

# Default component property values
REQUIRED_DEFAULTS = {
    'density_g_ml': 0.75, 'LHV_MJ_kg': 42.0, 'RON': 90.0, 'MON': 85.0,
    'cetane': 50.0, 'PMI': 20.0, 'TSI': 20.0, 'OC_ratio': 0.0, 'ring_count': 0,
    'vapor_pressure_kPa': 10.0, 'T10_C': 50.0, 'T50_C': 100.0, 'T90_C': 150.0,
    'cost_eur_L': 1.0, 'max_vol_frac': 1.0, 'flags': '',
    'waste_credit_eur_L': 0.0
}

# Spec penalty multipliers
PENALTY_MULTIPLIER_T10 = 0.1
PENALTY_MULTIPLIER_T90 = 0.1
PENALTY_MULTIPLIER_RVP = 0.1
PENALTY_MULTIPLIER_RON = 0.1
PENALTY_MULTIPLIER_CETANE = 0.1
PENALTY_MULTIPLIER_DENSITY = 10.0
PENALTY_MULTIPLIER_PMI = 0.1

# Theme weight multipliers
THEME_WEIGHT_MULTIPLIER = {
    'high_ron': 1.3,
    'low_pmi': 1.5,
    'high_lhv': 1.3,
    'low_cost': 1.5
}

# Greedy selection defaults
DEFAULT_GREEDY_INIT_FRAC = 0.1
MIN_COMPONENT_THRESHOLD = 1e-6


