"""
Fuel Engine - Deterministic fuel blend optimization engine.

A modular, clean, and efficient implementation for fuel blend optimization
from fixed component libraries.
"""

__version__ = "1.0.0"

from fuel_engine.optimization import greedy_then_refine, project_simplex_with_caps
from fuel_engine.data_loader import load_components, load_baseline
from fuel_engine.properties import blend_props, spec_penalties, score
from fuel_engine.io import write_top_n_csv, write_portfolio_csv, write_run_json

# GPU-accelerated batch processing (optional)
try:
    from fuel_engine.properties_gpu import blend_props_batch, blend_props_batch_gpu_available
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    blend_props_batch = None
    blend_props_batch_gpu_available = None

# Molecular descriptors (optional)
try:
    from fuel_engine.descriptors import (
        calculate_basic_descriptors,
        calculate_mordred_descriptors,
        add_descriptors_to_components,
        get_descriptor_availability
    )
    DESCRIPTORS_AVAILABLE = True
except ImportError:
    DESCRIPTORS_AVAILABLE = False
    calculate_basic_descriptors = None
    calculate_mordred_descriptors = None
    add_descriptors_to_components = None
    get_descriptor_availability = None

__all__ = [
    'greedy_then_refine',
    'project_simplex_with_caps',
    'load_components',
    'load_baseline',
    'blend_props',
    'spec_penalties',
    'score',
    'write_top_n_csv',
    'write_portfolio_csv',
    'write_run_json',
    'GPU_AVAILABLE',
    'DESCRIPTORS_AVAILABLE',
]

# Add GPU functions if available
if GPU_AVAILABLE:
    __all__.extend(['blend_props_batch', 'blend_props_batch_gpu_available'])

# Add descriptor functions if available
if DESCRIPTORS_AVAILABLE:
    __all__.extend([
        'calculate_basic_descriptors',
        'calculate_mordred_descriptors',
        'add_descriptors_to_components',
        'get_descriptor_availability'
    ])

