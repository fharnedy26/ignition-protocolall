"""
Molecular Generator Module - Genetic algorithm for novel fuel molecule discovery.

This module provides functionality for generating novel fuel molecules using
genetic algorithms with fragment-based assembly and evolution.
"""

__version__ = "1.0.0"

from molecular_generator.fragments import (
    load_fragments,
    clean_fragment_list,
    generate_molecule_from_fragments,
    merge_fragments,
    is_fuel_relevant_fragment,
    load_fragments_streaming,
    preprocess_gdb11,
    sample_fragments_reservoir,
    build_fragment_cache,
)

from molecular_generator.fitness import (
    compute_fuel_fitness,
    compute_fuel_fitness_v3,
    compute_home_synthesis_fitness,
    compute_stability_test_fitness,
    calculate_molecular_properties,
    estimate_synthesis_feasibility,
    estimate_oxidation_susceptibility,
    estimate_phase_separation_risk,
    estimate_ph_change,
)

from molecular_generator.generator import (
    evolve_population,
    evolve_with_fragments,
    evolve_home_synthesizable_molecules,
    mutate_smiles,
    is_valid_fuel,
)

from molecular_generator.molecule_to_component import (
    molecule_to_component,
    molecules_to_component_csv,
    estimate_fuel_properties,
)

from molecular_generator.synthesis import (
    suggest_synthesis_route,
    format_synthesis_instructions,
)

from molecular_generator.home_synthesis_fragments import (
    HOME_SYNTHESIS_FRAGMENTS,
    is_home_synthesizable,
    get_fragments_by_category,
)

__all__ = [
    # Fragment handling
    'load_fragments',
    'clean_fragment_list',
    'generate_molecule_from_fragments',
    'merge_fragments',
    'is_fuel_relevant_fragment',
    'load_fragments_streaming',
    'preprocess_gdb11',
    'sample_fragments_reservoir',
    'build_fragment_cache',
    # Fitness functions
    'compute_fuel_fitness',
    'compute_fuel_fitness_v3',
    'compute_home_synthesis_fitness',
    'compute_stability_test_fitness',
    'calculate_molecular_properties',
    'estimate_synthesis_feasibility',
    'estimate_oxidation_susceptibility',
    'estimate_phase_separation_risk',
    'estimate_ph_change',
    # Generation
    'evolve_population',
    'evolve_with_fragments',
    'evolve_home_synthesizable_molecules',
    'mutate_smiles',
    'is_valid_fuel',
    # Component conversion
    'molecule_to_component',
    'molecules_to_component_csv',
    'estimate_fuel_properties',
    # Synthesis
    'suggest_synthesis_route',
    'format_synthesis_instructions',
    # Home synthesis
    'HOME_SYNTHESIS_FRAGMENTS',
    'is_home_synthesizable',
    'get_fragments_by_category',
]

