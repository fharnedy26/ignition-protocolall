# Ignition Protocol: Fuel Blend Optimization Engine

A deterministic, fast, lab-ready engine for optimizing fuel blends from component libraries. Combines **molecular generation** (genetic algorithms for novel fuel discovery) with **blend optimization** (deterministic optimization of fuel mixtures) to create an end-to-end fuel design system.

## Overview

The Ignition Protocol engine is a production-ready fuel blend optimization system with two core components:

1. **Blend Optimizer Engine** (`engine.py`): Deterministic, fast optimization engine that produces lab-ready fuel blend recipes
2. **Molecular Generator**: Genetic algorithm-based system that evolves novel fuel molecules from fragment libraries (GDB-11)

The engine enables a complete workflow: **Generate → Convert → Optimize → Validate**

## Engine Features

### Blend Optimization Engine
- **Deterministic**: Same seed produces identical results - perfect for reproducibility
- **Fast**: 0.1-0.3s for typical component libraries
- **Lab-ready**: Volume fractions ready for immediate mixing
- **Fuel-type specific**: Gasoline, diesel, jet fuel specifications built-in
- **Delta mode**: Small edits around existing blends with L1 budget control
- **Novel components**: Optional micro-additives with configurable caps
- **Multi-objective**: Optimizes RON, cetane, LHV, PMI, cost simultaneously
- **Production-ready**: Handles component constraints, spec penalties, and waste credits

### Molecular Generation (Optional)
- **Fragment-based evolution**: Supports GDB-11 database (gdb11.tgz) with up to 11 atoms
- **GDB-11 integration**: Memory-efficient streaming, fuel-relevance filtering (C/H/O only)
- **Genetic algorithms**: Crossover, mutation, and selection for molecule evolution
- **Fitness functions**: Multi-objective optimization (energy, oxygen balance, volatility, H-bonding)
- **100% chemical validity**: RDKit validation ensures all molecules are chemically valid

## Project Structure

```
syste-26-main/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── engine.py                    # Blend optimizer CLI
├── fuel_engine/                 # Blend optimization module
│   ├── __init__.py
│   ├── constants.py
│   ├── data_loader.py
│   ├── descriptors.py
│   ├── io.py
│   ├── optimization.py
│   ├── properties.py
│   └── properties_gpu.py
├── molecular_generator/         # Molecular generation module
│   ├── __init__.py
│   ├── generator.py            # Core genetic algorithm
│   ├── fitness.py              # Fitness functions
│   ├── fragments.py            # Fragment handling
│   └── molecule_to_component.py # Bridge: convert molecules to components
├── notebooks/                   # Jupyter notebooks
│   ├── unified_workflow.ipynb  # End-to-end workflow demo
│   ├── molecular-fuel-generator.ipynb
│   └── ... (other versions)
├── data/                        # Data files
│   ├── components.csv
│   ├── components_novel.csv
│   ├── baseline_blend.csv
│   └── fragments/              # Fragment libraries
│       ├── gdb11_fuel_filtered.smi  # Preprocessed GDB-11 (C/H/O only)
│       ├── high_entropy_fragments.txt
│       ├── cleaned_fragments.txt
│       └── mergeable_fragments.txt
├── scripts/                     # Utility scripts
│   ├── clean_outputs.py
│   └── preprocess_gdb11.py     # GDB-11 preprocessing utility
├── runs/                        # Optimization run outputs
└── exports/                     # Molecular generation exports
```

## Installation

```bash
# Install dependencies
pip install -r requirements.txt
```

### Key Dependencies

- **Core**: numpy, pandas, scipy
- **Molecular**: rdkit-pypi, mordred
- **Optimization**: scikit-optimize, numba
- **GPU** (optional): jax, jax[cuda12]
- **Visualization**: matplotlib, Pillow
- **Parallel**: joblib

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the Engine

The engine is ready to use immediately with the provided component library:

```bash
# Basic gasoline optimization
python engine.py --components components.csv --fuel-type gasoline --K 5 --top 10 --seed 42

# Diesel optimization
python engine.py --components components.csv --fuel-type diesel --K 4 --top 15 --seed 123

# Jet fuel optimization
python engine.py --components components.csv --fuel-type jet --K 5 --top 20 --seed 456
```

### 3. View Results

Results are saved to `runs/<timestamp>_<seed>/`:
- `Top-N.csv`: Optimized blends with properties and lab recipes
- `RUN.json`: Complete run metadata

### 4. Advanced Usage

```bash
# With novel components
python engine.py --components components.csv --novel components_novel.csv --allow-novel 1 --max-add-novel 0.1 --fuel-type gasoline --K 5 --top 10

# Delta mode (small edits around baseline)
python engine.py --components components.csv --baseline baseline_blend.csv --fuel-type gasoline --K 5 --top 15 --delta-budget 0.2

# Using presets
python engine.py --components components.csv --preset baseline --fuel-type gasoline --top 5
```

### 5. Molecular Generation (Optional)

Generate novel fuel molecules using genetic algorithms:

```python
from molecular_generator import load_fragments, evolve_with_fragments, molecules_to_component_csv

# Load fragments and evolve
fragments = load_fragments("data/fragments/high_entropy_fragments.txt", max_fragments=10000)
evolution_df, history = evolve_with_fragments(fragments, generations=100, population_size=50)

# Convert to components and use in engine
top_molecules = evolution_df.nlargest(10, 'Fitness')['Best_SMILES'].tolist()
molecules_to_component_csv(top_molecules, "components_generated.csv", cost_eur_L=2.0, max_vol_frac=0.1, is_novel=1)
```

See `notebooks/unified_workflow.ipynb` for a complete end-to-end example.

## Engine Usage Examples

### Basic Optimization

```bash
# Gasoline optimization
python engine.py --components data/components.csv --fuel-type gasoline --K 3 --top 10 --seed 42

# Diesel optimization  
python engine.py --components data/components.csv --fuel-type diesel --K 4 --top 15 --seed 123

# Jet fuel optimization
python engine.py --components data/components.csv --fuel-type jet --K 5 --top 20 --seed 456

# Delta mode (small edits around baseline)
python engine.py --components data/components.csv --baseline data/baseline_blend.csv --fuel-type gasoline --K 5 --top 15 --delta-budget 0.2

# With novel components
python engine.py --components data/components.csv --novel data/components_novel.csv --allow-novel 1 --max-add-novel 0.1 --fuel-type gasoline --K 5 --top 10
```

### Presets

```bash
# Baseline preset (strict specs, no novel)
python engine.py --components data/components.csv --preset baseline --fuel-type gasoline --top 5

# Delta preset (small changes around baseline)
python engine.py --components data/components.csv --baseline data/baseline_blend.csv --preset delta --fuel-type gasoline --top 5

# Novel preset (allow up to 3% novel components)
python engine.py --components data/components.csv --novel data/components_novel.csv --preset novel3 --fuel-type gasoline --top 5
```

## Input Files

### Data File Locations
Data files can be located in either the project root directory or the `data/` subdirectory. The engine will accept paths relative to the current working directory.

**Example locations:**
- `components.csv` (root) or `data/components.csv`
- `baseline_blend.csv` (root) or `data/baseline_blend.csv`
- `components_novel.csv` (root) or `data/components_novel.csv`

### Component Library (`components.csv` or `data/components.csv`)
Required columns:
- `id`: Unique component identifier
- `name`: Component name
- `density_g_ml`: Density in g/mL
- `LHV_MJ_kg`: Lower heating value in MJ/kg
- `RON`, `MON`: Octane numbers
- `cetane`: Cetane number (for diesel)
- `PMI`, `TSI`: Sooting indices
- `cost_eur_L`: Cost per liter
- `max_vol_frac`: Maximum volume fraction allowed
- `is_novel`: 0 or 1 (novel component flag)

**Note**: Example files (`components.csv`, `components_novel.csv`, `baseline_blend.csv`) are provided in the project root. You can use these as templates or move them to the `data/` directory.

### Baseline Blend (`baseline_blend.csv` or `data/baseline_blend.csv`)
Optional CSV with columns: `id`, `vol_frac`

### Novel Components (`components_novel.csv` or `data/components_novel.csv`)
Same format as component library, with `is_novel=1`

### Fragment Libraries (`data/fragments/`)
Text files with one SMILES string per line. Supports:
- Plain text files (.txt, .smi)
- Gzipped files (.gz)
- Tar.gz archives (.tgz) - automatically extracts and processes

## GDB-11 Integration

The system now supports GDB-11 database (gdb11.tgz), which contains molecules with up to 11 atoms and is better suited for fuel applications than larger databases.

### Preprocessing GDB-11

Preprocess GDB-11 to extract only fuel-relevant fragments (C/H/O only):

```bash
# Preprocess gdb11.tgz to extract fuel-relevant fragments
python scripts/preprocess_gdb11.py \
    --input gdb11.tgz \
    --output data/fragments/gdb11_fuel_filtered.smi \
    --progress-interval 100000
```

This creates a filtered file containing only molecules with C, H, O atoms, which can be loaded much faster than processing the archive each time.

### Using GDB-11 in Code

```python
from molecular_generator import (
    load_fragments,
    build_fragment_cache,
    preprocess_gdb11,
    evolve_with_fragments
)

# Option 1: Use preprocessed file (fastest)
fragments = load_fragments("data/fragments/gdb11_fuel_filtered.smi", max_fragments=10000)

# Option 2: Stream directly from archive (slower, no preprocessing needed)
from molecular_generator.fragments import is_fuel_relevant_fragment
fragments = build_fragment_cache(
    "gdb11.tgz",
    cache_size=100000,
    filter_func=is_fuel_relevant_fragment
)

# Option 3: Pass file path directly to evolve_with_fragments
evolution_df, history = evolve_with_fragments(
    fragments="gdb11.tgz",  # Automatically builds cache with fuel filter
    generations=100,
    population_size=50
)
```

### GDB-11 Features

- **Memory-efficient streaming**: Processes large archives without loading everything into memory
- **Fuel-relevance filtering**: Automatically filters for C/H/O-only molecules
- **Tar.gz support**: Handles compressed archives seamlessly
- **Reservoir sampling**: Uniform random sampling from large datasets
- **Cache building**: Builds fixed-size caches for efficient random access

## Output Files

### Blend Optimization Outputs

- `runs/<timestamp>_<seed>/Top-N.csv`: Optimized blends with:
  - Properties: RON, MON, cetane, LHV, PMI, TSI, etc.
  - Component fractions
  - Lab recipes (mL and g per 1000mL)
  - Cost analysis
  - Novel component usage

- `runs/<timestamp>_<seed>/RUN.json`: Complete run metadata

### Molecular Generation Outputs

- `exports/top_molecules.csv`: Top evolved molecules with properties
- `exports/fitness_history.csv`: Fitness evolution over generations
- `exports/fitness_progress.png`: Fitness curve visualization

## Command Line Options

```bash
python engine.py --help
```

Key options:
- `--components`: Component CSV file (required)
- `--novel`: Novel components CSV file
- `--baseline`: Baseline blend CSV file
- `--fuel-type`: gasoline | diesel | jet
- `--K`: Maximum components in blend
- `--top`: Number of top results
- `--preset`: baseline | delta | novel3
- `--allow-novel`: 0 or 1
- `--max-add-novel`: Maximum novel fraction (e.g., 0.1)
- `--delta-budget`: L1 deviation budget for delta mode
- `--seed`: Random seed for reproducibility

## Advanced Features

### GPU Acceleration

The blend optimizer supports GPU acceleration via JAX:

```bash
# Install GPU support
pip install jax[cuda12]
```

### Molecular Descriptors

Calculate advanced molecular descriptors:

```python
from fuel_engine.descriptors import add_descriptors_to_components

# Add descriptors to component library
components = add_descriptors_to_components(components_df)
```

### Custom Fitness Functions

Modify fitness functions in `molecular_generator/fitness.py` to optimize for different criteria.

## Scientific Background

### Molecular Generation
- Uses genetic algorithms to evolve molecules from fragment libraries
- Fitness considers: molecular weight, logP, O/C ratio, H-bonding, ring count
- Preserves diversity through Tanimoto similarity metrics
- Validates all molecules using RDKit

### Blend Optimization
- Deterministic greedy-then-refine algorithm
- Multi-objective scoring: RON, cetane, LHV, PMI, cost
- Fuel-type specific constraints (gasoline, diesel, jet)
- Handles component caps, novel limits, and delta budgets

## Testing

```bash
# Run self-test
python engine.py --selftest
```

## Cleanup & Organization

```bash
# Move CSVs to data/ directory
python scripts/clean_outputs.py --move-csvs 1 --root .

# Prune old runs
python engine.py --components data/components.csv --clean-runs 14 ...
```

## Troubleshooting

### Common Issues

#### "Components file not found"
- **Solution**: Ensure `components.csv` exists in the project root or `data/` directory
- Check the file path is correct (use absolute path if needed)
- Verify file permissions allow reading

#### "ModuleNotFoundError" or import errors
- **Solution**: Install dependencies: `pip install -r requirements.txt`
- For GPU support: `pip install jax[cuda12]` (optional)
- Verify Python version: Python 3.8+ required

#### "Invalid components file" error
- **Solution**: Check CSV format matches required columns (see Input Files section)
- Ensure CSV is properly formatted (no missing commas, proper escaping)
- Verify required columns exist: `id`, `density_g_ml`, `LHV_MJ_kg`, `RON`, `MON`, etc.

#### GPU acceleration not working
- **Solution**: GPU support is optional - engine works on CPU
- For JAX GPU: Install CUDA 12.x and `jax[cuda12]`
- Engine automatically falls back to CPU if GPU unavailable

#### Large data files (gdb11.tgz) causing issues
- **Note**: Large data files should not be committed to Git (see `.gitignore`)
- Preprocess GDB-11 to create filtered fragments (see GDB-11 Integration section)
- Use preprocessed files for faster loading

### Getting Help

1. Run self-test: `python engine.py --selftest`
2. Check error messages - they often indicate the specific issue
3. Verify file paths and formats match documentation
4. Ensure all dependencies are installed

## Contributing

This project combines two complementary systems. When adding features:

- **Molecular generation**: Add to `molecular_generator/`
- **Blend optimization**: Add to `fuel_engine/`
- **Bridge functionality**: Add to `molecular_generator/molecule_to_component.py`

## License

Developed for educational and research purposes.

## References

- **RDKit**: Open-source cheminformatics toolkit
- **GDB-11**: Chemical space database (University of Bern) - molecules with up to 11 atoms, optimized for fuel applications
- **JAX**: GPU-accelerated numerical computing
- **scikit-optimize**: Bayesian optimization

---

*Ignition Protocol: A deterministic, production-ready engine for fuel blend optimization and molecular discovery.*
