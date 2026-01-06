# Plan: Home-Synthesizable Fuel Molecule Generation

## Overview
Modify the SYSTE-26 engine to generate fuel molecules that can be synthesized at home, enabling real-world testing and validation of novel fuel concepts.

## Goals
1. Generate molecules optimized for both fuel performance AND home synthesis feasibility
2. Provide synthesis route predictions for generated molecules
3. Filter molecules by synthesis complexity
4. Add stability testing prediction capabilities

---

## Implementation Plan

### Phase 1: Synthesis Feasibility Assessment

#### 1.1 Add Synthesis Feasibility Functions (`molecular_generator/fitness.py`)
- [x] `estimate_synthesis_feasibility(smiles)` - Score 0-100 for home synthesis ease
- [x] `count_functional_groups(mol)` - Helper to count reactive groups
- [x] `has_ester_group(mol)` - Detect ester functionality
- [x] `has_ether_group(mol)` - Detect ether functionality
- [x] `has_carboxylic_acid(mol)` - Detect acid groups
- [x] `has_alcohol_group(mol)` - Detect alcohol groups

**Scoring Criteria:**
- +30 points: Contains ester group (easy esterification)
- +20 points: Contains ether group (moderate Williamson synthesis)
- +20 points: Molecular weight < 200 (smaller = easier)
- -30 points: >2 rings (complex structures)
- -20 points: >3 functional groups (multi-step needed)
- -15 points: >2 double bonds (harder to control)

#### 1.2 Create Home Synthesis Fragment Library (`molecular_generator/home_synthesis_fragments.py`)
- [x] Define `HOME_SYNTHESIS_FRAGMENTS` list with common building blocks:
  - Simple alcohols (ethanol, propanol, butanol, isopropanol)
  - Simple carboxylic acids (acetic, propionic, butyric)
  - Simple aldehydes (acetaldehyde, propionaldehyde)
  - Simple ketones (acetone, butanone)
  - Simple ethers (dimethyl ether, ethyl methyl ether)
  - Alkyl chains (ethane through pentane)
- [x] `is_home_synthesizable(smiles)` - Check if molecule uses only home-accessible fragments

---

### Phase 2: Synthesis Route Prediction

#### 2.1 Create Synthesis Module (`molecular_generator/synthesis.py`)
- [x] `suggest_synthesis_route(smiles)` - Main function to suggest synthesis
- [x] `identify_reaction_type(mol)` - Detect what reaction type is needed
- [x] `get_starting_materials(mol, reaction_type)` - Suggest starting materials
- [x] `get_reaction_conditions(reaction_type)` - Provide reaction conditions
- [x] `estimate_difficulty(reaction_type, complexity)` - Rate difficulty (easy/moderate/hard)

**Reaction Types Supported:**
1. **Esterification** (Easy)
   - Starting: Alcohol + Carboxylic Acid
   - Conditions: Acid catalyst (H2SO4), heat, 1-2 hours
   - Yield: 70-90%

2. **Williamson Ether Synthesis** (Moderate)
   - Starting: Alcohol + Alkyl Halide
   - Conditions: Base (NaOH), heat, 2-4 hours
   - Yield: 60-80%

3. **Simple Reduction** (Moderate)
   - Starting: Aldehyde/Ketone
   - Conditions: NaBH4 or LiAlH4, room temp
   - Yield: 70-85%

4. **Oxidation** (Moderate-Hard)
   - Starting: Alcohol
   - Conditions: KMnO4 or CrO3, controlled conditions
   - Yield: 50-70%

---

### Phase 3: Enhanced Fitness Functions

#### 3.1 Update Fitness Functions (`molecular_generator/fitness.py`)
- [x] `compute_home_synthesis_fitness(smiles)` - Combined fitness:
  - 50% weight: Synthesis feasibility
  - 30% weight: Fuel properties
  - 20% weight: Simplicity (inverse of complexity)

- [x] `compute_stability_test_fitness(smiles, target_oxidation=None)` - For stability testing:
  - Targets specific oxidation susceptibility levels
  - Can generate high-risk (for testing) or low-risk (for controls)

#### 3.2 Add Stability Prediction Functions
- [x] `estimate_oxidation_susceptibility(smiles)` - Predict oxidation risk
- [x] `estimate_phase_separation_risk(components, vol_fractions)` - Predict phase separation
- [x] `estimate_ph_change(smiles)` - Predict pH changes

---

### Phase 4: Generator Updates

#### 4.1 Update Generator (`molecular_generator/generator.py`)
- [x] `evolve_home_synthesizable_molecules()` - New evolution function:
  - Uses constrained fragment library
  - Applies synthesis-aware fitness
  - Filters by minimum feasibility threshold
  - Returns synthesis metadata

**Parameters:**
- `min_feasibility`: Minimum feasibility score (default: 60.0)
- `synthesis_weight`: Weight of synthesis in fitness (default: 0.5)
- `fuel_weight`: Weight of fuel properties (default: 0.3)

---

### Phase 5: Component Conversion Updates

#### 5.1 Update Molecule-to-Component (`molecular_generator/molecule_to_component.py`)
- [x] Add synthesis metadata to component dictionary:
  - `synthesis_feasibility`: Feasibility score
  - `synthesis_complexity`: Complexity rating (1-10)
  - `synthesis_route`: Suggested synthesis route
  - `synthesis_difficulty`: easy/moderate/hard
  - `oxidation_susceptibility`: Oxidation risk score
  - `stability_notes`: Stability-related notes

---

### Phase 6: Notebook Updates

#### 6.1 Update Unified Workflow Notebook (`notebooks/unified_workflow.ipynb`)
- [x] Add cell for home-synthesis mode
- [x] Display synthesis routes for top molecules
- [x] Filter and rank by synthesis feasibility
- [x] Export synthesis instructions
- [x] Add stability testing candidate selection

---

## File Structure

```
molecular_generator/
├── fitness.py                    [MODIFY] Add synthesis & stability functions
├── generator.py                  [MODIFY] Add home-synthesis evolution
├── molecule_to_component.py      [MODIFY] Add synthesis metadata
├── home_synthesis_fragments.py   [NEW] Constrained fragment library
└── synthesis.py                  [NEW] Synthesis route prediction

notebooks/
└── unified_workflow.ipynb        [MODIFY] Add synthesis-aware workflow
```

---

## Usage Example

```python
from molecular_generator import evolve_home_synthesizable_molecules
from molecular_generator.home_synthesis_fragments import HOME_SYNTHESIS_FRAGMENTS

# Generate home-synthesizable molecules
df, history = evolve_home_synthesizable_molecules(
    fragments=HOME_SYNTHESIS_FRAGMENTS,
    generations=50,
    min_feasibility=60.0,
    synthesis_weight=0.5,
    fuel_weight=0.3
)

# Each molecule will have:
# - Synthesis feasibility score
# - Suggested synthesis route
# - Starting materials
# - Reaction conditions
# - Difficulty rating
```

---

## Testing Strategy

1. **Unit Tests:**
   - Test synthesis feasibility scoring on known molecules
   - Test route prediction for esters, ethers, alcohols
   - Verify fragment library constraints

2. **Integration Tests:**
   - Generate molecules and verify synthesis metadata
   - Check that feasibility filtering works
   - Validate synthesis routes are reasonable

3. **Validation:**
   - Compare predicted vs. actual synthesis difficulty
   - Verify molecules can be made from suggested routes
   - Test stability predictions against known data

---

## Safety Considerations

⚠️ **Important Notes:**
- All synthesis suggestions assume proper safety equipment
- Users must understand chemical hazards
- Proper ventilation and waste disposal required
- Some reactions may require permits or professional supervision
- Engine will flag high-risk molecules but user discretion required

---

## Future Enhancements

1. **Multi-step Synthesis Planning:**
   - Suggest 2-3 step syntheses for complex molecules
   - Identify intermediate compounds

2. **Cost Estimation:**
   - Estimate cost of starting materials
   - Calculate cost per liter of fuel

3. **Yield Optimization:**
   - Suggest conditions for maximum yield
   - Identify side products

4. **Safety Scoring:**
   - Rate safety of each synthesis route
   - Flag hazardous reactions

5. **Equipment Requirements:**
   - List required equipment (glassware, heating, etc.)
   - Estimate setup time

---

## Implementation Status

- [x] Plan created
- [ ] Phase 1: Synthesis feasibility functions
- [ ] Phase 2: Synthesis route prediction
- [ ] Phase 3: Enhanced fitness functions
- [ ] Phase 4: Generator updates
- [ ] Phase 5: Component conversion updates
- [ ] Phase 6: Notebook updates
- [ ] Testing and validation




















