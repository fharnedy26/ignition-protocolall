# Quick Improvements Guide

## Immediate Actions (5-30 minutes each)

### 1. Reduce Novel Penalty (5 minutes)

**File**: `fuel_engine/constants.py`  
**Change**: Line 23
```python
NOVEL_PEN_PER_VOL = 0.02  # Reduced from 0.1 (was 10%, now 2%)
```

**Why**: Novel molecules are penalized too heavily, so optimizer avoids them even when they're better.

**Test**: Re-run engine and check if novel molecules appear in blends.

---

### 2. Run Multiple Evolutions (30 minutes)

**Command**:
```bash
python run_multiple_evolutions.py
```

**What it does**:
- Runs 5 independent evolutions with different seeds
- Combines results to get more unique molecules
- Creates `exports/gdb11_multiple_runs_components.csv` with 10-20 unique molecules

**Then use in engine**:
```bash
python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

---

### 3. Increase Diversity Rate (10 minutes)

**File**: `molecular_generator/generator.py`  
**Change**: Line 298
```python
base_diversity_rate = 0.25  # Increased from 0.15 (was 15%, now 25%)
```

**Why**: More diversity injection = more unique molecules per run.

**Test**: Re-run `run_maximum_gdb11.py` and check diversity ratio.

---

### 4. Increase Fragment Cache (5 minutes)

**File**: `run_maximum_gdb11.py`  
**Change**: Line 35
```python
fragments = build_fragment_cache(
    fragment_file="data/fragments/gdb11_fuel_filtered.smi",
    cache_size=500000,  # Increased from 200000
    filter_func=is_fuel_relevant_fragment
)
```

**Why**: Larger fragment pool = more diverse molecules.

---

## Expected Results After Quick Fixes

### Before:
- 2 unique molecules
- Novel molecules not used in blends
- 7.4% diversity

### After:
- 10-20 unique molecules (from multiple runs)
- Novel molecules may appear in blends (lower penalty)
- 15-25% diversity (higher injection rate)

---

## Next Steps After Quick Fixes

1. **Evaluate Results**:
   - Check if novel molecules appear in blends
   - Check diversity improvement
   - Check blend score improvement

2. **If Still Not Good Enough**:
   - See `IMPROVEMENT_ROADMAP.md` for Phase 2 improvements
   - Consider multi-objective optimization
   - Consider advanced diversity techniques

3. **If Good Enough**:
   - Focus on QOL improvements
   - Add visualization
   - Improve reporting

---

## Testing Checklist

After making changes, verify:

- [ ] Novel molecules appear in top blends (check blend CSV)
- [ ] Diversity ratio >15% in evolution output
- [ ] 10+ unique molecules in component CSV
- [ ] Blend score improves or stays competitive
- [ ] No errors in evolution or optimization

---

## Quick Command Reference

```bash
# 1. Run multiple evolutions (get more molecules)
python run_maximum_gdb11.py

# 2. Or run multiple independent runs
python run_multiple_evolutions.py

# 3. Use novel components in blend optimizer
python engine.py --components components.csv \
  --novel exports/gdb11_multiple_runs_components.csv \
  --allow-novel 1 --max-add-novel 0.15 \
  --fuel-type gasoline --K 5 --top 15

# 4. Check results
python compare_to_fuels.py
```















