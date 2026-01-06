# Quick Start Guide—Testing Improvements

## Step 1: Generate More Novel Molecules

First, run multiple evolutions to obtain 10–20 unique molecules:

```powershell
python run_multiple_evolutions.py
```

This will:
- Run 10 independent evolutions
- Generate 20–30 unique molecules
- Create `exports/gdb11_multiple_runs_components.csv`

**Time**: ~10–15 minutes

---

## Step 2: Utilise Novel Molecules in Blend

Once additional novel molecules are available, utilise them:

```powershell
python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

---

## Step 3: Examine Results

```powershell
python check_blend.py
```

This will display:
- Which components were utilised
- Whether novel molecules appear
- Available novel molecules

---

## Current Status

**Issue**: Only 1 novel molecule available, and it is not competitive enough:
- Novel: RON 101.9, PMI 19.2
- Best Blend: RON 110.8, PMI 14.4

**Solution**: Run multiple evolutions to obtain more diverse, potentially superior molecules.

---

## PowerShell Command Reference

**Single line (easiest)**:
```powershell
python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

**Multi-line (using backticks)**:
```powershell
python engine.py --components components.csv `
  --novel exports/gdb11_multiple_runs_components.csv `
  --allow-novel 1 --max-add-novel 0.15 `
  --fuel-type gasoline --K 5 --top 15
```

**Note**: PowerShell uses backticks (`) for line continuation, NOT backslashes (\).















