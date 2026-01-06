# PowerShell Command Reference

## Important: PowerShell Line Continuation

In PowerShell, use **backticks** (`) for line continuation, NOT backslashes (\).

## Correct Commands

### Single Line (Recommended)
```powershell
python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

### Multi-Line (Using Backticks)
```powershell
python engine.py --components components.csv `
  --novel exports/gdb11_multiple_runs_components.csv `
  --allow-novel 1 --max-add-novel 0.15 `
  --fuel-type gasoline --K 5 --top 15
```

## Common Commands

### 1. Run Multiple Evolutions
```powershell
python run_multiple_evolutions.py
```

### 2. Run Maximum Evolution
```powershell
python run_maximum_gdb11.py
```

### 3. Use Novel Components in Blend
```powershell
python engine.py --components components.csv --novel exports/gdb11_multiple_runs_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

### 4. Use Maximum Components
```powershell
python engine.py --components components.csv --novel exports/gdb11_maximum_components.csv --allow-novel 1 --max-add-novel 0.15 --fuel-type gasoline --K 5 --top 15
```

## Notes

- PowerShell uses backticks (`) not backslashes (\)
- No spaces after the backtick
- Or just use single-line commands (easier!)















