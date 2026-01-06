Run terminal command to execute the engine
// Log the current command to be run and check if it works
import subprocess
import sys
from pathlib import Path

# Default command with components file
command = "python engine.py --components data/components.csv --fuel-type gasoline --K 5 --top 10"

# Check if components file exists (try data/ first, then root)
components_file = Path("data/components.csv")
if not components_file.exists():
    components_file = Path("components.csv")
if not components_file.exists():
    print(f"[WARN] Components file '{components_file}' not found.")
    print("Falling back to selftest mode...")
    command = "python engine.py --selftest"
else:
    print(f"[INFO] Using components file: {components_file}")

print(f"\n[INFO] About to run command: {command}\n")
print("=" * 60)

try:
    result = subprocess.run(
        command, 
        shell=True, 
        capture_output=True, 
        text=True, 
        check=True,
        timeout=300  # 5 minute timeout
    )
    print("=" * 60)
    print("\n[SUCCESS] Command executed successfully!\n")
    print("Command output:")
    print("-" * 60)
    print(result.stdout)
    if result.stderr:
        print("\n[WARN] Stderr output:")
        print("-" * 60)
        print(result.stderr)
    print("\n[SUCCESS] Engine run completed successfully.")
    
except subprocess.CalledProcessError as e:
    print("=" * 60)
    print("\n[ERROR] Command failed with return code:", e.returncode)
    print("-" * 60)
    if e.stdout:
        print("Stdout output:")
        print(e.stdout)
    if e.stderr:
        print("\nStderr output:")
        print(e.stderr)
    if e.output:
        print("\nError output:")
        print(e.output)
    print("\n[ERROR] Engine run failed.")
    sys.exit(e.returncode)
    
except subprocess.TimeoutExpired:
    print("=" * 60)
    print("\n[ERROR] Command timed out after 5 minutes.")
    sys.exit(1)
    
except Exception as e:
    print("=" * 60)
    print(f"\n[ERROR] Unexpected error: {type(e).__name__}: {e}")
    sys.exit(1)

print("\n" + "=" * 60)
print("\nYou can also run the engine directly in your terminal:")
print(f"\n    {command}\n")
print("Or with custom arguments:")
print("    python engine.py --components data/components.csv --fuel-type gasoline --K 5 --top 10 --restarts 2\n")
