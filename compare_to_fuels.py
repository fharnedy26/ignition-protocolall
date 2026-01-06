#!/usr/bin/env python3
"""
Compare generated molecules to common fuels.
"""

from molecular_generator.fitness import compute_fuel_fitness_v3
from molecular_generator.molecule_to_component import estimate_fuel_properties

# Common fuels with their SMILES
common_fuels = {
    'Isooctane': ('CC(C)CC(C)(C)C', 'RON 100, LHV 44.3 MJ/kg'),
    'Ethanol': ('CCO', 'RON 108, LHV 26.8 MJ/kg'),
    'MTBE': ('CC(C)OC(C)(C)C', 'RON 118, LHV 35.2 MJ/kg'),
    'Toluene': ('Cc1ccccc1', 'RON 120, LHV 40.6 MJ/kg'),
    'Hexane': ('CCCCCC', 'RON 25, LHV 44.7 MJ/kg'),
    'Octane': ('CCCCCCCC', 'RON 0, LHV 44.4 MJ/kg'),
}

# Your generated molecule
your_molecule = 'C=C(C=O)CC#CC=O.O=CC1C=CC=C1'

print("=" * 70)
print("Fitness Score Comparison: Generated vs Common Fuels")
print("=" * 70)

print("\nCommon Fuels:")
print(f"{'Fuel':<15} {'Fitness':<10} {'Notes'}")
print("-" * 70)

fuel_scores = []
for name, (smiles, notes) in common_fuels.items():
    fitness = compute_fuel_fitness_v3(smiles)
    fuel_scores.append((name, fitness))
    print(f"{name:<15} {fitness:<10.2f} {notes}")

print("\n" + "-" * 70)
print("\nYour Generated Molecule:")
your_fitness = compute_fuel_fitness_v3(your_molecule)
props = estimate_fuel_properties(your_molecule)

print(f"  SMILES: {your_molecule}")
print(f"  Fitness Score: {your_fitness:.2f}")
print(f"\n  Properties:")
print(f"    RON: {props.get('RON', 0):.1f}")
print(f"    MON: {props.get('MON', 0):.1f}")
print(f"    LHV: {props.get('LHV_MJ_kg', 0):.2f} MJ/kg")
print(f"    PMI: {props.get('PMI', 0):.1f}")
print(f"    Density: {props.get('density_g_ml', 0):.3f} g/mL")

# Ranking
all_scores = fuel_scores + [('Your Molecule', your_fitness)]
all_scores.sort(key=lambda x: x[1], reverse=True)

print("\n" + "=" * 70)
print("Ranking (by Fitness Score):")
print("=" * 70)
for i, (name, score) in enumerate(all_scores, 1):
    marker = " <-- YOURS" if name == "Your Molecule" else ""
    print(f"{i}. {name:<20} {score:.2f}{marker}")

print("\n" + "=" * 70)
print("Analysis:")
print("=" * 70)

if your_fitness > max([s for _, s in fuel_scores]):
    print("[EXCELLENT] Your molecule has HIGHER fitness than all common fuels!")
elif your_fitness > sum([s for _, s in fuel_scores]) / len(fuel_scores):
    avg = sum([s for _, s in fuel_scores]) / len(fuel_scores)
    print(f"[GOOD] Your molecule ({your_fitness:.2f}) is ABOVE average ({avg:.2f})")
    print(f"  Better than: {', '.join([n for n, s in fuel_scores if s < your_fitness])}")
else:
    print("Your molecule is competitive with common fuels")

print(f"\nFitness Score Meaning:")
print(f"  - Maximum possible: ~50 (perfect match to all targets)")
print(f"  - Your score: {your_fitness:.2f} ({your_fitness/50*100:.1f}% of maximum)")
print(f"  - Targets: MW~150, logP~1.5, O/C~0.2, low H-donors, low rings")

print(f"\nPerformance Notes:")
print(f"  - RON {props.get('RON', 0):.1f}: {'Excellent' if props.get('RON', 0) > 100 else 'Good' if props.get('RON', 0) > 90 else 'Moderate'} for gasoline")
print(f"  - LHV {props.get('LHV_MJ_kg', 0):.2f} MJ/kg: {'Good' if props.get('LHV_MJ_kg', 0) > 38 else 'Moderate'} energy density")
print(f"  - PMI {props.get('PMI', 0):.1f}: {'Low sooting' if props.get('PMI', 0) < 15 else 'Moderate sooting' if props.get('PMI', 0) < 20 else 'Higher sooting'}")

print("=" * 70)

