"""
Synthesis route prediction for home-synthesizable molecules.

Provides functions to suggest synthesis routes, starting materials,
and reaction conditions for molecules that can be made at home.
"""

from typing import Dict, List, Optional, Any
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from molecular_generator.fitness import (
    has_ester_group,
    has_ether_group,
    has_alcohol_group,
    has_carboxylic_acid,
    estimate_synthesis_feasibility
)


def identify_reaction_type(mol) -> Optional[str]:
    """
    Identify the primary reaction type needed to synthesize the molecule.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Reaction type string or None
    """
    if has_ester_group(mol):
        return 'esterification'
    elif has_ether_group(mol):
        return 'williamson_ether'
    elif has_alcohol_group(mol):
        # Could be reduction or hydration
        return 'reduction'
    elif has_carboxylic_acid(mol):
        return 'oxidation'
    else:
        return None


def get_starting_materials(mol, reaction_type: str) -> List[str]:
    """
    Suggest starting materials for a given reaction type.
    
    Args:
        mol: RDKit molecule object
        reaction_type: Type of reaction needed
        
    Returns:
        List of suggested starting material names (not SMILES, as we'd need
        to decompose the molecule which is complex)
    """
    if reaction_type == 'esterification':
        return ['alcohol', 'carboxylic_acid']
    elif reaction_type == 'williamson_ether':
        return ['alcohol', 'alkyl_halide']
    elif reaction_type == 'reduction':
        return ['aldehyde', 'ketone']
    elif reaction_type == 'oxidation':
        return ['alcohol']
    else:
        return ['unknown_starting_materials']


def get_reaction_conditions(reaction_type: str) -> Dict[str, str]:
    """
    Get reaction conditions for a given reaction type.
    
    Args:
        reaction_type: Type of reaction
        
    Returns:
        Dictionary with reaction conditions
    """
    conditions = {
        'esterification': {
            'catalyst': 'H2SO4 (concentrated)',
            'temperature': '60-80°C',
            'time': '1-2 hours',
            'equipment': 'Round-bottom flask, condenser, heating mantle',
            'notes': 'Use excess alcohol or acid to drive equilibrium'
        },
        'williamson_ether': {
            'base': 'NaOH or KOH',
            'solvent': 'Water or ethanol',
            'temperature': '60-80°C',
            'time': '2-4 hours',
            'equipment': 'Round-bottom flask, condenser, heating mantle',
            'notes': 'Alkyl halide should be primary for best yield'
        },
        'reduction': {
            'reagent': 'NaBH4 (sodium borohydride) or LiAlH4',
            'solvent': 'Ethanol or THF',
            'temperature': 'Room temperature (NaBH4) or 0°C (LiAlH4)',
            'time': '30 minutes - 2 hours',
            'equipment': 'Round-bottom flask, ice bath (for LiAlH4)',
            'notes': 'NaBH4 is safer and easier to handle than LiAlH4'
        },
        'oxidation': {
            'reagent': 'KMnO4 or CrO3 (Jones reagent)',
            'solvent': 'Water or acetone',
            'temperature': 'Room temperature to 60°C',
            'time': '1-3 hours',
            'equipment': 'Round-bottom flask, condenser',
            'notes': 'Requires careful control to avoid over-oxidation'
        }
    }
    
    return conditions.get(reaction_type, {
        'notes': 'Complex multi-step synthesis required'
    })


def estimate_difficulty(reaction_type: Optional[str], complexity: float) -> str:
    """
    Estimate synthesis difficulty.
    
    Args:
        reaction_type: Type of reaction
        complexity: Complexity score from synthesis feasibility (1-10)
        
    Returns:
        Difficulty rating: 'easy', 'moderate', or 'hard'
    """
    if complexity <= 3.0:
        return 'easy'
    elif complexity <= 6.0:
        return 'moderate'
    else:
        return 'hard'


def estimate_yield(reaction_type: Optional[str], complexity: float) -> str:
    """
    Estimate expected yield range.
    
    Args:
        reaction_type: Type of reaction
        complexity: Complexity score
        
    Returns:
        Yield estimate as string
    """
    if reaction_type == 'esterification':
        if complexity <= 3.0:
            return '70-90%'
        else:
            return '50-70%'
    elif reaction_type == 'williamson_ether':
        if complexity <= 4.0:
            return '60-80%'
        else:
            return '40-60%'
    elif reaction_type == 'reduction':
        return '70-85%'
    elif reaction_type == 'oxidation':
        return '50-70%'
    else:
        if complexity <= 3.0:
            return '60-80%'
        elif complexity <= 6.0:
            return '40-60%'
        else:
            return '20-40%'


def suggest_synthesis_route(smiles: str) -> Dict[str, Any]:
    """
    Suggest a simple synthesis route for home synthesis.
    
    Args:
        smiles: SMILES string of the target molecule
        
    Returns:
        Dictionary with synthesis route information:
        - reaction_type: Type of reaction needed
        - starting_materials: List of starting material names
        - conditions: Dictionary with reaction conditions
        - difficulty: 'easy', 'moderate', or 'hard'
        - yield_estimate: Expected yield range
        - notes: Additional notes
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {
                'reaction_type': 'unknown',
                'difficulty': 'hard',
                'note': 'Invalid molecule structure'
            }
        
        # Get synthesis feasibility
        synth = estimate_synthesis_feasibility(smiles)
        
        # Identify reaction type
        reaction_type = identify_reaction_type(mol)
        
        if not reaction_type:
            return {
                'reaction_type': 'unknown',
                'difficulty': 'hard',
                'note': 'Multi-step synthesis required - too complex for simple home synthesis'
            }
        
        # Get starting materials
        starting_materials = get_starting_materials(mol, reaction_type)
        
        # Get reaction conditions
        conditions = get_reaction_conditions(reaction_type)
        
        # Estimate difficulty
        difficulty = estimate_difficulty(reaction_type, synth['complexity'])
        
        # Estimate yield
        yield_estimate = estimate_yield(reaction_type, synth['complexity'])
        
        # Build response
        route = {
            'reaction_type': reaction_type,
            'starting_materials': starting_materials,
            'conditions': conditions,
            'difficulty': difficulty,
            'yield_estimate': yield_estimate,
            'complexity_score': synth['complexity'],
            'feasibility_score': synth['feasibility']
        }
        
        # Add notes based on complexity
        if synth['complexity'] > 6.0:
            route['note'] = 'This molecule may require professional laboratory equipment and expertise'
        elif synth['complexity'] > 4.0:
            route['note'] = 'Moderate complexity - ensure proper safety equipment and ventilation'
        else:
            route['note'] = 'Relatively straightforward synthesis with proper precautions'
        
        return route
        
    except Exception as e:
        return {
            'reaction_type': 'unknown',
            'difficulty': 'hard',
            'note': f'Error analyzing molecule: {str(e)}'
        }


def format_synthesis_instructions(route: Dict[str, Any]) -> str:
    """
    Format synthesis route as human-readable instructions.
    
    Args:
        route: Synthesis route dictionary from suggest_synthesis_route()
        
    Returns:
        Formatted string with synthesis instructions
    """
    if route['reaction_type'] == 'unknown':
        return f"⚠️ {route.get('note', 'Unknown synthesis route')}"
    
    lines = [
        f"**Synthesis Route: {route['reaction_type'].replace('_', ' ').title()}**",
        "",
        f"**Difficulty:** {route['difficulty'].upper()}",
        f"**Expected Yield:** {route['yield_estimate']}",
        f"**Feasibility Score:** {route['feasibility_score']:.1f}/100",
        "",
        "**Starting Materials:**",
    ]
    
    for material in route['starting_materials']:
        lines.append(f"  - {material.replace('_', ' ').title()}")
    
    lines.append("")
    lines.append("**Reaction Conditions:**")
    
    conditions = route['conditions']
    for key, value in conditions.items():
        lines.append(f"  - {key.replace('_', ' ').title()}: {value}")
    
    if 'note' in route:
        lines.append("")
        lines.append(f"**Note:** {route['note']}")
    
    return "\n".join(lines)




















