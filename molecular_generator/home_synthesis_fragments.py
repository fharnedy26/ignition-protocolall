"""
Home-synthesizable molecular fragments library.

Provides a curated list of molecular building blocks that are:
1. Commercially available or easily accessible
2. Can be used in simple home synthesis reactions
3. Safe to handle with proper precautions
"""

# Simple alcohols (common starting materials)
SIMPLE_ALCOHOLS = [
    'CCO',      # Ethanol
    'CCCO',     # Propanol (1-propanol)
    'CCCCO',    # Butanol (1-butanol)
    'CC(C)O',   # Isopropanol (2-propanol)
    'CC(C)(C)O', # tert-Butanol
    'CCCCCO',   # Pentanol (1-pentanol)
    'CC(C)CO',  # Isobutanol
]

# Simple carboxylic acids
SIMPLE_ACIDS = [
    'CC(=O)O',      # Acetic acid
    'CCC(=O)O',     # Propionic acid
    'CCCC(=O)O',    # Butyric acid
    'CCCCC(=O)O',   # Valeric acid (pentanoic acid)
    'OC(=O)C',      # Formic acid
]

# Simple aldehydes
SIMPLE_ALDEHYDES = [
    'CC=O',     # Acetaldehyde
    'CCC=O',    # Propionaldehyde
    'CCCC=O',   # Butyraldehyde
]

# Simple ketones
SIMPLE_KETONES = [
    'CC(=O)C',      # Acetone
    'CC(=O)CC',     # Butanone (MEK)
    'CC(=O)CCC',    # Pentanone
]

# Simple ethers (can be synthesized)
SIMPLE_ETHERS = [
    'COC',      # Dimethyl ether
    'CCOC',     # Ethyl methyl ether
    'CCOCC',    # Diethyl ether
]

# Alkyl chains (for building)
ALKYL_CHAINS = [
    'CC',       # Ethane
    'CCC',      # Propane
    'CCCC',     # Butane
    'CCCCC',    # Pentane
    'CCCCCC',   # Hexane
]

# Esters (common and can be made via esterification)
COMMON_ESTERS = [
    'CC(=O)OCC',        # Ethyl acetate
    'CC(=O)OCCCC',      # Butyl acetate
    'CCC(=O)OCC',       # Ethyl propionate
    'CCCC(=O)OCC',      # Ethyl butyrate
]

# All home-synthesizable fragments
HOME_SYNTHESIS_FRAGMENTS = (
    SIMPLE_ALCOHOLS +
    SIMPLE_ACIDS +
    SIMPLE_ALDEHYDES +
    SIMPLE_KETONES +
    SIMPLE_ETHERS +
    ALKYL_CHAINS +
    COMMON_ESTERS
)

# Fragment categories for targeted generation
FRAGMENT_CATEGORIES = {
    'alcohols': SIMPLE_ALCOHOLS,
    'acids': SIMPLE_ACIDS,
    'aldehydes': SIMPLE_ALDEHYDES,
    'ketones': SIMPLE_KETONES,
    'ethers': SIMPLE_ETHERS,
    'alkyls': ALKYL_CHAINS,
    'esters': COMMON_ESTERS,
}


def is_home_synthesizable(smiles: str) -> bool:
    """
    Check if a molecule can be made from home-accessible fragments.
    
    This is a simple check - a molecule is considered synthesizable if:
    1. It's in the fragment library (already available)
    2. It can be made via 1-2 simple reactions from fragments
    
    Args:
        smiles: SMILES string to check
        
    Returns:
        True if molecule is likely synthesizable at home
    """
    from rdkit import Chem
    
    # If it's already in our fragment library, it's available
    if smiles in HOME_SYNTHESIS_FRAGMENTS:
        return True
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        
        # Check if it's a simple combination of known fragments
        # This is a heuristic - molecules with simple structures are more likely
        # to be synthesizable
        
        # Count functional groups
        ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O')))
        ether_count = len(mol.GetSubstructMatches(Chem.MolFromSmiles('COC')))
        alcohol_count = len(mol.GetSubstructMatches(Chem.MolFromSmiles('CO')))
        
        # Simple molecules with common functional groups are likely synthesizable
        if ester_count > 0 or ether_count > 0 or alcohol_count > 0:
            # Check molecular weight (smaller = easier)
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            if mw < 200:  # Reasonable size for home synthesis
                return True
        
        return False
    except:
        return False


def get_fragments_by_category(category: str) -> list:
    """
    Get fragments from a specific category.
    
    Args:
        category: One of 'alcohols', 'acids', 'aldehydes', 'ketones', 
                  'ethers', 'alkyls', 'esters'
    
    Returns:
        List of SMILES strings for that category
    """
    return FRAGMENT_CATEGORIES.get(category, [])




















