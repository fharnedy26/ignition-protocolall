"""
Fragment handling for molecular generation.

Provides functions for loading, cleaning, and merging molecular fragments
from databases like GDB-11.
"""

import random
import gzip
import tarfile
from pathlib import Path
from typing import List, Optional, Iterator, Callable
from rdkit import Chem

# Silence RDKit warnings (handle different RDKit versions)
try:
    from rdkit.Chem import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        # RDLogger not available in this RDKit version, warnings will show
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)


def load_fragments(file_path: str, max_fragments: Optional[int] = None) -> List[str]:
    """
    Load molecular fragments from a file.
    
    Args:
        file_path: Path to fragment file (supports .txt, .smi, .gz)
        max_fragments: Maximum number of fragments to load (None for all)
    
    Returns:
        List of SMILES strings representing fragments
    """
    fragments = []
    path = Path(file_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Fragment file not found: {file_path}")
    
    # Handle gzipped files
    if path.suffix == '.gz' or '.gz' in path.name:
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(path, mode) as f:
        for i, line in enumerate(f):
            if max_fragments and i >= max_fragments:
                break
            line = line.strip()
            if line and not line.startswith('#'):
                # Handle space-separated files (like GDB-17)
                smiles = line.split()[0] if line.split() else line
                fragments.append(smiles)
    
    return fragments


def clean_fragment_list(raw_fragments: List[str], require_carbon: bool = True) -> List[str]:
    """
    Clean a list of fragments, keeping only valid molecules.
    
    Args:
        raw_fragments: List of raw SMILES strings
        require_carbon: If True, only keep fragments with at least one carbon atom
    
    Returns:
        List of cleaned, canonical SMILES strings
    """
    cleaned = []
    seen = set()
    
    for smi in raw_fragments:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            
            atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
            if require_carbon and "C" not in atoms:
                continue
            
            Chem.SanitizeMol(mol)
            canonical = Chem.MolToSmiles(mol)
            
            # Deduplicate
            if canonical not in seen:
                seen.add(canonical)
                cleaned.append(canonical)
        except:
            continue
    
    return cleaned


def merge_fragments(frag1_smiles: str, frag2_smiles: str) -> Optional[str]:
    """
    Merge two fragments into a single molecule.
    
    Args:
        frag1_smiles: SMILES string of first fragment
        frag2_smiles: SMILES string of second fragment
    
    Returns:
        Merged SMILES string, or None if merge fails
    """
    try:
        mol1 = Chem.MolFromSmiles(frag1_smiles)
        mol2 = Chem.MolFromSmiles(frag2_smiles)
        
        if not mol1 or not mol2:
            return None
        
        # Combine molecules
        combo = Chem.CombineMols(mol1, mol2)
        if combo:
            # Try to sanitize and return canonical SMILES
            try:
                Chem.SanitizeMol(combo)
                return Chem.MolToSmiles(combo)
            except:
                return Chem.MolToSmiles(combo)
        return None
    except:
        return None


def generate_molecule_from_fragments(
    fragments: List[str],
    n_frags: int = 3,
    max_attempts: int = 100
) -> Optional[str]:
    """
    Generate a molecule by randomly combining fragments.
    
    Args:
        fragments: List of fragment SMILES strings
        n_frags: Number of fragments to combine
        max_attempts: Maximum attempts to generate valid molecule
    
    Returns:
        Generated SMILES string, or None if generation fails
    """
    if not fragments:
        return None
    
    for _ in range(max_attempts):
        selected = random.sample(fragments, k=min(n_frags, len(fragments)))
        smiles = "".join(selected)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                return Chem.MolToSmiles(mol)
            except:
                continue
    
    return None


def is_fuel_relevant_fragment(smiles: str) -> bool:
    """
    Check if fragment contains only C, H, O atoms (fuel-relevant).
    
    Args:
        smiles: SMILES string to check
    
    Returns:
        True if fragment contains only C, H, O and has at least one carbon atom
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
        allowed = {'C', 'H', 'O'}
        
        # Must contain only C, H, O and at least one carbon
        return atoms.issubset(allowed) and 'C' in atoms
    except:
        return False


def load_fragments_streaming(
    file_path: str,
    max_fragments: Optional[int] = None,
    filter_func: Optional[Callable[[str], bool]] = None
) -> Iterator[str]:
    """
    Stream fragments from file instead of loading all into memory.
    
    Supports .txt, .smi, .gz, and .tgz (tar.gz) files.
    For .tgz files, extracts and streams from the first SMILES file found.
    
    Args:
        file_path: Path to fragment file
        max_fragments: Maximum fragments to yield (None for all)
        filter_func: Optional function to filter fragments (returns True to keep)
    
    Yields:
        SMILES strings one at a time
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Fragment file not found: {file_path}")
    
    count = 0
    
    # Handle tar.gz archives
    if path.suffix == '.tgz' or path.name.endswith('.tgz'):
        with tarfile.open(path, 'r:gz') as tar:
            # Find ALL .smi or .txt files in archive (GDB-11 has multiple files)
            smi_files = []
            for member in tar.getmembers():
                if member.isfile() and (member.name.endswith('.smi') or member.name.endswith('.txt')):
                    smi_files.append(member)
            
            if not smi_files:
                raise ValueError(f"No SMILES file found in archive: {file_path}")
            
            # Stream from all SMILES files in the archive
            for smi_file in smi_files:
                if max_fragments and count >= max_fragments:
                    break
                
                f = tar.extractfile(smi_file)
                if f is None:
                    continue  # Skip if can't extract
                
                try:
                    for line in f:
                        if max_fragments and count >= max_fragments:
                            break
                        line = line.decode('utf-8').strip()
                        if line and not line.startswith('#'):
                            smiles = line.split()[0] if line.split() else line
                            if filter_func is None or filter_func(smiles):
                                yield smiles
                                count += 1
                finally:
                    f.close()
    
    # Handle gzipped files
    elif path.suffix == '.gz' or '.gz' in path.name:
        with gzip.open(path, 'rt') as f:
            for line in f:
                if max_fragments and count >= max_fragments:
                    break
                line = line.strip()
                if line and not line.startswith('#'):
                    smiles = line.split()[0] if line.split() else line
                    if filter_func is None or filter_func(smiles):
                        yield smiles
                        count += 1
    
    # Handle plain text files
    else:
        with open(path, 'r') as f:
            for line in f:
                if max_fragments and count >= max_fragments:
                    break
                line = line.strip()
                if line and not line.startswith('#'):
                    smiles = line.split()[0] if line.split() else line
                    if filter_func is None or filter_func(smiles):
                        yield smiles
                        count += 1


def preprocess_gdb11(
    input_file: str,
    output_file: str,
    max_fragments: Optional[int] = None,
    progress_interval: int = 100000
) -> int:
    """
    Pre-process GDB-11 file to extract only fuel-relevant fragments (C/H/O only).
    
    Handles tar.gz archives (.tgz), gzipped files (.gz), and plain text files.
    For .tgz files, extracts and processes the first SMILES file found.
    
    Args:
        input_file: Path to GDB-11 file (supports .tgz, .gz, .smi, .txt)
        output_file: Path to write filtered fragments (plain text, one SMILES per line)
        max_fragments: Maximum fragments to process (None for all)
        progress_interval: Show progress every N fragments
    
    Returns:
        Number of fragments written
    """
    count = 0
    path = Path(input_file)
    
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    with open(output_file, 'w') as outfile:
        # Handle tar.gz archives
        if path.suffix == '.tgz' or path.name.endswith('.tgz'):
            with tarfile.open(path, 'r:gz') as tar:
                # Find ALL .smi or .txt files in archive (GDB-11 has multiple files)
                smi_files = []
                for member in tar.getmembers():
                    if member.isfile() and (member.name.endswith('.smi') or member.name.endswith('.txt')):
                        smi_files.append(member)
                
                if not smi_files:
                    raise ValueError(f"No SMILES file found in archive: {input_file}")
                
                # Process all SMILES files in the archive
                for smi_file in smi_files:
                    if max_fragments and count >= max_fragments:
                        break
                    
                    f = tar.extractfile(smi_file)
                    if f is None:
                        continue  # Skip if can't extract
                    
                    try:
                        for line in f:
                            if max_fragments and count >= max_fragments:
                                break
                            line = line.decode('utf-8').strip()
                            if line and not line.startswith('#'):
                                smiles = line.split()[0] if line.split() else line
                                if is_fuel_relevant_fragment(smiles):
                                    outfile.write(smiles + '\n')
                                    count += 1
                                    if count % progress_interval == 0:
                                        print(f"Processed {count:,} fragments...")
                    finally:
                        f.close()
        
        # Handle gzipped files
        elif path.suffix == '.gz' or '.gz' in path.name:
            with gzip.open(path, 'rt') as infile:
                for line in infile:
                    if max_fragments and count >= max_fragments:
                        break
                    line = line.strip()
                    if line and not line.startswith('#'):
                        smiles = line.split()[0] if line.split() else line
                        if is_fuel_relevant_fragment(smiles):
                            outfile.write(smiles + '\n')
                            count += 1
                            if count % progress_interval == 0:
                                print(f"Processed {count:,} fragments...")
        
        # Handle plain text files
        else:
            with open(path, 'r') as infile:
                for line in infile:
                    if max_fragments and count >= max_fragments:
                        break
                    line = line.strip()
                    if line and not line.startswith('#'):
                        smiles = line.split()[0] if line.split() else line
                        if is_fuel_relevant_fragment(smiles):
                            outfile.write(smiles + '\n')
                            count += 1
                            if count % progress_interval == 0:
                                print(f"Processed {count:,} fragments...")
    
    return count


def sample_fragments_reservoir(
    fragment_file: str,
    n_samples: int,
    filter_func: Optional[Callable[[str], bool]] = None
) -> List[str]:
    """
    Sample fragments from large file using reservoir sampling.
    
    Useful for GDB-11 where loading all fragments is impractical.
    Ensures uniform random sampling from the entire file.
    
    Args:
        fragment_file: Path to fragment file (supports .tgz, .gz, .smi, .txt)
        n_samples: Number of fragments to sample
        filter_func: Optional function to filter fragments (returns True to keep)
    
    Returns:
        List of sampled fragment SMILES strings
    """
    samples = []
    count = 0
    
    for smiles in load_fragments_streaming(fragment_file, filter_func=filter_func):
        count += 1
        if len(samples) < n_samples:
            samples.append(smiles)
        else:
            # Reservoir sampling: replace with probability n_samples/count
            j = random.randint(0, count - 1)
            if j < n_samples:
                samples[j] = smiles
    
    return samples


def build_fragment_cache(
    fragment_file: str,
    cache_size: int = 100000,
    filter_func: Optional[Callable[[str], bool]] = None
) -> List[str]:
    """
    Build a fixed-size cache of fragments for random sampling.
    
    Uses reservoir sampling to ensure uniform distribution from large files.
    This is useful when working with GDB-11 or other large fragment databases.
    
    Args:
        fragment_file: Path to fragment file (supports .tgz, .gz, .smi, .txt)
        cache_size: Maximum number of fragments to cache
        filter_func: Optional function to filter fragments (returns True to keep)
    
    Returns:
        List of cached fragment SMILES strings
    """
    return sample_fragments_reservoir(fragment_file, n_samples=cache_size, filter_func=filter_func)

