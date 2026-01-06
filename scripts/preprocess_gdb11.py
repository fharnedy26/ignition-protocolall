#!/usr/bin/env python3
"""
Preprocess GDB-11 database to extract fuel-relevant fragments.

Filters GDB-11 database (gdb11.tgz) to extract only molecules containing
C, H, O atoms (fuel-relevant). Supports tar.gz archives, gzipped files,
and plain text files.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import molecular_generator
sys.path.insert(0, str(Path(__file__).parent.parent))

from molecular_generator.fragments import preprocess_gdb11


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Preprocess GDB-11 database to extract fuel-relevant fragments (C/H/O only)"
    )
    
    parser.add_argument(
        '--input',
        required=True,
        help='Path to GDB-11 file (supports .tgz, .gz, .smi, .txt)'
    )
    
    parser.add_argument(
        '--output',
        required=True,
        help='Path to write filtered fragments (plain text, one SMILES per line)'
    )
    
    parser.add_argument(
        '--max-fragments',
        type=int,
        default=None,
        help='Maximum fragments to process (None for all)'
    )
    
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=100000,
        help='Progress update frequency (default: 100000)'
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {args.input}")
        return 1
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"[INFO] Preprocessing GDB-11 file: {args.input}")
    print(f"[INFO] Output file: {args.output}")
    if args.max_fragments:
        print(f"[INFO] Maximum fragments: {args.max_fragments:,}")
    print(f"[INFO] Progress interval: {args.progress_interval:,}")
    print()
    
    try:
        count = preprocess_gdb11(
            input_file=args.input,
            output_file=args.output,
            max_fragments=args.max_fragments,
            progress_interval=args.progress_interval
        )
        
        print()
        print(f"[SUCCESS] Preprocessing complete!")
        print(f"[RESULT] Wrote {count:,} fuel-relevant fragments to {args.output}")
        return 0
    
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        return 1
    except ValueError as e:
        print(f"[ERROR] {e}")
        return 1
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())


