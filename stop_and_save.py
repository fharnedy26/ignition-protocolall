#!/usr/bin/env python3
"""
Stop current evolution run and save results.

If evolution has plateaued, you can stop it and save what you have.
"""

import signal
import sys
from pathlib import Path

def signal_handler(sig, frame):
    print("\n\n[STOPPED] Evolution interrupted by user")
    print("Results will be saved when the script completes naturally.")
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

print("=" * 70)
print("Early Stopping Helper")
print("=" * 70)
print("\nIf your evolution has plateaued, you can:")
print("1. Press Ctrl+C in the terminal running the evolution")
print("2. The script will save results when it completes")
print("\nOr wait for the script to finish - it will save all results.")
print("=" * 70)

