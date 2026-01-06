#!/usr/bin/env python3
"""
Utility to prune old run outputs and optionally move CSV data files into a data/ folder.

Usage examples:
  python scripts/clean_outputs.py --runs "runs" --older-than-days 14
  python scripts/clean_outputs.py --move-csvs 1 --root "."
"""
import argparse
import shutil
from pathlib import Path
import time


def prune_runs(runs_dir: Path, older_than_days: int) -> int:
    if not runs_dir.exists():
        return 0
    cutoff = time.time() - older_than_days * 86400.0
    removed = 0
    for child in runs_dir.iterdir():
        if child.is_dir():
            try:
                if child.stat().st_mtime < cutoff:
                    shutil.rmtree(child, ignore_errors=True)
                    removed += 1
            except Exception:
                pass
    return removed


def move_csvs_to_data(root: Path) -> int:
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    for csv_name in ("components.csv", "components_novel.csv", "baseline_blend.csv"):
        src = root / csv_name
        if src.exists() and src.is_file():
            dst = data_dir / csv_name
            try:
                shutil.move(str(src), str(dst))
                moved += 1
            except Exception:
                pass
    return moved


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", default="runs")
    p.add_argument("--older-than-days", type=int, default=14)
    p.add_argument("--move-csvs", type=int, default=0, choices=[0,1])
    p.add_argument("--root", default=".")
    args = p.parse_args()

    removed = prune_runs(Path(args.runs), int(args.older_than_days))
    if removed:
        print(f"Removed {removed} old run folder(s)")

    if int(args.move_csvs) == 1:
        moved = move_csvs_to_data(Path(args.root))
        if moved:
            print(f"Moved {moved} CSV file(s) to data/")


if __name__ == "__main__":
    main()

