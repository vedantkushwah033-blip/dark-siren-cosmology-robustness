"""Extend fixed-truth coverage diagnostics across the prior-boundary transition.

The simulation generator and inference implementation are reused unchanged.
This is a synthetic diagnostic, not an observational H0 analysis.
"""
from __future__ import annotations
import argparse, csv, importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("prior_support_source", HERE / "run_prior_support_sensitivity.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load prior-support experiment")
prior = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prior)

TRUTH_GRID = (52.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 88.0)
PRIOR_GRID = ((50.0, 90.0), (40.0, 100.0))

def write_csv(path, rows):
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        check = list(reader)
        if reader.fieldnames != fields or len(check) != len(rows):
            raise RuntimeError(f"CSV round-trip validation failed: {path}")
    return fields

def run(replicates=300, seed=20261025, grid_size=4001):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates per truth")
    old_truths, old_priors = prior.TRUTHS, prior.PRIORS
    try:
        prior.TRUTHS = TRUTH_GRID
        prior.PRIORS = PRIOR_GRID
        rows, summaries = prior.run_audit(replicates=replicates, seed=seed, grid_size=grid_size)
    finally:
        prior.TRUTHS, prior.PRIORS = old_truths, old_priors
    expected_rows = len(TRUTH_GRID) * len(PRIOR_GRID) * replicates
    expected_summaries = len(TRUTH_GRID) * len(PRIOR_GRID) * len(prior.METHODS)
    if len(rows) != expected_rows or len(summaries) != expected_summaries:
        raise RuntimeError("Coverage-transition experiment is incomplete")
    return rows, summaries

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/coverage-transition")
    args = parser.parse_args()
    rows, summaries = run(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir)
    write_csv(out / "per_run.csv", rows)
    write_csv(out / "summary.csv", summaries)
    print(f"PASS: {len(rows)} observation/prior rows; {len(summaries)} summaries")
    for row in summaries:
        print(row)

if __name__ == "__main__":
    main()
