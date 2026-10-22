"""Reproduce the original repository's complete-catalogue coverage baseline.

Run from the robustness repository after the original project has been cloned
into --original-dir. The original code is imported and called directly; this
script does not reimplement its likelihood.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/original-audit")
    args = parser.parse_args()

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").exists():
        raise FileNotFoundError(f"Cannot find original src/experiments.py under {original}")
    sys.path.insert(0, str(original))

    from src.experiments import run_repeated
    from src.metrics_v2 import summarize_runs

    rows = run_repeated(
        n_runs=args.replicates,
        completeness=1.0,
        missingness="complete",
        host_generation="uniform",
        inference_weighting="uniform",
        seed=args.seed,
    )
    summary = summarize_runs(rows, true_h0=70.0)
    # Save actual number of completed rows, as the original simulation can
    # skip a replicate if the observed catalogue is empty.
    summary["requested_replicates"] = args.replicates
    summary["seed"] = args.seed
    summary["scenario"] = "complete catalogue; uniform host generation and inference weighting"
    summary["original_repository"] = "vedantkushwah033-blip/bayesian-dark-siren-cosmology"
    summary["original_commit"] = "a726eae6018c021b79dd2a7fb245a2dd63d6d1c4"

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "per_run.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
