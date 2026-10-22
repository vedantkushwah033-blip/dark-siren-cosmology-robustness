"""Reproduce both original baseline calibration controls from a pinned source tree."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/original-audit")
    args = parser.parse_args()

    if args.replicates < 2:
        parser.error("--replicates must be at least 2 to calculate sample pull SD")
    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").exists():
        raise FileNotFoundError(f"Cannot find original src/experiments.py under {original}")
    sys.path.insert(0, str(original))

    from src.experiments import run_repeated
    from src.metrics_v2 import summarize_runs
    from run_baseline_diagnostics import run_known_host

    # Unknown-host, complete-catalogue baseline, matching run_calibration.py.
    mixture_rows = run_repeated(
        n_runs=args.replicates,
        completeness=1.0,
        missingness="complete",
        seed=args.seed,
    )
    mixture_summary = summarize_runs(mixture_rows, true_h0=70.0)

    # Known-host control, matching run_baseline_diagnostics.py.
    known_rows = run_known_host(n_runs=args.replicates, H0_true=70.0, seed=args.seed)
    known_summary = summarize_runs(known_rows, true_h0=70.0)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(output_dir / "per_run_unknown_host.csv", mixture_rows)
    write_rows(output_dir / "per_run_known_host.csv", known_rows)

    metadata = {
        "requested_replicates": args.replicates,
        "seed": args.seed,
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "unknown_host_scenario": "complete catalogue; uniform host generation and inference weighting",
        "known_host_scenario": "true host supplied as sole candidate",
        "unknown_host": mixture_summary,
        "known_host": known_summary,
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
