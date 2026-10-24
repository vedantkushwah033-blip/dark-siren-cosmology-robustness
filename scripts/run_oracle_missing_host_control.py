"""Oracle-control experiment for missing-host contributions.

Compares observed-catalogue-only inference with two missing-host mixture terms.
The oracle variant uses hidden simulated population information and is a
simulation validation control, not an implementable real-survey estimator.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from statistics import median

SCENARIOS = (
    ("random70_observed_only", "random", False, "inference"),
    ("random70_missing_term_inference_weights", "random", True, "inference"),
    ("random70_missing_term_oracle_weights", "random", True, "oracle_true"),
    ("redshift70_observed_only", "redshift", False, "inference"),
    ("redshift70_missing_term_inference_weights", "redshift", True, "inference"),
    ("redshift70_missing_term_oracle_weights", "redshift", True, "oracle_true"),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--output-dir", default="results/oracle-missing-host")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise SystemExit(f"Could not find original simulator under {original}")
    sys.path.insert(0, str(original))
    from src.experiments import run_repeated
    from src.metrics_v2 import summarize_runs

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summaries, per_run = [], []

    for i, (name, missingness, include_term, missing_weighting) in enumerate(SCENARIOS):
        rows = run_repeated(
            n_runs=args.replicates,
            n_galaxies=500,
            H0_true=70.0,
            completeness=0.7,
            missingness=missingness,
            host_generation="luminosity",
            inference_weighting="uniform",
            include_missing_host_term=include_term,
            missing_host_weighting=missing_weighting,
            seed=args.seed + i * 1_000_000,
        )
        metrics = summarize_runs(rows, true_h0=70.0)
        # Match the main factorial's explicit definition: median across runs
        # of each run's posterior-median error, not the mean of those errors.
        metrics["median_bias"] = float(median(float(row["median"]) - 70.0 for row in rows))
        summaries.append({
            "scenario": name,
            "missingness": missingness,
            "target_completeness": 0.7,
            "host_generation": "luminosity",
            "inference_weighting": "uniform",
            "include_missing_host_term": include_term,
            "missing_host_weighting": missing_weighting,
            "requested_replicates": args.replicates,
            "seed": args.seed + i * 1_000_000,
            **metrics,
        })
        for row in rows:
            per_run.append({"scenario": name, **row})

    for filename, rows in (
        ("scenario_summary.csv", summaries),
        ("per_run_results.csv", per_run),
    ):
        with (out / filename).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    metadata = {
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0": 70.0,
        "n_galaxies": 500,
        "requested_replicates_per_scenario": args.replicates,
        "base_seed": args.seed,
        "scenario_count": len(SCENARIOS),
        "warning": "Synthetic validation only. Oracle hidden-host weighting uses information unavailable to a real survey.",
        "scenarios": summaries,
    }
    (out / "scenario_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
