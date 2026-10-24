"""Targeted factorial sensitivity study using the pinned original simulator.

Separates catalogue completeness/missingness from host-generation versus
inference-weighting mismatch. This is synthetic-model analysis, not LVK data.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

SCENARIOS = (
    ("complete_uniform_matched", 1.0, "complete", "uniform", "uniform"),
    ("complete_luminosity_matched", 1.0, "complete", "luminosity", "luminosity"),
    ("complete_luminosity_host_uniform_inference", 1.0, "complete", "luminosity", "uniform"),
    ("random70_uniform_matched", 0.7, "random", "uniform", "uniform"),
    ("random70_luminosity_matched", 0.7, "random", "luminosity", "luminosity"),
    ("random70_luminosity_host_uniform_inference", 0.7, "random", "luminosity", "uniform"),
    ("redshift70_uniform_matched", 0.7, "redshift", "uniform", "uniform"),
    ("redshift70_luminosity_matched", 0.7, "redshift", "luminosity", "luminosity"),
    ("redshift70_luminosity_host_uniform_inference", 0.7, "redshift", "luminosity", "uniform"),
)



def median_bias(posterior_medians: list[float], true_h0: float) -> float:
    """Return the median replicate-level posterior-median error."""
    if not posterior_medians:
        raise ValueError("posterior_medians must not be empty")
    true_value = float(true_h0)
    values = [float(value) for value in posterior_medians]
    if not math.isfinite(true_value):
        raise ValueError("true_h0 must be finite")
    if not all(math.isfinite(value) for value in values):
        raise ValueError("posterior_medians must contain only finite values")
    errors = sorted(value - true_value for value in values)
    n_errors = len(errors)
    if n_errors % 2:
        return errors[n_errors // 2]
    return (errors[n_errors // 2 - 1] + errors[n_errors // 2]) / 2.0
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/factorial-sensitivity")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise FileNotFoundError(f"Original source tree not found under {original}")
    sys.path.insert(0, str(original))
    from src.experiments import run_repeated
    from src.metrics_v2 import summarize_runs

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    per_run = []
    for i, (name, completeness, missingness, host_generation, inference_weighting) in enumerate(SCENARIOS):
        rows = run_repeated(
            n_runs=args.replicates,
            n_galaxies=500,
            H0_true=70.0,
            completeness=completeness,
            missingness=missingness,
            host_generation=host_generation,
            inference_weighting=inference_weighting,
            seed=args.seed + i * 1_000_000,
        )
        metrics = summarize_runs(rows, true_h0=70.0)
        # The pinned legacy summarizer labels the mean replicate error as
        # median_bias. Recompute the median of replicate-level posterior-median
        # errors explicitly so the scenario summary matches the metric name.
        metrics["median_bias"] = median_bias([float(row["median"]) for row in rows], 70.0)
        summaries.append({
            "scenario": name,
            "target_completeness": completeness,
            "missingness": missingness,
            "host_generation": host_generation,
            "inference_weighting": inference_weighting,
            "requested_replicates": args.replicates,
            "seed": args.seed + i * 1_000_000,
            **metrics,
        })
        for row in rows:
            # Preserve the full prespecified design in every raw record.
            # The independent auditor must not infer settings from names.
            per_run.append({
                **row,
                # Keep the design labels authoritative even if the legacy
                # simulator starts emitting fields with the same names.
                "scenario": name,
                "missingness": missingness,
                "host_generation": host_generation,
                "inference_weighting": inference_weighting,
            })

    with (out / "scenario_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    with (out / "per_run_results.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(per_run[0]))
        writer.writeheader()
        writer.writerows(per_run)
    metadata = {
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0": 70.0,
        "n_galaxies": 500,
        "requested_replicates_per_scenario": args.replicates,
        "base_seed": args.seed,
        "scenario_count": len(SCENARIOS),
        "interpretation_warning": "Synthetic model experiment; not observational inference. Scenarios use distinct seed blocks.",
        "scenarios": summaries,
    }
    (out / "scenario_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
