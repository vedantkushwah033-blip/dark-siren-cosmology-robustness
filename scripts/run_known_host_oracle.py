"""Known-true-host oracle calibration control for the synthetic dark-siren model.

This intentionally gives inference the actual simulated host redshift. It is an
upper-information calibration reference, not a usable dark-siren estimator.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from statistics import median

import numpy as np

def median_bias_from_rows(rows, true_h0=70.0):
    """Median across runs of posterior-median errors."""
    if not rows:
        raise ValueError("rows must not be empty")
    return float(median(float(row["median"]) - true_h0 for row in rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/known-host-oracle")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise SystemExit(f"Could not find original simulator under {original}")
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0, posterior_summary
    from src.metrics_v2 import summarize_runs

    grid = np.linspace(50.0, 90.0, 1601)
    rows = []
    for r in range(args.replicates):
        s = args.seed + r
        cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=70.0, seed=s)
        host = choose_host(cat, weighting="luminosity", seed=s + 100000)
        gw = simulate_gw_event(cat["true_distance_mpc"][host], seed=s + 200000)
        # Oracle: use only the true host's known redshift in the likelihood.
        posterior = posterior_h0(
            grid,
            np.asarray([cat["redshift"][host]], dtype=float),
            gw["observed_distance_mpc"],
            gw["sigma_mpc"],
            weights=np.ones(1),
        )
        summary = posterior_summary(grid, posterior)
        rows.append({
            "run": r,
            "scenario": "known_true_host_oracle",
            "seed": s,
            "H0_true": 70.0,
            "host_generation": "luminosity",
            "host_redshift": float(cat["redshift"][host]),
            "observed_distance_mpc": float(gw["observed_distance_mpc"]),
            "distance_sigma_mpc": float(gw["sigma_mpc"]),
            "median": float(summary["median"]),
            "mean": float(summary["mean"]),
            "posterior_sd": float(summary["posterior_sd"]),
            "lower_68": float(summary["lower_68"]),
            "upper_68": float(summary["upper_68"]),
            "lower_95": float(summary["lower_95"]),
            "upper_95": float(summary["upper_95"]),
            "width_68": float(summary["upper_68"] - summary["lower_68"]),
            "width_95": float(summary["upper_95"] - summary["lower_95"]),
            "pull_median": float((summary["median"] - 70.0) / max(summary["posterior_sd"], 1e-12)),
            "pull_mean": float((summary["mean"] - 70.0) / max(summary["posterior_sd"], 1e-12)),
            "host_observed": True,
            "realized_completeness": 1.0,
        })
    metrics = summarize_runs(rows, true_h0=70.0)
    metrics["median_bias"] = median_bias_from_rows(rows, true_h0=70.0)
    summary_row = {
        "scenario": "known_true_host_oracle",
        "host_generation": "luminosity",
        "inference_weighting": "known true host redshift",
        "uses_hidden_population": True,
        "requested_replicates": args.replicates,
        "seed": args.seed,
        **metrics,
    }
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "per_run_results.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (out / "scenario_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_row))
        writer.writeheader()
        writer.writerow(summary_row)
    metadata = {
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0": 70.0,
        "n_galaxies": 500,
        "requested_replicates_per_scenario": args.replicates,
        "base_seed": args.seed,
        "scenario_count": 1,
        "interpretation_warning": (
            "Simulation-only calibration reference. The true host redshift is supplied to inference; "
            "this information is not available for a real dark siren with an unknown host."
        ),
        "scenarios": [summary_row],
    }
    (out / "scenario_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
