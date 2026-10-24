"""Full-population oracle benchmark for synthetic missing-host experiments.

This simulation-only benchmark uses the complete mock galaxy population and its
true host-generation weights in the posterior. It is an upper-information
reference, not a practical estimator for an incomplete real survey.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from statistics import median


def median_bias_from_rows(rows, true_h0=70.0):
    """Median across runs of posterior-median errors, consistently defined."""
    if not rows:
        raise ValueError("rows must not be empty")
    return float(median(float(row["median"]) - true_h0 for row in rows))


SCENARIOS = (("random70_full_population_oracle", "random"),
             ("redshift70_full_population_oracle", "redshift"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/full-population-oracle")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    import numpy as np

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise SystemExit(f"Could not find original simulator under {original}")
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0, posterior_summary
    from src.experiments import _mask, _weights
    from src.metrics_v2 import summarize_runs

    grid = np.linspace(50.0, 90.0, 1601)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summaries, per_run = [], []

    # Use the same mock universe and GW noise for both mask mechanisms at
    # each replicate index. The oracle posterior uses the full population, so
    # posterior summaries must be identical across masks; only mask diagnostics
    # (realized completeness and whether the host is observed) may differ.
    for name, missingness in SCENARIOS:
        scenario_seed = args.seed
        rows = []
        for r in range(args.replicates):
            s = scenario_seed + r
            cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=70.0, seed=s)
            true_weights = _weights(cat, "luminosity")
            host = choose_host(cat, weighting="luminosity", seed=s + 100000)
            gw = simulate_gw_event(cat["true_distance_mpc"][host], seed=s + 200000)
            # Generate the incomplete observed mask for an explicit paired
            # selection context, but oracle inference deliberately uses every
            # galaxy and its true host-generation weights.
            mask = _mask(cat, 0.7, missingness, s + 300000)
            posterior = posterior_h0(
                grid, cat["redshift"], gw["observed_distance_mpc"],
                gw["sigma_mpc"], weights=true_weights,
            )
            summary = posterior_summary(grid, posterior)
            rows.append({
                "run": r, "seed": s, "H0_true": 70.0, "scenario": name,
                "host_redshift": float(cat["redshift"][host]),
                "observed_distance_mpc": float(gw["observed_distance_mpc"]),
                "distance_sigma_mpc": float(gw["sigma_mpc"]),
                "missingness": missingness, "target_completeness": 0.7,
                "realized_completeness": float(mask.mean()),
                "host_observed": bool(mask[host]),
                "median": float(summary["median"]), "mean": float(summary["mean"]),
                "posterior_sd": float(summary["posterior_sd"]),
                "lower_68": float(summary["lower_68"]),
                "upper_68": float(summary["upper_68"]),
                "lower_95": float(summary["lower_95"]),
                "upper_95": float(summary["upper_95"]),
                "width_68": float(summary["upper_68"] - summary["lower_68"]),
                "width_95": float(summary["upper_95"] - summary["lower_95"]),
                "pull_median": float((summary["median"] - 70.0) / max(summary["posterior_sd"], 1e-12)),
                "pull_mean": float((summary["mean"] - 70.0) / max(summary["posterior_sd"], 1e-12)),
            })
        metrics = summarize_runs(rows, true_h0=70.0)
        # Consistent definition: median across runs of posterior-median errors.
        metrics["median_bias"] = median_bias_from_rows(rows, true_h0=70.0)
        summaries.append({
            "scenario": name, "missingness": missingness,
            "target_completeness": 0.7, "host_generation": "luminosity",
            "inference_weighting": "true luminosity weights over full population",
            "uses_hidden_population": True, "requested_replicates": args.replicates,
            "seed": scenario_seed, **metrics,
        })
        per_run.extend(rows)

    for filename, rows in (("scenario_summary.csv", summaries),
                           ("per_run_results.csv", per_run)):
        with (out / filename).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    metadata = {
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0": 70.0, "n_galaxies": 500,
        "requested_replicates_per_scenario": args.replicates,
        "base_seed": args.seed, "scenario_count": len(SCENARIOS),
        "interpretation_warning": (
            "Simulation-only upper-information benchmark: inference uses all "
            "galaxies and true luminosity host weights, including galaxies hidden "
            "by the simulated catalogue mask. Not implementable from an incomplete real catalogue."
        ),
        "scenarios": summaries,
    }
    (out / "scenario_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
