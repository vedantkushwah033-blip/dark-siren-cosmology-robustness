"""Paired synthetic flat-Lambda-CDM versus low-redshift model stress test.

This is a controlled Gaussian-distance toy experiment, not a detector-level or
observational dark-siren analysis.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from dark_siren_robustness.cosmological_inference import (
    distance_grid_for_model,
    posterior_h0_from_distance_grid,
)
from dark_siren_robustness.hardened_inference import posterior_summary

TRUTHS = (60.0, 70.0, 80.0)
REDSHIFTS = (0.02, 0.05, 0.10, 0.20, 0.40, 0.60)
OMEGA_M = 0.3
PRIOR = (40.0, 100.0)


def run(output_dir: Path, replicates: int = 200, seed: int = 20261024, grid_size: int = 2001) -> dict:
    if replicates < 50:
        raise ValueError("Use at least 50 replicates per truth/redshift cell.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    h0_grid = np.linspace(25.0, 125.0, grid_size)
    raw: list[dict] = []

    for z in REDSHIFTS:
        lcdm_distances = distance_grid_for_model(
            h0_grid, [z], "flat_lcdm", omega_m=OMEGA_M, integration_steps=4000
        )
        lowz_distances = distance_grid_for_model(h0_grid, [z], "low_z")
        for truth in TRUTHS:
            true_distance = float(
                distance_grid_for_model([truth], [z], "flat_lcdm", omega_m=OMEGA_M, integration_steps=8000)[0, 0]
            )
            sigma = max(20.0, 0.05 * true_distance)
            for replicate in range(replicates):
                observed = float(rng.normal(true_distance, sigma))
                if observed <= 0:
                    # Keep the Gaussian observation model explicit while ensuring
                    # the likelihood's positive-distance input contract.
                    observed = float(abs(observed) + np.finfo(float).eps)
                for model, distances in (("matched_flat_lcdm", lcdm_distances), ("mismatched_low_z", lowz_distances)):
                    posterior = posterior_h0_from_distance_grid(
                        h0_grid,
                        distances,
                        observed,
                        sigma,
                        host_weights=[1.0],
                        prior_min=PRIOR[0],
                        prior_max=PRIOR[1],
                    )
                    summary = posterior_summary(h0_grid, posterior)
                    median = float(summary["median"])
                    raw.append({
                        "true_H0": truth,
                        "redshift": z,
                        "replicate": replicate,
                        "observed_distance_mpc": observed,
                        "distance_sigma_mpc": sigma,
                        "inference_model": model,
                        "posterior_median_H0": median,
                        "error_H0": median - truth,
                        "absolute_error_H0": abs(median - truth),
                        "covered_68": int(summary["lower_68"] <= truth <= summary["upper_68"]),
                        "covered_95": int(summary["lower_95"] <= truth <= summary["upper_95"]),
                        "width_68": float(summary["upper_68"] - summary["lower_68"]),
                        "width_95": float(summary["upper_95"] - summary["lower_95"]),
                    })

    raw_path = output_dir / "per_run.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(raw[0]))
        writer.writeheader()
        writer.writerows(raw)

    summaries: list[dict] = []
    for z in REDSHIFTS:
        for truth in TRUTHS:
            for model in ("matched_flat_lcdm", "mismatched_low_z"):
                rows = [r for r in raw if r["redshift"] == z and r["true_H0"] == truth and r["inference_model"] == model]
                errors = np.array([r["error_H0"] for r in rows], dtype=float)
                summaries.append({
                    "redshift": z,
                    "true_H0": truth,
                    "inference_model": model,
                    "n": len(rows),
                    "bias_H0": float(errors.mean()),
                    "mae_H0": float(np.abs(errors).mean()),
                    "rmse_H0": float(np.sqrt(np.mean(errors ** 2))),
                    "coverage_68": float(np.mean([r["covered_68"] for r in rows])),
                    "coverage_95": float(np.mean([r["covered_95"] for r in rows])),
                    "mean_width_68": float(np.mean([r["width_68"] for r in rows])),
                    "mean_width_95": float(np.mean([r["width_95"] for r in rows])),
                })
    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    manifest = {
        "experiment": "paired flat-Lambda-CDM generation with matched and low-z inference",
        "scope": "synthetic Gaussian-distance diagnostic; not observational dark-siren inference",
        "seed": seed,
        "replicates_per_truth_redshift_model": replicates,
        "truths_H0": list(TRUTHS),
        "redshifts": list(REDSHIFTS),
        "omega_m_fixed": OMEGA_M,
        "prior_H0": list(PRIOR),
        "grid_size": grid_size,
        "grid_range_H0": [25.0, 125.0],
        "distance_sigma_rule": "max(20 Mpc, 0.05 times true flat-Lambda-CDM luminosity distance)",
        "raw_rows": len(raw),
        "summary_rows": len(summaries),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    (output_dir / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/cosmological-model-robustness"))
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--grid-size", type=int, default=2001)
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir, args.replicates, args.seed, args.grid_size), indent=2))


if __name__ == "__main__":
    main()
