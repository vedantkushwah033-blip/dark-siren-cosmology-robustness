"""Coverage audit across injected H0 values, host separations, and weight misspecification.

Synthetic low-redshift model only. Each observation is generated from one latent
host using 70/30 probabilities. Its posterior is evaluated with either correct
or reversed host weights. The noise scale is fixed within each replicate and is
shared by the generator and likelihood.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_VALUES = (52.0, 55.0, 60.0, 70.0, 80.0, 85.0, 88.0)
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10)
TRUE_WEIGHTS = np.array([0.7, 0.3])
INFERENCE_WEIGHTS = {
    "correct": np.array([0.7, 0.3]),
    "reversed": np.array([0.3, 0.7]),
}
GRID = np.linspace(25.0, 125.0, 2001)


def run_coverage_audit(replicates=300, seed=20261029):
    if replicates < 200:
        raise ValueError("Use at least 200 replicates per design cell.")
    rng = np.random.default_rng(seed)
    rows, summaries = [], []
    for h0_true in H0_VALUES:
        for z2 in Z2_VALUES:
            redshifts = np.array([Z1, z2])
            for noise_fraction in NOISE_FRACTIONS:
                sigma = noise_fraction * C_KM_S * 0.03 / h0_true
                cell = []
                for replicate in range(replicates):
                    host = int(rng.choice(2, p=TRUE_WEIGHTS))
                    d_true = C_KM_S * redshifts[host] / h0_true
                    d_obs = float(rng.normal(d_true, sigma))
                    while d_obs <= 0:
                        d_obs = float(rng.normal(d_true, sigma))
                    for weighting, weights in INFERENCE_WEIGHTS.items():
                        density = posterior_h0(
                            GRID, redshifts, d_obs, sigma, weights=weights,
                            prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                        )
                        post = posterior_summary(GRID, density)
                        row = {
                            "H0_true": h0_true,
                            "z1": Z1,
                            "z2": z2,
                            "redshift_separation": z2 - Z1,
                            "noise_fraction": noise_fraction,
                            "replicate": replicate,
                            "true_host_index": host,
                            "observed_distance_mpc": d_obs,
                            "sigma_mpc": sigma,
                            "inference_weighting": weighting,
                            "posterior_median": post["median"],
                            "median_error": post["median"] - h0_true,
                            "covered_68": int(post["lower_68"] <= h0_true <= post["upper_68"]),
                            "covered_95": int(post["lower_95"] <= h0_true <= post["upper_95"]),
                            "width_68": post["upper_68"] - post["lower_68"],
                            "width_95": post["upper_95"] - post["lower_95"],
                            "lower_68": post["lower_68"],
                            "upper_68": post["upper_68"],
                            "lower_95": post["lower_95"],
                            "upper_95": post["upper_95"],
                        }
                        rows.append(row)
                        cell.append(row)
                for weighting in INFERENCE_WEIGHTS:
                    group = [r for r in cell if r["inference_weighting"] == weighting]
                    cov68 = float(np.mean([r["covered_68"] for r in group]))
                    cov95 = float(np.mean([r["covered_95"] for r in group]))
                    errors = np.array([r["median_error"] for r in group])
                    summaries.append({
                        "H0_true": h0_true,
                        "z1": Z1,
                        "z2": z2,
                        "redshift_separation": z2 - Z1,
                        "noise_fraction": noise_fraction,
                        "inference_weighting": weighting,
                        "replicates": replicates,
                        "mean_median_error": float(errors.mean()),
                        "median_rmse": float(np.sqrt(np.mean(errors**2))),
                        "coverage_68": cov68,
                        "coverage_68_mcse": float(np.sqrt(cov68 * (1 - cov68) / replicates)),
                        "coverage_95": cov95,
                        "coverage_95_mcse": float(np.sqrt(cov95 * (1 - cov95) / replicates)),
                        "mean_width_68": float(np.mean([r["width_68"] for r in group])),
                        "mean_width_95": float(np.mean([r["width_95"] for r in group])),
                    })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--output-dir", default="results/multi-h0-coverage-audit")
    args = parser.parse_args()
    rows, summaries = run_coverage_audit(args.replicates, args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Completed {len(rows)} posterior fits across {len(summaries)} design cells.")
    print("Selected aggregate summaries by H0 and inference weighting:")
    for h0 in H0_VALUES:
        for weighting in INFERENCE_WEIGHTS:
            group = [r for r in summaries if r["H0_true"] == h0 and r["inference_weighting"] == weighting]
            print({
                "H0_true": h0,
                "weighting": weighting,
                "mean_68_coverage_across_designs": float(np.mean([r["coverage_68"] for r in group])),
                "mean_95_coverage_across_designs": float(np.mean([r["coverage_95"] for r in group])),
                "mean_rmse_across_designs": float(np.mean([r["median_rmse"] for r in group])),
            })


if __name__ == "__main__":
    main()
