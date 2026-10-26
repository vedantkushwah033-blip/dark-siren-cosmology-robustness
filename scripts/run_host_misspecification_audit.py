"""Controlled stress tests for misspecified host weights and omitted hosts.

All observations are synthetic. Each replicate has a fixed injected H0; host
identity is latent and distance noise is generated with the same fixed sigma
used by the likelihood. This is a misspecification diagnostic, not real-data
validation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_TRUE = 70.0
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
OBSERVED_REDSHIFTS = np.array([0.02, 0.04])
HIDDEN_REDSHIFT = 0.03
SIGMA_MPC = 0.10 * C_KM_S * 0.03 / 70.0
GRID = np.linspace(25.0, 125.0, 4001)
WILSON_Z_95 = 1.959963984540054


def wilson_interval(successes: int, trials: int, z: float = WILSON_Z_95):
    """Return a two-sided Wilson score interval for a binomial proportion."""
    if trials <= 0 or successes < 0 or successes > trials:
        raise ValueError("Require 0 <= successes <= trials and trials > 0.")
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half_width = z * np.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials**2)) / denominator
    return max(0.0, float(center - half_width)), min(1.0, float(center + half_width))


def generate_dataset(replicates: int, seed: int):
    if replicates < 200:
        raise ValueError("Use at least 200 replicates per scenario.")
    rng = np.random.default_rng(seed)
    specs = {
        "correct_weights": (OBSERVED_REDSHIFTS, np.array([0.7, 0.3])),
        "reversed_weights": (OBSERVED_REDSHIFTS, np.array([0.7, 0.3])),
        "hidden_host_ignored": (np.array([0.02, 0.04, HIDDEN_REDSHIFT]), np.array([0.35, 0.35, 0.30])),
        "hidden_host_modelled": (np.array([0.02, 0.04, HIDDEN_REDSHIFT]), np.array([0.35, 0.35, 0.30])),
    }
    data = []
    for scenario, (z_support, true_weights) in specs.items():
        for i in range(replicates):
            host_idx = int(rng.choice(len(z_support), p=true_weights))
            z_true = float(z_support[host_idx])
            d_true = C_KM_S * z_true / H0_TRUE
            d_obs = float(rng.normal(d_true, SIGMA_MPC))
            while d_obs <= 0:
                d_obs = float(rng.normal(d_true, SIGMA_MPC))
            data.append({
                "scenario": scenario,
                "replicate": i,
                "true_H0": H0_TRUE,
                "true_host_index": host_idx,
                "true_host_redshift": z_true,
                "observed_distance_mpc": d_obs,
                "sigma_mpc": SIGMA_MPC,
            })
    return data


def infer_row(row):
    scenario = row["scenario"]
    d, sigma = row["observed_distance_mpc"], row["sigma_mpc"]
    if scenario == "correct_weights":
        density = posterior_h0(GRID, OBSERVED_REDSHIFTS, d, sigma, weights=[0.7, 0.3],
                               prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
    elif scenario == "reversed_weights":
        density = posterior_h0(GRID, OBSERVED_REDSHIFTS, d, sigma, weights=[0.3, 0.7],
                               prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
    elif scenario == "hidden_host_ignored":
        # Deliberately omit z=0.03; renormalize the two observed hosts.
        density = posterior_h0(GRID, OBSERVED_REDSHIFTS, d, sigma, weights=[0.5, 0.5],
                               prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
    elif scenario == "hidden_host_modelled":
        # The model assigns 30% prior weight to the hidden-host component.
        density = posterior_h0(GRID, OBSERVED_REDSHIFTS, d, sigma, weights=[0.5, 0.5],
                               missing_weight=0.30, missing_redshifts=[HIDDEN_REDSHIFT],
                               prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    return posterior_summary(GRID, density)


def run_stress_test(replicates=1000, seed=20261027):
    data = generate_dataset(replicates, seed)
    rows, summaries = [], []
    for row in data:
        post = infer_row(row)
        record = {
            **row,
            "posterior_median": post["median"],
            "median_error": post["median"] - row["true_H0"],
            "lower_68": post["lower_68"],
            "upper_68": post["upper_68"],
            "lower_95": post["lower_95"],
            "upper_95": post["upper_95"],
            "width_68": post["upper_68"] - post["lower_68"],
            "width_95": post["upper_95"] - post["lower_95"],
            "covered_68": int(post["lower_68"] <= row["true_H0"] <= post["upper_68"]),
            "covered_95": int(post["lower_95"] <= row["true_H0"] <= post["upper_95"]),
        }
        rows.append(record)
    for scenario in ("correct_weights", "reversed_weights", "hidden_host_ignored", "hidden_host_modelled"):
        subset = [r for r in rows if r["scenario"] == scenario]
        errors = np.array([r["median_error"] for r in subset])
        covered_68 = int(sum(r["covered_68"] for r in subset))
        covered_95 = int(sum(r["covered_95"] for r in subset))
        coverage_68_interval = wilson_interval(covered_68, len(subset))
        coverage_95_interval = wilson_interval(covered_95, len(subset))
        summaries.append({
            "scenario": scenario,
            "replicates": len(subset),
            "true_H0": H0_TRUE,
            "mean_median_error": float(errors.mean()),
            "median_median_error": float(np.median(errors)),
            "rmse_median": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": covered_68 / len(subset),
            "coverage_68_wilson_95_low": coverage_68_interval[0],
            "coverage_68_wilson_95_high": coverage_68_interval[1],
            "coverage_95": covered_95 / len(subset),
            "coverage_95_wilson_95_low": coverage_95_interval[0],
            "coverage_95_wilson_95_high": coverage_95_interval[1],
            "mean_width_68": float(np.mean([r["width_68"] for r in subset])),
            "mean_width_95": float(np.mean([r["width_95"] for r in subset])),
        })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--output-dir", default="results/host-misspecification-audit")
    args = parser.parse_args()
    rows, summaries = run_stress_test(args.replicates, args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Completed {len(rows)} synthetic fits across {len(summaries)} stress-test scenarios.")
    for summary in summaries:
        print(summary)


if __name__ == "__main__":
    main()
