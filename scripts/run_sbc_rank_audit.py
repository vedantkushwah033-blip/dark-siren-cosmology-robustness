"""Simulation-based calibration ranks for the matched single-host toy model.

For each prior, draw H0 from that same prior, generate one distance from the
specified Gaussian likelihood, infer H0 under the same prior, and evaluate the
posterior CDF at the injected H0. Under a correct prior-predictive Bayesian
calculation, these CDF values should be approximately Uniform(0, 1), up to
Monte Carlo and finite-grid error. This is not an astrophysical realism test.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import kstest, chisquare

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0

PRIORS = {"narrow_50_90": (50.0, 90.0), "broad_30_120": (30.0, 120.0)}
Z_FIXED = 0.03


def posterior_cdf_at(grid, posterior, truth):
    cdf = np.concatenate(([0.0], cumulative_trapezoid(posterior, grid)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Invalid posterior CDF normalization")
    cdf = cdf / cdf[-1]
    return float(np.clip(np.interp(truth, grid, cdf), 0.0, 1.0))


def run_sbc(replicates=1000, seed=20261025, grid_size=5001):
    if replicates < 200:
        raise ValueError("Use at least 200 prior-predictive replicates per prior.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    grid = np.linspace(25.0, 125.0, grid_size)
    rng = np.random.default_rng(seed)
    rows, summaries = [], []
    for prior_name, (prior_min, prior_max) in PRIORS.items():
        ranks = []
        for replicate in range(replicates):
            h0_true = float(rng.uniform(prior_min, prior_max))
            true_distance = C_KM_S * Z_FIXED / h0_true
            sigma = 0.10 * true_distance
            observed_distance = float(rng.normal(true_distance, sigma))
            while observed_distance <= 0:
                observed_distance = float(rng.normal(true_distance, sigma))
            posterior = posterior_h0(
                grid, [Z_FIXED], observed_distance, sigma,
                prior_min=prior_min, prior_max=prior_max,
            )
            rank = posterior_cdf_at(grid, posterior, h0_true)
            ranks.append(rank)
            rows.append({
                "prior": prior_name,
                "replicate": replicate,
                "true_H0": h0_true,
                "observed_distance_mpc": observed_distance,
                "sigma_mpc": sigma,
                "posterior_cdf_rank": rank,
            })
        ranks = np.asarray(ranks)
        ks = kstest(ranks, "uniform")
        counts, edges = np.histogram(ranks, bins=np.linspace(0.0, 1.0, 11))
        chi = chisquare(counts, f_exp=np.full(10, replicates / 10.0))
        summaries.append({
            "prior": prior_name,
            "prior_min": prior_min,
            "prior_max": prior_max,
            "replicates": replicates,
            "mean_rank": float(ranks.mean()),
            "rank_variance": float(ranks.var(ddof=1)),
            "ks_statistic": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
            "decile_chi_square": float(chi.statistic),
            "decile_chi_square_pvalue": float(chi.pvalue),
            **{f"decile_{i+1}_count": int(counts[i]) for i in range(10)},
        })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--grid-size", type=int, default=5001)
    parser.add_argument("--output-dir", default="results/sbc-rank-audit")
    args = parser.parse_args()
    rows, summaries = run_sbc(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Validated design: {len(PRIORS)} prior supports × {args.replicates} prior-predictive replicates = {len(rows)} fits.")
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
