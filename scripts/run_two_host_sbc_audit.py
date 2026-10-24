"""Prior-predictive SBC for a correctly specified two-host distance mixture.

H0 is drawn from the same uniform prior used for inference. A latent host is
drawn from the same catalogue weights used in the mixture likelihood, then a
distance is generated with fixed Gaussian uncertainty. Posterior CDF ranks
should be approximately Uniform(0,1) when the generator and likelihood match.
This is a toy statistical calibration test, not an observational analysis.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import chisquare, kstest

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0

PRIOR_MIN = 50.0
PRIOR_MAX = 90.0
HOST_REDSHIFTS = np.array([0.02, 0.04])
SCENARIOS = {
    "two_host_70_30": np.array([0.70, 0.30]),
    "two_host_equal": np.array([0.50, 0.50]),
}
Z_REFERENCE = 0.03
SIGMA_FIXED_MPC = 0.10 * C_KM_S * Z_REFERENCE / 70.0


def posterior_cdf_at(grid, density, value):
    """Evaluate posterior CDF on the nonzero support, excluding exterior zero bins."""
    grid = np.asarray(grid, dtype=float)
    density = np.asarray(density, dtype=float)
    if grid.ndim != 1 or density.ndim != 1 or grid.size != density.size or grid.size < 2:
        raise ValueError("grid and density must be same-length one-dimensional arrays.")
    if not np.all(np.isfinite(grid)) or not np.all(np.diff(grid) > 0):
        raise ValueError("grid must be finite and strictly increasing.")
    if not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ValueError("density must be finite and non-negative.")
    support = np.flatnonzero(density > 0)
    if support.size < 2:
        raise ValueError("density must have at least two positive support points.")
    lo, hi = int(support[0]), int(support[-1])
    support_grid = grid[lo:hi + 1]
    support_density = density[lo:hi + 1]
    cdf = np.concatenate(([0.0], cumulative_trapezoid(support_density, support_grid)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Invalid posterior CDF normalization.")
    cdf /= cdf[-1]
    return float(np.clip(np.interp(value, support_grid, cdf), 0.0, 1.0))


def generate_replicates(replicates, seed):
    if replicates < 200:
        raise ValueError("Use at least 200 replicates per host-weight scenario.")
    rng = np.random.default_rng(seed)
    data = []
    for scenario, weights in SCENARIOS.items():
        for replicate in range(replicates):
            h0_true = float(rng.uniform(PRIOR_MIN, PRIOR_MAX))
            host_index = int(rng.choice(len(HOST_REDSHIFTS), p=weights))
            z_host = float(HOST_REDSHIFTS[host_index])
            true_distance = C_KM_S * z_host / h0_true
            observed_distance = float(rng.normal(true_distance, SIGMA_FIXED_MPC))
            while observed_distance <= 0:
                observed_distance = float(rng.normal(true_distance, SIGMA_FIXED_MPC))
            data.append({
                "scenario": scenario,
                "replicate": replicate,
                "true_H0": h0_true,
                "host_index": host_index,
                "host_redshift": z_host,
                "observed_distance_mpc": observed_distance,
                "sigma_mpc": SIGMA_FIXED_MPC,
            })
    return data


def run_audit(replicates=1000, seed=20261026, grid_sizes=(1001, 4001)):
    data = generate_replicates(replicates, seed)
    all_rows, summaries = [], []
    for grid_size in grid_sizes:
        grid = np.linspace(25.0, 125.0, int(grid_size))
        for scenario, weights in SCENARIOS.items():
            subset = [row for row in data if row["scenario"] == scenario]
            ranks = []
            for row in subset:
                density = posterior_h0(
                    grid,
                    HOST_REDSHIFTS,
                    row["observed_distance_mpc"],
                    row["sigma_mpc"],
                    weights=weights,
                    prior_min=PRIOR_MIN,
                    prior_max=PRIOR_MAX,
                )
                rank = posterior_cdf_at(grid, density, row["true_H0"])
                ranks.append(rank)
                all_rows.append({
                    **row,
                    "grid_size": int(grid_size),
                    "posterior_cdf_rank": rank,
                })
            ranks = np.asarray(ranks)
            ks = kstest(ranks, "uniform")
            counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, 11))
            chi = chisquare(counts, f_exp=np.full(10, len(ranks) / 10.0))
            summaries.append({
                "scenario": scenario,
                "grid_size": int(grid_size),
                "replicates": len(ranks),
                "mean_rank": float(ranks.mean()),
                "rank_variance": float(ranks.var(ddof=1)),
                "ks_statistic": float(ks.statistic),
                "ks_pvalue": float(ks.pvalue),
                "decile_chi_square": float(chi.statistic),
                "decile_chi_square_pvalue": float(chi.pvalue),
                **{f"decile_{i+1}_count": int(counts[i]) for i in range(10)},
            })
    rank_map = {(r["scenario"], r["replicate"], r["grid_size"]): r["posterior_cdf_rank"] for r in all_rows}
    resolution_diffs = []
    if len(grid_sizes) >= 2:
        coarse, fine = int(grid_sizes[0]), int(grid_sizes[-1])
        for scenario in SCENARIOS:
            diffs = [
                abs(rank_map[(scenario, i, coarse)] - rank_map[(scenario, i, fine)])
                for i in range(replicates)
            ]
            resolution_diffs.append({
                "scenario": scenario,
                "coarse_grid_size": coarse,
                "fine_grid_size": fine,
                "median_abs_rank_difference": float(np.median(diffs)),
                "p95_abs_rank_difference": float(np.quantile(diffs, 0.95)),
                "max_abs_rank_difference": float(np.max(diffs)),
            })
    return all_rows, summaries, resolution_diffs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--grid-sizes", type=int, nargs="+", default=[1001, 4001])
    parser.add_argument("--output-dir", default="results/two-host-sbc-audit")
    args = parser.parse_args()
    rows, summaries, diffs = run_audit(args.replicates, args.seed, tuple(args.grid_sizes))
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries), ("grid_resolution.csv", diffs)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Validated design: {len(SCENARIOS)} host-weight scenarios × {args.replicates} draws × {len(args.grid_sizes)} grid resolutions = {len(rows)} posterior fits.")
    for row in summaries:
        print(row)
    for row in diffs:
        print({"grid_resolution": row})


if __name__ == "__main__":
    main()
