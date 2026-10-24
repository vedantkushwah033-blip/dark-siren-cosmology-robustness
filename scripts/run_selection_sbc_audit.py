"""Prior-predictive SBC for the toy selection-normalized two-host model.

H0 is drawn from the analysis prior. For each fixed H0, joint host/distance draws
are repeated until the toy detection rule d_obs < threshold is satisfied. The
posterior then uses the same host mixture divided by beta(H0). This conditions
the data on detection without reweighting the H0 prior itself.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import chisquare, kstest

from dark_siren_robustness.hardened_inference import C_KM_S
from dark_siren_robustness.selection import detection_efficiency_h0, posterior_h0_selected

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
SIGMA_MPC = 20.0
DETECTION_DISTANCE_MAX = 180.0


def posterior_cdf_at(grid, density, value):
    grid = np.asarray(grid, dtype=float)
    density = np.asarray(density, dtype=float)
    support = np.flatnonzero(density > 0)
    if support.size < 2:
        raise ValueError("Posterior must have at least two positive support points.")
    lo, hi = int(support[0]), int(support[-1])
    x, y = grid[lo:hi + 1], density[lo:hi + 1]
    increments = 0.5 * (y[1:] + y[:-1]) * np.diff(x)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Invalid posterior CDF normalization.")
    return float(np.clip(np.interp(value, x, cdf / cdf[-1]), 0.0, 1.0))


def draw_detected_observation(h0_true, rng, max_attempts=100000):
    """Draw host and noisy distance jointly conditional on toy detection."""
    for _ in range(max_attempts):
        host_index = int(rng.choice(len(HOST_REDSHIFTS), p=HOST_WEIGHTS))
        mean = C_KM_S * HOST_REDSHIFTS[host_index] / h0_true
        distance = float(rng.normal(mean, SIGMA_MPC))
        if 0.0 < distance < DETECTION_DISTANCE_MAX:
            return host_index, distance
    raise RuntimeError("Exceeded attempts while sampling a detected observation.")


def run_audit(replicates=2000, seed=20261027, grid_size=4001):
    if replicates < 500:
        raise ValueError("Use at least 500 prior-predictive replicates.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    rows = []
    for replicate in range(replicates):
        # H0 is drawn from the analysis prior, then data are generated conditional
        # on detection at that fixed H0. Do not reject H0 draws based on detection.
        h0_true = float(rng.uniform(PRIOR_MIN, PRIOR_MAX))
        host_index, observed_distance = draw_detected_observation(h0_true, rng)
        posterior = posterior_h0_selected(
            grid, HOST_REDSHIFTS, observed_distance, SIGMA_MPC,
            DETECTION_DISTANCE_MAX, weights=HOST_WEIGHTS,
            prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
        )
        uncorrected_posterior = posterior_h0(
            grid, HOST_REDSHIFTS, observed_distance, SIGMA_MPC,
            weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
        )
        rank = posterior_cdf_at(grid, posterior, h0_true)
        uncorrected_rank = posterior_cdf_at(grid, uncorrected_posterior, h0_true)
        rows.append({
            "replicate": replicate,
            "true_H0": h0_true,
            "host_index": host_index,
            "host_redshift": float(HOST_REDSHIFTS[host_index]),
            "observed_distance_mpc": observed_distance,
            "posterior_cdf_rank": rank,
            "uncorrected_posterior_cdf_rank": uncorrected_rank,
        })
    ranks = np.asarray([r["posterior_cdf_rank"] for r in rows])
    uncorrected_ranks = np.asarray([r["uncorrected_posterior_cdf_rank"] for r in rows])
    ks = kstest(ranks, "uniform")
    counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, 11))
    chi = chisquare(counts, f_exp=np.full(10, replicates / 10.0))
    uncorrected_ks = kstest(uncorrected_ranks, "uniform")
    uncorrected_counts, _ = np.histogram(uncorrected_ranks, bins=np.linspace(0.0, 1.0, 11))
    uncorrected_chi = chisquare(uncorrected_counts, f_exp=np.full(10, replicates / 10.0))
    beta_grid = detection_efficiency_h0(
        grid, HOST_REDSHIFTS, SIGMA_MPC, DETECTION_DISTANCE_MAX, HOST_WEIGHTS,
        detection_distance_min=0.0,
    )
    summary = {
        "replicates": replicates,
        "seed": seed,
        "grid_size": grid_size,
        "prior_min": PRIOR_MIN,
        "prior_max": PRIOR_MAX,
        "threshold_mpc": DETECTION_DISTANCE_MAX,
        "sigma_mpc": SIGMA_MPC,
        "mean_rank": float(ranks.mean()),
        "rank_variance": float(ranks.var(ddof=1)),
        "ks_statistic": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
        "decile_chi_square": float(chi.statistic),
        "decile_chi_square_pvalue": float(chi.pvalue),
        "uncorrected_mean_rank": float(uncorrected_ranks.mean()),
        "uncorrected_rank_variance": float(uncorrected_ranks.var(ddof=1)),
        "uncorrected_ks_statistic": float(uncorrected_ks.statistic),
        "uncorrected_ks_pvalue": float(uncorrected_ks.pvalue),
        "uncorrected_decile_chi_square": float(uncorrected_chi.statistic),
        "uncorrected_decile_chi_square_pvalue": float(uncorrected_chi.pvalue),
        "beta_min_on_prior_grid": float(beta_grid[(grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)].min()),
        "beta_max_on_prior_grid": float(beta_grid[(grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)].max()),
        **{f"decile_{i+1}_count": int(counts[i]) for i in range(10)},
        **{f"uncorrected_decile_{i+1}_count": int(uncorrected_counts[i]) for i in range(10)},
    }
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/selection-sbc-audit")
    args = parser.parse_args()
    rows, summary = run_audit(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "per_run.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)
    print(f"PASS: {len(rows)} prior-predictive selected-data SBC replicates.")
    print(summary)


if __name__ == "__main__":
    main()
