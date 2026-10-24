"""Prior-predictive SBC for a one-host selected-data toy model.

This isolates selection normalization from multi-host mixture ambiguity. It is
not a realistic gravitational-wave selection model.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.stats import chisquare, kstest
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0
from dark_siren_robustness.selection import detection_efficiency_h0, posterior_h0_selected

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
HOST_REDSHIFTS = np.array([0.03])
HOST_WEIGHTS = np.array([1.0])
SIGMA_MPC = 20.0
DETECTION_MIN, DETECTION_MAX = 0.0, 180.0

def posterior_cdf_at(grid, density, value):
    support = np.flatnonzero(np.asarray(density) > 0)
    if support.size < 2:
        raise ValueError("Posterior must have at least two positive support points.")
    lo, hi = int(support[0]), int(support[-1])
    x = np.asarray(grid)[lo:hi+1]
    y = np.asarray(density)[lo:hi+1]
    increments = 0.5*(y[1:]+y[:-1])*np.diff(x)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed.")
    return float(np.clip(np.interp(value, x, cdf/cdf[-1]), 0.0, 1.0))

def draw_detected_observation(h0_true, rng, max_attempts=100000):
    mean = C_KM_S * HOST_REDSHIFTS[0] / h0_true
    for _ in range(max_attempts):
        distance = float(rng.normal(mean, SIGMA_MPC))
        if DETECTION_MIN < distance < DETECTION_MAX:
            return distance
    raise RuntimeError("Exceeded attempts while sampling a detected observation.")

def run_audit(replicates=2000, seed=20261031, grid_size=4001):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    rows = []
    for rep in range(replicates):
        truth = float(rng.uniform(PRIOR_MIN, PRIOR_MAX))
        distance = draw_detected_observation(truth, rng)
        corrected = posterior_h0_selected(
            grid, HOST_REDSHIFTS, distance, SIGMA_MPC, DETECTION_MAX,
            weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
            detection_distance_min=DETECTION_MIN)
        uncorrected = posterior_h0(
            grid, HOST_REDSHIFTS, distance, SIGMA_MPC, weights=HOST_WEIGHTS,
            prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
        rows.append({"replicate": rep, "true_H0": truth,
            "observed_distance_mpc": distance,
            "corrected_rank": posterior_cdf_at(grid, corrected, truth),
            "uncorrected_rank": posterior_cdf_at(grid, uncorrected, truth)})
    ranks = {name: np.asarray([r[f"{name}_rank"] for r in rows]) for name in ("corrected", "uncorrected")}
    summaries = []
    for name, values in ranks.items():
        ks = kstest(values, "uniform")
        counts, _ = np.histogram(values, bins=np.linspace(0.0, 1.0, 11))
        chi = chisquare(counts, f_exp=np.full(10, replicates/10.0))
        summaries.append({
            "method": name, "replicates": replicates, "seed": seed,
            "grid_size": grid_size, "prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX,
            "host_redshift": float(HOST_REDSHIFTS[0]), "sigma_mpc": SIGMA_MPC,
            "detection_min_mpc": DETECTION_MIN, "detection_max_mpc": DETECTION_MAX,
            "mean_rank": float(values.mean()), "rank_variance": float(values.var(ddof=1)),
            "ks_statistic": float(ks.statistic), "ks_pvalue": float(ks.pvalue),
            "decile_chi_square": float(chi.statistic), "decile_chi_square_pvalue": float(chi.pvalue),
            **{f"decile_{i+1}_count": int(counts[i]) for i in range(10)},
            "beta_min_on_prior_grid": float(detection_efficiency_h0(
                grid, HOST_REDSHIFTS, SIGMA_MPC, DETECTION_MAX, HOST_WEIGHTS,
                detection_distance_min=DETECTION_MIN)[(grid>=PRIOR_MIN)&(grid<=PRIOR_MAX)].min()),
            "beta_max_on_prior_grid": float(detection_efficiency_h0(
                grid, HOST_REDSHIFTS, SIGMA_MPC, DETECTION_MAX, HOST_WEIGHTS,
                detection_distance_min=DETECTION_MIN)[(grid>=PRIOR_MIN)&(grid<=PRIOR_MAX)].max()),
        })
    return rows, summaries

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/single-host-selection-sbc")
    args = parser.parse_args()
    rows, summaries = run_audit(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out/name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    print(f"PASS: {len(rows)} single-host prior-predictive SBC replicates.")
    for row in summaries: print(row)
if __name__ == "__main__":
    main()
