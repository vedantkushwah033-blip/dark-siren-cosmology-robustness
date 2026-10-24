"""Map fixed-truth coverage across the prior for the toy selected-data model.

Synthetic experiment only. Credible intervals need not have nominal pointwise
frequentist coverage, so this script measures that behavior instead of assuming it.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import detection_efficiency_h0, posterior_h0_selected

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
TRUTHS = (52.0, 55.0, 60.0, 70.0, 80.0, 85.0, 88.0)
HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
SIGMA_MPC = 20.0
DETECTION_MIN, DETECTION_MAX = 0.0, 180.0

def draw_detected_observation(h0_true, rng, max_attempts=100000):
    """Sample the pre-detection host mixture and condition the noisy distance on detection."""
    for _ in range(max_attempts):
        host = int(rng.choice(len(HOST_REDSHIFTS), p=HOST_WEIGHTS))
        mean = C_KM_S * HOST_REDSHIFTS[host] / h0_true
        distance = float(rng.normal(mean, SIGMA_MPC))
        if DETECTION_MIN < distance < DETECTION_MAX:
            return host, distance
    raise RuntimeError("Exceeded attempts while sampling a detected observation.")

def wilson_interval(successes, n, z=1.959963984540054):
    if n <= 0:
        raise ValueError("n must be positive")
    p = successes / n
    denom = 1.0 + z*z/n
    center = (p + z*z/(2*n)) / denom
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    low = 0.0 if successes == 0 else max(0.0, float(center-half))
    high = 1.0 if successes == n else min(1.0, float(center+half))
    return low, high

def summarize_cell(rows, method):
    truth = float(rows[0]["true_H0"])
    median = np.asarray([r[f"{method}_median"] for r in rows], dtype=float)
    lo68 = np.asarray([r[f"{method}_lower_68"] for r in rows], dtype=float)
    hi68 = np.asarray([r[f"{method}_upper_68"] for r in rows], dtype=float)
    lo95 = np.asarray([r[f"{method}_lower_95"] for r in rows], dtype=float)
    hi95 = np.asarray([r[f"{method}_upper_95"] for r in rows], dtype=float)
    error = median - truth
    hit68 = (lo68 <= truth) & (truth <= hi68)
    hit95 = (lo95 <= truth) & (truth <= hi95)
    w68 = wilson_interval(int(hit68.sum()), len(rows))
    w95 = wilson_interval(int(hit95.sum()), len(rows))
    return {
        "n": len(rows), "bias": float(error.mean()),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "coverage_68": float(hit68.mean()), "coverage_68_wilson_low": w68[0],
        "coverage_68_wilson_high": w68[1],
        "coverage_95": float(hit95.mean()), "coverage_95_wilson_low": w95[0],
        "coverage_95_wilson_high": w95[1],
        "mean_width_68": float(np.mean(hi68-lo68)),
        "mean_width_95": float(np.mean(hi95-lo95)),
    }

def run_audit(replicates_per_truth=500, seed=20261030, grid_size=4001):
    if replicates_per_truth < 50:
        raise ValueError("Use at least 50 replicates per truth.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    beta = detection_efficiency_h0(grid, HOST_REDSHIFTS, SIGMA_MPC,
        DETECTION_MAX, weights=HOST_WEIGHTS, detection_distance_min=DETECTION_MIN)
    prior_mask = (grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)
    rows = []
    summaries = []
    for truth in TRUTHS:
        cell = []
        for rep in range(replicates_per_truth):
            host, distance = draw_detected_observation(truth, rng)
            base = posterior_h0(grid, HOST_REDSHIFTS, distance, SIGMA_MPC,
                weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
            corrected = posterior_h0_selected(grid, HOST_REDSHIFTS, distance,
                SIGMA_MPC, DETECTION_MAX, weights=HOST_WEIGHTS,
                prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                detection_distance_min=DETECTION_MIN)
            bsum, csum = posterior_summary(grid, base), posterior_summary(grid, corrected)
            row = {"true_H0": truth, "replicate": rep, "host_index": host,
                   "observed_distance_mpc": distance,
                   "beta_at_truth_interpolated": float(np.interp(truth, grid, beta))}
            for label, sm in (("uncorrected", bsum), ("corrected", csum)):
                for key in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    row[f"{label}_{key}"] = float(sm[key])
            row["paired_abs_error_change_corrected_minus_uncorrected"] = (
                abs(row["corrected_median"]-truth)-abs(row["uncorrected_median"]-truth))
            rows.append(row); cell.append(row)
        for method in ("uncorrected", "corrected"):
            summaries.append({"true_H0": truth, "method": method,
                **summarize_cell(cell, method), "seed": seed,
                "replicates_per_truth": replicates_per_truth, "grid_size": grid_size,
                "prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX,
                "sigma_mpc": SIGMA_MPC, "detection_min_mpc": DETECTION_MIN,
                "detection_max_mpc": DETECTION_MAX,
                "beta_min_on_prior_grid": float(beta[prior_mask].min()),
                "beta_max_on_prior_grid": float(beta[prior_mask].max())})
        change = np.asarray([r["paired_abs_error_change_corrected_minus_uncorrected"] for r in cell])
        summaries.append({"true_H0": truth,
            "method": "paired_corrected_minus_uncorrected_absolute_error",
            "n": len(cell), "mean_paired_change": float(change.mean()),
            "paired_change_mcse": float(change.std(ddof=1)/np.sqrt(len(change))),
            "seed": seed, "replicates_per_truth": replicates_per_truth,
            "grid_size": grid_size, "prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX})
    return rows, summaries

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates-per-truth", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/fixed-truth-coverage-map")
    args = parser.parse_args()
    rows, summaries = run_audit(args.replicates_per_truth, args.seed, args.grid_size)
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out/name).open("w", newline="", encoding="utf-8") as f:
            fields = list(dict.fromkeys(k for row in data for k in row))
            writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(data)
    print(f"PASS: {len(rows)} paired observations; {len(summaries)} fixed-truth summaries.")
    for row in summaries: print(row)
if __name__ == "__main__":
    main()
