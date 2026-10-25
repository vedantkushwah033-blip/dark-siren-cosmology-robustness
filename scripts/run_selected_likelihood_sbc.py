"""Run a reproducible prior-predictive PIT diagnostic for the toy selected likelihood.

This is a calibration smoke study under a matched synthetic model, not an
astrophysical validation. It writes replicate-level provenance and a compact
summary so results can be independently re-analysed.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import kstest

from dark_siren_robustness.selection import posterior_h0_selected

SEED = 20261025
N_REPLICATES = 1000
GRID_MIN, GRID_MAX, GRID_SIZE = 40.0, 100.0, 2401
REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.7, 0.3])
DISTANCE_SIGMA = 15.0
DETECTION_MIN, DETECTION_MAX = 50.0, 260.0


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054):
    """Return a two-sided Wilson score interval for a binomial proportion."""
    if trials <= 0:
        raise ValueError("trials must be positive")
    p = successes / trials
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denom
    return float(center - half), float(center + half)


def run_calibration(n_replicates: int = N_REPLICATES, seed: int = SEED):
    if n_replicates < 100:
        raise ValueError("Use at least 100 replicates for this diagnostic.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_SIZE)
    rows = []

    for replicate in range(n_replicates):
        # Draw truth once. Condition data generation on this truth by
        # rejection-sampling host and distance only.
        truth = float(rng.uniform(GRID_MIN, GRID_MAX))
        attempts = 0
        while True:
            attempts += 1
            host = int(rng.choice(len(REDSHIFTS), p=HOST_WEIGHTS))
            mean = 299792.458 * REDSHIFTS[host] / truth
            distance = float(rng.normal(mean, DISTANCE_SIGMA))
            if DETECTION_MIN < distance < DETECTION_MAX:
                break
            if attempts > 1_000_000:
                raise RuntimeError("Selection rejection sampler exceeded safety limit.")

        posterior = posterior_h0_selected(
            grid,
            REDSHIFTS,
            distance,
            DISTANCE_SIGMA,
            DETECTION_MAX,
            weights=HOST_WEIGHTS,
            prior_min=GRID_MIN,
            prior_max=GRID_MAX,
            detection_distance_min=DETECTION_MIN,
        )
        increments = 0.5 * (posterior[1:] + posterior[:-1]) * np.diff(grid)
        cdf = np.concatenate(([0.0], np.cumsum(increments)))
        cdf /= cdf[-1]
        pit = float(np.interp(truth, grid, cdf))
        rows.append({
            "replicate": replicate,
            "seed": seed,
            "truth_h0": truth,
            "host_index": host,
            "host_redshift": float(REDSHIFTS[host]),
            "observed_distance_mpc": distance,
            "rejection_attempts": attempts,
            "pit": pit,
        })

    pits = np.asarray([row["pit"] for row in rows], dtype=float)
    counts, edges = np.histogram(pits, bins=np.linspace(0.0, 1.0, 11))
    bins = []
    for i, count in enumerate(counts):
        lo, hi = wilson_interval(int(count), n_replicates)
        bins.append({
            "left": float(edges[i]),
            "right": float(edges[i + 1]),
            "count": int(count),
            "frequency": float(count / n_replicates),
            "wilson_95_low": lo,
            "wilson_95_high": hi,
        })

    ks = kstest(pits, "uniform")
    summary = {
        "status": "synthetic matched-model calibration diagnostic",
        "interpretation_limit": (
            "A pass or non-significant KS result is not proof of exact calibration, "
            "model correctness, or astrophysical validity."
        ),
        "n_replicates": n_replicates,
        "seed": seed,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "h0_prior": {"distribution": "uniform", "min": GRID_MIN, "max": GRID_MAX},
        "host_redshifts": REDSHIFTS.tolist(),
        "host_weights": HOST_WEIGHTS.tolist(),
        "distance_sigma_mpc": DISTANCE_SIGMA,
        "detection_interval_mpc_open": [DETECTION_MIN, DETECTION_MAX],
        "pit_mean": float(np.mean(pits)),
        "pit_mean_monte_carlo_se": float(np.std(pits, ddof=1) / np.sqrt(n_replicates)),
        "pit_variance": float(np.var(pits, ddof=1)),
        "ks_statistic": float(ks.statistic),
        "ks_p_value_descriptive_only": float(ks.pvalue),
        "max_abs_decile_frequency_deviation": float(np.max(np.abs(counts / n_replicates - 0.1))),
        "decile_bins_with_wilson_95_intervals": bins,
    }
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=N_REPLICATES)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--output-dir", default="artifacts/selected_likelihood_sbc")
    args = parser.parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows, summary = run_calibration(args.replicates, args.seed)
    raw_path = outdir / "pit_replicates.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = outdir / "pit_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"replicates={summary['n_replicates']} seed={summary['seed']}")
    print(f"mean PIT={summary['pit_mean']:.4f} (MC SE {summary['pit_mean_monte_carlo_se']:.4f})")
    print(f"KS statistic={summary['ks_statistic']:.4f} (descriptive p={summary['ks_p_value_descriptive_only']:.4g})")
    print(f"max absolute decile-frequency deviation={summary['max_abs_decile_frequency_deviation']:.4f}")
    print(f"wrote {raw_path} and {summary_path}")


if __name__ == "__main__":
    main()
