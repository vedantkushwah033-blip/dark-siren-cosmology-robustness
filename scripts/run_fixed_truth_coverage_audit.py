"""Audit fixed-truth credible-interval coverage across H0 and host knowledge.

This is a synthetic diagnostic. It reuses the original project's event simulator,
then computes central posterior intervals directly from the returned posterior
density to make interval construction explicit and auditable.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

H0_VALUES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
GRID_MIN, GRID_MAX, GRID_SIZE = 50.0, 90.0, 1601


def normalized_cdf(grid: np.ndarray, posterior: np.ndarray) -> np.ndarray:
    if np.any(~np.isfinite(posterior)) or np.any(posterior < 0):
        raise ValueError("Posterior contains invalid density values")
    areas = (posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0
    cdf = np.concatenate(([0.0], np.cumsum(areas)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed")
    return cdf / cdf[-1]


def central_interval(grid: np.ndarray, cdf: np.ndarray, level: float) -> tuple[float, float]:
    alpha = (1.0 - level) / 2.0
    low = float(np.interp(alpha, cdf, grid))
    high = float(np.interp(1.0 - alpha, cdf, grid))
    if low > high:
        raise ValueError("Interval endpoints reversed")
    return low, high


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--output-dir", default="results/fixed-truth-coverage-audit")
    args = parser.parse_args()
    if args.replicates < 50:
        parser.error("--replicates must be at least 50")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "bayesian_inference.py").is_file():
        raise FileNotFoundError(f"Original simulator source not found: {original}")
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0

    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_SIZE)
    rows: list[dict] = []
    for h_index, h0_true in enumerate(H0_VALUES):
        for replicate in range(args.replicates):
            seed = args.seed + h_index * 1_000_000 + replicate * 10
            cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=h0_true, seed=seed)
            host = choose_host(cat, weighting="uniform", seed=seed + 1)
            gw = simulate_gw_event(cat["true_distance_mpc"][host], seed=seed + 2)
            scenarios = {
                "known_host": (np.asarray([cat["redshift"][host]], dtype=float), np.asarray([1.0])),
                "unknown_host_complete_catalogue": (
                    np.asarray(cat["redshift"], dtype=float),
                    np.ones(len(cat["redshift"]), dtype=float),
                ),
            }
            for scenario, (redshifts, weights) in scenarios.items():
                posterior = np.asarray(posterior_h0(
                    grid, redshifts, gw["observed_distance_mpc"], gw["sigma_mpc"],
                    weights=weights, prior_min=GRID_MIN, prior_max=GRID_MAX,
                ), dtype=float)
                cdf = normalized_cdf(grid, posterior)
                lo68, hi68 = central_interval(grid, cdf, 0.68)
                lo95, hi95 = central_interval(grid, cdf, 0.95)
                median = float(np.interp(0.5, cdf, grid))
                rows.append({
                    "true_H0": h0_true,
                    "replicate": replicate,
                    "scenario": scenario,
                    "posterior_median": median,
                    "median_error": median - h0_true,
                    "lower_68": lo68,
                    "upper_68": hi68,
                    "width_68": hi68 - lo68,
                    "covers_68": int(lo68 <= h0_true <= hi68),
                    "lower_95": lo95,
                    "upper_95": hi95,
                    "width_95": hi95 - lo95,
                    "covers_95": int(lo95 <= h0_true <= hi95),
                })

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "per_run.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for h0_true in H0_VALUES:
        for scenario in ("known_host", "unknown_host_complete_catalogue"):
            subset = [r for r in rows if r["true_H0"] == h0_true and r["scenario"] == scenario]
            n = len(subset)
            c68 = float(np.mean([r["covers_68"] for r in subset]))
            c95 = float(np.mean([r["covers_95"] for r in subset]))
            errors = np.asarray([r["median_error"] for r in subset], dtype=float)
            summary.append({
                "true_H0": h0_true,
                "scenario": scenario,
                "replicates": n,
                "mean_signed_error": float(errors.mean()),
                "rmse_median": float(np.sqrt(np.mean(errors**2))),
                "coverage_68": c68,
                "coverage_68_mc_se": float(np.sqrt(c68 * (1-c68) / n)),
                "coverage_95": c95,
                "coverage_95_mc_se": float(np.sqrt(c95 * (1-c95) / n)),
                "mean_width_68": float(np.mean([r["width_68"] for r in subset])),
                "mean_width_95": float(np.mean([r["width_95"] for r in subset])),
            })
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f"Validated design: {len(H0_VALUES)} fixed truths × 2 scenarios × {args.replicates} replicates = {len(rows)} posterior fits.")
    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
