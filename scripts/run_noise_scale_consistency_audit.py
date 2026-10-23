"""Sensitivity experiment for H0-dependent versus fixed distance-noise scales.

This isolated known-host experiment reproduces the toy simulator's stated
noise rule sigma_true = fractional_sigma * d_true. It compares:
(1) treating the reported sigma as fixed across the H0 grid, and
(2) using sigma(H0,z) = fractional_sigma * c*z/H0 in the likelihood, including
    its Gaussian normalization.

It is a diagnostic of one-host noise-model consistency, not an LVK analysis.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

C_KM_S = 299792.458
H0_VALUES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)


def normalized_cdf(grid: np.ndarray, density: np.ndarray) -> np.ndarray:
    if np.any(~np.isfinite(density)) or np.any(density < 0):
        raise ValueError("density must be finite and nonnegative")
    areas = 0.5 * (density[1:] + density[:-1]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(areas)))
    if cdf[-1] <= 0 or not np.isfinite(cdf[-1]):
        raise ValueError("density has invalid normalization")
    return cdf / cdf[-1]


def posterior_for_event(
    grid: np.ndarray,
    z: float,
    observed_distance: float,
    reported_sigma: float,
    fractional_sigma: float,
    mode: str,
) -> np.ndarray:
    predicted_distance = C_KM_S * z / grid
    if mode == "fixed_reported_sigma":
        sigma = np.full_like(grid, reported_sigma)
    elif mode == "sigma_depends_on_h0":
        sigma = fractional_sigma * predicted_distance
    else:
        raise ValueError(f"unknown mode: {mode}")
    log_likelihood = -0.5 * ((observed_distance - predicted_distance) / sigma) ** 2
    if mode == "sigma_depends_on_h0":
        log_likelihood -= np.log(sigma)
    log_likelihood -= np.max(log_likelihood)
    likelihood = np.exp(log_likelihood)
    return likelihood / np.trapezoid(likelihood, grid)


def interval(grid: np.ndarray, cdf: np.ndarray, level: float) -> tuple[float, float]:
    tail = (1.0 - level) / 2.0
    return float(np.interp(tail, cdf, grid)), float(np.interp(1.0-tail, cdf, grid))


def run(replicates: int, seed: int, output_dir: Path, fractional_sigma: float) -> list[dict]:
    grid = np.linspace(30.0, 120.0, 1801)
    rng = np.random.default_rng(seed)
    per_run: list[dict] = []
    for h0_true in H0_VALUES:
        z = 0.035
        true_distance = C_KM_S * z / h0_true
        true_sigma = fractional_sigma * true_distance
        for rep in range(replicates):
            observed_distance = float(rng.normal(true_distance, true_sigma))
            for mode in ("fixed_reported_sigma", "sigma_depends_on_h0"):
                posterior = posterior_for_event(
                    grid, z, observed_distance, true_sigma, fractional_sigma, mode
                )
                cdf = normalized_cdf(grid, posterior)
                median = float(np.interp(0.5, cdf, grid))
                lo68, hi68 = interval(grid, cdf, 0.68)
                lo95, hi95 = interval(grid, cdf, 0.95)
                per_run.append({
                    "true_H0": h0_true,
                    "replicate": rep,
                    "mode": mode,
                    "observed_distance_mpc": observed_distance,
                    "reported_sigma_mpc": true_sigma,
                    "posterior_median": median,
                    "median_error": median - h0_true,
                    "covers_68": int(lo68 <= h0_true <= hi68),
                    "covers_95": int(lo95 <= h0_true <= hi95),
                    "width_68": hi68 - lo68,
                    "width_95": hi95 - lo95,
                })

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "per_run.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(per_run[0]))
        writer.writeheader()
        writer.writerows(per_run)

    summary = []
    for h0_true in H0_VALUES:
        for mode in ("fixed_reported_sigma", "sigma_depends_on_h0"):
            subset = [r for r in per_run if r["true_H0"] == h0_true and r["mode"] == mode]
            errors = np.asarray([r["median_error"] for r in subset])
            n = len(subset)
            c68 = float(np.mean([r["covers_68"] for r in subset]))
            c95 = float(np.mean([r["covers_95"] for r in subset]))
            summary.append({
                "true_H0": h0_true,
                "mode": mode,
                "replicates": n,
                "mean_signed_error": float(errors.mean()),
                "rmse_median": float(np.sqrt(np.mean(errors**2))),
                "coverage_68": c68,
                "coverage_68_mc_se": float(np.sqrt(c68*(1-c68)/n)),
                "coverage_95": c95,
                "coverage_95_mc_se": float(np.sqrt(c95*(1-c95)/n)),
                "mean_width_68": float(np.mean([r["width_68"] for r in subset])),
                "mean_width_95": float(np.mean([r["width_95"] for r in subset])),
            })
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f"Completed {len(per_run)} posterior fits: {len(H0_VALUES)} truths × {replicates} replicates × 2 likelihood modes.")
    for row in summary:
        print(row)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--fractional-sigma", type=float, default=0.15)
    parser.add_argument("--output-dir", default="results/noise-scale-consistency")
    args = parser.parse_args()
    if args.replicates < 100:
        parser.error("--replicates must be at least 100")
    if not 0 < args.fractional_sigma < 1:
        parser.error("--fractional-sigma must be between 0 and 1")
    run(args.replicates, args.seed, Path(args.output_dir), args.fractional_sigma)


if __name__ == "__main__":
    main()
