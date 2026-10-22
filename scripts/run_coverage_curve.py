"""Estimate coverage curves for equal-tail H0 credible intervals.

Uses the pinned original simulation and inference functions. This is a
diagnostic extension: it calculates multiple posterior quantiles directly
rather than relying only on the original hard-coded 68% and 95% summaries.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np


LEVELS = (0.50, 0.68, 0.80, 0.90, 0.95)


def equal_tail_interval(grid: np.ndarray, posterior: np.ndarray, level: float) -> tuple[float, float]:
    """Match the original posterior_summary trapezoidal CDF convention."""
    cdf = np.concatenate(
        [[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2)]
    )
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed.")
    cdf /= cdf[-1]
    alpha = (1.0 - level) / 2.0
    return float(np.interp(alpha, cdf, grid)), float(np.interp(1.0 - alpha, cdf, grid))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", default="results/original-audit/coverage_curve.csv")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0

    grid = np.linspace(50.0, 90.0, 1601)
    scenarios = {
        "unknown_host": {level: {"covered": 0, "widths": []} for level in LEVELS},
        "known_host": {level: {"covered": 0, "widths": []} for level in LEVELS},
    }

    for r in range(args.replicates):
        s = args.seed + r
        cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=70.0, seed=s)
        host = choose_host(cat, weighting="uniform", seed=s + 100000)
        gw = simulate_gw_event(cat["true_distance_mpc"][host], seed=s + 200000)

        posterior_unknown = posterior_h0(
            grid, cat["redshift"], gw["observed_distance_mpc"], gw["sigma_mpc"],
            weights=np.ones(len(cat["redshift"])),
        )
        posterior_known = posterior_h0(
            grid, np.array([cat["redshift"][host]]),
            gw["observed_distance_mpc"], gw["sigma_mpc"], weights=np.array([1.0]),
        )

        for scenario, posterior in (
            ("unknown_host", posterior_unknown),
            ("known_host", posterior_known),
        ):
            for level in LEVELS:
                low, high = equal_tail_interval(grid, posterior, level)
                stats = scenarios[scenario][level]
                stats["covered"] += int(low <= 70.0 <= high)
                stats["widths"].append(high - low)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        fields = ["scenario", "nominal_level", "replicates", "coverage", "coverage_mc_se", "mean_interval_width"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for scenario, by_level in scenarios.items():
            for level, stats in by_level.items():
                coverage = stats["covered"] / args.replicates
                writer.writerow({
                    "scenario": scenario,
                    "nominal_level": level,
                    "replicates": args.replicates,
                    "coverage": coverage,
                    "coverage_mc_se": float(np.sqrt(coverage * (1.0 - coverage) / args.replicates)),
                    "mean_interval_width": float(np.mean(stats["widths"])),
                })
    print(f"Saved coverage curve to {output}")


if __name__ == "__main__":
    main()
