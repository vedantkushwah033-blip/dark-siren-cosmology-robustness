"""Run a controlled repeated-sampling check for the toy single-host model.

Example: python scripts/run_single_host_calibration.py --replicates 500 --seed 20261022

This is a diagnostic for the simplified low-redshift model only. It is not an
LVK analysis and does not assess catalogue incompleteness.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dark_siren_robustness.calibration import equal_tail_interval, summarize_interval_performance  # noqa: E402
from dark_siren_robustness.inference import posterior_h0_low_z  # noqa: E402


def run(replicates: int, seed: int, distance_sigma_mpc: float) -> dict[str, float]:
    """Simulate distance observations and assess nominal 68% interval coverage."""
    if replicates < 1:
        raise ValueError("replicates must be >= 1")
    if not np.isfinite(distance_sigma_mpc) or distance_sigma_mpc <= 0:
        raise ValueError("distance_sigma_mpc must be finite and positive")

    rng = np.random.default_rng(seed)
    true_h0, redshift = 70.0, 0.02
    true_distance = 299_792.458 * redshift / true_h0
    grid = np.linspace(50.0, 90.0, 4001)
    estimates, lows, highs = [], [], []

    for _ in range(replicates):
        observed = rng.normal(true_distance, distance_sigma_mpc)
        while observed <= 0:
            observed = rng.normal(true_distance, distance_sigma_mpc)
        posterior = posterior_h0_low_z(
            grid, [redshift], float(observed), distance_sigma_mpc,
            prior_min=50.0, prior_max=90.0,
        )
        low, high = equal_tail_interval(grid, posterior, level=0.68)
        estimates.append(float(grid[np.argmax(posterior)]))
        lows.append(low)
        highs.append(high)

    metrics = summarize_interval_performance(
        np.full(replicates, true_h0), estimates, lows, highs
    )
    metrics.update({
        "true_h0": true_h0,
        "redshift": redshift,
        "distance_sigma_mpc": distance_sigma_mpc,
        "seed": float(seed),
        "nominal_interval_level": 0.68,
    })
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261022)
    parser.add_argument("--distance-sigma-mpc", type=float, default=10.0)
    parser.add_argument("--output", default="results/tables/single_host_calibration.csv")
    args = parser.parse_args()
    result = run(args.replicates, args.seed, args.distance_sigma_mpc)

    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)

    for key, value in result.items():
        print(f"{key}: {value}")
    print(f"Saved: {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
