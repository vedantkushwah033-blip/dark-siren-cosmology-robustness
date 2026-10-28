"""Audit whether a fixed posterior grid resolves narrow single-host likelihoods.

This deterministic stress test compares the production implementation against a
much finer grid. It is a numerical diagnostic for the simplified model, not a
realistic gravitational-wave simulation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import (
    C_KM_S,
    posterior_h0,
    posterior_summary,
)

COARSE_POINTS = 1601
FINE_POINTS = 400001
SIGMAS_MPC = (2.0, 0.2, 0.05, 0.02, 0.01)
QUANTILES = ("lower_68", "upper_68", "lower_95", "upper_95", "median")


def _summary(sigma_mpc: float, points: int) -> dict[str, float]:
    grid = np.linspace(50.0, 90.0, points)
    z = 0.05
    distance = C_KM_S * z / 70.0
    density = posterior_h0(
        grid, [z], distance, sigma_mpc,
        prior_min=50.0, prior_max=90.0,
    )
    return posterior_summary(grid, density)


def run_audit(
    *, coarse_points: int = COARSE_POINTS,
    fine_points: int = FINE_POINTS,
    quantile_tolerance: float = 0.001,
    sd_tolerance: float = 0.001,
) -> dict:
    """Compare coarse summaries with fine-grid reference summaries."""
    if coarse_points < 3 or fine_points <= coarse_points:
        raise ValueError("fine_points must be greater than coarse_points >= 3")
    if quantile_tolerance <= 0 or sd_tolerance <= 0:
        raise ValueError("tolerances must be positive")
    rows = []
    for sigma in SIGMAS_MPC:
        coarse = _summary(sigma, coarse_points)
        fine = _summary(sigma, fine_points)
        quantile_diffs = {
            key: abs(coarse[key] - fine[key]) for key in QUANTILES
        }
        max_quantile_diff = max(quantile_diffs.values())
        sd_diff = abs(coarse["posterior_sd"] - fine["posterior_sd"])
        rows.append({
            "sigma_distance_mpc": sigma,
            "coarse_points": coarse_points,
            "fine_points": fine_points,
            "coarse_spacing_h0": 40.0 / (coarse_points - 1),
            "fine_spacing_h0": 40.0 / (fine_points - 1),
            "coarse_posterior_sd": coarse["posterior_sd"],
            "fine_posterior_sd": fine["posterior_sd"],
            "absolute_sd_difference": sd_diff,
            "absolute_quantile_differences": quantile_diffs,
            "maximum_absolute_quantile_difference": max_quantile_diff,
            "status": (
                "UNDERRESOLVED"
                if max_quantile_diff > quantile_tolerance or sd_diff > sd_tolerance
                else "WITHIN_TOLERANCE"
            ),
        })
    return {
        "scope": "deterministic single-host synthetic numerical stress test",
        "not_real_gw_data": True,
        "quantile_tolerance_h0": quantile_tolerance,
        "posterior_sd_tolerance_h0": sd_tolerance,
        "rows": rows,
        "underresolved_case_count": sum(r["status"] == "UNDERRESOLVED" for r in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--coarse-points", type=int, default=COARSE_POINTS)
    parser.add_argument("--fine-points", type=int, default=FINE_POINTS)
    parser.add_argument("--quantile-tolerance", type=float, default=0.001)
    parser.add_argument("--sd-tolerance", type=float, default=0.001)
    args = parser.parse_args()
    report = run_audit(
        coarse_points=args.coarse_points,
        fine_points=args.fine_points,
        quantile_tolerance=args.quantile_tolerance,
        sd_tolerance=args.sd_tolerance,
    )
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
