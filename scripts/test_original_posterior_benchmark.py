"""Integration benchmark for the pinned original posterior_h0 implementation.

This test imports the actual inference function from the pinned source checkout,
rather than reconstructing the likelihood independently.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    args = parser.parse_args()
    original = Path(args.original_dir).resolve()
    sys.path.insert(0, str(original))

    from src.bayesian_inference import posterior_h0

    c_km_s = 299792.458
    h0_true = 72.0
    z = 0.035
    observed_distance = c_km_s * z / h0_true
    sigma_mpc = 0.01 * observed_distance
    grid = np.linspace(30.0, 120.0, 18001)

    posterior = np.asarray(
        posterior_h0(
            grid,
            np.asarray([z]),
            observed_distance,
            sigma_mpc,
            weights=np.asarray([1.0]),
            prior_min=30.0,
            prior_max=120.0,
        ),
        dtype=float,
    )
    if posterior.shape != grid.shape:
        raise SystemExit(f"Posterior shape {posterior.shape} != grid shape {grid.shape}")
    if not np.all(np.isfinite(posterior)) or np.any(posterior < 0):
        raise SystemExit("Original posterior_h0 returned invalid density values")
    area = float(np.trapezoid(posterior, grid))
    if not np.isclose(area, 1.0, rtol=0.0, atol=2e-3):
        raise SystemExit(f"Posterior density does not integrate to 1: area={area}")

    # Integrate the actual returned density and compute its median.
    increments = (posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    median = float(np.interp(0.5, cdf, grid))
    if abs(median - h0_true) > 0.5:
        raise SystemExit(
            f"Single-host benchmark failed: injected H0={h0_true}, posterior median={median}"
        )
    print(
        "PASS: pinned posterior_h0 returned a finite normalized density; "
        f"single-host injected H0={h0_true:.2f}, posterior median={median:.4f}."
    )

if __name__ == "__main__":
    main()
