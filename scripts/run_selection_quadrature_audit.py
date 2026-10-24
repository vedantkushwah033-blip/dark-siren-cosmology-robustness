"""Independent adaptive-quadrature check of toy selection-conditioned posterior normalization.

Synthetic numerical audit only; this is not an observational inference pipeline.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.integrate import quad, trapezoid
from scipy.stats import norm
from dark_siren_robustness.hardened_inference import C_KM_S
from dark_siren_robustness.selection import detection_efficiency_h0

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
REDSHIFTS = np.array([0.02, 0.04])
WEIGHTS = np.array([0.70, 0.30])
SIGMA_MPC = 20.0
DETECTION_MIN, DETECTION_MAX = 0.0, 180.0
DISTANCES = (50.0, 90.0, 140.0)

def unnormalized_selected_likelihood(h0: float, distance: float) -> float:
    means = C_KM_S * REDSHIFTS / h0
    base = float(np.sum(WEIGHTS * norm.pdf(distance, loc=means, scale=SIGMA_MPC)))
    beta = float(detection_efficiency_h0(
        np.array([h0, h0 + 1e-8]), REDSHIFTS, SIGMA_MPC, DETECTION_MAX,
        weights=WEIGHTS, detection_distance_min=DETECTION_MIN
    )[0])
    return base / beta

def run_audit(grid_size: int = 4001) -> list[dict]:
    if grid_size < 1001 or (grid_size - 1) % 10 != 0:
        raise ValueError("grid_size must be >=1001 and chosen so the prior endpoints lie on the grid.")
    grid = np.linspace(25.0, 125.0, grid_size)
    support = (grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)
    rows = []
    for distance in DISTANCES:
        if not DETECTION_MIN < distance < DETECTION_MAX:
            raise ValueError("Audit distance must lie strictly inside the detection interval.")
        grid_values = np.array([unnormalized_selected_likelihood(float(h), distance) for h in grid[support]])
        grid_z = float(trapezoid(grid_values, grid[support]))
        quad_z, quad_error = quad(
            lambda h: unnormalized_selected_likelihood(float(h), distance),
            PRIOR_MIN, PRIOR_MAX, epsabs=1e-11, epsrel=1e-10, limit=300
        )
        relative_error = abs(grid_z - quad_z) / abs(quad_z)
        rows.append({
            "observed_distance_mpc": distance, "grid_size": grid_size,
            "prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX,
            "grid_normalization_constant": grid_z,
            "adaptive_quadrature_normalization_constant": float(quad_z),
            "quadrature_estimated_absolute_error": float(quad_error),
            "relative_absolute_difference": float(relative_error),
            "pass_relative_tolerance_1e-4": bool(relative_error < 1e-4),
        })
    return rows

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/selection-quadrature-audit")
    args = parser.parse_args()
    rows = run_audit(args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "quadrature_comparison.csv"
    required = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=required)
        writer.writeheader()
        writer.writerows(rows)
    # Read-back schema and row-count check protects against malformed artifacts.
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != required:
            raise RuntimeError("CSV artifact schema does not match required columns.")
        read_rows = list(reader)
    if len(read_rows) != len(DISTANCES):
        raise RuntimeError("CSV artifact row count is incorrect.")
    if not all(row["pass_relative_tolerance_1e-4"] for row in rows):
        raise RuntimeError("Grid normalization differs from adaptive quadrature beyond tolerance.")
    print(f"PASS: {len(rows)} normalization checks; artifact schema validated at {path}")
    for row in rows:
        print(row)

if __name__ == "__main__":
    main()
