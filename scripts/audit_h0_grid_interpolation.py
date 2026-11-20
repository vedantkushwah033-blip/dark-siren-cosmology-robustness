#!/usr/bin/env python3
"""Audit H0-grid quadrature using archived likelihood curves.

Interpolation-only diagnostic: this does not recompute likelihoods at new H0
values. It measures changes in posterior integration when the archived
281-point log-likelihood curves are integrated on finer grids using PCHIP.
"""
from __future__ import annotations
import argparse, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid, trapezoid

GAMMAS = ("gamma_m1", "gamma_0", "gamma_p1")
ARMS = ("A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw", "D_map_plugin_ipw", "E_uncertainty_integrated_ipw")
DEGREES, ORDERS, VALIDATION = (1, 2), (5, 7, 9), (100, 500)
REFINEMENTS = (561, 1121, 2241)

def summarize(grid: np.ndarray, loglike: np.ndarray):
    shifted = np.asarray(loglike, float) - np.max(loglike)
    raw = np.exp(shifted)
    area = float(trapezoid(raw, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("non-finite posterior normalization")
    density = raw / area
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    return area, density, np.interp([.025, .16, .5, .84, .975], cdf, grid)

def audit(root: Path):
    rows = []
    for gamma in GAMMAS:
        for degree in DEGREES:
            for order in ORDERS:
                path = root / gamma / f"degree{degree}_order{order}" / "results" / "log_likelihoods.npz"
                if not path.is_file():
                    raise FileNotFoundError(path)
                with np.load(path) as data:
                    x = np.asarray(data["H0_grid"], float)
                    if x.shape != (281,) or not np.allclose(x, np.linspace(30, 110, 281), rtol=0, atol=1e-12):
                        raise ValueError(f"unexpected archived H0 grid in {path}")
                    for n in VALIDATION:
                        for arm in ARMS:
                            y = np.asarray(data[f"{arm}_n{n}"], float)
                            if y.shape != x.shape or not np.all(np.isfinite(y)):
                                raise ValueError(f"invalid likelihood curve in {path}")
                            area0, _, q0 = summarize(x, y)
                            interp = PchipInterpolator(x, y - np.max(y), extrapolate=False)
                            for m in REFINEMENTS:
                                g = np.linspace(x[0], x[-1], m)
                                raw = np.exp(interp(g))
                                area = float(trapezoid(raw, g))
                                if not np.isfinite(area) or area <= 0:
                                    raise ValueError("invalid refined normalization")
                                density = raw / area
                                cdf = np.r_[0.0, cumulative_trapezoid(density, g)]
                                cdf /= cdf[-1]
                                q = np.interp([.025, .16, .5, .84, .975], cdf, g)
                                density0_fine = raw / area0
                                rows.append({
                                    "gamma": gamma, "degree": degree, "quadrature_order": order,
                                    "validation_size": n, "arm": arm, "grid_points": m,
                                    "posterior_l1_vs_281": float(trapezoid(np.abs(density-density0_fine), g)),
                                    "normalization_relative_change": abs(area-area0)/area,
                                    "median_shift": float(q[2]-q0[2]),
                                    "q025_shift": float(q[0]-q0[0]),
                                    "q975_shift": float(q[4]-q0[4]),
                                    "median": float(q[2]), "q025": float(q[0]), "q975": float(q[4]),
                                })
    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_zip", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("h0_grid_interpolation_audit"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.artifact_zip) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError(f"corrupt zip member: {bad}")
        archive.extractall(args.output_dir / "extracted")
    result = audit(args.output_dir / "extracted")
    result.to_csv(args.output_dir / "h0_grid_interpolation_cells.csv", index=False)
    summary = result.groupby(["arm", "grid_points"]).agg(
        max_posterior_l1=("posterior_l1_vs_281", "max"),
        max_normalization_relative_change=("normalization_relative_change", "max"),
        max_abs_median_shift=("median_shift", lambda s: float(np.max(np.abs(s)))),
        max_abs_q025_shift=("q025_shift", lambda s: float(np.max(np.abs(s)))),
        max_abs_q975_shift=("q975_shift", lambda s: float(np.max(np.abs(s)))),
    ).reset_index()
    summary.to_csv(args.output_dir / "h0_grid_interpolation_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"PASS: {len(result)} posterior curves; interpolation-only diagnostic.")

if __name__ == "__main__":
    main()
