"""Independently reconstruct catalogue-only H0 posteriors on expanded support.

Uses saved host-density arrays and detected-event inputs from a completed
five-arm artifact. It recomputes the selected-event likelihood, verifies the
archived [50, 90] curve, and compares expanded [30, 110] posterior summaries
on 561- and 1121-point grids. This is a diagnostic, not a new model fit.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import ndtr

C_KM_S = 299792.458


def _trap_weights(x: np.ndarray) -> np.ndarray:
    dx = np.diff(x)
    w = np.empty_like(x)
    w[0], w[-1] = dx[0] / 2, dx[-1] / 2
    w[1:-1] = (dx[:-1] + dx[1:]) / 2
    return w


def _catalogue_log_likelihood(h0: np.ndarray, z: np.ndarray, density: np.ndarray,
                              distances: np.ndarray, sigmas: np.ndarray,
                              threshold: float, chunk: int = 128) -> np.ndarray:
    qw = _trap_weights(z)
    result = np.zeros(h0.size)
    for start in range(0, h0.size, chunk):
        stop = min(start + chunk, h0.size)
        means = C_KM_S * z[None, :] / h0[start:stop, None]
        chunk_ll = np.zeros(stop - start)
        for distance, sigma in zip(distances, sigmas):
            x = (distance - means) / sigma
            numerator_kernel = np.exp(-0.5 * x * x) / (np.sqrt(2 * np.pi) * sigma)
            numerator = (numerator_kernel * qw[None, :]) @ density
            efficiency = (ndtr((means - threshold) / sigma) * qw[None, :]) @ density
            chunk_ll += np.log(numerator) - np.log(efficiency)
        result[start:stop] = chunk_ll
    return result


def _posterior_summary(grid: np.ndarray, logl: np.ndarray) -> dict[str, float]:
    unnorm = np.exp(logl - np.max(logl))
    density = unnorm / trapezoid(unnorm, grid)
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf = np.maximum.accumulate(cdf / cdf[-1])
    probs = np.array([0.025, 0.05, 0.16, 0.50, 0.84, 0.95, 0.975])
    quantiles = np.interp(probs, cdf, grid)

    def mass(lo: float, hi: float) -> float:
        mask = (grid >= lo) & (grid <= hi)
        return float(trapezoid(density[mask], grid[mask]))

    return {
        "q025": float(quantiles[0]), "q05": float(quantiles[1]),
        "q16": float(quantiles[2]), "median": float(quantiles[3]),
        "q84": float(quantiles[4]), "q95": float(quantiles[5]),
        "q975": float(quantiles[6]),
        "mass_30_40": mass(30, 40), "mass_30_45": mass(30, 45),
        "mass_30_50": mass(30, 50), "mass_30_55": mass(30, 55),
        "mass_105_110": mass(105, 110),
    }


def audit(root: Path, output: Path) -> dict:
    root = root.resolve()
    cfg = json.loads((root / "config.json").read_text())
    reps = int(cfg["replicates"])
    gammas = [float(x) for x in cfg["gammas"]]
    sizes = [int(x) for x in cfg["validation_sizes"]]
    if 500 not in sizes:
        raise AssertionError("artifact does not include validation size 500")
    h0_min, h0_max = 30.0, 110.0
    rows = []
    max_reconstruction_error = 0.0
    for rep in range(reps):
        for gamma in gammas:
            suffix = f"g{gamma:+g}".replace("+", "p").replace("-", "m")
            child = root / "replicates" / f"rep{rep:03d}_{suffix}"
            archive = np.load(child / "results" / "log_likelihoods.npz")
            z = archive["redshift_grid"]
            old_grid = archive["H0_grid"]
            event_frame = pd.read_csv(child / "observable" / "gw_events_observed.csv")
            distances = event_frame.distance_mpc.to_numpy(float)
            sigmas = np.full(distances.size, float(cfg["sigma_mpc"]))
            density = archive["density_A_catalogue_only_n500"]
            old_ll = _catalogue_log_likelihood(
                old_grid, z, density, distances, sigmas,
                float(cfg["detection_threshold_mpc"]),
            )
            stored_ll = archive["A_catalogue_only_n500"]
            err = float(np.max(np.abs(old_ll - stored_ll)))
            max_reconstruction_error = max(max_reconstruction_error, err)
            if err > 1e-9:
                raise AssertionError(f"archived likelihood mismatch rep={rep}, gamma={gamma}: {err}")
            summaries = {}
            for npoints in (561, 1121):
                grid = np.linspace(h0_min, h0_max, npoints)
                ll = _catalogue_log_likelihood(
                    grid, z, density, distances, sigmas,
                    float(cfg["detection_threshold_mpc"]),
                )
                summaries[npoints] = _posterior_summary(grid, ll)
            row = {
                "replicate": rep, "gamma": gamma, "original_prior_min": 50.0,
                "expanded_prior_min": h0_min, "expanded_prior_max": h0_max,
                "max_original_grid_loglike_abs_diff": err,
            }
            for npoints in (561, 1121):
                for key, value in summaries[npoints].items():
                    row[f"{key}_{npoints}"] = value
            row["median_grid_abs_diff"] = abs(summaries[561]["median"] - summaries[1121]["median"])
            row["q025_grid_abs_diff"] = abs(summaries[561]["q025"] - summaries[1121]["q025"])
            row["mass_30_50_grid_abs_diff"] = abs(summaries[561]["mass_30_50"] - summaries[1121]["mass_30_50"])
            rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(output, index=False)
    summary = frame.groupby("gamma").agg(
        replicates=("replicate", "count"),
        mean_median_561=("median_561", "mean"),
        mean_q025_561=("q025_561", "mean"),
        mean_mass_30_50_561=("mass_30_50_561", "mean"),
        mean_mass_30_55_561=("mass_30_55_561", "mean"),
        max_median_grid_abs_diff=("median_grid_abs_diff", "max"),
        max_q025_grid_abs_diff=("q025_grid_abs_diff", "max"),
        max_mass_grid_abs_diff=("mass_30_50_grid_abs_diff", "max"),
    ).reset_index()
    report = {
        "status": "PASS" if max_reconstruction_error <= 1e-9
            and frame.median_grid_abs_diff.max() <= 0.001
            and frame.q025_grid_abs_diff.max() <= 0.01
            and frame.mass_30_50_grid_abs_diff.max() <= 0.001 else "FAIL",
        "replicate_gamma_rows": len(frame),
        "max_original_grid_loglike_abs_diff": max_reconstruction_error,
        "max_median_grid_abs_diff": float(frame.median_grid_abs_diff.max()),
        "max_q025_grid_abs_diff": float(frame.q025_grid_abs_diff.max()),
        "max_mass_30_50_grid_abs_diff": float(frame.mass_30_50_grid_abs_diff.max()),
        "expanded_prior": [h0_min, h0_max],
        "summary_by_gamma": summary.to_dict(orient="records"),
        "interpretation": "Catalogue-only posterior boundary diagnostic only; no estimator calibration or observational inference.",
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--output", type=Path, default=Path("catalogue_only_prior_boundary.csv"))
    args = p.parse_args()
    print(json.dumps(audit(args.artifact_dir, args.output), indent=2))


if __name__ == "__main__":
    main()
