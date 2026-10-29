"""Run a synthetic hidden-host-fraction × event-selection-cutoff experiment.

This is a simplified low-redshift Gaussian-distance toy, not real LVK inference.
The hidden-host fraction is defined before event selection. Each data set is
shared across all f_model arms to support paired comparisons.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm, truncnorm

from dark_siren_robustness.selection_conditioned_hidden_host import (
    C_KM_S,
    selection_conditioned_hidden_host_log_likelihood,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary

F_TRUE = 0.30
F_MODELS = (0.10, 0.30, 0.50)
D_MAX_LEVELS = (500.0, 1000.0, 2000.0)
H0_TRUE = 70.0
SIGMA_MPC = 140.0
D_MIN = 0.0


def stable_seed(seed: int, configuration: dict) -> int:
    payload = json.dumps([int(seed), configuration], sort_keys=True, separators=(",", ":"))
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "big") % (2**32 - 1)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_gaussian_density(z: np.ndarray, mean: float, sd: float) -> np.ndarray:
    density = np.exp(-0.5 * ((z - mean) / sd) ** 2)
    area = float(trapezoid(density, z))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("Could not normalize component redshift density")
    return density / area


def sample_from_grid_density(z: np.ndarray, density: np.ndarray, rng: np.random.Generator) -> float:
    cdf = np.r_[0.0, cumulative_trapezoid(density, z)]
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ArithmeticError("Cannot sample a non-positive redshift density")
    cdf /= cdf[-1]
    # Remove duplicate CDF values from numerically negligible tails before interpolation.
    unique_cdf, idx = np.unique(cdf, return_index=True)
    if unique_cdf.size < 2:
        raise ArithmeticError("Redshift CDF is degenerate")
    return float(np.interp(rng.random(), unique_cdf, z[idx]))


def generate_selected_dataset(
    *, z: np.ndarray, g_catalogue: np.ndarray, g_missing: np.ndarray,
    h0_true: float, f_true: float, sigma: float, d_min: float, d_max: float,
    n_events: int, seed: int,
) -> list[dict]:
    """Sample from the pre-selection mixture conditional on observed-distance selection."""
    rng = np.random.default_rng(seed)
    means = C_KM_S * z / h0_true
    acceptance = norm.cdf((d_max - means) / sigma) - norm.cdf((d_min - means) / sigma)
    alpha_c = float(trapezoid(g_catalogue * acceptance, z))
    alpha_m = float(trapezoid(g_missing * acceptance, z))
    beta = (1.0 - f_true) * alpha_c + f_true * alpha_m
    if not np.isfinite(beta) or beta <= 0:
        raise ArithmeticError("True mixture has non-positive selection efficiency")
    p_missing_detected = f_true * alpha_m / beta
    rows = []
    for event_index in range(n_events):
        is_missing = bool(rng.random() < p_missing_detected)
        density = g_missing if is_missing else g_catalogue
        selected_density = density * acceptance
        selected_density /= float(trapezoid(selected_density, z))
        z_true = sample_from_grid_density(z, selected_density, rng)
        mu = C_KM_S * z_true / h0_true
        a = (d_min - mu) / sigma
        b = (d_max - mu) / sigma
        distance = float(truncnorm.rvs(a, b, loc=mu, scale=sigma, random_state=rng))
        if not d_min <= distance <= d_max:
            raise ArithmeticError("Truncated-normal generator returned an unselected distance")
        rows.append({
            "event_index": event_index,
            "z_true_for_audit_only": z_true,
            "component_for_audit_only": "missing" if is_missing else "catalogue",
            "distance_mpc": distance,
            "d_min_true_mpc": d_min,
            "d_max_true_mpc": d_max,
        })
    return rows


def wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid Wilson interval counts")
    z = 1.959963984540054
    p = successes / trials
    denom = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / denom
    half = z*np.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / denom
    return max(0.0, float(center-half)), min(1.0, float(center+half))


def run_experiment(
    *, replicates: int, n_events: int, seed: int,
    h0_points: int, z_points: int, d_max_levels: tuple[float, ...] = D_MAX_LEVELS,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if replicates < 1 or n_events < 1:
        raise ValueError("replicates and n_events must be positive")
    if h0_points < 31 or z_points < 101:
        raise ValueError("h0_points must be >=31 and z_points >=101")
    if len(d_max_levels) < 1 or any((not np.isfinite(x) or x <= D_MIN) for x in d_max_levels):
        raise ValueError("d_max_levels must be finite and greater than d_min")
    z = np.linspace(0.01, 0.40, z_points)
    h0_grid = np.linspace(50.0, 90.0, h0_points)
    g_catalogue = normalized_gaussian_density(z, 0.12, 0.05)
    g_missing = normalized_gaussian_density(z, 0.25, 0.07)
    event_rows: list[dict] = []
    fit_rows: list[dict] = []
    started = time.perf_counter()

    for d_max in d_max_levels:
        for rep in range(replicates):
            configuration = {"d_max_true_mpc": float(d_max), "replicate": int(rep)}
            data_seed = stable_seed(seed, configuration)
            events = generate_selected_dataset(
                z=z, g_catalogue=g_catalogue, g_missing=g_missing,
                h0_true=H0_TRUE, f_true=F_TRUE, sigma=SIGMA_MPC,
                d_min=D_MIN, d_max=float(d_max), n_events=n_events, seed=data_seed,
            )
            for row in events:
                event_rows.append({
                    "replicate": rep, "data_seed": data_seed,
                    "d_max_true_mpc": float(d_max), **row,
                })
            distances = np.asarray([row["distance_mpc"] for row in events], dtype=float)
            for f_model in F_MODELS:
                fit = selection_conditioned_hidden_host_log_likelihood(
                    z, g_catalogue, g_missing, h0_grid, distances,
                    sigma=SIGMA_MPC, f_missing=f_model, d_min=D_MIN, d_max=float(d_max),
                )
                posterior = posterior_summary(h0_grid, fit["log_likelihood"])
                density = np.asarray(json.loads(posterior["posterior_density_grid_json"]), dtype=float)
                cdf = np.r_[0.0, cumulative_trapezoid(density, h0_grid)]
                cdf /= cdf[-1]
                ci95_low, median, ci95_high = (
                    float(np.interp(q, cdf, h0_grid)) for q in (0.025, 0.5, 0.975)
                )
                fit_rows.append({
                    "replicate": rep, "data_seed": data_seed,
                    "d_max_true_mpc": float(d_max), "d_max_model_mpc": float(d_max),
                    "f_true_preselection": F_TRUE, "f_model_preselection": f_model,
                    "h0_true": H0_TRUE, "sigma_mpc": SIGMA_MPC, "n_events": n_events,
                    "posterior_mean": posterior["posterior_mean"],
                    "posterior_median": median,
                    "ci68_low": posterior["ci68_low"], "ci68_high": posterior["ci68_high"],
                    "ci95_low": ci95_low, "ci95_high": ci95_high,
                    "ci68_covered": bool(posterior["ci68_low"] <= H0_TRUE <= posterior["ci68_high"]),
                    "ci95_covered": bool(ci95_low <= H0_TRUE <= ci95_high),
                    "ci68_width": posterior["ci68_high"] - posterior["ci68_low"],
                    "posterior_integral": posterior["posterior_integral"],
                    "posterior_density_grid_json": posterior["posterior_density_grid_json"],
                    "log_likelihood_grid_json": json.dumps(fit["log_likelihood"].tolist(), separators=(",", ":")),
                    "beta_grid_json": json.dumps(fit["beta"].tolist(), separators=(",", ":")),
                    "h0_grid_json": json.dumps(h0_grid.tolist(), separators=(",", ":")),
                })

    events_frame = pd.DataFrame(event_rows)
    fits_frame = pd.DataFrame(fit_rows)
    aggregates = []
    for d_max in d_max_levels:
        for f_model in F_MODELS:
            sub = fits_frame[(fits_frame.d_max_true_mpc == d_max) & (fits_frame.f_model_preselection == f_model)]
            errors = sub.posterior_median.to_numpy() - H0_TRUE
            cov68 = int(sub.ci68_covered.sum())
            cov95 = int(sub.ci95_covered.sum())
            lo68, hi68 = wilson(cov68, len(sub))
            lo95, hi95 = wilson(cov95, len(sub))
            aggregates.append({
                "d_max_true_mpc": float(d_max), "f_model_preselection": float(f_model),
                "n_replicates": int(len(sub)), "bias_median": float(np.mean(errors)),
                "rmse_median": float(np.sqrt(np.mean(errors**2))),
                "coverage68_count": cov68, "coverage68_denominator": int(len(sub)),
                "coverage68": cov68 / len(sub), "coverage68_wilson95": [lo68, hi68],
                "coverage95_count": cov95, "coverage95_denominator": int(len(sub)),
                "coverage95": cov95 / len(sub), "coverage95_wilson95": [lo95, hi95],
                "mean_ci68_width": float(sub.ci68_width.mean()),
            })
    contrasts = []
    for d_max in d_max_levels:
        baseline = fits_frame[(fits_frame.d_max_true_mpc == d_max) & (fits_frame.f_model_preselection == F_TRUE)].set_index("replicate")
        for f_model in F_MODELS:
            if np.isclose(f_model, F_TRUE):
                continue
            other = fits_frame[(fits_frame.d_max_true_mpc == d_max) & (fits_frame.f_model_preselection == f_model)].set_index("replicate")
            delta = np.abs(other.posterior_median.to_numpy() - H0_TRUE) - np.abs(baseline.posterior_median.to_numpy() - H0_TRUE)
            contrasts.append({
                "d_max_true_mpc": float(d_max), "f_model_preselection": float(f_model),
                "reference_f_model": F_TRUE,
                "estimand": "abs_error_misspecified_f_minus_abs_error_correct_f",
                "positive_means": "correct_f_model_has_smaller_absolute_median_error",
                "mean_paired_contrast": float(np.mean(delta)),
                "standard_error": float(np.std(delta, ddof=1) / np.sqrt(len(delta))) if len(delta) > 1 else None,
                "n_pairs": int(len(delta)),
            })
    summary = {
        "status": "synthetic_hidden_host_fraction_selection_cutoff_pilot",
        "claims": "toy simulation only; not real LVK data or an observational H0 measurement",
        "seed": int(seed), "replicates_per_cutoff": int(replicates), "n_events_per_dataset": int(n_events),
        "h0_true": H0_TRUE, "f_true_preselection": F_TRUE, "f_model_levels": list(F_MODELS),
        "d_max_levels_mpc": list(map(float, d_max_levels)), "sigma_mpc": SIGMA_MPC,
        "distance_acceptance_region": "[0, d_max] Mpc; Gaussian distances truncated to accepted region",
        "redshift_support": [float(z[0]), float(z[-1]), int(z_points)],
        "h0_grid": [float(h0_grid[0]), float(h0_grid[-1]), int(h0_points)],
        "catalogue_component_density": "normalized Gaussian on redshift grid, mean 0.12, SD 0.05",
        "missing_component_density": "normalized Gaussian on redshift grid, mean 0.25, SD 0.07",
        "fraction_definition": "pre-selection host-population mixture fraction",
        "selection_model": "hard observed-distance cut; same true and assumed cutoff in this pilot",
        "aggregate_metrics": aggregates, "paired_contrasts": contrasts,
        "runtime_seconds": float(time.perf_counter() - started),
        "limitations": [
            "Synthetic Gaussian-distance low-redshift toy; no real GW strain or detector response.",
            "Component redshift densities are declared toy distributions, not empirical galaxy populations.",
            "No event-rate information, peculiar velocities, redshift uncertainty, or survey-level catalogue-generation model.",
            "The pilot varies the selection regime but does not yet misspecify the model cutoff.",
            "PPC/SBC diagnostics are not included in this first implementation phase.",
            "Results are conditional on the fixed H0 truth and declared prior; coverage is imprecise at small replicate counts.",
        ],
    }
    return events_frame, fits_frame, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--n-events", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--h0-points", type=int, default=61)
    parser.add_argument("--z-points", type=int, default=151)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--d-max-levels", type=float, nargs="+", default=list(D_MAX_LEVELS))
    args = parser.parse_args()
    started = time.perf_counter()
    events, fits, summary = run_experiment(
        replicates=args.replicates, n_events=args.n_events, seed=args.seed,
        h0_points=args.h0_points, z_points=args.z_points,
        d_max_levels=tuple(args.d_max_levels),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    events_path = args.output_dir / "events.csv"
    fits_path = args.output_dir / "fits.csv"
    summary_path = args.output_dir / "summary.json"
    events.to_csv(events_path, index=False)
    fits.to_csv(fits_path, index=False)
    summary["runtime_seconds_wall_including_serialization"] = time.perf_counter() - started
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        commit = __import__("subprocess").run(["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
    except Exception:
        commit = "unavailable"
    packages = {}
    for package in ("numpy", "scipy", "pandas"):
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = "unavailable"
    manifest = {
        "experiment": "hidden_fraction_selection_cutoff_v1",
        "git_commit": commit, "python": sys.version, "platform": platform.platform(),
        "packages": packages, "command": " ".join(sys.argv),
        "configuration": {"replicates": args.replicates, "n_events": args.n_events, "seed": args.seed,
                          "h0_points": args.h0_points, "z_points": args.z_points,
                          "d_max_levels_mpc": args.d_max_levels, "f_true_preselection": F_TRUE,
                          "f_model_levels": F_MODELS, "sigma_mpc": SIGMA_MPC},
        "files": {
            events_path.name: {"sha256": sha256(events_path), "bytes": events_path.stat().st_size},
            fits_path.name: {"sha256": sha256(fits_path), "bytes": fits_path.stat().st_size},
            summary_path.name: {"sha256": sha256(summary_path), "bytes": summary_path.stat().st_size},
        },
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
