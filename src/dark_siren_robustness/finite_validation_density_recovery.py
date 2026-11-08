"""Repository-native finite-validation IPW host-density recovery benchmark.

This module is deliberately separate from the H0 likelihood benchmark. It
scores reconstructed merger-weighted host-redshift densities against a declared
analytic synthetic target; it makes no posterior-calibration claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import scipy
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import t as student_t

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.estimated_inclusion import fit_inclusion_logistic

Z_MIN, Z_MAX = 0.01, 0.40
PRIMARY_GAMMA = (-1.2, 0.0, 1.2)
PRIMARY_RETENTION = 0.25
PRIMARY_VALIDATION_SIZES = (100, 300, 1000)
PRIMARY_DEGREES = (1, 2)
PRIMARY_MULTIPLIER = 0.15
SENSITIVITY_MULTIPLIERS = (0.40, 1.00)


def parent_intensity(z: np.ndarray) -> np.ndarray:
    """Unnormalised, declared two-overdensity parent intensity."""
    z = np.asarray(z, dtype=float)
    return (
        1.0
        + 2.2 * np.exp(-0.5 * ((z - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((z - 0.275) / 0.040) ** 2)
    )


def _parent_cdf_grid(points: int = 20001):
    grid = np.linspace(Z_MIN, Z_MAX, points)
    intensity = parent_intensity(grid)
    cdf = np.r_[0.0, cumulative_trapezoid(intensity, grid)]
    cdf /= cdf[-1]
    return grid, cdf


def _draw_parent_redshifts(rng: np.random.Generator, n: int) -> np.ndarray:
    grid, cdf = _parent_cdf_grid()
    return np.interp(rng.random(n), cdf, grid)


def _curvature_center() -> float:
    grid = np.linspace(Z_MIN, Z_MAX, 20001)
    x = (grid - Z_MIN) / (Z_MAX - Z_MIN)
    return float(trapezoid(parent_intensity(grid) * x**2, grid)
                 / trapezoid(parent_intensity(grid), grid))


def _calibrate_intercept(gamma: float, target_retention: float, center: float) -> float:
    grid = np.linspace(Z_MIN, Z_MAX, 20001)
    x = (grid - Z_MIN) / (Z_MAX - Z_MIN)
    intensity = parent_intensity(grid)
    denominator = trapezoid(intensity, grid)

    def mean_retention(alpha: float) -> float:
        probability = expit(alpha - x + gamma * (x**2 - center))
        return float(trapezoid(intensity * probability, grid) / denominator)

    return float(brentq(lambda a: mean_retention(a) - target_retention, -30.0, 30.0))


def _analytic_target(grid: np.ndarray) -> np.ndarray:
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    raw = parent_intensity(grid) * merger
    return raw / trapezoid(raw, grid)


def _weight_diagnostics(inclusion_probability: np.ndarray | None, n: int) -> dict:
    weights = (np.ones(n, dtype=float) if inclusion_probability is None
               else 1.0 / np.asarray(inclusion_probability, dtype=float))
    scaled = weights / np.max(weights)
    normalized = scaled / np.sum(scaled)
    entropy = -float(np.sum(normalized * np.log(normalized)))
    ordered = np.sort(normalized)[::-1]
    return {
        "kish_ess": float(np.sum(scaled)**2 / np.sum(scaled**2)),
        "entropy_ess": float(np.exp(entropy)),
        "max_normalized_weight": float(np.max(normalized)),
        "top_1pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(0.01*n)))])),
        "top_5pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(0.05*n)))])),
    }


def _density_iae(estimate: np.ndarray, truth: np.ndarray, grid: np.ndarray) -> float:
    return float(trapezoid(np.abs(estimate - truth), grid))


def _safe_density(grid, sample, probabilities, merger, multiplier):
    return estimate_host_density(
        grid, sample, inclusion_probability=probabilities,
        merger_weight=merger, scott_multiplier=multiplier,
    )


def _paired_contrast_summary(attempts: pd.DataFrame, *, n_replicates: int,
                             seed: int = 20261108, bootstrap_resamples: int = 2000) -> list[dict]:
    """Summarize replicate-paired estimated-IPW minus catalogue-only IAE."""
    results: list[dict] = []
    group_keys = ["gamma", "target_retention", "n_validation", "bandwidth_multiplier"]
    rng = np.random.default_rng(seed)
    for key, group in attempts.groupby(group_keys, sort=True):
        primary = group.loc[group.arm.isin(
            ["catalogue_only", "estimated_degree_1_ipw", "estimated_degree_2_ipw"]
        )]
        wide = primary.pivot(index="replicate", columns="arm", values="iae")
        if "catalogue_only" not in wide:
            continue
        for degree in (1, 2):
            arm = f"estimated_degree_{degree}_ipw"
            if arm not in wide:
                continue
            paired = wide[["catalogue_only", arm]].dropna()
            difference = (paired[arm] - paired["catalogue_only"]).to_numpy(dtype=float)
            n = int(difference.size)
            record = {
                "gamma": float(key[0]), "target_retention": float(key[1]),
                "n_validation": int(key[2]), "bandwidth_multiplier": float(key[3]),
                "degree": degree, "contrast": f"{arm}_minus_catalogue_only",
                "n_pairs": n, "n_replicates_planned": int(n_replicates),
                "mean_difference": float(np.mean(difference)) if n else None,
                "monte_carlo_se": float(np.std(difference, ddof=1)/np.sqrt(n))
                    if n > 1 else None,
                "bootstrap_95_interval": None,
                "paired_t_95_interval": None,
            }
            if n >= 2:
                draws = rng.integers(0, n, size=(bootstrap_resamples, n))
                boot_means = difference[draws].mean(axis=1)
                record["bootstrap_95_interval"] = [
                    float(np.quantile(boot_means, 0.025)),
                    float(np.quantile(boot_means, 0.975)),
                ]
                se = float(np.std(difference, ddof=1)/np.sqrt(n))
                critical = float(student_t.ppf(0.975, n-1))
                mean = float(np.mean(difference))
                record["paired_t_95_interval"] = [mean-critical*se, mean+critical*se]
            results.append(record)
    return results


def run_density_recovery_benchmark(
    *,
    n_replicates: int = 100,
    n_parent: int = 4000,
    seed: int = 20261108,
    gamma_values: tuple[float, ...] = PRIMARY_GAMMA,
    retention_values: tuple[float, ...] = (PRIMARY_RETENTION,),
    validation_sizes: tuple[int, ...] = PRIMARY_VALIDATION_SIZES,
    degrees: tuple[int, ...] = PRIMARY_DEGREES,
    bandwidth_multipliers: tuple[float, ...] = (PRIMARY_MULTIPLIER, *SENSITIVITY_MULTIPLIERS),
    grid_points: int = 801,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run paired density-recovery experiments and retain parent-level raw data.

    Validation sets are nested by size within each replicate and sampled
    independently of retention labels. Validation units are excluded from the
    analysis catalogue. True inclusion probabilities are never passed to the
    fitted logistic models.
    """
    if n_replicates < 1 or n_parent < 1100 or grid_points < 201:
        raise ValueError("require replicates >= 1, parent size >= 1100, grid >= 201")
    if any(n < 5 or n >= n_parent for n in validation_sizes):
        raise ValueError("validation sizes must be >=5 and smaller than parent size")
    if not retention_values or any(not 0.0 < r < 1.0 for r in retention_values):
        raise ValueError("retention_values must be nonempty and lie strictly between 0 and 1")
    if not gamma_values or any(not np.isfinite(g) for g in gamma_values):
        raise ValueError("gamma_values must be nonempty and finite")
    if tuple(sorted(set(validation_sizes))) != tuple(validation_sizes):
        raise ValueError("validation_sizes must be unique and increasing")
    if not degrees or any(d not in (1, 2) for d in degrees):
        raise ValueError("degrees must contain only 1 and/or 2")
    if any(b <= 0 or not np.isfinite(b) for b in bandwidth_multipliers):
        raise ValueError("bandwidth multipliers must be finite and positive")

    grid = np.linspace(Z_MIN, Z_MAX, grid_points)
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    truth = _analytic_target(grid)
    center = _curvature_center()
    attempt_rows: list[dict] = []
    unit_rows: list[dict] = []

    for rep in range(n_replicates):
        parent_rng = np.random.default_rng(seed + 1000003 * rep)
        parent_z = _draw_parent_redshifts(parent_rng, n_parent)
        object_ids = np.arange(n_parent, dtype=np.int64)
        split_rng = np.random.default_rng(seed + 2000003 * rep)
        validation_order = split_rng.permutation(n_parent)
        validation_rank = np.empty(n_parent, dtype=np.int64)
        validation_rank[validation_order] = np.arange(n_parent, dtype=np.int64)
        nested_validation = {
            n: np.sort(validation_order[:n]) for n in validation_sizes
        }
        base_retention_uniform = np.random.default_rng(seed + 3000017 * rep).random(n_parent)

        for gamma in gamma_values:
            for target_retention in retention_values:
                alpha = _calibrate_intercept(float(gamma), float(target_retention), center)
                x_parent = (parent_z - Z_MIN) / (Z_MAX - Z_MIN)
                true_probability = expit(
                    alpha - x_parent + gamma * (x_parent**2 - center)
                )
                retained = (base_retention_uniform < true_probability).astype(np.int8)
                for n_validation in validation_sizes:
                    validation_idx = nested_validation[n_validation]
                    is_validation = np.zeros(n_parent, dtype=bool)
                    is_validation[validation_idx] = True
                    analysis_mask = (~is_validation) & (retained == 1)
                    analysis_idx = np.flatnonzero(analysis_mask)
                    val_z = parent_z[validation_idx]
                    val_y = retained[validation_idx]
                    analysis_z = parent_z[analysis_idx]
                    analysis_pi = true_probability[analysis_idx]

                    # One raw parent row per object/scenario/split size makes the
                    # disjoint partition and oracle truth independently inspectable.
                    for idx in range(n_parent):
                        unit_rows.append({
                            "replicate": rep, "gamma": float(gamma),
                            "target_retention": float(target_retention),
                            "n_validation": int(n_validation),
                            "object_id": int(object_ids[idx]),
                            "redshift": float(parent_z[idx]),
                            "true_inclusion_probability": float(true_probability[idx]),
                            "retention_uniform": float(base_retention_uniform[idx]),
                            "validation_rank": int(validation_rank[idx]),
                            "retained": int(retained[idx]),
                            "role": ("VALIDATION" if is_validation[idx]
                                     else ("ANALYSIS_CATALOGUE" if retained[idx]
                                           else "NON_RETAINED_HOLDOUT")),
                        })

                    base_arms = {
                        "catalogue_only": (analysis_z, None),
                        "known_inclusion_ipw": (analysis_z, analysis_pi),
                        "parent_sample_kde": (parent_z, None),
                    }
                    fits = {
                        degree: fit_inclusion_logistic(val_z, val_y, degree=degree)
                        for degree in degrees
                    }
                    for multiplier in bandwidth_multipliers:
                        for arm, (sample, probabilities) in base_arms.items():
                            status, message, density = "SUCCESS", "", None
                            try:
                                density = _safe_density(
                                    grid, sample, probabilities, merger, multiplier
                                )
                            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                                status, message = "NON_EVALUABLE", str(exc)
                            weight_metrics = (_weight_diagnostics(probabilities, int(sample.size))
                                              if sample.size else {})
                            attempt_rows.append({
                                "replicate": rep, "seed": seed + 1000003 * rep,
                                "gamma": float(gamma), "target_retention": float(target_retention),
                                "n_validation": int(n_validation), "degree": 0,
                                "bandwidth_multiplier": float(multiplier), "arm": arm,
                                "status": status, "message": message,
                                "n_parent": n_parent, "n_validation_actual": n_validation,
                                "n_analysis_catalogue": int(analysis_z.size),
                                "iae": (_density_iae(density, truth, grid)
                                        if density is not None else np.nan),
                                "density_json": (json.dumps(density.tolist(), separators=(",", ":"))
                                                 if density is not None else ""),
                                **weight_metrics,
                                "ise": (float(trapezoid((density-truth)**2, grid))
                                        if density is not None else np.nan),
                                "max_abs_density_error": (float(np.max(np.abs(density-truth)))
                                                          if density is not None else np.nan),
                                "coefficients_json": "", "fit_status": "NOT_APPLICABLE",
                                "fit_message": "",
                            })
                        for degree, fit in fits.items():
                            arm = f"estimated_degree_{degree}_ipw"
                            status, message, density = fit.status, fit.message, None
                            if fit.status == "SUCCESS":
                                try:
                                    estimated_pi = fit.predict(analysis_z)
                                    density = _safe_density(
                                        grid, analysis_z, estimated_pi, merger, multiplier
                                    )
                                    status, message = "SUCCESS", ""
                                except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                                    status, message = "NON_EVALUABLE", str(exc)
                            weight_metrics = (_weight_diagnostics(estimated_pi, int(analysis_z.size))
                                              if density is not None else {})
                            attempt_rows.append({
                                "replicate": rep, "seed": seed + 1000003 * rep,
                                "gamma": float(gamma), "target_retention": float(target_retention),
                                "n_validation": int(n_validation), "degree": int(degree),
                                "bandwidth_multiplier": float(multiplier), "arm": arm,
                                "status": status, "message": message,
                                "n_parent": n_parent, "n_validation_actual": n_validation,
                                "n_analysis_catalogue": int(analysis_z.size),
                                "iae": (_density_iae(density, truth, grid)
                                        if density is not None else np.nan),
                                "density_json": (json.dumps(density.tolist(), separators=(",", ":"))
                                                 if density is not None else ""),
                                **weight_metrics,
                                "ise": (float(trapezoid((density-truth)**2, grid))
                                        if density is not None else np.nan),
                                "max_abs_density_error": (float(np.max(np.abs(density-truth)))
                                                          if density is not None else np.nan),
                                "coefficients_json": (
                                    json.dumps(fit.coefficients.tolist(), separators=(",", ":"))
                                    if fit.coefficients is not None else ""),
                                "fit_status": fit.status, "fit_message": fit.message,
                                "fit_center": fit.center, "fit_scale": fit.scale,
                                "fit_covariance_json": (
                                    json.dumps(fit.covariance.tolist(), separators=(",", ":"))
                                    if fit.covariance is not None else ""),
                                "fit_gradient_inf_norm": fit.gradient_inf_norm,
                                "fit_hessian_condition_number": fit.hessian_condition_number,
                            })

    attempts = pd.DataFrame(attempt_rows)
    units = pd.DataFrame(unit_rows)
    success = attempts.loc[attempts.status == "SUCCESS"]
    summary = {
        "artifact_schema_version": 1,
        "status": "synthetic_density_recovery_benchmark_no_H0_inference",
        "estimand": "analytic merger-weighted parent-host redshift density",
        "primary_metric": "integrated absolute error against analytic target",
        "dgp": {
            "parent_intensity": "1 + 2.2*N(0.105,0.025^2) + 1.6*N(0.275,0.040^2), unnormalised",
            "curvature_center": center,
            "inclusion_model": "expit(alpha - x + gamma*(x^2-c))",
            "intercept_calibration": "parent-intensity weighted quadrature to each declared retention target",
            "validation_design": "nested SRS from complete parent frame, independent of retention",
        },
        "n_replicates": n_replicates, "n_parent": n_parent, "seed_start": seed,
        "gamma_values": list(gamma_values), "retention_values": list(retention_values),
        "validation_sizes": list(validation_sizes), "degrees": list(degrees),
        "bandwidth_multipliers": list(bandwidth_multipliers),
        "grid": [Z_MIN, Z_MAX, grid_points],
        "analytic_target_density": truth.tolist(),
        "analytic_target_sha256": hashlib.sha256(
            json.dumps(truth.tolist(), separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "attempt_count": int(len(attempts)),
        "status_counts": {str(k): int(v) for k, v in attempts.status.value_counts().items()},
        "success_count": int((attempts.status == "SUCCESS").sum()),
        "conditional_success_mean_iae": {
            str(k): float(v.iae.mean())
            for k, v in success.groupby("arm", sort=True)
        },
        "paired_ipw_minus_catalogue_only_contrasts": _paired_contrast_summary(
            attempts, n_replicates=n_replicates, seed=seed + 9109
        ),
        "limitations": [
            "Synthetic low-redshift density recovery only; no H0 likelihood or posterior coverage.",
            "The 0.15 Scott multiplier is a frozen candidate from separate exploratory work, not universal.",
            "Parent-sample KDE is a finite-sample diagnostic, not the analytic target.",
            "Model-based inverse-information covariance is not a Bayesian posterior.",
        ],
    }
    return attempts, units, summary


def write_artifacts(attempts: pd.DataFrame, units: pd.DataFrame, summary: dict,
                    output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "attempts": output_dir / "attempt_metrics.csv",
        "units": output_dir / "parent_units.csv",
        "summary": output_dir / "summary.json",
        "config": output_dir / "config.json",
    }
    config_keys = (
        "artifact_schema_version", "status", "dgp", "n_replicates", "n_parent",
        "seed_start", "gamma_values", "retention_values", "validation_sizes",
        "degrees", "bandwidth_multipliers", "grid",
    )
    config = {key: summary[key] for key in config_keys}
    paths["config"].write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    attempts.to_csv(paths["attempts"], index=False)
    units.to_csv(paths["units"], index=False)
    paths["summary"].write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    hashes = {
        key: hashlib.sha256(path.read_bytes()).hexdigest()
        for key, path in paths.items()
    }
    manifest = {
        "artifact_schema_version": 1,
        "status": summary["status"],
        "files": {
            key: {"path": path.name, "bytes": path.stat().st_size, "sha256": hashes[key]}
            for key, path in paths.items()
        },
        "attempt_count": summary["attempt_count"],
        "success_count": summary["success_count"],
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "pandas": pd.__version__,
            "github_sha": os.environ.get("GITHUB_SHA"),
            "github_ref": os.environ.get("GITHUB_REF"),
        },
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--grid-points", type=int, default=801)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("artifacts/finite_validation_density_recovery"))
    args = parser.parse_args()
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=args.replicates, n_parent=args.n_parent,
        seed=args.seed, grid_points=args.grid_points,
    )
    manifest = write_artifacts(attempts, units, summary, args.output_dir)
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
