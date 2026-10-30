"""Independent structural audit for finite-parent SRSWOR benchmark outputs.

This audit verifies design bookkeeping and replay of validation samples. It does
not independently reconstruct the statistical estimator; that is a separate
validation requirement before confirmatory claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import minimize
from scipy.special import expit, log_expit, logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _quad_weights(grid: np.ndarray) -> np.ndarray:
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0], weights[-1] = delta[0] / 2.0, delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _independent_density_and_fit(grid, observed, validation_z, validation_y, rate):
    """Reconstruct bounded logistic fit, IP weights, and reflected weighted KDE."""
    midpoint = 0.5 * (grid[0] + grid[-1])
    half_range = 0.5 * (grid[-1] - grid[0])
    x = (validation_z - midpoint) / half_range
    y = validation_y.astype(float)

    def objective(theta):
        eta = theta[0] + theta[1] * x
        return float(-np.sum(y * log_expit(eta) + (1.0-y) * log_expit(-eta)))

    fit = minimize(objective, np.zeros(2), method="L-BFGS-B",
                   bounds=((-25.0, 25.0), (-50.0, 50.0)),
                   options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8})
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError("independent logistic refit failed")
    intercept, slope = map(float, fit.x)
    p = expit(intercept + slope * ((observed-midpoint)/half_range))
    inverse = 1.0/p
    scaled = inverse/np.max(inverse)
    ess = float(np.sum(scaled)**2/np.sum(scaled**2))
    mean = float(np.average(observed, weights=scaled))
    variance = float(np.average((observed-mean)**2, weights=scaled))
    if variance <= 0 or not np.isfinite(variance):
        raise ArithmeticError("independent weighted variance invalid")
    bandwidth = ess**(-0.2)*np.sqrt(variance)
    reflected = np.concatenate((observed, 2*grid[0]-observed, 2*grid[-1]-observed))
    kernel_weights = np.tile((scaled/np.sum(scaled))/3.0, 3)
    kernels = np.exp(-0.5*((grid[:,None]-reflected[None,:])/bandwidth)**2)
    raw_density = (kernels @ kernel_weights)/(np.sqrt(2*np.pi)*bandwidth)
    raw_density *= rate
    density = raw_density/trapezoid(raw_density, grid)
    boundary = bool(np.isclose(abs(intercept),25.0,atol=1e-5,rtol=0) or
                    np.isclose(abs(slope),50.0,atol=1e-5,rtol=0))
    return density, intercept, slope, ess, float(np.max(inverse)), float(np.quantile(inverse,0.99)), boundary


def _independent_posterior(z, density, h0, distances, sigma, threshold):
    q = _quad_weights(z)
    valid = (density > 0) & (q > 0)
    log_mass = np.log(density[valid]) + np.log(q[valid])
    means = C_KM_S*z[valid][None,:]/h0[:,None]
    terms = norm.logpdf(distances[None,:,None],loc=means[:,None,:],scale=sigma) + log_mass[None,None,:]
    joint = np.sum(logsumexp(terms,axis=2),axis=1)
    detect = norm.cdf((C_KM_S*z[None,:]/h0[:,None]-threshold)/sigma)
    beta = np.sum(density[None,:]*q[None,:]*detect,axis=1)
    selected = joint-distances.size*np.log(beta)
    post = np.exp(selected-np.max(selected))
    post /= trapezoid(post,h0)
    cdf = np.r_[0.0,cumulative_trapezoid(post,h0)]
    cdf /= cdf[-1]
    low, median, high = (float(np.interp(qv,cdf,h0)) for qv in (0.16,0.50,0.84))
    mean = float(trapezoid(h0*post,h0))
    return mean, median, low, high


def _close(actual, expected, label, *, rtol=3e-6, atol=3e-7):
    if not np.isclose(float(actual),float(expected),rtol=rtol,atol=atol):
        raise AssertionError(f"{label} fails independent numerical reconstruction: {actual} != {expected}")

from dark_siren_robustness.generative_selection import generate_synthetic_population


def audit(artifact_dir: Path) -> dict:
    records = pd.read_csv(artifact_dir / "records.csv")
    summary = json.loads((artifact_dir / "summary.json").read_text(encoding="utf-8"))
    required = {
        "scenario", "draw", "validation_size", "status", "validation_seed",
        "simulation_seed", "thinning_seed", "true_intercept", "true_slope",
        "true_quadratic", "sigma_mpc", "detection_threshold_mpc",
        "validation_indices", "parent_size", "parent_frame_sha256",
        "catalogue_sha256", "gw_distance_sha256",
    }
    missing = required.difference(records.columns)
    if missing:
        raise AssertionError(f"records missing columns: {sorted(missing)}")
    if summary.get("schema") != "finite_parent_srswor_smoke_v1":
        raise AssertionError("unsupported finite-parent artifact schema")
    if summary.get("status") != "SMOKE_ONLY_NOT_CONFIRMATORY":
        raise AssertionError("artifact must retain non-confirmatory status")
    attempted = len(records)
    successful = int(records["status"].eq("success").sum())
    failed = int(records["status"].eq("failed").sum())
    if attempted != int(summary["n_attempted"]) or successful != int(summary["n_successful"]):
        raise AssertionError("summary attempt/success counts disagree with records")
    if failed != int(summary["n_failed"]) or attempted != successful + failed:
        raise AssertionError("failed-attempt accounting does not reconcile")
    if not records["status"].isin(["success", "failed"]).all():
        raise AssertionError("unknown attempt status")

    if not {"h0_points", "z_points", "h0_grid_bounds", "redshift_grid_bounds"}.issubset(summary):
        raise AssertionError("summary omits inference grid specification")
    z_bounds = summary["redshift_grid_bounds"]
    h_bounds = summary["h0_grid_bounds"]
    z_grid = np.linspace(float(z_bounds[0]), float(z_bounds[1]), int(summary["z_points"]))
    h0_grid = np.linspace(float(h_bounds[0]), float(h_bounds[1]), int(summary["h0_points"]))
    zmin, zmax = float(z_grid[0]), float(z_grid[-1])
    rate = 1.0 + 1.5*(z_grid-zmin)/(zmax-zmin)
    target = 1.0 + 2.2*np.exp(-0.5*((z_grid-0.105)/0.025)**2) + 1.6*np.exp(-0.5*((z_grid-0.275)/0.040)**2)
    target *= rate
    target /= trapezoid(target,z_grid)
    sizes = [int(v) for v in summary["validation_sizes"]]
    n_scenarios = int(summary["scenarios"])
    draws = int(summary["draws_per_scenario_and_size"])
    expected = n_scenarios * draws * len(sizes)
    if attempted != expected:
        raise AssertionError(f"expected {expected} attempts, found {attempted}")

    hash_cols = ["parent_frame_sha256", "catalogue_sha256", "gw_distance_sha256"]
    for (scenario, size), group in records.groupby(["scenario", "validation_size"]):
        if len(group) != draws:
            raise AssertionError("scenario/size group has incorrect draw count")
        if group["draw"].astype(int).tolist() != list(range(draws)):
            raise AssertionError("draw identifiers are missing or out of order")
    for scenario, group in records.groupby("scenario"):
        for col in hash_cols:
            values = group[col].dropna().astype(str).unique()
            if len(values) != 1 or len(values[0]) != 64:
                raise AssertionError(f"scenario {scenario} has inconsistent {col}")
            if any(c not in "0123456789abcdef" for c in values[0]):
                raise AssertionError(f"scenario {scenario} has malformed {col}")
        if group["validation_seed"].dropna().nunique() != 1:
            raise AssertionError("validation seed must be fixed within scenario")

        seed_values = group["validation_seed"].dropna().unique()
        if len(seed_values) != 1:
            raise AssertionError("missing validation seed")
        first = group.iloc[0]
        simulation = generate_synthetic_population(
            seed=int(first["simulation_seed"]),
            n_parent=int(first["parent_size"]),
            n_detected=int(first["n_detected"]),
            h0_true=float(first["h0_true"]),
            sigma_mpc=float(first["sigma_mpc"]),
            detection_threshold_mpc=float(first["detection_threshold_mpc"]),
            z_min=0.01, z_max=0.40,
        )
        parent_z = simulation["parent_truth"]["redshift"].to_numpy(dtype=float)
        distances = simulation["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        midpoint, half_range = 0.205, 0.195
        x_parent = (parent_z - midpoint) / half_range
        p_inclusion = expit(
            float(first["true_intercept"])
            + float(first["true_slope"]) * x_parent
            + float(first["true_quadratic"]) * (x_parent*x_parent - 1.0/3.0)
        )
        retained = np.random.default_rng(int(first["thinning_seed"])).random(len(parent_z)) < p_inclusion
        observed_z = parent_z[retained]
        hashes = {
            "parent_frame_sha256": hashlib.sha256(np.ascontiguousarray(parent_z).view(np.uint8)).hexdigest(),
            "catalogue_sha256": hashlib.sha256(np.ascontiguousarray(observed_z).view(np.uint8)).hexdigest(),
            "gw_distance_sha256": hashlib.sha256(np.ascontiguousarray(distances).view(np.uint8)).hexdigest(),
        }
        for col, digest in hashes.items():
            if str(first[col]) != digest:
                raise AssertionError(f"independent regeneration does not match {col}")
        has_samples = group["validation_indices"].fillna("").astype(str).str.len().gt(0).any()
        if has_samples:
            rng = np.random.default_rng(int(seed_values[0]))
            for draw in range(draws):
                permutation = rng.permutation(int(first["parent_size"]))
                for size in sizes:
                    selected = group.loc[
                        group["draw"].eq(draw) & group["validation_size"].eq(size)
                    ]
                    if len(selected) != 1:
                        raise AssertionError("missing or duplicate scenario/draw/size record")
                    row = selected.iloc[0]
                    encoded = "" if pd.isna(row["validation_indices"]) else str(row["validation_indices"])
                    if not encoded:
                        if row["status"] == "success":
                            raise AssertionError("successful attempt is missing sample indices")
                        continue
                    indices = np.fromstring(encoded, sep=",", dtype=int)
                    expected_indices = permutation[:size]
                    if indices.size != size or np.unique(indices).size != size:
                        raise AssertionError("stored validation sample is not exact SRSWOR")
                    if np.any((indices < 0) | (indices >= int(row["parent_size"]))):
                        raise AssertionError("validation index outside parent frame")
                    if not np.array_equal(indices, expected_indices):
                        raise AssertionError("nested validation sample does not replay from declared seed")
                    if int(row["n_validation_included"]) != int(np.sum(retained[indices])):
                        raise AssertionError("validation inclusion count does not match regenerated parent frame")
                    if row["status"] == "success":
                        if np.unique(retained[indices]).size != 2:
                            raise AssertionError("successful logistic fit lacks both validation classes")
                        rebuilt = _independent_density_and_fit(
                            z_grid, observed_z, parent_z[indices], retained[indices].astype(int), rate
                        )
                        density, intercept, slope, ess, max_weight, p99, boundary = rebuilt
                        _close(row["fit_status"] == "boundary_limited", boundary, "fit status")
                        _close(row["fit_boundary_limited"], boundary, "fit boundary flag")
                        _close(row["weight_ess"], ess, "inverse-weight ESS")
                        _close(row["max_inverse_weight"], max_weight, "maximum inverse weight")
                        _close(row["p99_inverse_weight"], p99, "99th percentile inverse weight")
                        x_all = (parent_z-0.205)/0.195
                        p_true = expit(float(row["true_intercept"]) + float(row["true_slope"])*x_all +
                                       float(row["true_quadratic"])*(x_all*x_all-1.0/3.0))
                        p_fit = expit(intercept+slope*x_all)
                        _close(row["inclusion_probability_mae"], np.mean(np.abs(p_fit-p_true)),
                               "inclusion-probability MAE")
                        density_ise = trapezoid((density-target)**2,z_grid)
                        _close(row["density_ise"],density_ise,"host-density ISE")
                        post = _independent_posterior(
                            z_grid,density,h0_grid,distances,float(row["sigma_mpc"]),
                            float(row["detection_threshold_mpc"])
                        )
                        for column,value,label in zip(
                            ("posterior_mean","posterior_median","ci68_low","ci68_high"),
                            post,("posterior mean","posterior median","68% lower endpoint","68% upper endpoint")
                        ):
                            _close(row[column],value,label)
                        _close(row["ci68_width"],post[3]-post[2],"68% interval width")
                        _close(row["absolute_error"],abs(post[0]-float(row["h0_true"])),"absolute posterior-mean error")
                        if bool(row["covered_68"]) != (post[2] <= float(row["h0_true"]) <= post[3]):
                            raise AssertionError("68% coverage indicator fails independent posterior reconstruction")
    # Re-derive every reported validation-size summary from audited raw records.
    successful_records = records.loc[records["status"].eq("success")]
    reported = summary.get("metrics_by_validation_size", {})
    for size in sizes:
        group = successful_records.loc[successful_records["validation_size"].eq(size)]
        if not len(group):
            if str(size) in reported:
                raise AssertionError("summary reports metrics for a validation size with no successes")
            continue
        if str(size) not in reported:
            raise AssertionError("summary omits a validation-size metric block")
        metric = reported[str(size)]
        errors = group["posterior_mean"].to_numpy(dtype=float)-group["h0_true"].to_numpy(dtype=float)
        derived = {
            "posterior_mean_bias": float(np.mean(errors)),
            "posterior_mean_rmse": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": float(group["covered_68"].astype(float).mean()),
            "mean_ci68_width": float(group["ci68_width"].mean()),
            "mean_density_ise": float(group["density_ise"].mean()),
            "mean_inclusion_probability_mae": float(group["inclusion_probability_mae"].mean()),
            "mean_weight_ess": float(group["weight_ess"].mean()),
            "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
            "boundary_limited_fraction": float(group["fit_boundary_limited"].astype(float).mean()),
            "n_scenarios_with_successes": int(group["scenario"].nunique()),
            "scenario_balanced_bias_mean": float(
                group.assign(error=errors).groupby("scenario")["error"].mean().mean()
            ),
            "scenario_balanced_rmse_mean": float(
                group.assign(squared_error=errors**2).groupby("scenario")["squared_error"].mean().pow(0.5).mean()
            ),
            "scenario_balanced_coverage_mean": float(
                group.groupby("scenario")["covered_68"].mean().mean()
            ),
        }
        if group["scenario"].nunique() > 1:
            derived["scenario_coverage_sd"] = float(
                group.groupby("scenario")["covered_68"].mean().std(ddof=1)
            )
        else:
            derived["scenario_coverage_sd"] = None
        for name, value in derived.items():
            if value is None:
                if metric.get(name) is not None:
                    raise AssertionError(f"summary field {name} should be null for one scenario")
            else:
                _close(metric.get(name, np.nan), value, f"summary {size}/{name}")

    return {
        "status": "STRUCTURAL_AUDIT_PASS",
        "n_attempted": attempted,
        "n_successful": successful,
        "n_failed": failed,
        "n_scenarios": n_scenarios,
        "validation_sizes": sizes,
        "scope": "independent synthetic-input regeneration, SRSWOR bookkeeping, provenance, and seed replay; not estimator reconstruction",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    args = parser.parse_args()
    result = audit(args.artifact_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
