"""Estimate catalogue inclusion probabilities from an external validation sample.

This module is intentionally separate from the existing oracle-weighted KDE.
The validation data must include reference redshifts and catalogue-retention
indicators sampled independently of retention, or under a documented sampling
design. The incomplete catalogue alone does not identify the inclusion law.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.special import expit, log_expit

from dark_siren_robustness.finite_parent_uncertainty import paired_bootstrap_indices
from dark_siren_robustness.design_aware_srswor_bootstrap import (\n    design_aware_srswor_bootstrap_indices,\n)


def fit_logistic_inclusion_model(
    validation_redshift: ArrayLike,
    included: ArrayLike,
    *,
    z_min: float,
    z_max: float,
) -> dict[str, float | str | bool]:
    """Fit logit P(included=1 | z) = intercept + slope * standardized_z.

    Standardization is fixed by the declared domain: standardized_z is
    (z - (z_min + z_max)/2) / ((z_max - z_min)/2). It is not estimated from
    the validation sample, so fitted coefficients are comparable across
    replicates. The caller must ensure the validation sampling design supports
    this conditional model and must not treat fitted probabilities as known.

    Raises ValueError for malformed data or a validation sample containing only
    one outcome class. Raises RuntimeError if bounded maximum likelihood fails.
    """
    z = np.asarray(validation_redshift, dtype=float)
    y_raw = np.asarray(included)
    if z.ndim != 1 or z.size < 10 or not np.all(np.isfinite(z)):
        raise ValueError("validation_redshift must contain at least 10 finite values")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or not z_min < z_max:
        raise ValueError("z_min and z_max must be finite with z_min < z_max")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("validation_redshift must lie inside the declared domain")
    if y_raw.shape != z.shape or y_raw.ndim != 1:
        raise ValueError("included must be a one-dimensional array matching validation_redshift")
    if not np.all(np.isin(y_raw, [0, 1, False, True])):
        raise ValueError("included must contain only binary 0/1 values")
    y = y_raw.astype(float)
    if np.unique(y).size != 2:
        raise ValueError("validation sample must contain both included and not-included outcomes")

    midpoint = (float(z_min) + float(z_max)) / 2.0
    half_range = (float(z_max) - float(z_min)) / 2.0
    x = (z - midpoint) / half_range

    def negative_log_likelihood(theta: NDArray[np.float64]) -> float:
        eta = theta[0] + theta[1] * x
        value = -np.sum(y * log_expit(eta) + (1.0 - y) * log_expit(-eta))
        return float(value) if np.isfinite(value) else float("inf")

    result = minimize(
        negative_log_likelihood,
        x0=np.array([0.0, 0.0], dtype=float),
        method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0)),
        options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success or not np.all(np.isfinite(result.x)) or not np.isfinite(result.fun):
        raise RuntimeError(f"logistic inclusion fit failed: {result.message}")
    coefficient_bound_hit = bool(
        np.isclose(abs(float(result.x[0])), 25.0, rtol=0.0, atol=1e-5)
        or np.isclose(abs(float(result.x[1])), 50.0, rtol=0.0, atol=1e-5)
    )
    return {
        "model": "logistic_inclusion_v1",
        "fit_status": "boundary_limited" if coefficient_bound_hit else "converged",
        "coefficient_bound_hit": coefficient_bound_hit,
        "intercept": float(result.x[0]),
        "slope": float(result.x[1]),
        "z_min": float(z_min),
        "z_max": float(z_max),
        "n_validation": int(z.size),
        "n_included": int(np.sum(y)),
        "n_not_included": int(y.size - np.sum(y)),
        "negative_log_likelihood": float(result.fun),
    }


def predict_inclusion_probability(
    redshift: ArrayLike,
    model: Mapping[str, object],
) -> NDArray[np.float64]:
    """Predict inclusion probabilities using a fitted model, without clipping.

    No implicit probability floor or weight cap is applied. If the caller
    chooses to stabilize weights, that rule must be explicit and reported as
    part of the estimator rather than silently embedded in this function.
    """
    if model.get("model") != "logistic_inclusion_v1":
        raise ValueError("model must be a logistic_inclusion_v1 result")
    try:
        intercept = float(model["intercept"])
        slope = float(model["slope"])
        z_min = float(model["z_min"])
        z_max = float(model["z_max"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("model is missing valid coefficients or domain bounds") from exc
    if not all(np.isfinite(value) for value in (intercept, slope, z_min, z_max)) or z_min >= z_max:
        raise ValueError("model coefficients and domain bounds must be finite and valid")
    # These are the fitting bounds in fit_logistic_inclusion_model. Enforcing
    # the same model contract here prevents externally supplied coefficients
    # from producing subnormal probabilities and overflowing inverse weights.
    if abs(intercept) > 25.0 or abs(slope) > 50.0:
        raise ValueError("model coefficients exceed the bounded logistic-fit domain")
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)):
        raise ValueError("redshift must be a finite one-dimensional array")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("redshift values must lie inside the fitted model domain")
    midpoint = (z_min + z_max) / 2.0
    half_range = (z_max - z_min) / 2.0
    return np.asarray(expit(intercept + slope * ((z - midpoint) / half_range)), dtype=float)



def estimate_host_density_from_validation(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_included: ArrayLike,
    *,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> dict[str, object]:
    """Estimate host density using inclusion probabilities learned from audit data.

    The validation sample must represent the parent-frame population under a
    known sampling design and contain both catalogue-retention outcomes. The
    returned density uses only the observed catalogue and the audit records;
    the true inclusion probabilities and hidden parent objects are never inputs.

    merger_weight, when supplied, is a positive rate multiplier defined on
    redshift_grid (not an object-level weight). This wrapper intentionally
    applies no inverse-inclusion weight clipping. Its diagnostics make
    potentially unstable inverse weights visible; scientific use should
    preregister any positivity guard or cap rather than tuning it to H0 results.
    """
    from dark_siren_robustness.catalogue_density import estimate_host_density

    grid = np.asarray(redshift_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    validation_z = np.asarray(validation_redshift, dtype=float)
    if grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid)):
        raise ValueError("redshift_grid must be a finite 1D array with at least 21 points")
    if np.any(np.diff(grid) <= 0) or grid[0] <= 0:
        raise ValueError("redshift_grid must be strictly increasing and positive")
    if observed.ndim != 1 or observed.size < 5 or not np.all(np.isfinite(observed)):
        raise ValueError("observed_redshift must contain at least five finite values")
    if np.any((observed < grid[0]) | (observed > grid[-1])):
        raise ValueError("observed_redshift must lie inside redshift_grid bounds")
    model = fit_logistic_inclusion_model(
        validation_z,
        validation_included,
        z_min=float(grid[0]),
        z_max=float(grid[-1]),
    )
    probability = predict_inclusion_probability(observed, model)
    inverse_weights = 1.0 / probability
    scaled = inverse_weights / np.max(inverse_weights)
    effective_sample_size = float(np.sum(scaled) ** 2 / np.sum(scaled**2))
    density = estimate_host_density(
        grid,
        observed,
        inclusion_probability=probability,
        merger_weight=merger_weight,
        bandwidth_factor=bandwidth_factor,
    )
    return {
        "density": density,
        "model": model,
        "observed_inclusion_probability": probability,
        "inverse_weight_effective_sample_size": effective_sample_size,
        "inverse_weight_maximum": float(np.max(inverse_weights)),
        "inverse_weight_p95": float(np.quantile(inverse_weights, 0.95)),
        "inverse_weight_p99": float(np.quantile(inverse_weights, 0.99)),
        "n_observed_catalogue": int(observed.size),
        "n_validation": int(model["n_validation"]),
    }



def bootstrap_validation_inclusion_density(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_included: ArrayLike,
    *,
    n_bootstrap: int = 200,
    seed: int = 20261030,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> dict[str, object]:
    """Bootstrap validation-sample uncertainty in fitted inclusion and host density.

    The observed catalogue is held fixed. Each draw resamples the paired
    (validation redshift, retention indicator) records with replacement,
    refits the logistic inclusion model, and rebuilds the weighted host
    density. Thus the returned distribution captures finite-validation-sample
    variability conditional on this catalogue, the assumed logistic model,
    the declared validation sampling design, and the KDE bandwidth rule.

    This is a frequentist nonparametric bootstrap diagnostic, not a Bayesian
    posterior over the selection function and not a full uncertainty
    propagation for a real GW analysis. Resamples with only one retention
    class or otherwise invalid fits are recorded as failures, never silently
    replaced. A fixed seed makes the resampling reproducible.
    """
    if (
        isinstance(n_bootstrap, (bool, np.bool_))
        or not isinstance(n_bootstrap, (int, np.integer))
        or n_bootstrap < 2
    ):
        raise ValueError("n_bootstrap must be an integer of at least 2")
    if (
        isinstance(seed, (bool, np.bool_))
        or not isinstance(seed, (int, np.integer))
        or seed < 0
    ):
        raise ValueError("seed must be a nonnegative integer")
    validation_z = np.asarray(validation_redshift, dtype=float)
    validation_y = np.asarray(validation_included)
    if (
        validation_z.ndim != 1
        or validation_z.size < 10
        or not np.all(np.isfinite(validation_z))
        or validation_y.shape != validation_z.shape
        or validation_y.ndim != 1
        or not np.all(np.isin(validation_y, [0, 1, False, True]))
    ):
        raise ValueError("validation records must be matching finite redshifts and binary outcomes")
    validation_y = validation_y.astype(int)
    if np.unique(validation_y).size != 2:
        raise ValueError("original validation sample must contain both retention classes")

    grid = np.asarray(redshift_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    if (
        grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid))
        or np.any(np.diff(grid) <= 0) or grid[0] <= 0
    ):
        raise ValueError("redshift_grid must be finite, positive, and strictly increasing")
    if (
        observed.ndim != 1 or observed.size < 5 or not np.all(np.isfinite(observed))
        or np.any((observed < grid[0]) | (observed > grid[-1]))
        or np.unique(observed).size < 3
    ):
        raise ValueError(
            "observed_redshift must be finite, inside grid bounds, and contain >=3 distinct values"
        )
    if np.any((validation_z < grid[0]) | (validation_z > grid[-1])):
        raise ValueError("validation redshifts must lie inside redshift_grid bounds")
    if merger_weight is not None:
        rate = np.asarray(merger_weight, dtype=float)
        if rate.shape != grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
            raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    if bandwidth_factor is not None and (
        not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0
    ):
        raise ValueError("bandwidth_factor must be finite and positive")
    rng = np.random.default_rng(int(seed))
    density_draws: list[NDArray[np.float64]] = []
    intercept_draws: list[float] = []
    slope_draws: list[float] = []
    boundary_flags: list[bool] = []
    ess_draws: list[float] = []
    max_weight_draws: list[float] = []
    p99_weight_draws: list[float] = []
    failure_reasons: list[str] = []
    resample_indices: list[NDArray[np.int64]] = []
    success_mask = np.zeros(int(n_bootstrap), dtype=bool)

    for bootstrap_index in range(int(n_bootstrap)):
        indices = paired_bootstrap_indices(rng, validation_z.size)
        resample_indices.append(indices.astype(np.int64, copy=True))
        try:
            result = estimate_host_density_from_validation(
                grid,
                observed_redshift,
                validation_z[indices],
                validation_y[indices],
                merger_weight=merger_weight,
                bandwidth_factor=bandwidth_factor,
            )
        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError) as exc:
            failure_reasons.append(f"{type(exc).__name__}: {exc}")
            continue
        success_mask[bootstrap_index] = True
        model = result["model"]
        density_draws.append(np.asarray(result["density"], dtype=float))
        intercept_draws.append(float(model["intercept"]))
        slope_draws.append(float(model["slope"]))
        boundary_flags.append(bool(model["coefficient_bound_hit"]))
        ess_draws.append(float(result["inverse_weight_effective_sample_size"]))
        max_weight_draws.append(float(result["inverse_weight_maximum"]))
        p99_weight_draws.append(float(result["inverse_weight_p99"]))

    density_array = (
        np.stack(density_draws, axis=0)
        if density_draws
        else np.empty((0, grid.size), dtype=float)
    )
    return {
        "schema": "validation_inclusion_density_bootstrap_v1",
        "seed": int(seed),
        "n_attempted": int(n_bootstrap),
        "n_successful": int(len(density_draws)),
        "n_failed": int(len(failure_reasons)),
        "success_fraction": float(len(density_draws) / n_bootstrap),
        "density_draws": density_array,
        "resample_indices": np.stack(resample_indices, axis=0),
        "success_mask": success_mask,
        "intercept_draws": np.asarray(intercept_draws, dtype=float),
        "slope_draws": np.asarray(slope_draws, dtype=float),
        "coefficient_boundary_limited_draws": np.asarray(boundary_flags, dtype=bool),
        "inverse_weight_ess_draws": np.asarray(ess_draws, dtype=float),
        "inverse_weight_maximum_draws": np.asarray(max_weight_draws, dtype=float),
        "inverse_weight_p99_draws": np.asarray(p99_weight_draws, dtype=float),
        "failure_reasons": failure_reasons,
        "bootstrap_method": "paired validation-record resampling with replacement",
        "finite_population_correction_applied": False,
        "uncertainty_scope": (
            "approximate iid validation-sample bootstrap; observed catalogue fixed; "
            "finite-population correction not applied; logistic selection family assumed; "
            "not Bayesian or full GW uncertainty"
        ),
    }


def bootstrap_validation_inclusion_density_design_aware(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_included: ArrayLike,
    *,
    parent_population_size: int,
    n_bootstrap: int = 200,
    seed: int = 20261030,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> dict[str, object]:
    """Exploratory pseudo-population bootstrap for SRSWOR validation records.

    Each replicate first generates a size-N pseudo-population by multinomial
    replication of the observed paired validation records, then draws a
    size-n SRSWOR sample from that pseudo-population. The same indices are
    applied to redshift and retention labels. This follows the SRSWOR
    pseudo-population construction of Wang, Peng & Kim (2022) for finite-
    population mean inference, but its theoretical guarantees do not
    automatically transfer to this nonlinear logistic/KDE estimator.

    The observed catalogue remains fixed. The returned distribution is an
    empirical comparator only; no design-consistency claim is made here.
    """
    if (
        isinstance(parent_population_size, (bool, np.bool_))
        or not isinstance(parent_population_size, (int, np.integer))
        or parent_population_size < 10
    ):
        raise ValueError("parent_population_size must be an integer of at least 10")
    if (
        isinstance(n_bootstrap, (bool, np.bool_))
        or not isinstance(n_bootstrap, (int, np.integer))
        or n_bootstrap < 2
    ):
        raise ValueError("n_bootstrap must be an integer of at least 2")
    if (
        isinstance(seed, (bool, np.bool_))
        or not isinstance(seed, (int, np.integer))
        or seed < 0
    ):
        raise ValueError("seed must be a nonnegative integer")
    validation_z = np.asarray(validation_redshift, dtype=float)
    validation_y = np.asarray(validation_included)
    if (
        validation_z.ndim != 1 or validation_z.size < 10
        or not np.all(np.isfinite(validation_z))
        or validation_y.shape != validation_z.shape or validation_y.ndim != 1
        or not np.all(np.isin(validation_y, [0, 1, False, True]))
    ):
        raise ValueError("validation records must be matching finite redshifts and binary outcomes")
    validation_y = validation_y.astype(int)
    if np.unique(validation_y).size != 2:
        raise ValueError("original validation sample must contain both retention classes")
    if parent_population_size < validation_z.size:
        raise ValueError("parent_population_size must be >= validation sample size")

    grid = np.asarray(redshift_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    if (
        grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid))
        or np.any(np.diff(grid) <= 0) or grid[0] <= 0
    ):
        raise ValueError("redshift_grid must be finite, positive, and strictly increasing")
    if (
        observed.ndim != 1 or observed.size < 5 or not np.all(np.isfinite(observed))
        or np.any((observed < grid[0]) | (observed > grid[-1]))
        or np.unique(observed).size < 3
    ):
        raise ValueError("observed_redshift must be inside grid and contain >=3 distinct values")
    if np.any((validation_z < grid[0]) | (validation_z > grid[-1])):
        raise ValueError("validation redshifts must lie inside redshift_grid bounds")
    if merger_weight is not None:
        rate = np.asarray(merger_weight, dtype=float)
        if rate.shape != grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
            raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    if bandwidth_factor is not None and (
        not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0
    ):
        raise ValueError("bandwidth_factor must be finite and positive")

    rng = np.random.default_rng(int(seed))
    density_draws: list[NDArray[np.float64]] = []
    intercept_draws: list[float] = []
    slope_draws: list[float] = []
    boundary_flags: list[bool] = []
    ess_draws: list[float] = []
    max_weight_draws: list[float] = []
    p99_weight_draws: list[float] = []
    pseudo_population_indices: list[NDArray[np.int64]] = []
    bootstrap_sample_indices: list[NDArray[np.int64]] = []
    failure_reasons: list[str] = []
    success_mask = np.zeros(int(n_bootstrap), dtype=bool)

    for bootstrap_index in range(int(n_bootstrap)):
        pseudo_idx, sample_idx = design_aware_srswor_bootstrap_indices(
            rng, validation_z.size, int(parent_population_size), validation_z.size
        )
        pseudo_population_indices.append(pseudo_idx)
        bootstrap_sample_indices.append(sample_idx)
        try:
            result = estimate_host_density_from_validation(
                grid, observed, validation_z[sample_idx], validation_y[sample_idx],
                merger_weight=merger_weight, bandwidth_factor=bandwidth_factor,
            )
        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError) as exc:
            failure_reasons.append(f"{type(exc).__name__}: {exc}")
            continue
        success_mask[bootstrap_index] = True
        model = result["model"]
        density_draws.append(np.asarray(result["density"], dtype=float))
        intercept_draws.append(float(model["intercept"]))
        slope_draws.append(float(model["slope"]))
        boundary_flags.append(bool(model["coefficient_bound_hit"]))
        ess_draws.append(float(result["inverse_weight_effective_sample_size"]))
        max_weight_draws.append(float(result["inverse_weight_maximum"]))
        p99_weight_draws.append(float(result["inverse_weight_p99"]))

    return {
        "schema": "validation_inclusion_density_pseudopop_srswor_bootstrap_v1",
        "seed": int(seed),
        "parent_population_size": int(parent_population_size),
        "observed_validation_size": int(validation_z.size),
        "n_attempted": int(n_bootstrap),
        "n_successful": int(len(density_draws)),
        "n_failed": int(len(failure_reasons)),
        "success_fraction": float(len(density_draws) / n_bootstrap),
        "density_draws": (
            np.stack(density_draws) if density_draws
            else np.empty((0, grid.size), dtype=float)
        ),
        "pseudo_population_indices": np.stack(pseudo_population_indices),
        "bootstrap_sample_indices": np.stack(bootstrap_sample_indices),
        "success_mask": success_mask,
        "intercept_draws": np.asarray(intercept_draws, dtype=float),
        "slope_draws": np.asarray(slope_draws, dtype=float),
        "coefficient_boundary_limited_draws": np.asarray(boundary_flags, dtype=bool),
        "inverse_weight_ess_draws": np.asarray(ess_draws, dtype=float),
        "inverse_weight_maximum_draws": np.asarray(max_weight_draws, dtype=float),
        "inverse_weight_p99_draws": np.asarray(p99_weight_draws, dtype=float),
        "failure_reasons": failure_reasons,
        "bootstrap_method": (
            "multinomial pseudo-population of size N followed by SRSWOR of size n"
        ),
        "finite_population_correction_applied": False,
        "uncertainty_scope": (
            "exploratory design-aware pseudo-population bootstrap for SRSWOR validation; "
            "observed catalogue fixed; nonlinear-estimator validity not established; "
            "not a Bayesian posterior or full GW uncertainty"
        ),
    }
