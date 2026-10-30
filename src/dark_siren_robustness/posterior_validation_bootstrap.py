"""Posterior-level uncertainty propagation for finite-validation host weighting.

Exploratory synthetic toy pipeline: holds the observed catalogue and GW-distance
sample fixed, refits the inclusion model on every resample, rebuilds weighted
host density, and recomputes the selection-aware H0 posterior. Not a real GW
likelihood and not a proof of bootstrap validity.
"""
from __future__ import annotations

import json
import numpy as np
from scipy.integrate import trapezoid
from numpy.typing import ArrayLike

from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation
from dark_siren_robustness.design_aware_srswor_bootstrap import design_aware_srswor_bootstrap_indices
from dark_siren_robustness.finite_parent_uncertainty import (
    bootstrap_standard_error, finite_population_correction, paired_bootstrap_indices,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def _posterior_from_validation(z, h0, observed, validation_z, validation_y, distances,
                               *, sigma_mpc, detection_threshold_mpc, merger_weight,
                               bandwidth_factor):
    density_result = estimate_host_density_from_validation(
        z, observed, validation_z, validation_y,
        merger_weight=merger_weight, bandwidth_factor=bandwidth_factor,
    )
    likelihood = selection_aware_toy_likelihood(
        z, density_result["density"], h0, distances,
        sigma=sigma_mpc, detection_threshold=detection_threshold_mpc,
    )
    return density_result, posterior_summary(h0, likelihood["selected_log_likelihood"])


def bootstrap_validation_h0_posterior(
    redshift_grid: ArrayLike,
    h0_grid: ArrayLike,
    observed_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_included: ArrayLike,
    observed_distances_mpc: ArrayLike,
    *,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    parent_population_size: int,
    method: str = "paired_bootstrap",
    n_bootstrap: int = 200,
    seed: int = 20261030,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> dict[str, object]:
    """Propagate validation resampling through the full toy H0 posterior.

    paired_bootstrap resamples paired records with replacement. The
    design_aware_srswor method generates a multinomial pseudo-population of
    size N and then draws an SRSWOR sample of size n. For this latter method,
    each replicate's pseudo-population posterior is retained as its reference.

    Intervals summarize resampling variability of posterior means, not H0
    posterior credible intervals. The scalar FPC adjustment is diagnostic
    only; no validity theorem is claimed for this nonlinear pipeline.
    """
    if method not in {"paired_bootstrap", "design_aware_srswor"}:
        raise ValueError("method must be paired_bootstrap or design_aware_srswor")
    for name, value, minimum in (
        ("n_bootstrap", n_bootstrap, 2),
        ("parent_population_size", parent_population_size, 10),
    ):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
        if value < minimum:
            raise ValueError(f"{name} must be >= {minimum}")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")

    z = np.asarray(redshift_grid, dtype=float)
    h0 = np.asarray(h0_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    vz = np.asarray(validation_redshift, dtype=float)
    vy_raw = np.asarray(validation_included)
    distances = np.asarray(observed_distances_mpc, dtype=float)
    if z.ndim != 1 or z.size < 21 or not np.all(np.isfinite(z)) or np.any(z <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError("redshift_grid must be finite, positive, and strictly increasing")
    if h0.ndim != 1 or h0.size < 31 or not np.all(np.isfinite(h0)) or np.any(h0 <= 0) or np.any(np.diff(h0) <= 0):
        raise ValueError("h0_grid must be finite, positive, and strictly increasing")
    if observed.ndim != 1 or observed.size < 5 or not np.all(np.isfinite(observed)) or np.any((observed < z[0]) | (observed > z[-1])) or np.unique(observed).size < 3:
        raise ValueError("observed_redshift must contain at least five distinct in-grid values")
    if vz.ndim != 1 or vz.size < 10 or not np.all(np.isfinite(vz)) or np.any((vz < z[0]) | (vz > z[-1])):
        raise ValueError("validation_redshift must contain at least ten finite in-grid values")
    if vy_raw.shape != vz.shape or vy_raw.ndim != 1 or not np.all(np.isin(vy_raw, [0, 1, False, True])):
        raise ValueError("validation_included must be a matching binary vector")
    vy = vy_raw.astype(int)
    if np.unique(vy).size != 2:
        raise ValueError("original validation sample must contain both retention classes")
    if parent_population_size < vz.size:
        raise ValueError("parent_population_size must be >= validation sample size")
    if not np.isfinite(sigma_mpc) or sigma_mpc <= 0 or not np.isfinite(detection_threshold_mpc) or detection_threshold_mpc <= 0:
        raise ValueError("sigma_mpc and detection_threshold_mpc must be finite and positive")
    if distances.ndim != 1 or distances.size < 1 or not np.all(np.isfinite(distances)) or np.any(distances < detection_threshold_mpc):
        raise ValueError("observed distances must be finite and satisfy the detection threshold")
    rate = None if merger_weight is None else np.asarray(merger_weight, dtype=float)
    if rate is not None and (rate.shape != z.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0)):
        raise ValueError("merger_weight must be positive, finite, and match redshift_grid")
    if bandwidth_factor is not None and (not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0):
        raise ValueError("bandwidth_factor must be finite and positive")

    if method == "paired_bootstrap":
        _, full_reference = _posterior_from_validation(
            z, h0, observed, vz, vy, distances,
            sigma_mpc=float(sigma_mpc), detection_threshold_mpc=float(detection_threshold_mpc),
            merger_weight=rate, bandwidth_factor=bandwidth_factor,
        )
        fixed_reference_mean = float(full_reference["posterior_mean"])
    else:
        fixed_reference_mean = float("nan")

    rng = np.random.default_rng(int(seed))
    success = np.zeros(int(n_bootstrap), dtype=bool)
    attempt_means = np.full(int(n_bootstrap), np.nan)
    reference_means = np.full(int(n_bootstrap), np.nan)
    medians = np.full(int(n_bootstrap), np.nan)
    lows = np.full(int(n_bootstrap), np.nan)
    highs = np.full(int(n_bootstrap), np.nan)
    ess = np.full(int(n_bootstrap), np.nan)
    max_inverse = np.full(int(n_bootstrap), np.nan)
    boundary = np.zeros(int(n_bootstrap), dtype=bool)
    posterior_draws = []
    sample_indices = []
    pseudo_indices = []
    selected_positions = []
    failures = []

    for b in range(int(n_bootstrap)):
        try:
            if method == "paired_bootstrap":
                idx = paired_bootstrap_indices(rng, vz.size)
                sample_indices.append(idx.copy())
                sample_z, sample_y = vz[idx], vy[idx]
                reference_mean = fixed_reference_mean
            else:
                pseudo_idx, sample_idx, positions = design_aware_srswor_bootstrap_indices(
                    rng, vz.size, int(parent_population_size), vz.size
                )
                pseudo_indices.append(pseudo_idx.copy())
                sample_indices.append(sample_idx.copy())
                selected_positions.append(positions.copy())
                pseudo_z, pseudo_y = vz[pseudo_idx], vy[pseudo_idx]
                sample_z, sample_y = vz[sample_idx], vy[sample_idx]
                _, pseudo_summary = _posterior_from_validation(
                    z, h0, observed, pseudo_z, pseudo_y, distances,
                    sigma_mpc=float(sigma_mpc), detection_threshold_mpc=float(detection_threshold_mpc),
                    merger_weight=rate, bandwidth_factor=bandwidth_factor,
                )
                reference_mean = float(pseudo_summary["posterior_mean"])

            density_result, summary = _posterior_from_validation(
                z, h0, observed, sample_z, sample_y, distances,
                sigma_mpc=float(sigma_mpc), detection_threshold_mpc=float(detection_threshold_mpc),
                merger_weight=rate, bandwidth_factor=bandwidth_factor,
            )
            posterior_density = np.asarray(
                json.loads(str(summary["posterior_density_grid_json"])), dtype=float
            )
            if (
                posterior_density.shape != h0.shape
                or not np.all(np.isfinite(posterior_density))
                or np.any(posterior_density < 0)
            ):
                raise ArithmeticError("posterior density grid is malformed")
            if not np.isclose(trapezoid(posterior_density, h0), 1.0, rtol=1e-8, atol=1e-8):
                raise ArithmeticError("posterior density grid is not normalized")
            posterior_draws.append(posterior_density)
            attempt_means[b] = float(summary["posterior_mean"])
            reference_means[b] = reference_mean
            medians[b] = float(summary["posterior_median"])
            lows[b] = float(summary["ci68_low"])
            highs[b] = float(summary["ci68_high"])
            ess[b] = float(density_result["inverse_weight_effective_sample_size"])
            max_inverse[b] = float(density_result["inverse_weight_maximum"])
            boundary[b] = bool(density_result["model"]["coefficient_bound_hit"])
            success[b] = True
        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError) as exc:
            failures.append({
                "bootstrap_index": int(b),
                "exception_type": type(exc).__name__,
                "reason": f"{type(exc).__name__}: {exc}",
            })
            if method == "paired_bootstrap" and len(sample_indices) <= b:
                sample_indices.append(np.array([], dtype=np.int64))

    successful_means = attempt_means[success]
    fpc = finite_population_correction(vz.size, int(parent_population_size))
    if successful_means.size >= 2:
        raw_se = bootstrap_standard_error(successful_means)
        adjusted_se = raw_se * fpc
        percentile = np.quantile(successful_means, [0.025, 0.975]).tolist()
    else:
        raw_se = adjusted_se = float("nan")
        percentile = [float("nan"), float("nan")]

    return {
        "schema": "validation_h0_posterior_resampling_v1",
        "method": method,
        "seed": int(seed),
        "n_attempted": int(n_bootstrap),
        "n_successful": int(np.sum(success)),
        "n_failed": int(np.sum(~success)),
        "success_mask": success,
        "attempt_posterior_mean": attempt_means,
        "posterior_mean_draws": successful_means,
        "posterior_median_draws_by_attempt": medians,
        "ci68_low_draws_by_attempt": lows,
        "ci68_high_draws_by_attempt": highs,
        "reference_posterior_mean_by_attempt": reference_means,
        "fixed_full_validation_reference_mean": fixed_reference_mean,
        "posterior_density_draws": (
            np.stack(posterior_draws) if posterior_draws else np.empty((0, h0.size))
        ),
        "sample_indices": sample_indices,
        "pseudo_population_indices": pseudo_indices,
        "selected_pseudo_population_positions": selected_positions,
        "inverse_weight_ess_by_attempt": ess,
        "inverse_weight_maximum_by_attempt": max_inverse,
        "coefficient_boundary_limited_by_attempt": boundary,
        "failure_records": failures,
        "raw_bootstrap_standard_error": raw_se,
        "fpc_factor": fpc,
        "fpc_adjusted_standard_error_diagnostic_only": adjusted_se,
        "posterior_mean_percentile_interval_95": percentile,
        "uncertainty_interpretation": (
            "Intervals summarize resampling variability of posterior means, not H0 posterior credible intervals. "
            "FPC-adjusted SE is a diagnostic comparator only; no validity theorem is claimed for this nonlinear "
            "logistic/KDE/selection-aware posterior pipeline. The observed catalogue and GW distances are fixed."
        ),
    }
