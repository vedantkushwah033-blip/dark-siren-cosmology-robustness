"""Bounded KDE host-density estimators for controlled selection experiments.

Inverse-inclusion weighting assumes supplied probabilities are known. It does
not estimate the selection function and is not a real-survey correction.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde


def estimate_host_density(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    *,
    inclusion_probability: ArrayLike | None = None,
    sample_weights: ArrayLike | None = None,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
    scott_multiplier: float | None = None,
    absolute_bandwidth_sd: float | None = None,
) -> NDArray[np.float64]:
    """Estimate a normalized host-redshift density on a bounded grid.

    For non-oracle analyses, observed_redshift must come from the observable
    catalogue, never hidden parent truth. Probabilities must be strictly
    positive: zero-probability regions cannot be recovered by inverse weights.
    Reflection reduces, but does not eliminate, KDE boundary bias. The final
    density is truncated and renormalized on the specified finite domain.

    ``sample_weights`` supplies positive galaxy-level base weights; when
    ``inclusion_probability`` is also supplied, the raw weights are
    ``sample_weights / inclusion_probability``. This allows one common
    base-population factor to be propagated through all comparison arms.

    Bandwidth semantics are intentionally explicit:
    - ``bandwidth_factor`` is an absolute SciPy covariance factor, not a
      multiplier of Scott's rule.
    - ``scott_multiplier`` multiplies the weighted Scott bandwidth computed
      from the original (unreflected) sample, including weighted standard
      deviation and effective sample size.
    - absolute_bandwidth_sd sets a common physical Gaussian-kernel SD in
      redshift units, independent of each arm's effective sample size and spread.
    - If all three are omitted, the weighted Scott factor is used.
    - The three explicit options are mutually exclusive.
    """
    grid = np.asarray(redshift_grid, dtype=float)
    sample = np.asarray(observed_redshift, dtype=float)
    if grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid)):
        raise ValueError("redshift_grid must be a finite 1D array with at least 21 points")
    if np.any(np.diff(grid) <= 0) or grid[0] <= 0:
        raise ValueError("redshift_grid must be strictly increasing and positive")
    if sample.ndim != 1 or sample.size < 5 or not np.all(np.isfinite(sample)):
        raise ValueError("observed_redshift must contain at least five finite values")
    if np.any((sample < grid[0]) | (sample > grid[-1])):
        raise ValueError("all sample redshifts must lie inside redshift_grid bounds")
    if np.unique(sample).size < 3:
        raise ValueError("observed_redshift must contain at least three distinct values")

    if inclusion_probability is None:
        weights = np.ones(sample.size, dtype=float)
    else:
        probability = np.asarray(inclusion_probability, dtype=float)
        if probability.shape != sample.shape or not np.all(np.isfinite(probability)):
            raise ValueError("inclusion_probability must be finite and match observed_redshift")
        if np.any((probability <= 0) | (probability > 1)):
            raise ValueError("inclusion_probability must lie in (0, 1]")
        with np.errstate(over="ignore", divide="ignore"):
            weights = 1.0 / probability
        if not np.all(np.isfinite(weights)):
            raise ValueError("inverse inclusion weights overflow; positivity is numerically unusable")
    if sample_weights is not None:
        supplied_weights = np.asarray(sample_weights, dtype=float)
        if (
            supplied_weights.shape != sample.shape
            or not np.all(np.isfinite(supplied_weights))
            or np.any(supplied_weights <= 0)
        ):
            raise ValueError("sample_weights must be finite, strictly positive, and match observed_redshift")
        with np.errstate(over="ignore", invalid="ignore"):
            weights = weights * supplied_weights
        if not np.all(np.isfinite(weights)):
            raise ValueError("combined sample weights overflow; positivity is numerically unusable")
    explicit_bandwidths = sum(value is not None for value in (bandwidth_factor, scott_multiplier, absolute_bandwidth_sd))
    if explicit_bandwidths > 1:
        raise ValueError("bandwidth_factor, scott_multiplier, and absolute_bandwidth_sd are mutually exclusive")
    if bandwidth_factor is not None and (
        not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0
    ):
        raise ValueError("bandwidth_factor must be finite and positive when supplied")
    if scott_multiplier is not None and (
        not np.isfinite(scott_multiplier) or scott_multiplier <= 0
    ):
        raise ValueError("scott_multiplier must be finite and positive when supplied")
    if absolute_bandwidth_sd is not None and (
        not np.isfinite(absolute_bandwidth_sd) or absolute_bandwidth_sd <= 0
    ):
        raise ValueError("absolute_bandwidth_sd must be finite and positive when supplied")

    left = 2.0 * grid[0] - sample
    right = 2.0 * grid[-1] - sample
    augmented_sample = np.concatenate((sample, left, right))
    # A common scale does not change normalized weighted KDEs. Scale before
    # squaring so large but representable inverse probabilities do not overflow.
    scaled_weights = weights / np.max(weights)
    if np.any(scaled_weights == 0.0):
        raise ArithmeticError(
            "weight scaling underflowed positive mass; refusing to hide support loss"
        )
    augmented_weights = np.tile(scaled_weights / 3.0, 3)
    if absolute_bandwidth_sd is not None:
        # gaussian_kde scales the augmented weighted covariance by its scalar
        # covariance factor. Set that factor so the physical kernel SD in the
        # original redshift coordinate equals the requested common width.
        augmented_probability = augmented_weights / np.sum(augmented_weights)
        augmented_mean = float(augmented_probability @ augmented_sample)
        augmented_denominator = 1.0 - float(augmented_probability @ augmented_probability)
        augmented_variance = float(np.sum(augmented_probability * (augmented_sample - augmented_mean) ** 2) / augmented_denominator)
        if not np.isfinite(augmented_variance) or augmented_variance <= 0:
            raise ValueError("weighted reflected KDE covariance is undefined")
        bw_method = float(absolute_bandwidth_sd) / np.sqrt(augmented_variance)
    elif bandwidth_factor is None:
        # Splitting each original weight over its three reflected copies changes
        # scipy's augmented-sample neff by a factor of three. Preserve the
        # original weighted-sample Scott factor instead of shrinking bandwidth
        # merely because each observation was reflected twice.
        original_neff = float(
            np.sum(scaled_weights) ** 2 / np.sum(scaled_weights ** 2)
        )
        weighted_scott_factor = original_neff ** (-1.0 / 5.0)
        multiplier = 1.0 if scott_multiplier is None else float(scott_multiplier)

        # gaussian_kde estimates its data covariance from the augmented
        # reflected sample. Correct the covariance factor so the physical
        # kernel standard deviation is based on the original weighted sample:
        # h = multiplier * n_eff**(-1/5) * sd_original.
        original_probability = scaled_weights / np.sum(scaled_weights)
        original_mean = float(original_probability @ sample)
        original_denominator = 1.0 - float(original_probability @ original_probability)
        original_variance = float(
            np.sum(original_probability * (sample - original_mean) ** 2)
            / original_denominator
        )
        augmented_probability = augmented_weights / np.sum(augmented_weights)
        augmented_mean = float(augmented_probability @ augmented_sample)
        augmented_denominator = 1.0 - float(
            augmented_probability @ augmented_probability
        )
        augmented_variance = float(
            np.sum(
                augmented_probability * (augmented_sample - augmented_mean) ** 2
            ) / augmented_denominator
        )
        if (
            not np.isfinite(original_variance) or original_variance <= 0
            or not np.isfinite(augmented_variance) or augmented_variance <= 0
        ):
            raise ValueError("weighted KDE covariance is undefined")
        covariance_adjustment = np.sqrt(original_variance / augmented_variance)
        bw_method = lambda _kde: (
            multiplier * weighted_scott_factor * covariance_adjustment
        )
    else:
        # Preserve the legacy absolute-factor behavior for existing callers.
        bw_method = bandwidth_factor
    kde = gaussian_kde(
        augmented_sample, weights=augmented_weights, bw_method=bw_method
    )
    galaxy_density = np.asarray(kde(grid), dtype=float)

    if merger_weight is None:
        rate = np.ones_like(grid)
    else:
        rate = np.asarray(merger_weight, dtype=float)
        if rate.shape != grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
            raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    unnormalized = galaxy_density * rate
    area = float(trapezoid(unnormalized, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("Estimated host density has invalid normalization")
    density = unnormalized / area
    if not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ArithmeticError("Estimated host density contains invalid values")
    if not np.isclose(trapezoid(density, grid), 1.0, rtol=1e-10, atol=1e-10):
        raise ArithmeticError("Estimated host density failed its unit-integral check")
    return density


def reflected_weighted_kde_bandwidth(
    observed_redshift: ArrayLike,
    *,
    lower_bound: float,
    upper_bound: float,
    sample_weights: ArrayLike,
    scott_multiplier: float | None = None,
    absolute_bandwidth_sd: float | None = None,
) -> float:
    if (scott_multiplier is None) == (absolute_bandwidth_sd is None):
        raise ValueError("supply exactly one of scott_multiplier or absolute_bandwidth_sd")
    if absolute_bandwidth_sd is not None:
        if not np.isfinite(absolute_bandwidth_sd) or absolute_bandwidth_sd <= 0:
            raise ValueError("absolute_bandwidth_sd must be finite and positive")
        return float(absolute_bandwidth_sd)
    """Return the physical redshift-space kernel SD used by the reflected KDE.

    The covariance adjustment in :func:`estimate_host_density` makes the
    physical bandwidth equal to the original-sample weighted SD multiplied
    by the weighted Scott factor. This reports that quantity directly.
    """
    sample = np.asarray(observed_redshift, dtype=float)
    weights = np.asarray(sample_weights, dtype=float)
    if (sample.ndim != 1 or sample.size < 5 or not np.all(np.isfinite(sample))
            or np.unique(sample).size < 3):
        raise ValueError("observed_redshift must contain at least five finite values and three distinct values")
    if (weights.shape != sample.shape or not np.all(np.isfinite(weights))
            or np.any(weights <= 0)):
        raise ValueError("sample_weights must be finite, positive, and match observed_redshift")
    if (not np.isfinite(lower_bound) or not np.isfinite(upper_bound)
            or lower_bound <= 0 or upper_bound <= lower_bound
            or np.any((sample < lower_bound) | (sample > upper_bound))):
        raise ValueError("bounds must be finite, positive, ordered, and contain all samples")
    if scott_multiplier is not None and (not np.isfinite(scott_multiplier) or scott_multiplier <= 0):
        raise ValueError("scott_multiplier must be finite and positive")
    scaled = weights / np.max(weights)
    if np.any(scaled == 0):
        raise ArithmeticError("weight scaling underflowed positive mass")
    probability = scaled / np.sum(scaled)
    denominator = 1.0 - float(probability @ probability)
    mean = float(probability @ sample)
    variance = float(np.sum(probability * (sample - mean) ** 2) / denominator)
    neff = float(np.sum(scaled) ** 2 / np.sum(scaled ** 2))
    width = float(scott_multiplier) * neff ** (-1.0 / 5.0) * np.sqrt(variance)
    if not np.isfinite(width) or width <= 0:
        raise ArithmeticError("realized reflected-KDE bandwidth is invalid")
    return width
