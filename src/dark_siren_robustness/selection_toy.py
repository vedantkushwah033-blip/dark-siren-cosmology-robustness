"""Numerically explicit low-redshift selection-aware toy likelihood.

This is a controlled synthetic validation primitive, not a GW detector model.
The host density is supplied as a normalized density over a redshift grid; the
same density and Gaussian noise scale define event likelihoods and detection
efficiency. The function conditions on a fixed number of detected events and
does not infer an event rate or catalogue incompleteness.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _grid(values, name, *, positive=False, min_size=2):
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size < min_size:
        raise ValueError(f"{name} must be a one-dimensional array with at least {min_size} values.")
    if not np.all(np.isfinite(result)) or np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be finite and strictly increasing.")
    if positive and np.any(result <= 0):
        raise ValueError(f"{name} must be strictly positive.")
    return result


def _trapezoid_weights(grid):
    """Return weights q such that sum(q * f) is trapezoidal integration of f."""
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    if grid.size > 2:
        weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def normalized_host_density(redshift_grid, galaxy_intensity, merger_weight=None):
    """Normalize lambda(z) * rho(z) on a strictly increasing redshift grid."""
    z = _grid(redshift_grid, "redshift_grid", positive=True)
    intensity = np.asarray(galaxy_intensity, dtype=float)
    rate = np.ones_like(z) if merger_weight is None else np.asarray(merger_weight, dtype=float)
    if intensity.shape != z.shape or rate.shape != z.shape:
        raise ValueError("galaxy_intensity and merger_weight must match redshift_grid.")
    if (not np.all(np.isfinite(intensity)) or not np.all(np.isfinite(rate))
            or np.any(intensity <= 0) or np.any(rate <= 0)):
        raise ValueError("galaxy_intensity and merger_weight must be finite and positive.")
    unnormalized = intensity * rate
    area = float(trapezoid(unnormalized, z))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Host-density normalization failed.")
    density = unnormalized / area
    if not np.isclose(trapezoid(density, z), 1.0, rtol=1e-12, atol=1e-12):
        raise ArithmeticError("Normalized host density failed its unit-integral check.")
    return density


def selection_aware_toy_likelihood(
    redshift_grid,
    host_density,
    H0_grid,
    distances,
    *,
    sigma,
    detection_threshold,
):
    """Compute unconditional and selection-conditioned log likelihoods.

    Model:
        d | z,H0 ~ Normal(c*z/H0, sigma**2)
        detected iff d >= detection_threshold.

    Returns arrays on H0_grid: unconditional joint log likelihood,
    detection efficiency beta(H0), and log likelihood conditional on detecting
    the fixed number of supplied events. Every supplied distance must satisfy
    the detection rule. No event-rate or realistic detector response is modeled.
    """
    z = _grid(redshift_grid, "redshift_grid", positive=True)
    h0 = _grid(H0_grid, "H0_grid", positive=True)
    density = np.asarray(host_density, dtype=float)
    d = np.asarray(distances, dtype=float)
    sigma = float(sigma)
    threshold = float(detection_threshold)
    if density.shape != z.shape or not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ValueError("host_density must be finite, nonnegative, and match redshift_grid.")
    if not np.isclose(trapezoid(density, z), 1.0, rtol=1e-7, atol=1e-9):
        raise ValueError("host_density must integrate to one on redshift_grid.")
    if d.ndim != 1 or d.size < 1 or not np.all(np.isfinite(d)) or np.any(d <= 0):
        raise ValueError("distances must be a nonempty one-dimensional array of positive finite values.")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be finite and positive.")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_threshold must be finite and positive.")
    if np.any(d < threshold):
        raise ValueError("Every supplied event must satisfy d >= detection_threshold.")
    quad = _trapezoid_weights(z)
    positive_mass = (density > 0) & (quad > 0)
    log_mass = np.log(density[positive_mass]) + np.log(quad[positive_mass])
    means = C_KM_S * z[positive_mass][None, :] / h0[:, None]
    # Shape: (H0, event, redshift); logsumexp prevents underflow in mixture sums.
    log_terms = norm.logpdf(
        d[None, :, None],
        loc=means[:, None, :],
        scale=sigma,
    ) + log_mass[None, None, :]
    event_log_likelihood = logsumexp(log_terms, axis=2)
    joint_log_likelihood = np.sum(event_log_likelihood, axis=1)
    detection_probability = norm.cdf(
        (C_KM_S * z[None, :] / h0[:, None] - threshold) / sigma
    )
    beta = np.sum(density[None, :] * quad[None, :] * detection_probability, axis=1)
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1.0 + 1e-12):
        raise ArithmeticError("Detection efficiency must be finite and in (0, 1].")
    selected = joint_log_likelihood - d.size * np.log(beta)
    return {
        "H0_grid": h0.copy(),
        "unconditional_log_likelihood": joint_log_likelihood,
        "detection_efficiency": beta,
        "selected_log_likelihood": selected,
        "n_detected": int(d.size),
    }
