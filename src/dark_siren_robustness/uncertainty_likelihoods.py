"""Explicit distance-uncertainty likelihoods for the synthetic known-host study.

These functions separate three assumptions that must not be conflated:
(1) a realized, fixed measurement uncertainty; (2) a candidate-H0-dependent
fractional uncertainty, including the Gaussian normalization; and (3) a
fixed absolute uncertainty calibrated at a reference H0.

This is a one-host low-redshift model, not a realistic GW distance posterior.
"""
from __future__ import annotations

import math
import numpy as np

C_KM_S = 299_792.458


def _positive_finite(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0):
        raise ValueError(f"{name} must be finite and strictly positive")
    return arr


def candidate_distance_mpc(h0, host_redshift):
    """Return d_L = c z / H0 (Mpc) for positive candidate H0 and redshift."""
    h = _positive_finite("h0", h0)
    z = _positive_finite("host_redshift", host_redshift)
    return C_KM_S * z / h


def gaussian_log_likelihood(observed_distance_mpc, mean_mpc, sigma_mpc):
    """Normalized Gaussian log likelihood, including -log(sigma)."""
    observed = np.asarray(observed_distance_mpc, dtype=float)
    mean = np.asarray(mean_mpc, dtype=float)
    sigma = _positive_finite("sigma_mpc", sigma_mpc)
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(mean)):
        raise ValueError("observed distance and mean must be finite")
    return -0.5 * ((observed - mean) / sigma) ** 2 - np.log(sigma) - 0.5 * math.log(2.0 * math.pi)


def log_likelihood_fixed_sigma(observed_distance_mpc, h0, host_redshift, sigma_mpc):
    """Arm A: hold the observed event's realized sigma fixed over candidate H0."""
    mean = candidate_distance_mpc(h0, host_redshift)
    return gaussian_log_likelihood(observed_distance_mpc, mean, sigma_mpc)


def log_likelihood_fractional_sigma(
    observed_distance_mpc, h0, host_redshift, fractional_error
):
    """Arm B: sigma(H0)=fractional_error*c*z/H0, including its normalization."""
    f = _positive_finite("fractional_error", fractional_error)
    mean = candidate_distance_mpc(h0, host_redshift)
    sigma = f * mean
    return gaussian_log_likelihood(observed_distance_mpc, mean, sigma)


def log_likelihood_fixed_reference_sigma(
    observed_distance_mpc, h0, host_redshift, fractional_error, reference_h0=70.0
):
    """Arm C: use sigma_ref=f*c*z/reference_h0, fixed over candidate H0."""
    f = _positive_finite("fractional_error", fractional_error)
    reference = _positive_finite("reference_h0", reference_h0)
    sigma_ref = f * candidate_distance_mpc(reference, host_redshift)
    return log_likelihood_fixed_sigma(
        observed_distance_mpc, h0, host_redshift, sigma_ref
    )
