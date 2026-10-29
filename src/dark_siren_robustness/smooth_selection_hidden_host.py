"""Selection-conditioned hidden-host likelihood with a smooth detection curve.

This is a synthetic methodological primitive, not a detector response model.
The detected-event likelihood conditions on a fixed number of accepted events.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import expit, logsumexp
from scipy.special import roots_hermitenorm
from scipy.stats import norm

from .selection_conditioned_hidden_host import (
    C_KM_S,
    _trap_weights,
    _validate_density,
    _validate_grid,
)


def smooth_selection_probability(distance, *, d50_mpc: float, width_mpc: float):
    """Logistic acceptance probability; d50 is 50% point and width > 0."""
    if not np.isfinite(d50_mpc) or not np.isfinite(width_mpc) or width_mpc <= 0:
        raise ValueError("d50_mpc must be finite and width_mpc finite and positive")
    return expit(-(np.asarray(distance, dtype=float) - d50_mpc) / width_mpc)


def _component_selection_efficiency(
    z: np.ndarray,
    density: np.ndarray,
    h0: np.ndarray,
    sigma: float,
    d50_mpc: float,
    width_mpc: float,
    quadrature_order: int,
) -> np.ndarray:
    # If X ~ N(0,1), E[s(mu + sigma X)] is evaluated by Gauss-Hermite
    # quadrature for the standard-normal measure.
    nodes, weights = roots_hermitenorm(quadrature_order)
    weights = weights / np.sqrt(2.0 * np.pi)
    means = C_KM_S * z[None, :] / h0[:, None]
    sampled_d = means[:, :, None] + sigma * nodes[None, None, :]
    selected = smooth_selection_probability(
        sampled_d, d50_mpc=d50_mpc, width_mpc=width_mpc
    )
    acceptance = np.sum(selected * weights[None, None, :], axis=2)
    return np.asarray(trapezoid(acceptance * density[None, :], z, axis=1), dtype=float)


def selection_conditioned_smooth_hidden_host_log_likelihood(
    redshift_grid,
    catalogue_density,
    missing_density,
    H0_grid,
    distances,
    *,
    sigma: float,
    f_missing: float,
    d50_mpc: float,
    width_mpc: float,
    quadrature_order: int = 40,
) -> dict[str, np.ndarray]:
    """Evaluate a smooth-selection-conditioned mixture likelihood.

    Component densities are normalized on redshift_grid. The missing-host
    fraction is defined before selection. Selection probability is logistic in
    observed distance, and the likelihood is conditioned on a fixed number of
    detected events. This function assumes the same smooth selection curve is
    used for data conditioning and inference; a mismatch experiment must store
    true and assumed curves separately and generate data under the true curve.
    """
    z = _validate_grid(redshift_grid, "redshift_grid", positive=True)
    h0 = _validate_grid(H0_grid, "H0_grid", positive=True)
    g_cat = _validate_density(catalogue_density, z, "catalogue_density")
    g_miss = _validate_density(missing_density, z, "missing_density")
    d = np.asarray(distances, dtype=float)
    if d.ndim != 1 or d.size < 1 or not np.all(np.isfinite(d)):
        raise ValueError("distances must be a nonempty finite 1D array")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be finite and positive")
    if not np.isfinite(f_missing) or not 0.0 <= f_missing <= 1.0:
        raise ValueError("f_missing must be in [0, 1]")
    if not np.isfinite(d50_mpc) or not np.isfinite(width_mpc) or width_mpc <= 0:
        raise ValueError("d50_mpc must be finite and width_mpc finite and positive")
    if not isinstance(quadrature_order, int) or quadrature_order < 8:
        raise ValueError("quadrature_order must be an integer >= 8")

    alpha_cat = _component_selection_efficiency(
        z, g_cat, h0, sigma, d50_mpc, width_mpc, quadrature_order
    )
    alpha_miss = _component_selection_efficiency(
        z, g_miss, h0, sigma, d50_mpc, width_mpc, quadrature_order
    )
    beta = (1.0 - f_missing) * alpha_cat + f_missing * alpha_miss
    if (not np.all(np.isfinite(beta)) or np.any(beta <= 0)
            or np.any(beta > 1.0 + 1e-10)):
        raise ArithmeticError("mixture selection efficiency is invalid")

    means = C_KM_S * z[None, :] / h0[:, None]
    measure = _trap_weights(z)
    with np.errstate(divide="ignore"):
        log_cat_measure = np.log(measure) + np.log(g_cat)
        log_miss_measure = np.log(measure) + np.log(g_miss)
        log_cat_weight = np.log1p(-f_missing)
        log_miss_weight = np.log(f_missing)

    log_q_cat = np.empty((d.size, h0.size), dtype=float)
    log_q_miss = np.empty_like(log_q_cat)
    for i, observed in enumerate(d):
        log_pdf = norm.logpdf(observed, loc=means, scale=sigma)
        log_q_cat[i] = logsumexp(log_pdf + log_cat_measure[None, :], axis=1)
        log_q_miss[i] = logsumexp(log_pdf + log_miss_measure[None, :], axis=1)

    log_mixture_density = np.logaddexp(
        log_cat_weight + log_q_cat,
        log_miss_weight + log_q_miss,
    )
    # The event-specific selection probability is independent of H0 and cancels
    # in posterior ratios for fixed selection parameters, but is retained here
    # to represent the full density conditional on acceptance.
    p_selected = smooth_selection_probability(
        d, d50_mpc=d50_mpc, width_mpc=width_mpc
    )
    if np.any(p_selected <= 0) or not np.all(np.isfinite(p_selected)):
        raise ArithmeticError("invalid event selection probability")
    log_likelihood = (
        np.sum(log_mixture_density + np.log(p_selected)[:, None], axis=0)
        - d.size * np.log(beta)
    )
    if not np.all(np.isfinite(log_likelihood)):
        raise ArithmeticError("smooth-selection log likelihood is non-finite")
    return {
        "log_likelihood": log_likelihood,
        "alpha_catalogue": alpha_cat,
        "alpha_missing": alpha_miss,
        "beta": beta,
    }
