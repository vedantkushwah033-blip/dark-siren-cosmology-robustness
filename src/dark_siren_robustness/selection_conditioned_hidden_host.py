"""Selection-conditioned hidden-host mixture likelihood for a synthetic toy.

The missing-host fraction f is defined before event selection. This module
conditions on a fixed number of accepted observations in [d_min, d_max]; it
does not model event-rate information or a realistic detector response.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _trap_weights(grid: np.ndarray) -> np.ndarray:
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    if grid.size > 2:
        weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _validate_grid(values, name: str, *, positive: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size < 2 or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be a finite 1D array with at least two values")
    if np.any(np.diff(arr) <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    if positive and np.any(arr <= 0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


def _validate_density(values, grid: np.ndarray, name: str) -> np.ndarray:
    density = np.asarray(values, dtype=float)
    if density.shape != grid.shape or not np.all(np.isfinite(density)):
        raise ValueError(f"{name} must be finite and match redshift_grid")
    if np.any(density < 0):
        raise ValueError(f"{name} cannot contain negative values")
    integral = float(trapezoid(density, grid))
    if not np.isfinite(integral) or not np.isclose(integral, 1.0, rtol=1e-6, atol=1e-8):
        raise ValueError(f"{name} must integrate to one on redshift_grid; got {integral}")
    return density


def selection_conditioned_hidden_host_log_likelihood(
    redshift_grid,
    catalogue_density,
    missing_density,
    H0_grid,
    distances,
    *,
    sigma: float,
    f_missing: float,
    d_min: float = 0.0,
    d_max: float,
) -> dict[str, np.ndarray]:
    """Evaluate the fixed-count selection-conditioned mixture log likelihood.

    Parameters
    ----------
    redshift_grid:
        Strictly increasing redshift grid.
    catalogue_density, missing_density:
        Normalized redshift densities for catalogue and missing-host components.
    H0_grid:
        Strictly increasing positive trial H0 values in km/s/Mpc.
    distances:
        Observed distances, all of which must lie inside [d_min, d_max].
    sigma:
        Fixed Gaussian distance uncertainty in Mpc.
    f_missing:
        Missing-host mixture fraction *before event selection*, in [0, 1].
    d_min, d_max:
        Inclusive observed-distance acceptance region in Mpc.

    Returns
    -------
    Dictionary with the total log likelihood over H0 and component/mixture
    detection efficiencies. The mixture normalization is
    beta(H) = (1-f)*alpha_catalogue(H) + f*alpha_missing(H).
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
    if not np.isfinite(d_min) or not np.isfinite(d_max) or d_min >= d_max:
        raise ValueError("distance bounds must be finite and satisfy d_min < d_max")
    if np.any(d < d_min) or np.any(d > d_max):
        raise ValueError("every distance must lie in the accepted interval [d_min, d_max]")

    means = C_KM_S * z[None, :] / h0[:, None]
    acceptance = norm.cdf((d_max - means) / sigma) - norm.cdf((d_min - means) / sigma)
    alpha_cat = np.asarray(trapezoid(acceptance * g_cat[None, :], z, axis=1), dtype=float)
    alpha_miss = np.asarray(trapezoid(acceptance * g_miss[None, :], z, axis=1), dtype=float)
    if (not np.all(np.isfinite(alpha_cat)) or not np.all(np.isfinite(alpha_miss))
            or np.any(alpha_cat <= 0) or np.any(alpha_miss <= 0)):
        raise ArithmeticError("component selection efficiency is non-finite or non-positive")

    quadrature_measure = _trap_weights(z)
    with np.errstate(divide="ignore"):
        log_cat_measure = np.log(quadrature_measure) + np.log(g_cat)
        log_miss_measure = np.log(quadrature_measure) + np.log(g_miss)
        log_cat_weight = np.log1p(-f_missing)
        log_miss_weight = np.log(f_missing)

    log_q_cat = np.empty((d.size, h0.size), dtype=float)
    log_q_miss = np.empty_like(log_q_cat)
    for i, observed_distance in enumerate(d):
        log_pdf = norm.logpdf(observed_distance, loc=means, scale=sigma)
        log_q_cat[i] = logsumexp(log_pdf + log_cat_measure[None, :], axis=1)
        log_q_miss[i] = logsumexp(log_pdf + log_miss_measure[None, :], axis=1)

    log_numerator = np.logaddexp(
        log_cat_weight + log_q_cat,
        log_miss_weight + log_q_miss,
    )
    beta = (1.0 - f_missing) * alpha_cat + f_missing * alpha_miss
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0):
        raise ArithmeticError("mixture selection efficiency is non-finite or non-positive")
    log_likelihood = np.sum(log_numerator, axis=0) - d.size * np.log(beta)
    if not np.all(np.isfinite(log_likelihood)):
        raise ArithmeticError("selection-conditioned log likelihood is non-finite")

    return {
        "log_likelihood": log_likelihood,
        "alpha_catalogue": alpha_cat,
        "alpha_missing": alpha_miss,
        "beta": beta,
    }
