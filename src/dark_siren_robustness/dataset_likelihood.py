"""Minimal dataset-level composition for the controlled synthetic event model.

This module supports conditionally independent events with a shared H0 and a
uniform prior on a declared grid. It does not implement selection conditioning,
event-rate inference, or shared nuisance parameters.
"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458

def _vector(value, name, *, min_size=1, nonnegative=False, increasing=False):
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or result.size < min_size:
        raise ValueError(f"{name} must be a one-dimensional array with at least {min_size} value(s).")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values.")
    if nonnegative and np.any(result < 0):
        raise ValueError(f"{name} must be non-negative.")
    if increasing and np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be strictly increasing.")
    return result

def _weights(value, size, name):
    values = np.ones(size, dtype=float) if value is None else _vector(value, name, nonnegative=True)
    if values.size != size or not np.any(values > 0):
        raise ValueError(f"{name} must match its support and have positive total mass.")
    values = values / np.max(values)
    return values / values.sum()

def event_log_likelihood(
    H0_grid, *, redshifts, gw_distance, gw_sigma, weights=None,
    missing_weight=0.0, missing_redshifts=None, missing_weights=None,
    positive_truncated=False,
):
    """Return a prior-free event log-likelihood vector over candidate H0."""
    grid = _vector(H0_grid, "H0_grid", min_size=2, increasing=True)
    if np.any(grid <= 0):
        raise ValueError("H0_grid must be positive.")
    z = _vector(redshifts, "redshifts", nonnegative=True)
    distance, sigma, missing_weight = float(gw_distance), float(gw_sigma), float(missing_weight)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("gw_distance must be finite and positive.")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("gw_sigma must be finite and positive.")
    if not np.isfinite(missing_weight) or not 0 <= missing_weight < 1:
        raise ValueError("missing_weight must be in [0, 1).")
    catalog_weights = _weights(weights, z.size, "weights") * (1.0 - missing_weight)

    def mixture_log_density(host_z, host_weights):
        means = C_KM_S * host_z[None, :] / grid[:, None]
        terms = norm.logpdf(distance, loc=means, scale=sigma)
        if positive_truncated:
            terms -= norm.logcdf(means / sigma)
        log_weights = np.full(host_weights.shape, -np.inf)
        positive = host_weights > 0
        log_weights[positive] = np.log(host_weights[positive])
        return logsumexp(terms + log_weights[None, :], axis=1)

    catalog = mixture_log_density(z, catalog_weights)
    if missing_weight == 0:
        return catalog
    if missing_redshifts is None:
        raise ValueError("missing_redshifts are required when missing_weight > 0.")
    mz = _vector(missing_redshifts, "missing_redshifts", nonnegative=True)
    mw = _weights(missing_weights, mz.size, "missing_weights")
    missing = mixture_log_density(mz, mw) + np.log(missing_weight)
    return np.logaddexp(catalog, missing)

def dataset_posterior_h0(H0_grid, events: Sequence[Mapping], *, prior_min=50.0, prior_max=90.0):
    """Combine independent event likelihoods and apply one uniform H0 prior.

    Each event mapping contains keyword arguments accepted by event_log_likelihood,
    except H0_grid. Selection effects, rates, and shared nuisance parameters are
    deliberately unsupported by this minimal prototype.
    """
    grid = _vector(H0_grid, "H0_grid", min_size=2, increasing=True)
    lo, hi = float(prior_min), float(prior_max)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("prior_min and prior_max must be finite with prior_min < prior_max.")
    if len(events) == 0:
        raise ValueError("At least one event is required.")
    if grid[0] > lo or grid[-1] < hi:
        raise ValueError(
            "H0_grid must span the full prior interval; extrapolation beyond the grid is unsupported."
        )
    support = (grid >= lo) & (grid <= hi)
    if np.count_nonzero(support) < 2:
        raise ValueError("H0_grid must contain at least two points inside the prior interval.")
    joint_log_likelihood = np.zeros(grid.size, dtype=float)
    for index, event in enumerate(events):
        if not isinstance(event, Mapping):
            raise ValueError(f"events[{index}] must be a mapping.")
        joint_log_likelihood += event_log_likelihood(grid, **dict(event))
    supported = joint_log_likelihood[support]
    finite = np.isfinite(supported)
    if not np.any(finite):
        raise ValueError("Joint posterior normalization failed: no finite likelihood on prior support.")

    # Integrate over the full prior interval even when lo/hi fall between grid
    # points. Interpolate log likelihood at inserted boundaries (piecewise-linear
    # in log likelihood); otherwise trapezoidal normalization silently truncates
    # the prior support to the first/last interior grid points.
    interior_grid = grid[support]
    interior_log_likelihood = joint_log_likelihood[support]
    norm_grid = np.unique(np.concatenate(([lo], interior_grid, [hi])))
    norm_log_likelihood = np.interp(norm_grid, interior_grid, interior_log_likelihood)
    finite_norm = np.isfinite(norm_log_likelihood)
    if not np.any(finite_norm):
        raise ValueError("Joint posterior normalization failed: no finite likelihood on prior support.")
    log_shift = float(np.max(norm_log_likelihood[finite_norm]))
    norm_density = np.exp(norm_log_likelihood - log_shift)
    norm_density[~finite_norm] = 0.0
    area = float(trapezoid(norm_density, norm_grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Joint posterior normalization failed.")

    posterior = np.zeros_like(grid)
    posterior[support] = np.exp(joint_log_likelihood[support] - log_shift) / area
    return posterior
