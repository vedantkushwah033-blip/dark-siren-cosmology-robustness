"""Continuous hidden-host population model for synthetic dark-siren inference.

Observation model: u=log(d) | z,H0 ~ Normal(log(c*z/H0), sigma_u).
The hidden redshift is a truncated Normal on fixed support (z_min,z_max).
This module is deliberately separate from the existing positive-truncated
Gaussian-distance experiment family.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _logdiffexp(log_a: float, log_b: float) -> float:
    """log(exp(log_a)-exp(log_b)) for log_a >= log_b."""
    if log_b > log_a:
        raise ValueError("logdiffexp requires log_a >= log_b")
    if np.isneginf(log_b):
        return float(log_a)
    delta = log_b - log_a
    if delta >= 0:
        return -np.inf
    return float(log_a + np.log1p(-np.exp(delta)))


def _log_truncnorm_mass(a: float, b: float, mu: float, sd: float) -> float:
    lo, hi = (a - mu) / sd, (b - mu) / sd
    # Use survival probabilities in the positive tail to avoid subtracting
    # two CDF values both rounded close to one.
    if lo > 0:
        return _logdiffexp(float(log_ndtr(-lo)), float(log_ndtr(-hi)))
    return _logdiffexp(float(log_ndtr(hi)), float(log_ndtr(lo)))


def _validate_observations(u, sigma_u):
    values = np.asarray(u, dtype=float)
    if values.ndim != 1 or values.size < 1 or not np.all(np.isfinite(values)):
        raise ValueError("u must be a non-empty finite one-dimensional array")
    sigma = float(sigma_u)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma_u must be finite and strictly positive")
    return values, sigma


def _validate_components(listed_redshifts, listed_weights, z_bounds, order):
    z = np.asarray(listed_redshifts, dtype=float)
    w = np.asarray(listed_weights, dtype=float)
    if z.ndim != 1 or z.size < 1 or not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("listed_redshifts must be finite and strictly positive")
    if w.shape != z.shape or not np.all(np.isfinite(w)) or np.any(w < 0) or not np.any(w > 0):
        raise ValueError("listed_weights must be finite, nonnegative and match listed_redshifts")
    w = w / w.sum()
    bounds = np.asarray(z_bounds, dtype=float)
    if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] <= 0 or bounds[0] >= bounds[1]:
        raise ValueError("z_bounds must be two finite values with 0 < lower < upper")
    if np.any(z < bounds[0]) or np.any(z > bounds[1]):
        # Listed hosts need not lie inside the hidden component's support.
        pass
    q = int(order)
    if q < 8:
        raise ValueError("quadrature_order must be at least 8")
    return z, w, (float(bounds[0]), float(bounds[1])), q


def hidden_redshift_log_density(z, mu_h, sigma_h, z_bounds=(0.005, 0.10)):
    """Log density of a Normal(mu_h,sigma_h) truncated to fixed z_bounds."""
    values = np.asarray(z, dtype=float)
    mu, sd = float(mu_h), float(sigma_h)
    a, b = map(float, z_bounds)
    if not np.isfinite(mu) or not np.isfinite(sd) or sd <= 0:
        raise ValueError("mu_h must be finite and sigma_h finite and positive")
    if not (0 < a < b) or not np.all(np.isfinite([a, b])):
        raise ValueError("z_bounds must satisfy 0 < lower < upper")
    if np.any(values < a) or np.any(values > b) or not np.all(np.isfinite(values)):
        raise ValueError("z values must lie within the truncation support")
    log_mass = _log_truncnorm_mass(a, b, mu, sd)
    return norm.logpdf(values, loc=mu, scale=sd) - log_mass


def _legendre_rule(a, b, order):
    x, w = np.polynomial.legendre.leggauss(order)
    z = (b - a) * (x + 1.0) / 2.0 + a
    weights = w * (b - a) / 2.0
    return z, np.log(weights)


def hidden_component_log_density(
    u, H0, mu_h, sigma_h, sigma_u, z_bounds=(0.005, 0.10), quadrature_order=64
):
    """Compute log p(u_i | H0, mu_h, sigma_h) by Gauss-Legendre quadrature."""
    observations, su = _validate_observations(u, sigma_u)
    h, mu, sh = float(H0), float(mu_h), float(sigma_h)
    a, b = map(float, z_bounds)
    if not np.isfinite(h) or h <= 0:
        raise ValueError("H0 must be finite and positive")
    if not np.isfinite(mu) or not np.isfinite(sh) or sh <= 0:
        raise ValueError("mu_h and sigma_h must be finite; sigma_h must be positive")
    if int(quadrature_order) < 8:
        raise ValueError("quadrature_order must be at least 8")
    z, log_weights = _legendre_rule(a, b, int(quadrature_order))
    log_pz = hidden_redshift_log_density(z, mu, sh, (a, b))
    means = np.log(C_KM_S * z / h)
    log_pu_given_z = norm.logpdf(observations[:, None], loc=means[None, :], scale=su)
    return logsumexp(log_pu_given_z + log_pz[None, :] + log_weights[None, :], axis=1)


def log_likelihood(
    u, H0, hidden_fraction, mu_h, sigma_h, sigma_u,
    listed_redshifts=(0.02, 0.04), listed_weights=(0.7, 0.3),
    z_bounds=(0.005, 0.10), quadrature_order=64,
):
    """Joint log likelihood for independent events under listed/hidden mixture."""
    observations, su = _validate_observations(u, sigma_u)
    z, weights, bounds, order = _validate_components(
        listed_redshifts, listed_weights, z_bounds, quadrature_order
    )
    h, f, mu, sh = map(float, (H0, hidden_fraction, mu_h, sigma_h))
    if not np.isfinite(h) or h <= 0 or not np.isfinite(f) or not 0 <= f <= 1:
        raise ValueError("H0 must be positive and hidden_fraction must lie in [0,1]")
    if not np.isfinite(mu) or not np.isfinite(sh) or sh <= 0:
        raise ValueError("mu_h must be finite and sigma_h finite and positive")

    listed_means = np.log(C_KM_S * z / h)
    listed_terms = norm.logpdf(
        observations[:, None], loc=listed_means[None, :], scale=su
    ) + np.log(weights)[None, :]
    log_listed = logsumexp(listed_terms, axis=1)

    if f == 0.0:
        return float(np.sum(log_listed))
    log_hidden = hidden_component_log_density(
        observations, h, mu, sh, su, bounds, order
    )
    if f == 1.0:
        return float(np.sum(log_hidden))
    return float(np.sum(np.logaddexp(np.log1p(-f) + log_listed, np.log(f) + log_hidden)))


def _trapezoid_weights(grid):
    x = np.asarray(grid, dtype=float)
    w = np.empty_like(x)
    w[0] = (x[1] - x[0]) / 2
    w[-1] = (x[-1] - x[-2]) / 2
    if x.size > 2:
        w[1:-1] = (x[2:] - x[:-2]) / 2
    return w


def _validate_grid(grid, name, positive=False, lower=None, upper=None):
    x = np.asarray(grid, dtype=float)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
        raise ValueError(f"{name} must be a finite, strictly increasing 1D grid")
    if positive and np.any(x <= 0):
        raise ValueError(f"{name} must be strictly positive")
    if lower is not None and x[0] < lower:
        raise ValueError(f"{name} extends below its allowed lower bound")
    if upper is not None and x[-1] > upper:
        raise ValueError(f"{name} extends above its allowed upper bound")
    return x


def posterior_grid(
    u, sigma_u, H0_grid, fraction_grid, mu_grid, sigma_h_grid,
    listed_redshifts=(0.02, 0.04), listed_weights=(0.7, 0.3),
    z_bounds=(0.005, 0.10), quadrature_order=48,
):
    """Evaluate and normalize a 4D uniform-prior posterior on a Cartesian grid.

    Returns a density array with axes (H0, fraction, mu_h, sigma_h), normalized
    under the tensor-product trapezoidal integration rule.
    """
    observations, su = _validate_observations(u, sigma_u)
    hz = _validate_grid(H0_grid, "H0_grid", positive=True)
    fg = _validate_grid(fraction_grid, "fraction_grid", lower=0.0, upper=1.0)
    mg = _validate_grid(mu_grid, "mu_grid")
    sg = _validate_grid(sigma_h_grid, "sigma_h_grid", positive=True)
    z, weights, bounds, order = _validate_components(
        listed_redshifts, listed_weights, z_bounds, quadrature_order
    )
    if np.any(mg <= 0):
        raise ValueError("mu_grid must be strictly positive")

    wh, wf, wm, ws = (_trapezoid_weights(x) for x in (hz, fg, mg, sg))
    log_post = np.full((hz.size, fg.size, mg.size, sg.size), -np.inf, dtype=float)
    log_w_listed = np.log(weights)

    # The hidden component is independent of f. Compute it once per (mu,sigma_h)
    # pair and reuse it across the fraction grid.
    for im, mu in enumerate(mg):
        for isd, sh in enumerate(sg):
            hidden_by_h0 = []
            for h in hz:
                hidden_by_h0.append(hidden_component_log_density(
                    observations, h, mu, sh, su, bounds, order
                ))
            for ih, h in enumerate(hz):
                means = np.log(C_KM_S * z / h)
                log_listed = logsumexp(
                    norm.logpdf(observations[:, None], loc=means[None, :], scale=su)
                    + log_w_listed[None, :], axis=1
                )
                log_hidden = hidden_by_h0[ih]
                for jf, fraction in enumerate(fg):
                    if fraction == 0:
                        ll = np.sum(log_listed)
                    elif fraction == 1:
                        ll = np.sum(log_hidden)
                    else:
                        ll = np.sum(np.logaddexp(
                            np.log1p(-fraction) + log_listed,
                            np.log(fraction) + log_hidden,
                        ))
                    log_post[ih, jf, im, isd] = ll

    log_post -= np.max(log_post)
    density = np.exp(log_post)
    tensor_weights = wh[:, None, None, None] * wf[None, :, None, None] * wm[None, None, :, None] * ws[None, None, None, :]
    normalizer = float(np.sum(density * tensor_weights))
    if not np.isfinite(normalizer) or normalizer <= 0:
        raise ValueError("Posterior grid normalization failed")
    density /= normalizer
    return {
        "H0_grid": hz, "fraction_grid": fg, "mu_grid": mg, "sigma_h_grid": sg,
        "density": density,
    }


def _marginal_density(density, axes, axis):
    other = tuple(i for i in range(density.ndim) if i != axis)
    integrated = density
    # Integrate from highest axis down so axis indices remain valid.
    for current in sorted(other, reverse=True):
        integrated = trapezoid(integrated, x=axes[current], axis=current)
    return integrated


def _quantile_from_density(grid, density, probability):
    area = float(trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Marginal posterior has non-positive area")
    d = density / area
    cdf = np.concatenate(([0.0], np.cumsum((d[1:] + d[:-1]) * np.diff(grid) / 2)))
    cdf /= cdf[-1]
    return float(np.interp(probability, cdf, grid))


def posterior_marginals(posterior):
    """Return marginal posterior densities and median/central 68/95% intervals."""
    names = ("H0", "fraction", "mu_h", "sigma_h")
    keys = ("H0_grid", "fraction_grid", "mu_grid", "sigma_h_grid")
    axes = [np.asarray(posterior[k], dtype=float) for k in keys]
    density = np.asarray(posterior["density"], dtype=float)
    expected_shape = tuple(x.size for x in axes)
    if density.shape != expected_shape or not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ValueError("Posterior density shape/values do not match the grids")
    results = {}
    for i, (name, grid) in enumerate(zip(names, axes)):
        marginal = _marginal_density(density, axes, i)
        results[name] = {
            "grid": grid,
            "density": marginal,
            "median": _quantile_from_density(grid, marginal, 0.5),
            "lower_68": _quantile_from_density(grid, marginal, 0.16),
            "upper_68": _quantile_from_density(grid, marginal, 0.84),
            "lower_95": _quantile_from_density(grid, marginal, 0.025),
            "upper_95": _quantile_from_density(grid, marginal, 0.975),
        }
    return results
