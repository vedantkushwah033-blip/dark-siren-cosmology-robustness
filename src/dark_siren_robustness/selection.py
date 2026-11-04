"""Toy event-detection selection normalization for synthetic dark-siren experiments.

This module implements a deliberately simplified selection rule: an event is
"detected" when its noisy distance observation lies inside a fixed interval.
It is a statistical testbed for the beta(H0) normalization, not a model of a
real gravitational-wave detector or LVK search pipeline.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import erf, log_ndtr, ndtr


def _normal_interval_probability(lower_z, upper_z):
    """Stable standard-normal probability between ordered bounds.

    Use survival-function symmetry for intervals wholly above zero and a
    log-CDF difference for intervals wholly below zero. The latter avoids
    losing representable probability when both CDF values are tiny and close
    together. Intervals touching or crossing zero use an erf difference to
    avoid cancellation when both CDF values round to one half.
    """
    lower_z, upper_z = np.broadcast_arrays(
        np.asarray(lower_z, dtype=float), np.asarray(upper_z, dtype=float)
    )
    if not np.all(np.isfinite(lower_z)) or not np.all(np.isfinite(upper_z)):
        raise ValueError("normal interval bounds must be finite")
    if np.any(upper_z < lower_z):
        raise ValueError("upper normal bound must not be below lower bound")

    probability = np.empty(lower_z.shape, dtype=float)
    positive_tail = lower_z > 0.0
    negative_tail = upper_z < 0.0
    central = ~(positive_tail | negative_tail)

    if np.any(positive_tail):
        # Phi(b)-Phi(a) = Phi(-a)-Phi(-b). Evaluate the difference in log
        # space so close, far-positive-tail probabilities retain precision.
        log_lower_survival = log_ndtr(-lower_z[positive_tail])
        log_upper_survival = log_ndtr(-upper_z[positive_tail])
        probability[positive_tail] = np.exp(log_lower_survival) * (
            -np.expm1(log_upper_survival - log_lower_survival)
        )
    if np.any(negative_tail):
        # exp(log Phi(b)) * (1 - exp(log Phi(a)-log Phi(b))).
        # expm1 retains relative precision when the two log-CDFs are close.
        log_lower = log_ndtr(lower_z[negative_tail])
        log_upper = log_ndtr(upper_z[negative_tail])
        probability[negative_tail] = np.exp(log_upper) * (
            -np.expm1(log_lower - log_upper)
        )
    if np.any(central):
        # For intervals crossing zero, CDF values can both round to 0.5 when
        # the interval is extremely narrow. The erf identity subtracts values
        # with opposite signs and preserves the representable interval mass.
        probability[central] = 0.5 * (
            erf(upper_z[central] / np.sqrt(2.0))
            - erf(lower_z[central] / np.sqrt(2.0))
        )
    return probability

from .hardened_inference import C_KM_S, posterior_h0


def detection_efficiency_h0(H0_grid, redshifts, gw_sigma, detection_distance_max,
                            weights=None, detection_distance_min=0.0):
    """Return beta(H0) = P(detected | H0) for the toy finite distance interval.

    A datum is detected when detection_distance_min < d_obs <
    detection_distance_max and d_obs | (H0, host j) ~
    Normal(c*z_j/H0, gw_sigma**2).
    """
    grid = np.asarray(H0_grid, dtype=float)
    z = np.asarray(redshifts, dtype=float)
    sigma = float(gw_sigma)
    threshold = float(detection_distance_max)
    lower_threshold = float(detection_distance_min)
    if not np.isfinite(lower_threshold) or lower_threshold < 0 or lower_threshold >= threshold:
        raise ValueError("detection_distance_min must be finite, non-negative and below the upper threshold.")
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)) or np.any(grid <= 0):
        raise ValueError("H0_grid must be a finite, positive one-dimensional grid.")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("H0_grid must be strictly increasing.")
    if z.ndim != 1 or z.size < 1 or not np.all(np.isfinite(z)) or np.any(z < 0):
        raise ValueError("redshifts must be a finite, non-negative one-dimensional array.")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("gw_sigma must be finite and positive.")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_distance_max must be finite and positive.")
    w = np.ones(z.size, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size != z.size or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("weights must be finite, non-negative and match redshifts.")
    if w.sum() <= 0:
        raise ValueError("weights must have a positive sum.")
    w = w / w.sum()
    means = C_KM_S * z[None, :] / grid[:, None]
    # Detection is the interval lower_threshold < d_obs < threshold.
    # Include both boundaries in beta so it matches the generator exactly.
    host_efficiency = _normal_interval_probability(
        (lower_threshold - means) / sigma,
        (threshold - means) / sigma,
    )
    beta = host_efficiency @ w
    # Roundoff can push a probability a few ulps above one, but a zero
    # computed efficiency must not be replaced by machine tiny: that would
    # disguise numerical underflow / effectively impossible detection as valid.
    beta = np.minimum(beta, 1.0)
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1):
        raise ValueError(
            "Computed detection efficiency is not representable as a positive probability "
            "on this grid; revise the model/grid or use numerically stable tail calculations."
        )
    return beta


def posterior_h0_selected(H0_grid, redshifts, observed_distance, gw_sigma,
                          detection_distance_max, weights=None,
                          prior_min=50.0, prior_max=90.0,
                          detection_distance_min=0.0):
    """Return posterior density conditioned on the toy event being detected.

    The selected likelihood is proportional to the unselected host-mixture
    likelihood divided by beta(H0). The observed distance must satisfy the
    same finite interval used to define detection.
    """
    distance = float(observed_distance)
    threshold = float(detection_distance_max)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("observed_distance must be finite and positive.")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_distance_max must be finite and positive.")
    lower_threshold = float(detection_distance_min)
    if not np.isfinite(lower_threshold) or lower_threshold < 0 or lower_threshold >= threshold:
        raise ValueError("detection_distance_min must be finite, non-negative and below the upper threshold.")
    if not lower_threshold < distance < threshold:
        raise ValueError("observed_distance must lie strictly inside the detection interval.")
    grid = np.asarray(H0_grid, dtype=float)
    base = posterior_h0(
        grid, redshifts, distance, gw_sigma, weights=weights,
        prior_min=prior_min, prior_max=prior_max,
    )
    support = (grid >= prior_min) & (grid <= prior_max)
    if np.count_nonzero(support) < 2:
        raise ValueError("H0_grid must contain at least two points inside the prior interval.")
    # Only evaluate the selection denominator where the posterior is defined.
    # A grid point outside prior support can have numerically unrepresentable
    # detection efficiency even when every in-prior value is well behaved; such
    # an irrelevant point must not invalidate the posterior calculation.
    beta_support = detection_efficiency_h0(
        grid[support], redshifts, gw_sigma, threshold, weights=weights,
        detection_distance_min=detection_distance_min,
    )
    selected = np.zeros_like(base)
    selected[support] = base[support] / beta_support
    area = float(trapezoid(selected[support], grid[support]))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Selection-conditioned posterior normalization failed.")
    selected[support] /= area
    return selected


def posterior_h0_selected_with_missing_hosts(
    H0_grid, catalog_redshifts, observed_distance, gw_sigma,
    detection_distance_max, *, catalog_weights=None, missing_weight=0.0,
    missing_redshifts=None, missing_weights=None, prior_min=50.0,
    prior_max=90.0, detection_distance_min=0.0,
):
    """Infer H0 with event-selection correction and an explicit missing-host mixture.

    missing_weight is the prior population mass assigned to the unobserved host
    component, not a fraction to apply after normalizing the full likelihood.
    Catalogue and missing-component weights are each normalized internally,
    then scaled by their population masses before combining. Full-population
    weights are used in both the event likelihood and detection denominator.

    This finite synthetic mixture is not a practical survey correction:
    callers must supply a defensible model for missing-host redshifts and mass.
    """
    catalog_z = np.asarray(catalog_redshifts, dtype=float)
    if catalog_z.ndim != 1 or catalog_z.size == 0 or not np.all(np.isfinite(catalog_z)) or np.any(catalog_z < 0):
        raise ValueError("catalog_redshifts must be a non-empty finite non-negative vector.")
    fraction = float(missing_weight)
    if not np.isfinite(fraction) or not 0.0 <= fraction < 1.0:
        raise ValueError("missing_weight must be finite and in [0, 1).")

    cw = np.ones(catalog_z.size, dtype=float) if catalog_weights is None else np.asarray(catalog_weights, dtype=float)
    if cw.ndim != 1 or cw.size != catalog_z.size or not np.all(np.isfinite(cw)) or np.any(cw < 0) or cw.sum() <= 0:
        raise ValueError("catalog_weights must be finite, non-negative, match catalogue hosts, and sum positive.")
    cw = cw / cw.sum() * (1.0 - fraction)

    if fraction > 0.0:
        if missing_redshifts is None:
            raise ValueError("missing_redshifts are required when missing_weight > 0.")
        missing_z = np.asarray(missing_redshifts, dtype=float)
        if missing_z.ndim != 1 or missing_z.size == 0 or not np.all(np.isfinite(missing_z)) or np.any(missing_z < 0):
            raise ValueError("missing_redshifts must be a non-empty finite non-negative vector.")
        mw = np.ones(missing_z.size, dtype=float) if missing_weights is None else np.asarray(missing_weights, dtype=float)
        if mw.ndim != 1 or mw.size != missing_z.size or not np.all(np.isfinite(mw)) or np.any(mw < 0) or mw.sum() <= 0:
            raise ValueError("missing_weights must be finite, non-negative, match missing hosts, and sum positive.")
        mw = mw / mw.sum() * fraction
        full_z = np.concatenate((catalog_z, missing_z))
        full_weights = np.concatenate((cw, mw))
    else:
        if missing_weights is not None:
            raise ValueError("missing_weights were supplied but missing_weight is zero.")
        full_z, full_weights = catalog_z, cw

    # Both the event likelihood mixture and beta(H0) use the same full-population
    # weights; no catalogue-only renormalization is applied to the missing term.
    return posterior_h0_selected(
        H0_grid, full_z, observed_distance, gw_sigma, detection_distance_max,
        weights=full_weights, prior_min=prior_min, prior_max=prior_max,
        detection_distance_min=detection_distance_min,
    )
