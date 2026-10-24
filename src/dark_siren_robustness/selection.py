"""Toy event-detection selection normalization for synthetic dark-siren experiments.

This module implements a deliberately simplified selection rule: an event is
"detected" when its noisy distance observation lies inside a fixed interval.
It is a statistical testbed for the beta(H0) normalization, not a model of a
real gravitational-wave detector or LVK search pipeline.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import ndtr

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
    host_efficiency = ndtr((threshold - means) / sigma) - ndtr((lower_threshold - means) / sigma)
    beta = host_efficiency @ w
    # Floating-point summation can exceed one by a few ulps in the no-selection limit.
    beta = np.clip(beta, np.finfo(float).tiny, 1.0)
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1):
        raise ValueError("Computed detection efficiency is outside (0, 1].")
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
    beta = detection_efficiency_h0(
        grid, redshifts, gw_sigma, threshold, weights=weights,
        detection_distance_min=detection_distance_min,
    )
    support = (grid >= prior_min) & (grid <= prior_max)
    selected = np.zeros_like(base)
    selected[support] = base[support] / beta[support]
    area = float(trapezoid(selected[support], grid[support]))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Selection-conditioned posterior normalization failed.")
    selected[support] /= area
    return selected
