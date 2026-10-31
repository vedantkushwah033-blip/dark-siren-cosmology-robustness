"""Precision-based planning for independent finite-population frames.

This planner consumes an externally estimated standard deviation of frame-level
contrasts. It does not infer that quantity from the three-frame pilot and does
not choose a confirmatory design automatically.
"""
from __future__ import annotations

import math

from scipy.stats import chi2, t


def plan_independent_frames(
    pilot_frame_sd: float,
    target_half_width: float,
    *,
    confidence_level: float = 0.95,
    variance_confidence: float = 0.95,
    minimum_frames: int = 3,
    maximum_frames: int = 1_000_000,
) -> dict[str, float | int | str]:
    """Find frame count meeting a t-based half-width using an upper SD bound.

    The upper confidence bound for the frame-level SD is
    s * sqrt(df / chi2.ppf(1 - variance_confidence, df)), where df is
    minimum_frames - 1. This is a conservative planning plug-in, not a
    guarantee: it assumes independent approximately normal frame-level
    contrasts and that the pilot SD is representative of the target design.
    """
    for name, value in (
        ("pilot_frame_sd", pilot_frame_sd),
        ("target_half_width", target_half_width),
        ("confidence_level", confidence_level),
        ("variance_confidence", variance_confidence),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if pilot_frame_sd <= 0:
        raise ValueError("pilot_frame_sd must be positive")
    if target_half_width <= 0:
        raise ValueError("target_half_width must be positive")
    if not 0 < confidence_level < 1 or not 0 < variance_confidence < 1:
        raise ValueError("confidence levels must lie strictly between zero and one")
    if isinstance(minimum_frames, bool) or not isinstance(minimum_frames, int):
        raise TypeError("minimum_frames must be an integer")
    if isinstance(maximum_frames, bool) or not isinstance(maximum_frames, int):
        raise TypeError("maximum_frames must be an integer")
    if minimum_frames < 3 or maximum_frames < minimum_frames:
        raise ValueError("require 3 <= minimum_frames <= maximum_frames")

    df = minimum_frames - 1
    lower_tail = 1.0 - variance_confidence
    variance_quantile = float(chi2.ppf(lower_tail, df))
    if not math.isfinite(variance_quantile) or variance_quantile <= 0:
        raise ValueError("variance confidence bound is numerically undefined")
    planning_sd = pilot_frame_sd * math.sqrt(df / variance_quantile)
    alpha = 1.0 - confidence_level

    selected = None
    selected_half_width = None
    for n_frames in range(minimum_frames, maximum_frames + 1):
        half_width = float(t.ppf(1.0 - alpha / 2.0, n_frames - 1)
                           * planning_sd / math.sqrt(n_frames))
        if half_width <= target_half_width:
            selected = n_frames
            selected_half_width = half_width
            break
    if selected is None:
        raise ValueError("target precision not reached by maximum_frames")

    return {
        "status": "PLANNING_ONLY_NOT_CONFIRMATORY",
        "minimum_frames_for_target": selected,
        "pilot_frame_sd": float(pilot_frame_sd),
        "upper_planning_sd": float(planning_sd),
        "target_half_width": float(target_half_width),
        "planned_half_width": float(selected_half_width),
        "confidence_level": float(confidence_level),
        "variance_confidence": float(variance_confidence),
        "minimum_frames_for_sd_bound": int(minimum_frames),
        "assumption": "independent approximately normal frame-level contrasts",
    }
