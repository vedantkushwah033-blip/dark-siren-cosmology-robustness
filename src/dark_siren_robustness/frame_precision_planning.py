"""Precision-based planning for independent finite-population frames.

This planner consumes an externally estimated standard deviation of frame-level
contrasts. It does not infer that quantity from the three-frame pilot and does
not choose a confirmatory design automatically.
"""
from __future__ import annotations

import math
from numbers import Real

from scipy.stats import chi2, t


def plan_independent_frames(
    pilot_frame_sd: float,
    target_half_width: float,
    *,
    confidence_level: float = 0.95,
    variance_confidence: float = 0.95,
    minimum_frames: int = 3,
    maximum_frames: int = 1_000_000,
    pilot_frame_count: int | None = None,
) -> dict[str, float | int | str]:
    """Find frame count meeting a t-based half-width using an upper SD bound.

    The upper confidence bound for the frame-level SD is
    s * sqrt(df / chi2.ppf(1 - variance_confidence, df)), where df is
    pilot_frame_count - 1. If omitted, pilot_frame_count defaults to
    minimum_frames for backward compatibility. The search for the target
    frame count is separate from the pilot degrees of freedom. This is a
    conservative planning plug-in, not a guarantee: it assumes independent
    approximately normal frame-level contrasts and a representative pilot SD.
    """
    for name, value in (
        ("pilot_frame_sd", pilot_frame_sd),
        ("target_half_width", target_half_width),
        ("confidence_level", confidence_level),
        ("variance_confidence", variance_confidence),
    ):
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"{name} must be a real number, not a boolean")
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
    if pilot_frame_count is None:
        pilot_frame_count = minimum_frames
    if isinstance(pilot_frame_count, bool) or not isinstance(pilot_frame_count, int):
        raise TypeError("pilot_frame_count must be an integer")
    if pilot_frame_count < 3:
        raise ValueError("pilot_frame_count must be at least 3")

    df = pilot_frame_count - 1
    lower_tail = 1.0 - variance_confidence
    variance_quantile = float(chi2.ppf(lower_tail, df))
    if not math.isfinite(variance_quantile) or variance_quantile <= 0:
        raise ValueError("variance confidence bound is numerically undefined")
    planning_sd = pilot_frame_sd * math.sqrt(df / variance_quantile)
    alpha = 1.0 - confidence_level

    def half_width(n_frames: int) -> float:
        return float(
            t.ppf(1.0 - alpha / 2.0, n_frames - 1)
            * planning_sd / math.sqrt(n_frames)
        )

    if half_width(maximum_frames) > target_half_width:
        raise ValueError("target precision not reached by maximum_frames")
    low, high = minimum_frames, maximum_frames
    while low < high:
        middle = (low + high) // 2
        if half_width(middle) <= target_half_width:
            high = middle
        else:
            low = middle + 1
    selected = low
    selected_half_width = half_width(selected)

    return {
        "status": "PLANNING_ONLY_NOT_CONFIRMATORY",
        "minimum_frames_for_target": selected,
        "pilot_frame_sd": float(pilot_frame_sd),
        "upper_planning_sd": float(planning_sd),
        "target_half_width": float(target_half_width),
        "planned_half_width": float(selected_half_width),
        "confidence_level": float(confidence_level),
        "variance_confidence": float(variance_confidence),
        "minimum_frames_for_sd_bound": int(pilot_frame_count),
        "pilot_frame_count": int(pilot_frame_count),
        "assumption": "independent approximately normal frame-level contrasts",
    }

def plan_frame_count_sensitivity(
    pilot_frame_sd: float,
    target_half_width: float,
    *,
    pilot_frame_count: int,
    sd_multipliers: tuple[float, ...] = (0.75, 1.0, 1.25, 1.5),
    variance_confidences: tuple[float, ...] = (0.80, 0.90, 0.95, 0.975),
    confidence_level: float = 0.95,
    minimum_frames: int = 3,
    maximum_frames: int = 1_000_000,
) -> list[dict[str, float | int | str]]:
    """Evaluate frame-count planning across prespecified SD/confidence stresses.

    This is a sensitivity table, not a data-driven choice of the confirmatory
    design. Every row retains the planning-only label from the base planner.
    """
    if not isinstance(sd_multipliers, tuple) or not sd_multipliers:
        raise TypeError("sd_multipliers must be a non-empty tuple")
    if not isinstance(variance_confidences, tuple) or not variance_confidences:
        raise TypeError("variance_confidences must be a non-empty tuple")
    for multiplier in sd_multipliers:
        if isinstance(multiplier, bool) or not isinstance(multiplier, Real):
            raise TypeError("all sd_multipliers must be real numbers, not booleans")
        if not math.isfinite(multiplier) or multiplier <= 0:
            raise ValueError("all sd_multipliers must be finite and positive")
    for level in variance_confidences:
        if isinstance(level, bool) or not isinstance(level, Real):
            raise TypeError("all variance_confidences must be real numbers, not booleans")
        if not math.isfinite(level) or not 0 < level < 1:
            raise ValueError("all variance_confidences must lie strictly between zero and one")

    results: list[dict[str, float | int | str]] = []
    for level in variance_confidences:
        for multiplier in sd_multipliers:
            row = plan_independent_frames(
                pilot_frame_sd * multiplier,
                target_half_width,
                confidence_level=confidence_level,
                variance_confidence=level,
                minimum_frames=minimum_frames,
                maximum_frames=maximum_frames,
                pilot_frame_count=pilot_frame_count,
            )
            row["sd_multiplier"] = float(multiplier)
            results.append(row)
    return results
