from __future__ import annotations

import pytest

from dark_siren_robustness.frame_precision_planning import plan_independent_frames


def test_frame_planner_returns_precision_attainment_and_planning_label():
    result = plan_independent_frames(
        0.08, 0.05, minimum_frames=8, variance_confidence=0.80
    )
    assert result["status"] == "PLANNING_ONLY_NOT_CONFIRMATORY"
    assert result["minimum_frames_for_target"] >= 8
    assert result["planned_half_width"] <= 0.05
    assert result["upper_planning_sd"] > 0.08


@pytest.mark.parametrize(
    "kwargs",
    [
        {"pilot_frame_sd": 0.0, "target_half_width": 0.1},
        {"pilot_frame_sd": 0.1, "target_half_width": 0.0},
        {"pilot_frame_sd": float("nan"), "target_half_width": 0.1},
        {"pilot_frame_sd": 0.1, "target_half_width": 0.1, "minimum_frames": 2},
        {"pilot_frame_sd": 0.1, "target_half_width": 0.1, "maximum_frames": 2},
        {"pilot_frame_sd": 0.1, "target_half_width": 0.1, "confidence_level": 1.0},
    ],
)
def test_frame_planner_rejects_invalid_inputs(kwargs):
    with pytest.raises((TypeError, ValueError)):
        plan_independent_frames(**kwargs)


def test_frame_planner_fails_explicitly_when_cap_cannot_meet_precision():
    with pytest.raises(ValueError, match="maximum_frames"):
        plan_independent_frames(
            10.0, 1e-12, minimum_frames=3, maximum_frames=3
        )
