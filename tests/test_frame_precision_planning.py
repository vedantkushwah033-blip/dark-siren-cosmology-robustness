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


def test_frame_planner_separates_pilot_df_from_confirmatory_frame_search():
    result = plan_independent_frames(
        0.08, 0.05, minimum_frames=8, maximum_frames=1000,
        pilot_frame_count=12, variance_confidence=0.80,
    )
    assert result["pilot_frame_count"] == 12
    assert result["minimum_frames_for_sd_bound"] == 12
    assert result["minimum_frames_for_target"] >= 8
    assert result["planned_half_width"] <= 0.05


def test_frame_planner_rejects_invalid_pilot_frame_count():
    with pytest.raises(ValueError, match="pilot_frame_count"):
        plan_independent_frames(0.1, 0.05, pilot_frame_count=2)
    with pytest.raises(TypeError, match="pilot_frame_count"):
        plan_independent_frames(0.1, 0.05, pilot_frame_count=True)


def test_frame_planner_fails_explicitly_when_cap_cannot_meet_precision():
    with pytest.raises(ValueError, match="maximum_frames"):
        plan_independent_frames(
            10.0, 1e-12, minimum_frames=3, maximum_frames=3
        )

def test_frame_planner_returns_smallest_candidate_meeting_precision():
    result = plan_independent_frames(
        0.12,
        0.04,
        minimum_frames=8,
        maximum_frames=20_000,
        pilot_frame_count=15,
        variance_confidence=0.90,
    )
    n = result["minimum_frames_for_target"]
    assert result["planned_half_width"] <= result["target_half_width"]
    if n > 8:
        from scipy.stats import t

        previous_half_width = (
            t.ppf(0.975, n - 2)
            * result["upper_planning_sd"]
            / (n - 1) ** 0.5
        )
        assert previous_half_width > result["target_half_width"]


def test_more_pilot_frames_tighten_upper_sd_bound():
    short = plan_independent_frames(
        0.1, 0.05, minimum_frames=8, pilot_frame_count=5,
        variance_confidence=0.90,
    )
    longer = plan_independent_frames(
        0.1, 0.05, minimum_frames=8, pilot_frame_count=30,
        variance_confidence=0.90,
    )
    assert longer["upper_planning_sd"] < short["upper_planning_sd"]
