from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_inclusion_posterior_calibration.py"
SPEC = importlib.util.spec_from_file_location("inclusion_calibration_runner", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_posterior_pit_uses_normalized_trapezoid_cdf():
    frame = pd.DataFrame({
        "H0": [0.0, 1.0, 2.0],
        "density": [0.0, 1.0, 0.0],
    })
    assert MODULE._posterior_pit(frame, 1.0, "density") == pytest.approx(0.5)
    assert MODULE._posterior_pit(frame, 0.0, "density") == pytest.approx(0.0)
    assert MODULE._posterior_pit(frame, 2.0, "density") == pytest.approx(1.0)


def test_posterior_pit_rejects_truth_outside_grid_and_invalid_density():
    frame = pd.DataFrame({"H0": [50.0, 60.0, 70.0], "density": [0.1, 0.1, 0.1]})
    with pytest.raises(ValueError, match="outside posterior support"):
        MODULE._posterior_pit(frame, 49.0, "density")
    bad = frame.copy()
    bad["density"] = 0.0
    with pytest.raises(ValueError, match="invalid integrated mass"):
        MODULE._posterior_pit(bad, 60.0, "density")
    negative = frame.copy()
    negative["density"] = [0.1, -0.1, 0.1]
    with pytest.raises(ValueError, match="invalid posterior grid/density"):
        MODULE._posterior_pit(negative, 60.0, "density")


def test_posterior_pit_rejects_nonfinite_generating_h0():
    frame = pd.DataFrame({"H0": [50.0, 60.0, 70.0], "density": [0.1, 0.1, 0.1]})
    with pytest.raises(ValueError, match="must be finite"):
        MODULE._posterior_pit(frame, float("nan"), "density")
    with pytest.raises(ValueError, match="must be finite"):
        MODULE._posterior_pit(frame, float("inf"), "density")


def test_wilson_interval_contains_observed_proportion_and_is_bounded():
    low, high = MODULE._wilson_interval(68, 100)
    assert 0.0 <= low < 0.68 < high <= 1.0
    empty = MODULE._wilson_interval(0, 0)
    assert np.isnan(empty[0]) and np.isnan(empty[1])


def test_wilson_interval_rejects_invalid_counts_and_z_score():
    for successes, trials in ((-1, 10), (11, 10), (1, -1), (1.5, 10)):
        with pytest.raises(ValueError):
            MODULE._wilson_interval(successes, trials)
    with pytest.raises(ValueError, match="finite and positive"):
        MODULE._wilson_interval(1, 10, z=float("nan"))


def test_seed_streams_are_deterministic_and_separated():
    a = MODULE._seed(1234, 0, 1, 2, 0)
    b = MODULE._seed(1234, 0, 1, 2, 0)
    c = MODULE._seed(1234, 0, 1, 2, 1)
    assert a == b
    assert a != c


def test_summary_keeps_non_evaluable_and_failed_attempts_visible():
    frame = pd.DataFrame([
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "n_validation": 100,
         "status": "SUCCESS", "covered_68": 1, "covered_90": 1, "covered_95": 1,
         "width_68": 10.0, "width_90": 20.0, "width_95": 25.0,
         "h0_median": 71.0, "h0_true": 70.0,
         "coefficient_quadrature_ess": 20.0,
         "minimum_catalogue_kish_ess": 500.0,
         "maximum_catalogue_normalized_weight": 0.02},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "n_validation": 100,
         "status": "NON_EVALUABLE"},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "n_validation": 100,
         "status": "FAILED"},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "n_validation": 100,
         "status": "NOT_RUN"},
    ])
    summary = MODULE._summarize(frame, bins=10)
    cell = summary["FIXED_TRUTH|floor=0|n=100"]
    assert cell["attempted"] == 3
    assert cell["success"] == 1
    assert cell["non_evaluable"] == 1
    assert cell["failed"] == 1
    assert cell["not_run"] == 1
    assert cell["evaluable_fraction"] == pytest.approx(1 / 3)
    assert cell["coverage_68"] == pytest.approx(1.0)
    assert cell["mean_posterior_median_bias"] == pytest.approx(1.0)
    assert cell["median_minimum_catalogue_kish_ess"] == pytest.approx(500.0)


def test_conditional_h0_pit_summary_reports_interval_inclusion():
    frame = pd.DataFrame([
        {"experiment": "H0_PIT_FIXED_INCLUSION_TRUTH", "catalogue_floor": 0.1, "n_validation": 300,
         "status": "SUCCESS", "pit": 0.4,
         "h0_pit_interval_included_68": 1, "h0_pit_interval_included_90": 1, "h0_pit_interval_included_95": 1},
        {"experiment": "H0_PIT_FIXED_INCLUSION_TRUTH", "catalogue_floor": 0.1, "n_validation": 300,
         "status": "SUCCESS", "pit": 0.6,
         "h0_pit_interval_included_68": 0, "h0_pit_interval_included_90": 1, "h0_pit_interval_included_95": 1},
    ])
    cell = MODULE._summarize(frame, bins=5)["H0_PIT_FIXED_INCLUSION_TRUTH|floor=0.1|n=300"]
    assert cell["pit_mean"] == pytest.approx(0.5)
    assert cell["prior_predictive_interval_inclusion_68"] == pytest.approx(0.5)
    assert cell["prior_predictive_interval_inclusion_90"] == pytest.approx(1.0)


def test_conditional_h0_pit_summary_is_not_labelled_joint_sbc():
    frame = pd.DataFrame([
        {"experiment": "H0_PIT_FIXED_INCLUSION_TRUTH", "catalogue_floor": 0.0,
         "n_validation": 100, "status": "SUCCESS", "pit": 0.25,
         "h0_pit_interval_included_68": 1, "h0_pit_interval_included_90": 1,
         "h0_pit_interval_included_95": 1},
        {"experiment": "H0_PIT_FIXED_INCLUSION_TRUTH", "catalogue_floor": 0.0,
         "n_validation": 100, "status": "SUCCESS", "pit": 0.75,
         "h0_pit_interval_included_68": 0, "h0_pit_interval_included_90": 1,
         "h0_pit_interval_included_95": 1},
    ])
    summary = MODULE._summarize(frame, bins=2)
    cell = summary["H0_PIT_FIXED_INCLUSION_TRUTH|floor=0|n=100"]
    assert cell["pit_mean"] == pytest.approx(0.5)
    assert cell["pit_bin_counts"] == [1, 1]


def test_all_attempt_binary_bounds_treat_failed_and_non_evaluable_as_unknown():
    frame = pd.DataFrame([
        {"status": "SUCCESS", "covered_68": 1},
        {"status": "SUCCESS", "covered_68": 0},
        {"status": "NON_EVALUABLE"},
        {"status": "FAILED"},
        {"status": "NOT_RUN"},
    ])
    lower, upper = MODULE._attempt_level_binary_bounds(frame, "covered_68")
    assert lower == pytest.approx(1 / 4)
    assert upper == pytest.approx(3 / 4)

def test_paired_validation_size_contrasts_match_replicates_and_exclude_unevaluable_pairs():
    frame = pd.DataFrame([
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "replicate": 0,
         "n_validation": 100, "status": "SUCCESS", "h0_true": 70.0,
         "h0_median": 71.0, "h0_q16": 66.0, "h0_q84": 76.0,
         "h0_q05": 62.0, "h0_q95": 82.0, "h0_q025": 60.0, "h0_q975": 85.0,
         "covered_68": 1, "covered_90": 1, "covered_95": 1},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "replicate": 0,
         "n_validation": 300, "status": "SUCCESS", "h0_true": 70.0,
         "h0_median": 70.0, "h0_q16": 66.0, "h0_q84": 74.0,
         "h0_q05": 64.0, "h0_q95": 78.0, "h0_q025": 62.0, "h0_q975": 80.0,
         "covered_68": 0, "covered_90": 1, "covered_95": 1},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "replicate": 1,
         "n_validation": 100, "status": "SUCCESS", "h0_true": 70.0,
         "h0_median": 69.0, "h0_q16": 64.0, "h0_q84": 74.0,
         "h0_q05": 60.0, "h0_q95": 80.0, "h0_q025": 58.0, "h0_q975": 82.0,
         "covered_68": 1, "covered_90": 1, "covered_95": 1},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.0, "replicate": 1,
         "n_validation": 300, "status": "NON_EVALUABLE", "h0_true": 70.0},
    ])
    result = MODULE._paired_validation_size_contrasts(frame)
    cell = result["FIXED_TRUTH|floor=0|100->300"]
    assert cell["attempted_common_replicates"] == 2
    assert cell["evaluable_pairs"] == 1
    assert cell["excluded_pairs"] == 1
    assert cell["metrics"]["posterior_median"]["mean_delta"] == pytest.approx(-1.0)
    assert cell["metrics"]["signed_median_error"]["mean_delta"] == pytest.approx(-1.0)
    assert cell["metrics"]["absolute_median_error"]["mean_delta"] == pytest.approx(-1.0)
    assert cell["metrics"]["width_68"]["mean_delta"] == pytest.approx(-2.0)
    assert cell["metrics"]["coverage_68"]["mean_delta"] == pytest.approx(-1.0)

def test_paired_contrasts_tolerate_cells_with_no_evaluable_posterior_columns():
    frame = pd.DataFrame([
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.25, "replicate": 0,
         "n_validation": 100, "status": "FAILED"},
        {"experiment": "FIXED_TRUTH", "catalogue_floor": 0.25, "replicate": 0,
         "n_validation": 300, "status": "NON_EVALUABLE"},
    ])
    result = MODULE._paired_validation_size_contrasts(frame)
    cell = result["FIXED_TRUTH|floor=0.25|100->300"]
    assert cell["attempted_common_replicates"] == 1
    assert cell["evaluable_pairs"] == 0
    assert cell["excluded_pairs"] == 1
    assert cell["metrics"] == {}

def test_prior_slope_multiplier_grid_requires_positive_unique_baseline():
    assert MODULE._parse_prior_slope_multipliers("0.5,1,2") == (0.5, 1.0, 2.0)
    for invalid in ("", "0.5,2", "1,1", "0,1", "-1,1", "nan,1", "one,1"):
        with pytest.raises(ValueError):
            MODULE._parse_prior_slope_multipliers(invalid)

