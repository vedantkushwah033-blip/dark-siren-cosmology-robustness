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


def test_wilson_interval_contains_observed_proportion_and_is_bounded():
    low, high = MODULE._wilson_interval(68, 100)
    assert 0.0 <= low < 0.68 < high <= 1.0
    empty = MODULE._wilson_interval(0, 0)
    assert np.isnan(empty[0]) and np.isnan(empty[1])


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
    ])
    summary = MODULE._summarize(frame, bins=10)
    cell = summary["FIXED_TRUTH|floor=0|n=100"]
    assert cell["attempted"] == 3
    assert cell["success"] == 1
    assert cell["non_evaluable"] == 1
    assert cell["failed"] == 1
    assert cell["coverage_68"] == pytest.approx(1.0)
    assert cell["mean_posterior_median_bias"] == pytest.approx(1.0)
    assert cell["median_minimum_catalogue_kish_ess"] == pytest.approx(500.0)
