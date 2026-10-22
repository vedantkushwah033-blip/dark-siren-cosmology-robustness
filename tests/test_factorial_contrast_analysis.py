import importlib.util
import random
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "analyze_factorial_contrasts.py"
SPEC = importlib.util.spec_from_file_location("analyze_factorial_contrasts", MODULE_PATH)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def sample_row(median=71.0, mean=70.5, lower68=68.0, upper68=72.0,
               lower95=66.0, upper95=74.0, posterior_sd=2.0, pull_median=""):
    return {
        "median": str(median), "mean": str(mean),
        "lower_68": str(lower68), "upper_68": str(upper68),
        "lower_95": str(lower95), "upper_95": str(upper95),
        "posterior_sd": str(posterior_sd), "pull_median": str(pull_median),
    }


def test_coverage_and_width_metrics():
    row = sample_row()
    assert analysis.value(row, "coverage_68") == 1.0
    assert analysis.value(row, "coverage_95") == 1.0
    assert analysis.value(row, "width_68") == 4.0
    assert analysis.value(row, "width_95") == 8.0


def test_error_and_pull_metrics():
    row = sample_row(median=72.0, mean=71.0, posterior_sd=2.0)
    assert analysis.value(row, "median_error") == 2.0
    assert analysis.value(row, "squared_median_error") == 4.0
    assert analysis.value(row, "pull_median") == 1.0


def test_bootstrap_difference_is_b_minus_a_and_reproducible():
    a = [0.0, 1.0, 2.0]
    b = [2.0, 3.0, 4.0]
    first = analysis.bootstrap_mean_difference(a, b, random.Random(123), 500)
    second = analysis.bootstrap_mean_difference(a, b, random.Random(123), 500)
    assert first == second
    assert first["difference_B_minus_A"] == pytest.approx(2.0)
    assert first["n_A"] == 3
    assert first["n_B"] == 3
    assert len(first["bootstrap_95_interval"]) == 2


def test_percentile_interpolates():
    assert analysis.percentile([0.0, 10.0], 0.25) == pytest.approx(2.5)
