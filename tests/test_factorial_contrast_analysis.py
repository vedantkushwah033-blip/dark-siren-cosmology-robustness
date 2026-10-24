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


def test_coverage_monte_carlo_standard_error():
    assert analysis.coverage_mc_se(0.68, 500) == pytest.approx((0.68 * 0.32 / 500) ** 0.5)
    assert analysis.coverage_mc_se(0.0, 500) == 0.0
    with pytest.raises(ValueError):
        analysis.coverage_mc_se(0.5, 0)


def test_wilson_interval_is_bounded_and_handles_extremes():
    low = analysis.wilson_interval(0, 500)
    middle = analysis.wilson_interval(340, 500)
    high = analysis.wilson_interval(500, 500)
    assert low[0] == 0.0 and low[1] > 0.0
    assert middle[0] < 0.68 < middle[1]
    assert high[1] == 1.0 and high[0] < 1.0
    with pytest.raises(ValueError):
        analysis.wilson_interval(2, 0)
    with pytest.raises(ValueError):
        analysis.wilson_interval(501, 500)

def test_percentile_rejects_empty_nonfinite_and_invalid_quantile():
    with pytest.raises(ValueError, match="at least one"):
        analysis.percentile([], 0.5)
    with pytest.raises(ValueError, match="between 0 and 1"):
        analysis.percentile([1.0, 2.0], 1.1)
    with pytest.raises(ValueError, match="finite"):
        analysis.percentile([1.0, float("nan")], 0.5)


@pytest.mark.parametrize(
    "a,b,repetitions",
    [([], [1.0], 100), ([1.0], [], 100), ([1.0], [2.0], 0)],
)
def test_bootstrap_rejects_empty_groups_or_nonpositive_repetitions(a, b, repetitions):
    with pytest.raises(ValueError):
        analysis.bootstrap_mean_difference(a, b, random.Random(1), repetitions)


def test_bootstrap_rejects_nonfinite_values():
    with pytest.raises(ValueError, match="finite"):
        analysis.bootstrap_mean_difference([1.0], [float("inf")], random.Random(1), 100)
