import importlib.util
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_known_host_calibration_surface.py"
SPEC = importlib.util.spec_from_file_location("run_known_host_calibration_surface", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(median, lo68, hi68, lo95, hi95, sd=2.0):
    return {
        "median": median,
        "mean": median,
        "posterior_sd": sd,
        "lower_68": lo68,
        "upper_68": hi68,
        "lower_95": lo95,
        "upper_95": hi95,
    }


def test_wilson_interval_has_expected_symmetry_for_half_successes():
    lo, hi = MODULE.wilson_interval(50, 100)
    assert lo == pytest.approx(0.4038, abs=0.001)
    assert hi == pytest.approx(0.5962, abs=0.001)


def test_wilson_interval_rejects_invalid_counts():
    with pytest.raises(ValueError, match="require n > 0"):
        MODULE.wilson_interval(0, 0)
    with pytest.raises(ValueError, match="require n > 0"):
        MODULE.wilson_interval(4, 3)


def test_summary_recomputes_coverage_bias_width_and_edge_fraction():
    rows = [
        row(68.0, 64.0, 72.0, 58.0, 78.0),
        row(70.0, 66.0, 74.0, 60.0, 80.0),
        row(70.0, 66.0, 74.0, 60.0, 80.0),
        row(76.0, 72.0, 80.0, 68.0, 86.0),
    ]
    result = MODULE.summarize_cell(rows, true_h0=70.0, fractional_error=0.1)
    assert result["n"] == 4
    assert result["coverage_68"] == pytest.approx(0.75)
    assert result["coverage_95"] == pytest.approx(1.0)
    assert result["median_bias"] == pytest.approx(0.0)
    assert result["mean_median_error"] == pytest.approx(1.0)
    assert result["mean_width_68"] == pytest.approx(8.0)
    assert result["mean_width_95"] == pytest.approx(20.0)
    assert result["truth_on_prior_boundary"] is False


def test_summary_rejects_non_nested_intervals():
    # The 95% interval must extend at least as far as the 68% interval
    # on both sides; test the lower-bound and upper-bound violations.
    invalid_rows = [
        row(70.0, 65.0, 75.0, 66.0, 76.0),  # lower 95% bound is too high
        row(70.0, 65.0, 75.0, 64.0, 74.0),  # upper 95% bound is too low
    ]
    for rows in ([candidate] for candidate in invalid_rows):
        with pytest.raises(ValueError, match="nested"):
            MODULE.summarize_cell(rows, true_h0=70.0, fractional_error=0.1)


def test_summary_flags_truth_on_prior_boundary():
    rows = [row(50.0, 50.0, 51.0, 50.0, 52.0)]
    result = MODULE.summarize_cell(rows, true_h0=50.0, fractional_error=0.05)
    assert result["truth_on_prior_boundary"] is True


@pytest.mark.parametrize(
    "candidate",
    [
        row(49.9, 50.0, 51.0, 50.0, 52.0),
        row(70.0, 65.0, 75.0, 49.9, 76.0),
        row(70.0, 65.0, 75.0, 64.0, 90.1),
        row(70.0, 65.0, 75.0, 64.0, 76.0, sd=0.0),
        row(70.0, 65.0, 75.0, 64.0, 76.0, sd=-1.0),
    ],
)
def test_summary_rejects_out_of_support_summaries_and_nonpositive_sd(candidate):
    with pytest.raises(ValueError):
        MODULE.summarize_cell([candidate], true_h0=70.0, fractional_error=0.1)


@pytest.mark.parametrize("fractional_error", [0.0, -0.1, float("nan"), float("inf")])
def test_summary_rejects_invalid_fractional_error(fractional_error):
    with pytest.raises(ValueError, match="truth must be finite"):
        MODULE.summarize_cell(
            [row(70.0, 65.0, 75.0, 64.0, 76.0)],
            true_h0=70.0,
            fractional_error=fractional_error,
        )


def test_summary_rejects_nonfinite_truth():
    with pytest.raises(ValueError, match="truth must be finite"):
        MODULE.summarize_cell(
            [row(70.0, 65.0, 75.0, 64.0, 76.0)],
            true_h0=float("nan"),
            fractional_error=0.1,
        )
