import numpy as np
import pytest

from dark_siren_robustness.calibration import summarize_interval_performance


def test_perfect_estimates_and_coverage_are_summarized():
    truth = np.array([1.0, 2.0, 3.0])
    result = summarize_interval_performance(
        truth, truth, truth - 0.5, truth + 0.5
    )
    assert result["bias"] == pytest.approx(0.0)
    assert result["rmse"] == pytest.approx(0.0)
    assert result["mean_interval_width"] == pytest.approx(1.0)
    assert result["empirical_coverage"] == pytest.approx(1.0)
    assert result["coverage_mc_se"] == pytest.approx(0.0)


def test_missed_intervals_reduce_coverage():
    result = summarize_interval_performance(
        [1.0, 2.0, 3.0],
        [1.1, 2.2, 3.3],
        [0.9, 2.1, 2.0],
        [1.2, 2.3, 2.5],
    )
    # Only the first interval contains its corresponding true value.
    assert result["empirical_coverage"] == pytest.approx(1 / 3)
    assert result["bias"] == pytest.approx(0.2)


@pytest.mark.parametrize(
    "args",
    [
        ([], [], [], []),
        ([1.0], [1.0, 2.0], [0.0], [2.0]),
        ([1.0], [np.nan], [0.0], [2.0]),
        ([1.0], [1.0], [2.0], [1.0]),
        ([[1.0]], [[1.0]], [[0.0]], [[2.0]]),
    ],
)
def test_invalid_metric_inputs_raise(args):
    with pytest.raises(ValueError):
        summarize_interval_performance(*args)
