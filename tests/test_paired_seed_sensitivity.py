import math

import pytest

from scripts.run_paired_seed_sensitivity import mean_se


def test_mean_se_known_values():
    mean, se = mean_se([1.0, 2.0, 3.0])
    assert mean == pytest.approx(2.0)
    assert se == pytest.approx(math.sqrt(1.0 / 3.0))


def test_mean_se_rejects_empty_input():
    with pytest.raises(ValueError, match="empty"):
        mean_se([])


def test_mean_se_rejects_non_finite_values():
    with pytest.raises(ValueError, match="finite"):
        mean_se([1.0, float("nan")])


def test_mean_se_single_value_has_nan_standard_error():
    mean, se = mean_se([4.0])
    assert mean == pytest.approx(4.0)
    assert math.isnan(se)
