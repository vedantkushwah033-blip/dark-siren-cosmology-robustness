import numpy as np
import pytest

from dark_siren_robustness.completeness import distance_dependent_completeness


def test_completeness_is_bounded_and_decreases_with_distance():
    probabilities = distance_dependent_completeness(
        np.array([0.0, 100.0, 500.0]), scale_mpc=200.0, floor=0.05
    )
    assert np.all((probabilities >= 0.0) & (probabilities <= 1.0))
    assert np.all(np.diff(probabilities) < 0)


def test_completeness_floor_is_respected():
    probabilities = distance_dependent_completeness([1e8], scale_mpc=100.0, floor=0.2)
    assert probabilities[0] == pytest.approx(0.2)


def test_invalid_inputs_raise():
    with pytest.raises(ValueError):
        distance_dependent_completeness([-1.0], scale_mpc=100.0)
    with pytest.raises(ValueError):
        distance_dependent_completeness([10.0], scale_mpc=0.0)
