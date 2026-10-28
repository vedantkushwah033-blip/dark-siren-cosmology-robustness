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


def test_thinning_is_seeded_and_respects_probability_endpoints():
    from dark_siren_robustness.completeness import thin_catalogue

    probabilities = np.array([0.0, 1.0, 0.0, 1.0, 0.4, 0.7])
    first = thin_catalogue(probabilities, np.random.default_rng(314159))
    second = thin_catalogue(probabilities, np.random.default_rng(314159))
    assert np.array_equal(first, second)
    assert np.array_equal(first[:4], np.array([False, True, False, True]))


def test_thinning_empirical_rate_matches_known_inclusion_probability():
    from dark_siren_robustness.completeness import thin_catalogue

    n = 100_000
    probability = 0.37
    mask = thin_catalogue(
        np.full(n, probability), np.random.default_rng(20261028)
    )
    observed_rate = mask.mean()
    standard_error = np.sqrt(probability * (1.0 - probability) / n)
    assert abs(observed_rate - probability) < 5.0 * standard_error


def test_thinning_accepts_vector_probabilities_and_rejects_invalid_inputs():
    from dark_siren_robustness.completeness import thin_catalogue

    mask = thin_catalogue(
        [0.0, 1.0, 0.25, 0.75], np.random.default_rng(7)
    )
    assert mask.shape == (4,)
    assert mask[0] == np.False_
    assert mask[1] == np.True_
    for invalid in ([np.nan], [-0.1], [1.1], [[0.2, 0.3]], []):
        with pytest.raises(ValueError):
            thin_catalogue(invalid, np.random.default_rng(1))
    with pytest.raises(TypeError, match="Generator"):
        thin_catalogue([0.5], np.random.RandomState(1))

