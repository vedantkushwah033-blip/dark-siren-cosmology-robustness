import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.dataset_likelihood import (
    dataset_posterior_h0,
    event_log_likelihood,
)
from dark_siren_robustness.hardened_inference import posterior_h0


GRID = np.linspace(40.0, 100.0, 3001)
PRIOR = {"prior_min": 50.0, "prior_max": 90.0}


def event(distance, sigma, redshifts=(0.02, 0.04), **kwargs):
    return {
        "redshifts": list(redshifts),
        "gw_distance": distance,
        "gw_sigma": sigma,
        **kwargs,
    }


def normalize_on_prior(log_likelihood):
    support = (GRID >= PRIOR["prior_min"]) & (GRID <= PRIOR["prior_max"])
    density = np.zeros_like(GRID)
    density[support] = np.exp(log_likelihood[support] - np.max(log_likelihood[support]))
    density[support] /= trapezoid(density[support], GRID[support])
    return density


def test_single_event_matches_existing_event_level_posterior():
    one = event(175.0, 18.0, weights=[0.7, 0.3])
    expected = posterior_h0(
        GRID, one["redshifts"], one["gw_distance"], one["gw_sigma"],
        weights=one["weights"], **PRIOR
    )
    actual = dataset_posterior_h0(GRID, [one], **PRIOR)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)


def test_two_event_joint_log_likelihood_is_sum_of_event_log_likelihoods():
    events = [event(175.0, 18.0), event(210.0, 25.0, redshifts=(0.03, 0.05))]
    actual = np.zeros_like(GRID)
    for item in events:
        actual += event_log_likelihood(GRID, **item)
    expected = event_log_likelihood(GRID, **events[0]) + event_log_likelihood(GRID, **events[1])
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def test_event_order_does_not_change_dataset_posterior():
    first = event(175.0, 18.0)
    second = event(210.0, 25.0, redshifts=(0.03, 0.05))
    ab = dataset_posterior_h0(GRID, [first, second], **PRIOR)
    ba = dataset_posterior_h0(GRID, [second, first], **PRIOR)
    np.testing.assert_allclose(ab, ba, rtol=2e-12, atol=2e-12)


def test_duplicating_event_multiplies_its_likelihood_but_not_prior():
    one = event(175.0, 18.0)
    log_l = event_log_likelihood(GRID, **one)
    expected = normalize_on_prior(2.0 * log_l)
    actual = dataset_posterior_h0(GRID, [one, one], **PRIOR)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)


def test_zero_missing_mass_is_exact_catalogue_only_limit():
    base = event(175.0, 18.0, weights=[0.7, 0.3])
    with_zero = {
        **base,
        "missing_weight": 0.0,
        "missing_redshifts": [0.06, 0.08],
        "missing_weights": [0.25, 0.75],
    }
    np.testing.assert_allclose(
        event_log_likelihood(GRID, **base),
        event_log_likelihood(GRID, **with_zero),
        rtol=0,
        atol=1e-13,
    )


def test_positive_truncation_matches_event_level_posterior():
    item = event(175.0, 18.0, positive_truncated=True)
    expected = posterior_h0(
        GRID, item["redshifts"], item["gw_distance"], item["gw_sigma"],
        positive_truncated=True, **PRIOR
    )
    actual = dataset_posterior_h0(GRID, [item], **PRIOR)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)


def test_empty_dataset_rejected():
    with pytest.raises(ValueError, match="At least one event"):
        dataset_posterior_h0(GRID, [], **PRIOR)


def test_missing_component_requires_support():
    with pytest.raises(ValueError, match="missing_redshifts"):
        event_log_likelihood(GRID, **event(175.0, 18.0, missing_weight=0.2))


def test_missing_host_mixture_matches_existing_event_level_posterior():
    item = event(
        175.0, 18.0, weights=[0.7, 0.3],
        missing_weight=0.25,
        missing_redshifts=[0.06, 0.08],
        missing_weights=[0.25, 0.75],
    )
    expected = posterior_h0(
        GRID, item["redshifts"], item["gw_distance"], item["gw_sigma"],
        weights=item["weights"],
        missing_weight=item["missing_weight"],
        missing_redshifts=item["missing_redshifts"],
        missing_weights=item["missing_weights"],
        **PRIOR,
    )
    actual = dataset_posterior_h0(GRID, [item], **PRIOR)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)


def test_dataset_posterior_is_zero_outside_prior_and_normalized_inside():
    items = [event(175.0, 18.0), event(210.0, 25.0, redshifts=(0.03, 0.05))]
    posterior = dataset_posterior_h0(GRID, items, **PRIOR)
    support = (GRID >= PRIOR["prior_min"]) & (GRID <= PRIOR["prior_max"])
    assert np.all(posterior[~support] == 0.0)
    assert trapezoid(posterior[support], GRID[support]) == pytest.approx(1.0, abs=2e-12)


def test_rescaling_host_weights_does_not_change_posterior():
    base = event(175.0, 18.0, weights=[0.7, 0.3])
    scaled = event(175.0, 18.0, weights=[7.0, 3.0])
    a = dataset_posterior_h0(GRID, [base], **PRIOR)
    b = dataset_posterior_h0(GRID, [scaled], **PRIOR)
    np.testing.assert_allclose(a, b, rtol=2e-12, atol=2e-12)


def test_missing_weight_one_is_rejected_by_declared_parameterization():
    with pytest.raises(ValueError, match=r"missing_weight must be in \[0, 1\)"):
        event_log_likelihood(
            GRID, redshifts=[0.02], gw_distance=175.0, gw_sigma=18.0,
            missing_weight=1.0, missing_redshifts=[0.03],
        )


def test_multi_event_posterior_matches_direct_joint_likelihood_normalization():
    items = [
        event(175.0, 18.0, weights=[0.7, 0.3]),
        event(210.0, 25.0, redshifts=(0.03, 0.05), weights=[0.4, 0.6]),
        event(260.0, 32.0, redshifts=(0.025, 0.055)),
    ]
    joint = np.zeros_like(GRID)
    for item in items:
        joint += event_log_likelihood(GRID, **item)
    expected = normalize_on_prior(joint)
    actual = dataset_posterior_h0(GRID, items, **PRIOR)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)
