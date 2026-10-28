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


def test_prior_boundaries_between_grid_points_are_included_in_normalization():
    grid = np.linspace(40.0, 100.0, 301)
    lo, hi = 50.07, 89.93
    item = event(175.0, 18.0)
    posterior = dataset_posterior_h0(
        grid, [item], prior_min=lo, prior_max=hi
    )
    log_likelihood = event_log_likelihood(grid, **item)
    support = (grid >= lo) & (grid <= hi)
    interior_grid = grid[support]
    assert interior_grid.size >= 2

    # Independent reconstruction: interpolate log likelihood at both prior
    # endpoints using the full original grid, including the grid points just
    # outside prior support, then integrate across the exact prior interval.
    norm_grid = np.unique(np.concatenate(([lo], interior_grid, [hi])))
    norm_log_likelihood = np.interp(norm_grid, grid, log_likelihood)
    shift = np.max(norm_log_likelihood)
    norm_density = np.exp(norm_log_likelihood - shift)
    area = trapezoid(norm_density, norm_grid)
    norm_density /= area
    expected = np.zeros_like(grid)
    expected[support] = np.exp(log_likelihood[support] - shift) / area
    np.testing.assert_allclose(posterior, expected, rtol=2e-12, atol=2e-12)
    assert trapezoid(norm_density, norm_grid) == pytest.approx(1.0, abs=2e-10)
    assert posterior[grid < lo].sum() == 0.0
    assert posterior[grid > hi].sum() == 0.0

@pytest.mark.parametrize(
    ("grid", "prior_min", "prior_max"),
    [
        (np.linspace(60.0, 100.0, 201), 50.0, 90.0),
        (np.linspace(40.0, 80.0, 201), 50.0, 90.0),
    ],
)
def test_grid_must_span_full_prior_interval(grid, prior_min, prior_max):
    with pytest.raises(ValueError, match="must span the full prior interval"):
        dataset_posterior_h0(
            grid, [event(175.0, 18.0)],
            prior_min=prior_min, prior_max=prior_max,
        )


def test_selection_conditioning_subtracts_log_detection_efficiency():
    from dark_siren_robustness.dataset_likelihood import (
        selection_conditioned_log_likelihood,
    )

    joint = np.array([-3.0, -2.0, -4.0])
    efficiency = np.array([0.2, 0.5, 0.8])
    actual = selection_conditioned_log_likelihood(
        joint, efficiency, n_events=4
    )
    np.testing.assert_allclose(actual, joint - 4 * np.log(efficiency))


def test_constant_detection_efficiency_changes_only_likelihood_offset():
    from dark_siren_robustness.dataset_likelihood import (
        selection_conditioned_log_likelihood,
    )

    joint = np.array([-4.0, -2.0, -3.0])
    corrected = selection_conditioned_log_likelihood(
        joint, np.full(joint.size, 0.25), n_events=7
    )
    np.testing.assert_allclose(corrected - joint, np.full(3, -7 * np.log(0.25)))
    np.testing.assert_allclose(
        np.diff(corrected), np.diff(joint), rtol=0, atol=1e-13
    )


@pytest.mark.parametrize(
    "efficiency",
    [
        np.array([0.2, 0.0, 0.8]),
        np.array([0.2, -0.1, 0.8]),
        np.array([0.2, 1.1, 0.8]),
        np.array([0.2, np.nan, 0.8]),
    ],
)
def test_selection_conditioning_rejects_invalid_efficiencies(efficiency):
    from dark_siren_robustness.dataset_likelihood import (
        selection_conditioned_log_likelihood,
    )

    with pytest.raises(ValueError):
        selection_conditioned_log_likelihood(
            np.array([-3.0, -2.0, -4.0]), efficiency, n_events=3
        )


def test_selection_conditioning_rejects_shape_mismatch_and_invalid_event_count():
    from dark_siren_robustness.dataset_likelihood import (
        selection_conditioned_log_likelihood,
    )

    with pytest.raises(ValueError, match="same size"):
        selection_conditioned_log_likelihood(
            np.array([-3.0, -2.0, -4.0]), np.array([0.2, 0.8]), n_events=2
        )
    for invalid_count in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            selection_conditioned_log_likelihood(
                np.array([-3.0, -2.0, -4.0]),
                np.array([0.2, 0.5, 0.8]),
                n_events=invalid_count,
            )


def test_toy_detection_efficiency_matches_analytic_single_host_gaussian():
    from scipy.stats import norm
    from dark_siren_robustness.dataset_likelihood import event_detection_efficiency

    grid = np.array([60.0, 70.0, 80.0])
    z = 0.03
    sigma = 20.0
    threshold = 100.0
    actual = event_detection_efficiency(
        grid, redshifts=[z], gw_sigma=sigma,
        detection_distance_min=threshold,
    )
    means = 299_792.458 * z / grid
    expected = norm.sf((threshold - means) / sigma)
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-15)
    assert np.all(np.diff(actual) < 0)


def test_toy_detection_efficiency_uses_same_catalogue_missing_mixture():
    from scipy.stats import norm
    from dark_siren_robustness.dataset_likelihood import event_detection_efficiency

    grid = np.array([65.0, 75.0])
    catalog_z = np.array([0.02, 0.04])
    hidden_z = np.array([0.06])
    catalog_weights = np.array([0.25, 0.75])
    hidden_weight = 0.3
    sigma = 25.0
    threshold = 120.0
    actual = event_detection_efficiency(
        grid, redshifts=catalog_z, weights=catalog_weights,
        gw_sigma=sigma, detection_distance_min=threshold,
        missing_weight=hidden_weight, missing_redshifts=hidden_z,
        missing_weights=[1.0],
    )
    means_cat = 299_792.458 * catalog_z[None, :] / grid[:, None]
    means_hidden = 299_792.458 * hidden_z[None, :] / grid[:, None]
    expected = (1 - hidden_weight) * np.sum(
        catalog_weights[None, :] * norm.sf((threshold - means_cat) / sigma),
        axis=1,
    ) + hidden_weight * norm.sf((threshold - means_hidden[:, 0]) / sigma)
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-15)


def test_positive_truncated_detection_efficiency_conditions_on_positive_distance():
    from scipy.stats import norm
    from dark_siren_robustness.dataset_likelihood import event_detection_efficiency

    grid = np.array([60.0, 75.0])
    z = 0.01
    sigma = 35.0
    threshold = 80.0
    actual = event_detection_efficiency(
        grid, redshifts=[z], gw_sigma=sigma,
        detection_distance_min=threshold, positive_truncated=True,
    )
    means = 299_792.458 * z / grid
    expected = norm.sf((threshold - means) / sigma) / norm.sf((0.0 - means) / sigma)
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-15)


def test_toy_detection_efficiency_rejects_missing_support_and_invalid_threshold():
    from dark_siren_robustness.dataset_likelihood import event_detection_efficiency

    with pytest.raises(ValueError, match="missing_redshifts"):
        event_detection_efficiency(
            GRID, redshifts=[0.02], gw_sigma=20.0,
            detection_distance_min=100.0, missing_weight=0.2,
        )
    with pytest.raises(ValueError, match="non-negative"):
        event_detection_efficiency(
            GRID, redshifts=[0.02], gw_sigma=20.0,
            detection_distance_min=-1.0, positive_truncated=True,
        )
