import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import logsumexp

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import (\n    arm_kde_selected_log_likelihood,\n    matched_five_arm_kde_selected_log_likelihood,\n)


def _fixture():
    parent = np.linspace(0.0101, 0.3999, 600)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 180)
    retained = (np.arange(validation.size) % 3 != 0).astype(int)
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
    )
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
    ]
    return arms, events


def test_kde_five_arm_likelihoods_use_normalized_densities_and_common_grids():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 501)
    h0 = np.linspace(55, 85, 41)
    rate = np.ones_like(z)
    result = matched_five_arm_kde_selected_log_likelihood(
        h0, arms, events, detection_threshold_mpc=250,
        redshift_grid=z, merger_weight=rate, scott_multiplier=0.7,
    )
    assert result.arm_a_log_likelihood.shape == h0.shape
    assert result.arm_e_node_log_likelihood.shape == (arms.arm_e_node_weights.shape[0], h0.size)
    for density in [*result.arm_density.values(), *result.arm_e_node_density]:
        assert np.all(np.isfinite(density))
        assert np.all(density >= 0)
        assert trapezoid(density, z) == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.isfinite(result.arm_e_marginal_log_likelihood))
    positive = arms.arm_e_quadrature_weights > 0
    expected_e = logsumexp(
        result.arm_e_node_log_likelihood[positive]
        + np.log(arms.arm_e_quadrature_weights[positive, None]),
        axis=0,
    )
    np.testing.assert_allclose(result.arm_e_marginal_log_likelihood, expected_e)
    assert result.arm_e_node_log_detection_efficiency.shape == (
        arms.arm_e_node_weights.shape[0], len(events), h0.size
    )


def test_kde_likelihood_rejects_invalid_bandwidth_and_rate():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 101)
    with pytest.raises(ValueError, match="scott_multiplier"):
        matched_five_arm_kde_selected_log_likelihood(
            [60, 70], arms, events, detection_threshold_mpc=250,
            redshift_grid=z, merger_weight=np.ones_like(z), scott_multiplier=0,
        )



def test_kde_propagates_complete_base_and_inverse_probability_weights():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 501)
    factor = 0.7
    result = matched_five_arm_kde_selected_log_likelihood(
        np.linspace(55, 85, 31), arms, events,
        detection_threshold_mpc=250, redshift_grid=z,
        merger_weight=np.ones_like(z), scott_multiplier=factor,
    )
    expected_a = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_a_catalogue_only,
        scott_multiplier=factor,
    )
    expected_b = estimate_host_density(
        z, arms.parent_redshift, sample_weights=arms.arm_b_parent_oracle,
        scott_multiplier=factor,
    )
    expected_c = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_c_known_probability,
        scott_multiplier=factor,
    )
    expected_d = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_d_map_plugin,
        scott_multiplier=factor,
    )
    np.testing.assert_allclose(result.arm_density["A_catalogue_only"], expected_a)
    np.testing.assert_allclose(result.arm_density["B_parent_oracle"], expected_b)
    np.testing.assert_allclose(result.arm_density["C_known_probability_ipw"], expected_c)
    np.testing.assert_allclose(result.arm_density["D_map_plugin_ipw"], expected_d)
    expected_e0 = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_e_node_weights[0],
        scott_multiplier=factor,
    )
    np.testing.assert_allclose(result.arm_e_node_density[0], expected_e0)




def test_kde_likelihood_applies_external_host_rate_to_inclusion_only_arm_weights():
    parent = np.linspace(0.0101, 0.3999, 600)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 180)
    retained = (np.arange(validation.size) % 3 != 0).astype(int)
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
    )

    z = np.linspace(0.01, 0.40, 501)
    rate = 1.0 + 1.5 * (z - 0.01) / 0.39
    h0 = np.linspace(55, 85, 31)
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
    ]
    result = matched_five_arm_kde_selected_log_likelihood(
        h0, arms, events, detection_threshold_mpc=250,
        redshift_grid=z, merger_weight=rate, absolute_bandwidth_sd=0.02,
    )

    expected_a = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_a_catalogue_only,
        merger_weight=rate, absolute_bandwidth_sd=0.02,
    )
    expected_b = estimate_host_density(
        z, arms.parent_redshift, sample_weights=arms.arm_b_parent_oracle,
        merger_weight=rate, absolute_bandwidth_sd=0.02,
    )
    expected_c = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_c_known_probability,
        merger_weight=rate, absolute_bandwidth_sd=0.02,
    )
    expected_d = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_d_map_plugin,
        merger_weight=rate, absolute_bandwidth_sd=0.02,
    )
    np.testing.assert_allclose(result.arm_density["A_catalogue_only"], expected_a)
    np.testing.assert_allclose(result.arm_density["B_parent_oracle"], expected_b)
    np.testing.assert_allclose(result.arm_density["C_known_probability_ipw"], expected_c)
    np.testing.assert_allclose(result.arm_density["D_map_plugin_ipw"], expected_d)
    expected_e = np.stack([
        estimate_host_density(
            z, arms.catalogue_redshift, sample_weights=weights,
            merger_weight=rate, absolute_bandwidth_sd=0.02,
        )
        for weights in arms.arm_e_node_weights
    ])
    np.testing.assert_allclose(result.arm_e_node_density, expected_e)
    assert np.all(np.isfinite(result.arm_e_marginal_log_likelihood))


def test_isolated_arm_likelihoods_reconstruct_all_matched_arm_outputs():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 501)
    h0 = np.linspace(55.0, 85.0, 41)
    rate = np.ones_like(z)
    common = dict(
        H0_grid=h0, events=events, detection_threshold_mpc=250.0,
        redshift_grid=z, merger_weight=rate, scott_multiplier=0.7,
    )
    matched = matched_five_arm_kde_selected_log_likelihood(
        h0, arms, events, detection_threshold_mpc=250.0,
        redshift_grid=z, merger_weight=rate, scott_multiplier=0.7,
    )
    one_dimensional = [
        (arms.catalogue_redshift, arms.arm_a_catalogue_only, matched.arm_a_log_likelihood),
        (arms.parent_redshift, arms.arm_b_parent_oracle, matched.arm_b_log_likelihood),
        (arms.catalogue_redshift, arms.arm_c_known_probability, matched.arm_c_log_likelihood),
        (arms.catalogue_redshift, arms.arm_d_map_plugin, matched.arm_d_log_likelihood),
    ]
    for sample, weights, expected in one_dimensional:
        result = arm_kde_selected_log_likelihood(
            sample_redshift=sample, raw_weights=weights, **common,
        )
        assert result.component_density.shape == (1, z.size)
        np.testing.assert_allclose(result.component_log_likelihood[0], expected, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(result.marginal_log_likelihood, expected, rtol=1e-11, atol=1e-11)

    integrated = arm_kde_selected_log_likelihood(
        sample_redshift=arms.catalogue_redshift,
        raw_weights=arms.arm_e_node_weights,
        mixture_weights=arms.arm_e_quadrature_weights,
        **common,
    )
    np.testing.assert_allclose(
        integrated.component_log_likelihood, matched.arm_e_node_log_likelihood,
        rtol=1e-11, atol=1e-11,
    )
    np.testing.assert_allclose(
        integrated.marginal_log_likelihood, matched.arm_e_marginal_log_likelihood,
        rtol=1e-11, atol=1e-11,
    )


def test_isolated_arm_failure_does_not_poison_other_arm_calculation():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 201)
    common = dict(
        H0_grid=np.linspace(55.0, 85.0, 21), events=events,
        detection_threshold_mpc=250.0, redshift_grid=z,
        merger_weight=np.ones_like(z), scott_multiplier=0.7,
    )
    valid = arm_kde_selected_log_likelihood(
        sample_redshift=arms.catalogue_redshift,
        raw_weights=arms.arm_a_catalogue_only, **common,
    )
    assert np.all(np.isfinite(valid.marginal_log_likelihood))
    with pytest.raises(ValueError, match="raw_weights"):
        arm_kde_selected_log_likelihood(
            sample_redshift=arms.catalogue_redshift,
            raw_weights=np.full(arms.catalogue_redshift.size, np.nan), **common,
        )
    replay = arm_kde_selected_log_likelihood(
        sample_redshift=arms.catalogue_redshift,
        raw_weights=arms.arm_a_catalogue_only, **common,
    )
    np.testing.assert_array_equal(valid.marginal_log_likelihood, replay.marginal_log_likelihood)


def test_isolated_arm_requires_explicit_mixture_weights_for_node_matrix():
    arms, events = _fixture()
    z = np.linspace(0.01, 0.40, 201)
    with pytest.raises(ValueError, match="mixture_weights"):
        arm_kde_selected_log_likelihood(
            np.linspace(55.0, 85.0, 21), arms.catalogue_redshift,
            arms.arm_e_node_weights, events, detection_threshold_mpc=250.0,
            redshift_grid=z, merger_weight=np.ones_like(z), scott_multiplier=0.7,
        )
