import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import matched_five_arm_kde_selected_log_likelihood


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
        bandwidth_factor=factor,
    )
    expected_b = estimate_host_density(
        z, arms.parent_redshift, sample_weights=arms.arm_b_parent_oracle,
        bandwidth_factor=factor,
    )
    expected_c = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_c_known_probability,
        bandwidth_factor=factor,
    )
    expected_d = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_d_map_plugin,
        bandwidth_factor=factor,
    )
    np.testing.assert_allclose(result.arm_density["A_catalogue_only"], expected_a)
    np.testing.assert_allclose(result.arm_density["B_parent_oracle"], expected_b)
    np.testing.assert_allclose(result.arm_density["C_known_probability_ipw"], expected_c)
    np.testing.assert_allclose(result.arm_density["D_map_plugin_ipw"], expected_d)
    expected_e0 = estimate_host_density(
        z, arms.catalogue_redshift, sample_weights=arms.arm_e_node_weights[0],
        bandwidth_factor=factor,
    )
    np.testing.assert_allclose(result.arm_e_node_density[0], expected_e0)
