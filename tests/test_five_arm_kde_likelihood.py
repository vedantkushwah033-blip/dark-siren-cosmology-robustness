import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import matched_five_arm_kde_selected_log_likelihood


def _fixture():
    parent = np.linspace(0.0101, 0.3999, 600)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 180)
    retained = (np.arange(validation.size) % 3 != 0).astype(int)
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
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
    rate = 1 + 1.5 * (z - 0.01) / 0.39
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
