import numpy as np
import pytest

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_likelihood import matched_five_arm_selected_log_likelihood
from dark_siren_robustness.shared_inclusion_likelihood import (
    shared_inclusion_selected_log_likelihood,
)


def _arms():
    parent = np.linspace(0.01, 0.40, 300)
    catalogue = parent[::3]
    validation = np.linspace(0.01, 0.40, 120)
    retained = (np.arange(validation.size) % 3 != 0).astype(np.int8)
    return build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=5,
    )


def _events():
    return [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
        {"gw_distance": 420.0, "gw_sigma": 140.0},
    ]


def test_five_arm_likelihoods_share_grid_and_e_uses_shared_node_marginalization():
    arms = _arms()
    grid = np.linspace(55.0, 85.0, 31)
    events = _events()
    result = matched_five_arm_selected_log_likelihood(
        grid, arms, events, detection_threshold_mpc=250.0
    )
    for value in (
        result.arm_a_log_likelihood,
        result.arm_b_log_likelihood,
        result.arm_c_log_likelihood,
        result.arm_d_log_likelihood,
        result.arm_e_marginal_log_likelihood,
    ):
        assert value.shape == grid.shape
        assert np.all(np.isfinite(value))
    assert result.arm_e_node_log_likelihood.shape == (
        arms.arm_e_node_weights.shape[0], grid.size
    )
    assert result.arm_e_node_log_detection_efficiency.shape == (
        arms.arm_e_node_weights.shape[0], len(events), grid.size
    )

    reference = shared_inclusion_selected_log_likelihood(
        grid,
        arms.catalogue_redshift,
        events,
        arms.inclusion_posterior,
        detection_threshold_mpc=250.0,
    )
    assert np.allclose(
        result.arm_e_marginal_log_likelihood,
        reference.marginal_log_likelihood,
        rtol=1e-11,
        atol=1e-11,
    )


def test_five_arm_likelihood_rejects_event_below_detection_threshold():
    with pytest.raises(ValueError, match="meet the detection threshold"):
        matched_five_arm_selected_log_likelihood(
            [60.0, 70.0],
            _arms(),
            [{"gw_distance": 200.0, "gw_sigma": 140.0}],
            detection_threshold_mpc=250.0,
        )


def test_five_arm_likelihood_rejects_invalid_h0_grid():
    with pytest.raises(ValueError, match="strictly increasing"):
        matched_five_arm_selected_log_likelihood(
            [70.0, 60.0], _arms(), _events(),
            detection_threshold_mpc=250.0,
        )
