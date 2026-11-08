import numpy as np
import pytest

from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature
from dark_siren_robustness.shared_inclusion_likelihood import shared_inclusion_selected_log_likelihood


def _inputs():
    validation_z = np.linspace(0.01, 0.40, 120)
    x = (validation_z - validation_z.mean()) / validation_z.std()
    probability = 1.0 / (1.0 + np.exp(-(-0.3 + 0.6 * x)))
    retained = (np.arange(validation_z.size) % 3 != 0).astype(int)
    posterior = fit_inclusion_posterior_quadrature(
        validation_z, retained, degree=1, order=5
    )
    catalogue_z = np.linspace(0.02, 0.39, 30)
    h0 = np.linspace(55.0, 85.0, 41)
    events = [
        {"gw_distance": 700.0, "gw_sigma": 120.0},
        {"gw_distance": 1050.0, "gw_sigma": 120.0},
    ]
    return posterior, catalogue_z, h0, events


def test_shared_selected_likelihood_uses_node_specific_weights_and_marginalizes_jointly():
    posterior, catalogue_z, h0, events = _inputs()
    result = shared_inclusion_selected_log_likelihood(
        h0, catalogue_z, events, posterior,
        detection_threshold_mpc=250.0,
        population_base_weights=1.0 + 1.5 * (catalogue_z - 0.01) / 0.39,
    )
    n_nodes = posterior.nodes.shape[0]
    assert result.node_log_selected_likelihood.shape == (n_nodes, h0.size)
    assert result.node_log_detection_efficiency.shape == (n_nodes, h0.size)
    assert result.marginal_log_likelihood.shape == h0.shape
    assert result.node_catalogue_weights.shape == (catalogue_z.size, n_nodes)
    assert np.allclose(result.node_catalogue_weights.sum(axis=0), 1.0)
    assert np.all(np.isfinite(result.marginal_log_likelihood))
    assert np.all((result.node_kish_ess >= 1.0) & (result.node_kish_ess <= catalogue_z.size))
    assert np.all((result.node_maximum_normalized_weight > 0) & (result.node_maximum_normalized_weight <= 1))
    assert np.all(result.node_log_detection_efficiency <= 0.0)


def test_shared_selected_likelihood_rejects_heterogeneous_noise_for_this_selection_model():
    posterior, catalogue_z, h0, events = _inputs()
    events[1] = {"gw_distance": 1050.0, "gw_sigma": 150.0}
    with pytest.raises(ValueError, match="common gw_sigma"):
        shared_inclusion_selected_log_likelihood(
            h0, catalogue_z, events, posterior, detection_threshold_mpc=250.0
        )


def test_shared_selected_likelihood_rejects_events_below_detection_threshold():
    posterior, catalogue_z, h0, events = _inputs()
    events[0] = {"gw_distance": 200.0, "gw_sigma": 120.0}
    with pytest.raises(ValueError, match="detection threshold"):
        shared_inclusion_selected_log_likelihood(
            h0, catalogue_z, events, posterior, detection_threshold_mpc=250.0
        )
