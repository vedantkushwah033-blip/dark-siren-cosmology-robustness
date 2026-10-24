import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.cosmological_inference import (
    distance_grid_for_model,
    posterior_h0_from_distance_grid,
)


def test_flat_lcdm_distance_grid_matches_cosmology_reference():
    h0 = np.array([60.0, 70.0, 80.0])
    z = np.array([0.02, 0.05, 0.1])
    result = distance_grid_for_model(h0, z, "flat_lcdm", integration_steps=5000)
    assert result.shape == (3, 3)
    assert np.all(np.diff(result, axis=0) < 0)
    assert np.all(np.diff(result, axis=1) > 0)
    assert result[1] == pytest.approx(
        distance_grid_for_model([70.0, 71.0], z, "flat_lcdm", integration_steps=5000)[0],
        rel=1e-12,
    )


def test_forward_distance_model_accepts_single_truth_value():
    distances = distance_grid_for_model([70.0], [0.1], "flat_lcdm")
    assert distances.shape == (1, 1)
    assert np.isfinite(distances[0, 0])
    assert distances[0, 0] > 0


def test_low_z_model_reproduces_linear_distance_relation():
    h0 = np.array([60.0, 70.0, 80.0])
    z = np.array([0.02, 0.04])
    got = distance_grid_for_model(h0, z, "low_z")
    expected = 299792.458 * z[None, :] / h0[:, None]
    assert got == pytest.approx(expected, rel=1e-12)


def test_matched_flat_lcdm_posterior_recovers_synthetic_h0():
    h0 = np.linspace(40.0, 100.0, 12001)
    truth = 70.0
    z = np.array([0.04])
    weights = np.array([1.0])
    distances = distance_grid_for_model(h0, z, "flat_lcdm", integration_steps=4000)
    # Use a distance generated at the true H0 from the same host-mixture model.
    at_truth = distance_grid_for_model([truth], z, "flat_lcdm", integration_steps=4000)[0]
    observed = float(np.dot(weights, at_truth))
    posterior = posterior_h0_from_distance_grid(
        h0, distances, observed, 8.0, host_weights=weights, prior_min=40.0, prior_max=100.0
    )
    median = h0[np.searchsorted(np.cumsum(posterior) / np.sum(posterior), 0.5)]
    assert trapezoid(posterior, h0) == pytest.approx(1.0, abs=1e-8)
    assert median == pytest.approx(truth, abs=0.5)


def test_wrong_low_z_model_differs_at_moderate_redshift():
    h0 = np.linspace(40.0, 100.0, 6001)
    z = np.array([0.6])
    weights = np.array([1.0])
    truth = 70.0
    true_distances = distance_grid_for_model([truth], z, "flat_lcdm", integration_steps=6000)[0]
    observed = float(np.dot(weights, true_distances))
    low_z = distance_grid_for_model(h0, z, "low_z")
    lcdm = distance_grid_for_model(h0, z, "flat_lcdm", integration_steps=6000)
    p_low = posterior_h0_from_distance_grid(h0, low_z, observed, 5.0, host_weights=weights)
    p_lcdm = posterior_h0_from_distance_grid(h0, lcdm, observed, 5.0, host_weights=weights)
    median_low = h0[np.searchsorted(np.cumsum(p_low) / np.sum(p_low), 0.5)]
    median_lcdm = h0[np.searchsorted(np.cumsum(p_lcdm) / np.sum(p_lcdm), 0.5)]
    assert abs(median_low - median_lcdm) > 0.5
    assert median_lcdm == pytest.approx(truth, abs=0.5)


@pytest.mark.parametrize(
    "h0,z,model",
    [
        ([60.0, 70.0], [-0.1], "flat_lcdm"),
        ([60.0, 70.0], [0.1], "bad_model"),
    ],
)
def test_invalid_distance_grid_inputs_are_rejected(h0, z, model):
    with pytest.raises(ValueError):
        distance_grid_for_model(h0, z, model)


def test_gaussian_distance_estimator_can_be_non_positive():
    h0 = np.linspace(40.0, 100.0, 1201)
    distances = distance_grid_for_model(h0, [0.02], "flat_lcdm")
    posterior = posterior_h0_from_distance_grid(
        h0, distances, observed_distance_mpc=-2.0, distance_sigma_mpc=20.0
    )
    assert np.all(np.isfinite(posterior))
    assert trapezoid(posterior, h0) == pytest.approx(1.0, abs=1e-8)


def test_posterior_rejects_bad_weights_and_distance_shapes():
    with pytest.raises(ValueError, match="host_weights"):
        posterior_h0_from_distance_grid([50, 70, 90], np.ones((3, 2)), 100, 10, host_weights=[0, 0])
    with pytest.raises(ValueError, match="shape"):
        posterior_h0_from_distance_grid([50, 70, 90], np.ones(3), 100, 10)
