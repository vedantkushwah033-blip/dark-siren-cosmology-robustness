import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.inference import posterior_h0_low_z


def test_single_host_recovers_known_h0_in_narrow_distance_likelihood():
    true_h0 = 70.0
    z = 0.02
    observed_distance = 299_792.458 * z / true_h0
    grid = np.linspace(50.0, 90.0, 8001)
    posterior = posterior_h0_low_z(
        grid,
        [z],
        observed_distance,
        distance_sigma_mpc=0.5,
    )
    assert grid[np.argmax(posterior)] == pytest.approx(true_h0, abs=0.02)
    assert trapezoid(posterior, grid) == pytest.approx(1.0, abs=1e-8)


def test_host_weights_are_normalization_invariant():
    grid = np.linspace(50.0, 90.0, 2001)
    args = (grid, [0.015, 0.025], 100.0, 15.0)
    a = posterior_h0_low_z(*args, host_weights=[1.0, 3.0])
    b = posterior_h0_low_z(*args, host_weights=[10.0, 30.0])
    assert np.allclose(a, b)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"h0_grid": [70.0, 60.0], "redshifts": [0.02],
         "observed_distance_mpc": 80.0, "distance_sigma_mpc": 5.0},
        {"h0_grid": [60.0, 70.0], "redshifts": [-0.02],
         "observed_distance_mpc": 80.0, "distance_sigma_mpc": 5.0},
        {"h0_grid": [60.0, 70.0], "redshifts": [0.02],
         "observed_distance_mpc": 80.0, "distance_sigma_mpc": 0.0},
        {"h0_grid": [60.0, 70.0], "redshifts": [0.02, 0.03],
         "observed_distance_mpc": 80.0, "distance_sigma_mpc": 5.0,
         "host_weights": [1.0]},
        {"h0_grid": [60.0, 70.0], "redshifts": [0.02],
         "observed_distance_mpc": 80.0, "distance_sigma_mpc": 5.0,
         "host_weights": [0.0]},
    ],
)
def test_invalid_inputs_raise(kwargs):
    with pytest.raises(ValueError):
        posterior_h0_low_z(**kwargs)


def test_prior_excludes_grid_points_outside_bounds():
    grid = np.linspace(40.0, 100.0, 601)
    posterior = posterior_h0_low_z(
        grid, [0.02], 299_792.458 * 0.02 / 70.0, 2.0,
        prior_min=60.0, prior_max=80.0,
    )
    assert np.all(posterior[(grid < 60.0) | (grid > 80.0)] == 0.0)
