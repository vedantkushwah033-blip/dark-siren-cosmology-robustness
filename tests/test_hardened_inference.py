import numpy as np
import pytest

from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary


def test_valid_posterior_is_finite_and_normalized():
    grid = np.linspace(50, 90, 1601)
    posterior = posterior_h0(grid, [0.02, 0.03, 0.04], 150.0, 20.0)
    assert np.all(np.isfinite(posterior))
    assert np.all(posterior >= 0)
    assert np.trapezoid(posterior, grid) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("grid", [[70.0], [70.0, 69.0], [70.0, np.nan], [[60.0, 70.0]]])
def test_invalid_h0_grid_rejected(grid):
    with pytest.raises(ValueError):
        posterior_h0(grid, [0.02], 100.0, 10.0)


@pytest.mark.parametrize("sigma", [0.0, -1.0, np.nan, np.inf])
def test_invalid_distance_uncertainty_rejected(sigma):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], [0.02], 100.0, sigma)


@pytest.mark.parametrize("z", [[np.nan], [np.inf], [-0.1], [[0.1]]])
def test_invalid_redshift_inputs_rejected(z):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], z, 100.0, 10.0)


@pytest.mark.parametrize("weights", [[1, np.nan], [1, np.inf], [-1, 2], [0, 0], [[1, 2]]])
def test_invalid_catalogue_weights_rejected(weights):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], [0.02, 0.03], 100.0, 10.0, weights=weights)


def test_summary_rejects_shape_mismatch_and_zero_area():
    with pytest.raises(ValueError, match="same length"):
        posterior_summary([50, 60, 70], [0.1, 0.2])
    with pytest.raises(ValueError, match="positive integrated area"):
        posterior_summary([50, 60, 70], [0, 0, 0])


def test_summary_normalizes_valid_density_and_returns_ordered_intervals():
    grid = np.linspace(50, 90, 401)
    density = np.exp(-0.5 * ((grid - 70) / 4) ** 2)
    summary = posterior_summary(grid, density)
    assert 50 < summary["lower_95"] < summary["lower_68"] < summary["median"]
    assert summary["median"] < summary["upper_68"] < summary["upper_95"] < 90


def test_missing_host_component_requires_support():
    with pytest.raises(ValueError, match="missing_redshifts"):
        posterior_h0([50, 70, 90], [0.02], 100.0, 10.0, missing_weight=0.2)
