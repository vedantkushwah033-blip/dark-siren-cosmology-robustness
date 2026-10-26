import numpy as np
import pytest

from dark_siren_robustness.cosmological_inference import posterior_h0_from_distance_grid


def test_cosmological_inference_normalizes_extreme_finite_host_weights():
    """Equivalent host mixtures must agree even when the raw weight sum overflows."""
    h0 = np.linspace(40.0, 100.0, 601)
    distances = np.column_stack((3000.0 / h0, 5200.0 / h0))
    common = dict(
        h0_grid=h0,
        predicted_distance_mpc=distances,
        observed_distance_mpc=48.0,
        distance_sigma_mpc=8.0,
        prior_min=40.0,
        prior_max=100.0,
    )

    ordinary = posterior_h0_from_distance_grid(
        **common, host_weights=[2.0, 1.0]
    )
    extreme = posterior_h0_from_distance_grid(
        **common, host_weights=[1.0e308, 5.0e307]
    )

    np.testing.assert_allclose(extreme, ordinary, rtol=2e-13, atol=2e-13)


def test_cosmological_inference_rejects_all_zero_host_weights():
    h0 = np.linspace(40.0, 100.0, 11)
    distances = np.column_stack((3000.0 / h0, 5200.0 / h0))

    with pytest.raises(ValueError, match="host_weights must be finite"):
        posterior_h0_from_distance_grid(
            h0,
            distances,
            observed_distance_mpc=48.0,
            distance_sigma_mpc=8.0,
            host_weights=[0.0, 0.0],
        )


def test_low_z_inference_normalizes_extreme_finite_host_weights():
    """The baseline mixture should accept finite weights whose raw sum overflows."""
    from dark_siren_robustness.inference import posterior_h0_low_z

    h0 = np.linspace(50.0, 90.0, 401)
    common = dict(
        h0_grid=h0,
        redshifts=[0.02, 0.04],
        observed_distance_mpc=100.0,
        distance_sigma_mpc=15.0,
        prior_min=50.0,
        prior_max=90.0,
    )
    ordinary = posterior_h0_low_z(**common, host_weights=[2.0, 1.0])
    extreme = posterior_h0_low_z(**common, host_weights=[1.0e308, 5.0e307])

    np.testing.assert_allclose(extreme, ordinary, rtol=2e-13, atol=2e-13)


def test_low_z_posterior_normalization_is_independent_of_out_of_prior_grid_points():
    """Adding zero-density grid points outside the prior must not change the posterior."""
    from dark_siren_robustness.inference import posterior_h0_low_z

    common = dict(
        redshifts=[0.02, 0.04],
        observed_distance_mpc=100.0,
        distance_sigma_mpc=15.0,
        prior_min=50.0,
        prior_max=90.0,
        host_weights=[2.0, 1.0],
    )
    in_prior_grid = np.linspace(50.0, 90.0, 401)
    wider_grid = np.linspace(40.0, 100.0, 601)

    narrow = posterior_h0_low_z(h0_grid=in_prior_grid, **common)
    wide = posterior_h0_low_z(h0_grid=wider_grid, **common)
    wide_support = (wider_grid >= 50.0) & (wider_grid <= 90.0)

    np.testing.assert_allclose(wide[wide_support], narrow, rtol=2e-13, atol=2e-13)
    assert np.all(wide[~wide_support] == 0.0)
