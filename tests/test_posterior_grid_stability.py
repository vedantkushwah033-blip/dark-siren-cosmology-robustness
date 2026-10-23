import numpy as np

from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary


def _posterior(grid, prior_min=50.0, prior_max=90.0):
    return posterior_h0(
        grid,
        redshifts=[0.012, 0.018, 0.026, 0.037, 0.051],
        gw_distance=165.0,
        gw_sigma=22.0,
        weights=[1.0, 1.5, 0.7, 2.0, 0.8],
        missing_weight=0.12,
        missing_redshifts=[0.009, 0.021, 0.044],
        missing_weights=[0.3, 0.5, 0.2],
        prior_min=prior_min,
        prior_max=prior_max,
    )


def test_hardened_posterior_grid_refinement_is_stable_for_fixed_inputs():
    coarse_grid = np.linspace(50.0, 90.0, 1601)
    fine_grid = np.linspace(50.0, 90.0, 3201)
    coarse = posterior_summary(coarse_grid, _posterior(coarse_grid))
    fine = posterior_summary(fine_grid, _posterior(fine_grid))
    assert abs(coarse["median"] - fine["median"]) < 0.02
    for key in ("lower_68", "upper_68", "lower_95", "upper_95"):
        assert abs(coarse[key] - fine[key]) < 0.05


def test_wider_grid_with_same_prior_matches_baseline():
    baseline_grid = np.linspace(50.0, 90.0, 1601)
    wider_grid = np.linspace(30.0, 110.0, 6401)
    baseline = posterior_summary(baseline_grid, _posterior(baseline_grid))
    wider = posterior_summary(wider_grid, _posterior(wider_grid))
    assert abs(baseline["median"] - wider["median"]) < 0.03
    for key in ("lower_68", "upper_68", "lower_95", "upper_95"):
        assert abs(baseline[key] - wider[key]) < 0.08
