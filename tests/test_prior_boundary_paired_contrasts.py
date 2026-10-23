import numpy as np

from scripts.analyze_prior_boundary_sensitivity import (
    ci_mean_difference,
    ci_rmse_difference,
)


def test_paired_mean_difference_bootstrap_constant_delta():
    delta = np.full(20, 2.5)
    mean, low, high = ci_mean_difference(delta, np.random.default_rng(7), 200)
    assert mean == 2.5
    assert low == 2.5
    assert high == 2.5


def test_paired_rmse_difference_sign_and_interval():
    baseline = np.array([-1.0, 1.0, -2.0, 2.0])
    alternative = np.array([-3.0, 3.0, -4.0, 4.0])
    mean, low, high = ci_rmse_difference(
        alternative, baseline, np.random.default_rng(9), 500
    )
    assert mean > 0
    assert low > 0
    assert high > 0
