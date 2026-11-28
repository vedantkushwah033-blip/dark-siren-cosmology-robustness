import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.hidden_population import (
    C_KM_S,
    hidden_component_log_density,
    hidden_redshift_log_density,
    log_likelihood,
    posterior_grid,
    posterior_marginals,
)


def test_truncated_normal_redshift_density_normalizes():
    z = np.linspace(0.005, 0.10, 20001)
    density = np.exp(hidden_redshift_log_density(z, 0.075, 0.012))
    assert trapezoid(density, z) == pytest.approx(1.0, abs=2e-7)


def test_hidden_component_quadrature_converges():
    u = np.log(np.array([100.0, 180.0, 260.0]))
    low = hidden_component_log_density(u, 70.0, 0.075, 0.012, 0.1, quadrature_order=48)
    high = hidden_component_log_density(u, 70.0, 0.075, 0.012, 0.1, quadrature_order=96)
    np.testing.assert_allclose(low, high, rtol=0, atol=2e-8)


def test_log_likelihood_zero_and_unit_fraction_limits():
    u = np.log(np.array([90.0, 130.0, 210.0]))
    listed = log_likelihood(u, 70, 0, 0.075, 0.012, 0.1)
    explicit_listed = np.sum(norm.logpdf(
        u[:, None],
        loc=np.log(C_KM_S * np.array([0.02, 0.04])[None, :] / 70),
        scale=0.1,
    ) + np.log([0.7, 0.3])[None, :])
    assert np.isfinite(listed)
    assert log_likelihood(u, 70, 0, 0.075, 0.012, 0.1) == pytest.approx(listed)
    hidden = log_likelihood(u, 70, 1, 0.075, 0.012, 0.1)
    assert np.isfinite(hidden)
    assert not np.isclose(listed, hidden)
    assert np.isfinite(explicit_listed)


def test_one_event_mixture_likelihood_is_between_component_limits():
    u = np.log(np.array([140.0]))
    listed_only = log_likelihood(
        u, 70, 0.0, 0.075, 0.012, 0.1,
        listed_redshifts=(0.02, 0.04), listed_weights=(0.7, 0.3),
    )
    hidden_only = log_likelihood(
        u, 70, 1.0, 0.075, 0.012, 0.1,
        listed_redshifts=(0.02, 0.04), listed_weights=(0.7, 0.3),
    )
    mixed = log_likelihood(
        u, 70, 0.37, 0.075, 0.012, 0.1,
        listed_redshifts=(0.02, 0.04), listed_weights=(0.7, 0.3),
    )
    assert min(listed_only, hidden_only) <= mixed <= max(listed_only, hidden_only)


def test_posterior_grid_normalizes_and_returns_ordered_marginals():
    rng = np.random.default_rng(20261128)
    true_h0, true_f, true_mu, true_sh, su = 70.0, 0.3, 0.075, 0.012, 0.1
    z = np.array([0.02, 0.04, true_mu])
    p = np.array([(1-true_f)*0.7, (1-true_f)*0.3, true_f])
    hosts = rng.choice(3, size=18, p=p)
    u = rng.normal(np.log(C_KM_S*z[hosts]/true_h0), su)
    posterior = posterior_grid(
        u, su,
        np.linspace(55, 85, 9), np.linspace(0, 0.7, 8),
        np.linspace(0.045, 0.095, 8), np.linspace(0.004, 0.025, 7),
        quadrature_order=32,
    )
    axes = [posterior[k] for k in ("H0_grid", "fraction_grid", "mu_grid", "sigma_h_grid")]
    w = [np.r_[np.diff(x)[0]/2, (x[2:]-x[:-2])/2, np.diff(x)[-1]/2] for x in axes]
    volume = w[0][:,None,None,None]*w[1][None,:,None,None]*w[2][None,None,:,None]*w[3][None,None,None,:]
    assert np.sum(posterior["density"]*volume) == pytest.approx(1.0, abs=1e-10)
    summaries = posterior_marginals(posterior)
    for name in ("H0", "fraction", "mu_h", "sigma_h"):
        s = summaries[name]
        assert s["lower_95"] <= s["lower_68"] <= s["median"] <= s["upper_68"] <= s["upper_95"]


@pytest.mark.parametrize("fraction", [-0.1, 1.1, np.nan])
def test_invalid_hidden_fraction_rejected(fraction):
    with pytest.raises(ValueError):
        log_likelihood(np.log([100, 120]), 70, fraction, 0.075, 0.012, 0.1)


def test_invalid_scales_rejected():
    with pytest.raises(ValueError):
        log_likelihood(np.log([100, 120]), 70, 0.3, 0.075, 0, 0.1)
    with pytest.raises(ValueError):
        hidden_redshift_log_density([0.02], 0.075, 0.012, (0.10, 0.005))


def test_known_hidden_scale_can_be_conditioned_on_as_singleton_axis():
    rng = np.random.default_rng(17)
    u = rng.normal(np.log(C_KM_S * np.array([0.02, 0.04, 0.075])[np.array([0, 1, 2, 2, 1, 0])] / 70), 0.1)
    posterior = posterior_grid(
        u, 0.1,
        np.linspace(55, 85, 7), np.linspace(0, 0.7, 6),
        np.linspace(0.045, 0.095, 6), [0.012],
        quadrature_order=24,
    )
    marginal = posterior_marginals(posterior)
    assert marginal["sigma_h"]["median"] == pytest.approx(0.012)
    assert marginal["sigma_h"]["lower_95"] == pytest.approx(0.012)
    assert marginal["sigma_h"]["upper_95"] == pytest.approx(0.012)
    assert 55 <= marginal["H0"]["median"] <= 85
