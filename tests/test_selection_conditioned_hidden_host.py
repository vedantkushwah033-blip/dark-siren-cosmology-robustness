from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.selection_conditioned_hidden_host import (
    C_KM_S,
    selection_conditioned_hidden_host_log_likelihood,
)


def _inputs():
    z = np.linspace(0.01, 0.08, 501)
    g1 = np.exp(-0.5 * ((z - 0.030) / 0.008) ** 2)
    g2 = np.exp(-0.5 * ((z - 0.055) / 0.010) ** 2)
    g1 /= trapezoid(g1, z)
    g2 /= trapezoid(g2, z)
    h0 = np.linspace(55.0, 85.0, 61)
    distances = np.array([180.0, 240.0, 310.0])
    return z, g1, g2, h0, distances


def _run(f=0.3, **overrides):
    z, g1, g2, h0, d = _inputs()
    args = dict(
        redshift_grid=z, catalogue_density=g1, missing_density=g2,
        H0_grid=h0, distances=d, sigma=100.0, f_missing=f,
        d_min=0.0, d_max=500.0,
    )
    args.update(overrides)
    return selection_conditioned_hidden_host_log_likelihood(**args)


def test_mixture_selection_efficiency_is_component_weighted_and_bounded():
    result = _run(f=0.3)
    expected = 0.7 * result["alpha_catalogue"] + 0.3 * result["alpha_missing"]
    np.testing.assert_allclose(result["beta"], expected, rtol=1e-13, atol=1e-14)
    assert np.all((result["alpha_catalogue"] > 0) & (result["alpha_catalogue"] <= 1))
    assert np.all((result["alpha_missing"] > 0) & (result["alpha_missing"] <= 1))
    assert np.all((result["beta"] > 0) & (result["beta"] <= 1))


def test_f_zero_and_one_reduce_to_component_likelihoods():
    z, g1, g2, h0, d = _inputs()
    result0 = _run(f=0.0)
    result1 = _run(f=1.0)
    for f, density, result in ((0.0, g1, result0), (1.0, g2, result1)):
        direct = []
        means = C_KM_S * z[None, :] / h0[:, None]
        alpha = trapezoid(
            (norm.cdf((500.0 - means) / 100.0) - norm.cdf((0.0 - means) / 100.0))
            * density[None, :], z, axis=1,
        )
        for observed in d:
            q = trapezoid(norm.pdf(observed, loc=means, scale=100.0) * density[None, :], z, axis=1)
            direct.append(np.log(q))
        direct_ll = np.sum(direct, axis=0) - len(d) * np.log(alpha)
        np.testing.assert_allclose(result["log_likelihood"], direct_ll, rtol=2e-11, atol=2e-11)


def test_identical_component_densities_make_likelihood_invariant_to_f():
    z, g1, _, h0, d = _inputs()
    a = _run(f=0.1, missing_density=g1)
    b = _run(f=0.5, missing_density=g1)
    c = _run(f=0.9, missing_density=g1)
    np.testing.assert_allclose(a["log_likelihood"], b["log_likelihood"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(b["log_likelihood"], c["log_likelihood"], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    "overrides",
    [
        {"f_missing": -0.1},
        {"f_missing": 1.1},
        {"sigma": 0.0},
        {"d_min": 10.0, "d_max": 10.0},
        {"d_max": 200.0},
    ],
)
def test_invalid_mixture_parameters_or_unselected_data_fail_closed(overrides):
    with pytest.raises(ValueError):
        _run(**overrides)


def test_density_normalization_is_required():
    z, g1, g2, h0, d = _inputs()
    with pytest.raises(ValueError, match="integrate to one"):
        selection_conditioned_hidden_host_log_likelihood(
            z, g1 * 2, g2, h0, d, sigma=100.0, f_missing=0.3, d_max=500.0
        )


def test_selection_cutoff_changes_selection_normalization():
    a = _run(f=0.3, d_max=500.0)
    b = _run(f=0.3, d_max=1000.0)
    assert not np.allclose(a["beta"], b["beta"])
