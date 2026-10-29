from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.smooth_selection_hidden_host import (
    selection_conditioned_smooth_hidden_host_log_likelihood,
    smooth_selection_probability,
)


def _inputs():
    z = np.linspace(0.01, 0.40, 301)
    gc = np.exp(-0.5 * ((z - 0.12) / 0.05) ** 2)
    gm = np.exp(-0.5 * ((z - 0.25) / 0.07) ** 2)
    gc /= trapezoid(gc, z)
    gm /= trapezoid(gm, z)
    h0 = np.linspace(50.0, 90.0, 61)
    return z, gc, gm, h0, np.array([350.0, 620.0, 840.0])


def _run(f=0.3, **overrides):
    z, gc, gm, h0, distances = _inputs()
    args = dict(
        redshift_grid=z, catalogue_density=gc, missing_density=gm,
        H0_grid=h0, distances=distances, sigma=140.0, f_missing=f,
        d50_mpc=900.0, width_mpc=120.0, quadrature_order=32,
    )
    args.update(overrides)
    return selection_conditioned_smooth_hidden_host_log_likelihood(**args)


def test_logistic_selection_is_bounded_and_half_at_midpoint():
    p = smooth_selection_probability(
        np.array([0.0, 900.0, 3000.0]), d50_mpc=900.0, width_mpc=120.0
    )
    assert np.all((p > 0.0) & (p < 1.0))
    assert p[1] == pytest.approx(0.5)
    assert p[0] > p[1] > p[2]


def test_selection_efficiencies_are_finite_and_bounded():
    result = _run()
    for key in ("alpha_catalogue", "alpha_missing", "beta"):
        assert np.all(np.isfinite(result[key]))
        assert np.all((result[key] > 0.0) & (result[key] <= 1.0))


def test_identical_components_make_likelihood_invariant_to_fraction():
    _, gc, _, _, _ = _inputs()
    a = _run(f=0.1, missing_density=gc)
    b = _run(f=0.5, missing_density=gc)
    c = _run(f=0.9, missing_density=gc)
    np.testing.assert_allclose(a["log_likelihood"], b["log_likelihood"], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(b["log_likelihood"], c["log_likelihood"], rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize(
    "kwargs",
    [{"f_missing": -0.1}, {"f_missing": 1.1}, {"sigma": 0.0},
     {"width_mpc": 0.0}, {"quadrature_order": 3}],
)
def test_invalid_parameters_fail_closed(kwargs):
    with pytest.raises(ValueError):
        _run(**kwargs)
