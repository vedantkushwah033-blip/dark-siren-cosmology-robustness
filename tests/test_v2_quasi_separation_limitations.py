"""Adversarial V1/V2 comparison for quadratic near-separation."""
import numpy as np
import pytest

from dark_siren_robustness.estimated_inclusion import (
    fit_quadratic_logistic_inclusion_model,
    predict_inclusion_probability,
)
from dark_siren_robustness.estimated_inclusion_v2 import (
    fit_quadratic_logistic_inclusion_model_v2,
    predict_inclusion_probability_v2,
)


def _near_separated_quadratic_fixture():
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 200)
    included = (z > 0.205).astype(int)
    # One discordant observation converts complete separation into a
    # deterministic near-separation stress case.
    included[99] = 1 - included[99]
    return z, included, z_min, z_max


def test_v1_and_v2_quadratic_fits_agree_on_near_separated_bounded_solution():
    z, included, z_min, z_max = _near_separated_quadratic_fixture()
    v1 = fit_quadratic_logistic_inclusion_model(
        z, included, z_min=z_min, z_max=z_max
    )
    v2 = fit_quadratic_logistic_inclusion_model_v2(
        z, included, z_min=z_min, z_max=z_max
    )

    assert v1["fit_status"] == "boundary_limited"
    assert v2["fit_status"] == "boundary_limited"
    assert v1["coefficient_bound_hit"] and v2["coefficient_bound_hit"]
    assert v1["projected_gradient_inf_norm"] <= 1e-5
    assert v2["projected_gradient_inf_norm"] <= 1e-7

    # Compare the common estimand (fitted probability curve), not just
    # coefficient vectors: rank-deficient designs can make coefficients
    # non-unique even when predictions are identified.
    grid = np.linspace(z_min, z_max, 301)
    p1 = predict_inclusion_probability(grid, v1)
    p2 = predict_inclusion_probability_v2(grid, v2)
    assert np.allclose(p1, p2, rtol=2e-4, atol=2e-5)
    assert v1["negative_log_likelihood"] == pytest.approx(
        v2["negative_log_likelihood"], rel=1e-8, abs=1e-8
    )
