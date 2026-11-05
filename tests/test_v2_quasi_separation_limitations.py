"""Characterization tests for V2 near-separation limitations.

These tests deliberately preserve a known limitation rather than treating
non-convergence as a successful estimator outcome. The quadratic fixture is
near-separated and V1 reaches its bounded optimum, while the current V2
damped-Newton implementation exhausts its iteration budget. This is a
release-gating characterization, not a claim that V2 is generally inferior.
"""
import numpy as np
import pytest

from dark_siren_robustness.estimated_inclusion import (
    fit_quadratic_logistic_inclusion_model,
)
from dark_siren_robustness.estimated_inclusion_v2 import (
    fit_quadratic_logistic_inclusion_model_v2,
)


def _near_separated_quadratic_fixture():
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 200)
    included = (z > 0.205).astype(int)
    # A single discordant observation turns complete separation into a
    # near-separation stress case while preserving both outcome classes.
    included[99] = 1 - included[99]
    return z, included, z_min, z_max


def test_v1_quadratic_near_separation_reaches_bounded_kkt_solution():
    z, included, z_min, z_max = _near_separated_quadratic_fixture()
    fit = fit_quadratic_logistic_inclusion_model(
        z, included, z_min=z_min, z_max=z_max
    )
    assert fit["fit_status"] == "boundary_limited"
    assert fit["coefficient_bound_hit"]
    assert fit["projected_gradient_inf_norm"] <= 1e-5
    assert np.isfinite(fit["negative_log_likelihood"])


def test_v2_quadratic_near_separation_nonconvergence_is_explicitly_reported():
    z, included, z_min, z_max = _near_separated_quadratic_fixture()
    # Current V2 reaches its 500-iteration limit with projected KKT residual
    # above 1e-7. It must raise rather than label this fit as converged.
    with pytest.raises(RuntimeError, match="did not converge"):
        fit_quadratic_logistic_inclusion_model_v2(
            z, included, z_min=z_min, z_max=z_max
        )
