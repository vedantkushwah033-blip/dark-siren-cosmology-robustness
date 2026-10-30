from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid

from dark_siren_robustness.posterior_validation_bootstrap import (
    bootstrap_validation_h0_posterior,
)
from dark_siren_robustness.posterior_validation_bootstrap_audit import (
    audit_posterior_validation_bootstrap,
)


@pytest.mark.parametrize("method", ["paired_bootstrap", "design_aware_srswor"])
def test_posterior_validation_bootstrap_replays_full_toy_pipeline(method):
    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(50.0, 90.0, 61)
    observed_z = np.linspace(0.025, 0.385, 45)
    validation_z = np.linspace(0.01, 0.40, 120)
    # Balanced outcomes with no deterministic separation in redshift.
    validation_y = (np.arange(validation_z.size) % 2).astype(int)
    distances = np.array([280.0, 330.0, 410.0, 520.0, 610.0])
    merger = 1.0 + 0.5 * (z_grid - z_grid[0]) / (z_grid[-1] - z_grid[0])

    first = bootstrap_validation_h0_posterior(
        z_grid, h0_grid, observed_z, validation_z, validation_y, distances,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
        parent_population_size=180, method=method, n_bootstrap=4,
        seed=9127, merger_weight=merger,
    )
    second = bootstrap_validation_h0_posterior(
        z_grid, h0_grid, observed_z, validation_z, validation_y, distances,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
        parent_population_size=180, method=method, n_bootstrap=4,
        seed=9127, merger_weight=merger,
    )
    assert first["schema"] == "validation_h0_posterior_resampling_v1"
    assert first["n_attempted"] == 4
    assert first["n_successful"] + first["n_failed"] == 4
    assert first["posterior_density_draws"].shape == (first["n_successful"], h0_grid.size)
    assert first["attempt_posterior_mean"].shape == (4,)
    assert len(first["sample_indices"]) == 4
    np.testing.assert_array_equal(first["success_mask"], second["success_mask"])
    np.testing.assert_allclose(
        first["attempt_posterior_mean"], second["attempt_posterior_mean"],
        equal_nan=True,
    )
    for density in first["posterior_density_draws"]:
        assert trapezoid(density, h0_grid) == pytest.approx(1.0, abs=1e-8)
    assert "not H0 posterior credible intervals" in first["uncertainty_interpretation"]
    audit = audit_posterior_validation_bootstrap(
        first, z_grid, h0_grid, observed_z, validation_z, validation_y, distances,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
        parent_population_size=180, merger_weight=merger,
    )
    assert audit["status"] == "PASS"
    assert audit["independently_reconstructed_posteriors"] == first["n_successful"]

    if method == "paired_bootstrap":
        assert first["pseudo_population_indices"] == []
        assert np.isfinite(first["fixed_full_validation_reference_mean"])
        assert np.all(first["reference_posterior_mean_by_attempt"][first["success_mask"]]
                      == first["fixed_full_validation_reference_mean"])
    else:
        assert len(first["pseudo_population_indices"]) == 4
        assert len(first["selected_pseudo_population_positions"]) == 4
        for pseudo, indices, positions in zip(
            first["pseudo_population_indices"],
            first["sample_indices"],
            first["selected_pseudo_population_positions"],
        ):
            assert pseudo.shape == (180,)
            assert indices.shape == (validation_z.size,)
            assert positions.shape == (validation_z.size,)
            assert np.unique(positions).size == validation_z.size
            np.testing.assert_array_equal(indices, pseudo[positions])
        assert np.all(np.isfinite(first["reference_posterior_mean_by_attempt"][first["success_mask"]]))


def test_posterior_validation_bootstrap_rejects_invalid_method_and_detection_sample():
    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(50.0, 90.0, 61)
    observed_z = np.linspace(0.02, 0.39, 20)
    validation_z = np.linspace(0.01, 0.40, 40)
    validation_y = (np.arange(40) % 2).astype(int)
    with pytest.raises(ValueError, match="method"):
        bootstrap_validation_h0_posterior(
            z_grid, h0_grid, observed_z, validation_z, validation_y, [300.0],
            sigma_mpc=140.0, detection_threshold_mpc=250.0,
            parent_population_size=80, method="unknown",
        )
    with pytest.raises(ValueError, match="detection threshold"):
        bootstrap_validation_h0_posterior(
            z_grid, h0_grid, observed_z, validation_z, validation_y, [200.0],
            sigma_mpc=140.0, detection_threshold_mpc=250.0,
            parent_population_size=80,
        )



def test_independent_posterior_bootstrap_auditor_rejects_tampered_posterior_mean():
    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(50.0, 90.0, 61)
    observed_z = np.linspace(0.025, 0.385, 45)
    validation_z = np.linspace(0.01, 0.40, 120)
    validation_y = (np.arange(validation_z.size) % 2).astype(int)
    distances = np.array([280.0, 330.0, 410.0, 520.0, 610.0])
    result = bootstrap_validation_h0_posterior(
        z_grid, h0_grid, observed_z, validation_z, validation_y, distances,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
        parent_population_size=180, method="paired_bootstrap",
        n_bootstrap=3, seed=713,
    )
    tampered = dict(result)
    tampered_means = np.asarray(result["attempt_posterior_mean"]).copy()
    tampered_means[0] += 1.0
    tampered["attempt_posterior_mean"] = tampered_means
    with pytest.raises(AssertionError, match="posterior mean"):
        audit_posterior_validation_bootstrap(
            tampered, z_grid, h0_grid, observed_z, validation_z, validation_y, distances,
            sigma_mpc=140.0, detection_threshold_mpc=250.0,
            parent_population_size=180,
        )
