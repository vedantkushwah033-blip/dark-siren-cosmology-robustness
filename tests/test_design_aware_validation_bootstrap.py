from __future__ import annotations

import numpy as np

from dark_siren_robustness.estimated_inclusion import (
    bootstrap_validation_inclusion_density_design_aware,
)


def test_design_aware_validation_bootstrap_records_both_sampling_layers():
    grid = np.linspace(0.01, 0.4, 51)
    validation_z = np.linspace(0.01, 0.4, 30)
    validation_included = (np.arange(30) % 2).astype(int)
    observed_z = np.linspace(0.015, 0.395, 12)
    result = bootstrap_validation_inclusion_density_design_aware(
        grid,
        observed_z,
        validation_z,
        validation_included,
        parent_population_size=60,
        n_bootstrap=4,
        seed=923,
    )
    assert result["schema"] == "validation_inclusion_density_pseudopop_srswor_bootstrap_v1"
    assert result["pseudo_population_indices"].shape == (4, 60)
    assert result["bootstrap_sample_indices"].shape == (4, 30)
    assert result["selected_pseudo_population_positions"].shape == (4, 30)
    assert all(
        np.unique(row).size == row.size
        for row in result["selected_pseudo_population_positions"]
    )
    assert result["pseudo_population_density_draws"].shape[1:] == (grid.size,)
    assert result["pseudo_population_intercept_draws"].shape == (result["n_successful"],)
    assert result["pseudo_population_slope_draws"].shape == (result["n_successful"],)
    assert result["success_mask"].shape == (4,)
    assert result["n_attempted"] == 4
    assert result["n_successful"] + result["n_failed"] == 4
    for pseudo, sample in zip(
        result["pseudo_population_indices"], result["bootstrap_sample_indices"]
    ):
        available = np.bincount(pseudo, minlength=30)
        selected = np.bincount(sample, minlength=30)
        assert np.all(selected <= available)
