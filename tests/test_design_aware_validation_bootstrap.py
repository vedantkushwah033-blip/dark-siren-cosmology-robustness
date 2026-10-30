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
    assert len(result["failure_records"]) == result["n_failed"]
    assert result["failure_records"] == []
    for pseudo, sample in zip(
        result["pseudo_population_indices"], result["bootstrap_sample_indices"]
    ):
        available = np.bincount(pseudo, minlength=30)
        selected = np.bincount(sample, minlength=30)
        assert np.all(selected <= available)



def test_design_aware_validation_bootstrap_indexes_failed_attempts(monkeypatch):
    import dark_siren_robustness.estimated_inclusion as module

    def forced_failure(*args, **kwargs):
        raise RuntimeError("controlled test failure")

    monkeypatch.setattr(module, "estimate_host_density_from_validation", forced_failure)
    grid = np.linspace(0.01, 0.4, 51)
    validation_z = np.linspace(0.01, 0.4, 30)
    validation_included = (np.arange(30) % 2).astype(int)
    observed_z = np.linspace(0.015, 0.395, 12)
    result = module.bootstrap_validation_inclusion_density_design_aware(
        grid,
        observed_z,
        validation_z,
        validation_included,
        parent_population_size=60,
        n_bootstrap=3,
        seed=924,
    )
    assert result["n_attempted"] == 3
    assert result["n_successful"] == 0
    assert result["n_failed"] == 3
    assert result["success_mask"].tolist() == [False, False, False]
    assert [record["bootstrap_index"] for record in result["failure_records"]] == [0, 1, 2]
    assert all(record["exception_type"] == "RuntimeError" for record in result["failure_records"])
    assert all("controlled test failure" in record["reason"] for record in result["failure_records"])
