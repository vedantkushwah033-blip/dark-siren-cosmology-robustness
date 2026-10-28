from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.selection_misspecification import run_selection_misspecification


def test_selection_misspecification_is_paired_and_reproducible():
    kwargs = dict(
        n_replicates=1, n_parent=1200, n_detected=8, seed=55129,
        h0_points=51, z_points=201,
    )
    first, summary = run_selection_misspecification(**kwargs)
    second, _ = run_selection_misspecification(**kwargs)
    assert list(first.arm) == list(second.arm)
    assert set(first.arm) == {
        "analytic_truth", "oracle", "catalogue_only", "correct_weighting",
        "shifted_transition", "shallow_slope",
    }
    assert len(first) == 6
    # The analytic-truth diagnostic must recover the declared target density
    # exactly on the grid, while the parent-sample KDE oracle need not.
    target = first.loc[first.arm == "analytic_truth", "density_ise"].iloc[0]
    assert target < 1e-20
    np.testing.assert_allclose(first.posterior_mean, second.posterior_mean)
    assert first.seed.nunique() == 1
    assert first.replicate.nunique() == 1
    assert first.posterior_integral.between(0.999999, 1.000001).all()
    assert first.density_ise.ge(0).all()
    assert first.weight_ess_fraction.between(0, 1 + 1e-9).all()
    assert summary["status"] == "synthetic_selection_misspecification_pilot_not_real_H0_measurement"
    assert summary["true_catalogue_selection"]["midpoint"] == 0.22
    assert summary["assumed_models"]["shifted_transition"]["midpoint"] == 0.27
    assert all(0 <= item["coverage_wilson_95_low"] <= item["coverage_wilson_95_high"] <= 1
               for item in summary["aggregate_metrics"])
    assert len(summary["paired_contrasts_vs_catalogue_only"]) == 5
    assert summary["paired_contrasts_vs_catalogue_only"][0]["squared_error_difference"]["standard_error"] is None


def test_selection_misspecification_rejects_invalid_grids():
    with pytest.raises(ValueError, match="grid"):
        run_selection_misspecification(n_replicates=1, h0_points=10)
