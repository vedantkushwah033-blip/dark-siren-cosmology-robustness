"""Focused tests for paired inclusion-degree quadrature sensitivity utilities."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_degree_matched_five_arm_quadrature_sensitivity import (
    _posterior, _compare_inputs, _assert_expected_config,
)


def test_independent_grid_posterior_normalizes_and_recovers_gaussian_quantiles():
    grid = np.linspace(-8.0, 8.0, 16001)
    loglike = -0.5 * ((grid - 1.25) / 0.8) ** 2
    density, quantiles = _posterior(grid, loglike)
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-10)
    assert quantiles["q500"] == pytest.approx(1.25, abs=1e-3)
    assert quantiles["q160"] == pytest.approx(1.25 - 0.9944579 * 0.8, abs=2e-3)
    assert quantiles["q840"] == pytest.approx(1.25 + 0.9944579 * 0.8, abs=2e-3)


def test_paired_input_auditor_accepts_identical_and_rejects_changed_inputs(tmp_path: Path):
    left = tmp_path / "left"
    right = tmp_path / "right"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
            "observable/validation_n100.csv",
        ):
            pd.DataFrame({"id": [1, 2], "value": [0.1, 0.2]}).to_csv(root / rel, index=False)
    assert _compare_inputs(left, right) == 6
    pd.DataFrame({"id": [1, 2], "value": [0.1, 0.3]}).to_csv(
        right / "observable/validation_n100.csv", index=False
    )
    with pytest.raises(AssertionError):
        _compare_inputs(left, right)


def test_posterior_rejects_nonfinite_likelihood():
    grid = np.linspace(0.0, 1.0, 11)
    with pytest.raises((AssertionError, ValueError)):
        _posterior(grid, np.full_like(grid, np.nan))



def test_frozen_child_config_gate_rejects_changed_simulation_setting():
    expected = {
        "seed": 20261120,
        "n_parent": 500,
        "n_detected": 4,
        "h0_prior_support": [30.0, 110.0],
    }
    _assert_expected_config(expected.copy(), expected, "matched")
    changed = expected | {"n_detected": 30}
    with pytest.raises(AssertionError, match="n_detected"):
        _assert_expected_config(changed, expected, "mismatched")
