from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_h0_known_inclusion_oracle_sensitivity import (
    _paired, assert_cross_gamma_pairing, make_h0_grid, selected_posterior_quantiles,
)


def inputs():
    z = np.array([0.04, 0.08, 0.12, 0.18, 0.25])
    weights = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    distances = np.array([450.0, 620.0])
    h0 = np.linspace(30.0, 150.0, 241)
    return z, weights, distances, h0


def test_oracle_posterior_quantiles_are_ordered_and_finite():
    z, weights, distances, h0 = inputs()
    q = selected_posterior_quantiles(z, weights, distances, h0, 100.0, 300.0)
    values = [q[k] for k in ("q025", "q050", "q160", "q500", "q840", "q950", "q975")]
    assert np.isfinite(values).all()
    assert np.all(np.diff(values) > 0)
    assert h0[0] <= values[0] < values[-1] <= h0[-1]


def test_oracle_posterior_is_invariant_to_common_weight_scaling():
    z, weights, distances, h0 = inputs()
    a = selected_posterior_quantiles(z, weights, distances, h0, 100.0, 300.0)
    b = selected_posterior_quantiles(z, 37.0 * weights, distances, h0, 100.0, 300.0)
    assert a == pytest.approx(b, abs=1e-12)


def test_oracle_rejects_events_outside_the_declared_detection_region():
    z, weights, distances, h0 = inputs()
    with pytest.raises(ValueError, match="detected-event"):
        selected_posterior_quantiles(z, weights, np.array([250.0]), h0, 100.0, 300.0)



def test_oracle_grid_can_be_refined_beyond_archived_config():
    grid = make_h0_grid(20.0, 300.0, configured_points=141, requested_points=1121)
    assert grid.size == 1121
    assert grid[0] == pytest.approx(20.0)
    assert grid[-1] == pytest.approx(300.0)
    assert np.diff(grid).mean() == pytest.approx(0.25)
    assert make_h0_grid(20.0, 300.0, 141).size == 141


def test_oracle_grid_rejects_invalid_resolution():
    with pytest.raises(ValueError, match="invalid H0 grid"):
        make_h0_grid(20.0, 300.0, configured_points=141, requested_points=1)



def test_cross_gamma_pairing_allows_only_expected_selection_changes():
    parent_a = __import__("pandas").DataFrame({
        "object_id": [1, 2], "redshift": [0.1, 0.2],
        "inclusion_probability": [0.2, 0.3],
    })
    parent_b = parent_a.copy()
    parent_b["inclusion_probability"] = [0.4, 0.5]
    events = __import__("pandas").DataFrame({
        "event_id": [1], "distance_mpc": [500.0],
    })
    validation_a = __import__("pandas").DataFrame({
        "validation_id": [1, 2], "redshift": [0.1, 0.2], "retained": [0, 1],
    })
    validation_b = validation_a.copy()
    validation_b["retained"] = [1, 0]
    assert_cross_gamma_pairing(
        parent_a, parent_b, events, events, validation_a, validation_b
    )


def test_cross_gamma_pairing_rejects_changed_event_or_validation_stream():
    import pandas as pd

    parent = pd.DataFrame({"object_id": [1, 2], "redshift": [0.1, 0.2]})
    events = pd.DataFrame({"event_id": [1], "distance_mpc": [500.0]})
    validation = pd.DataFrame({
        "validation_id": [1, 2], "redshift": [0.1, 0.2], "retained": [0, 1],
    })
    changed_events = events.copy()
    changed_events.loc[0, "distance_mpc"] += 1.0
    with pytest.raises(AssertionError, match="events differ"):
        assert_cross_gamma_pairing(parent, parent, events, changed_events, validation, validation)
    changed_validation = validation.copy()
    changed_validation.loc[0, "redshift"] += 1e-6
    with pytest.raises(AssertionError, match="validation IDs/redshifts"):
        assert_cross_gamma_pairing(parent, parent, events, events, validation, changed_validation)


def test_paired_oracle_contrast_excludes_non_evaluable_cells():
    result = _paired([1.0, float("nan"), 3.0])
    assert result["n"] == 2
    assert result["mean"] == pytest.approx(2.0)
    assert result["mcse"] == pytest.approx(1.0)
    assert result["ci95_low"] < result["mean"] < result["ci95_high"]
