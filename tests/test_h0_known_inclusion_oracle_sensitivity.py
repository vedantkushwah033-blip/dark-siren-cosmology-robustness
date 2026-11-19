from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_h0_known_inclusion_oracle_sensitivity import (
    assert_complete_planned_ledger, assert_cross_gamma_pairing,
    failure_inclusive_coverage_bounds, make_h0_grid,
    selected_posterior_quantiles, terminal_status_counts,
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



def test_complete_planned_ledger_accepts_exact_complete_factorial():
    import pandas as pd

    planned = pd.DataFrame([
        {"cell_id": "r0_gm1_n100", "replicate": 0, "gamma": -1.0,
         "n_validation": 100, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
        {"cell_id": "r0_gm1_n500", "replicate": 0, "gamma": -1.0,
         "n_validation": 500, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
        {"cell_id": "r0_gp0_n100", "replicate": 0, "gamma": 0.0,
         "n_validation": 100, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
        {"cell_id": "r0_gp0_n500", "replicate": 0, "gamma": 0.0,
         "n_validation": 500, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
    ])
    ledger = planned.copy()
    ledger["status"] = "SUCCESS"
    assert_complete_planned_ledger(planned, ledger)


def test_complete_planned_ledger_rejects_missing_cell_even_if_remaining_succeed():
    import pandas as pd

    planned = pd.DataFrame([
        {"cell_id": "a", "replicate": 0, "gamma": -1.0,
         "n_validation": 100, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
        {"cell_id": "b", "replicate": 0, "gamma": -1.0,
         "n_validation": 500, "analysis_seed": 10, "validation_seed": 20,
         "h0_true": 70.0},
    ])
    with pytest.raises(AssertionError, match="planned-cell inventory"):
        assert_complete_planned_ledger(planned, planned.iloc[:1].copy())


def test_complete_planned_ledger_rejects_changed_pairing_seed():
    import pandas as pd

    planned = pd.DataFrame([{
        "cell_id": "a", "replicate": 0, "gamma": -1.0,
        "n_validation": 100, "analysis_seed": 10, "validation_seed": 20,
        "h0_true": 70.0,
    }])
    ledger = planned.copy()
    ledger.loc[0, "validation_seed"] = 21
    with pytest.raises(AssertionError, match="planned-cell inventory"):
        assert_complete_planned_ledger(planned, ledger)


def test_complete_planned_ledger_rejects_duplicate_rows():
    import pandas as pd

    planned = pd.DataFrame([{
        "cell_id": "a", "replicate": 0, "gamma": -1.0,
        "n_validation": 100, "analysis_seed": 10, "validation_seed": 20,
        "h0_true": 70.0,
    }])
    with pytest.raises(AssertionError, match="duplicate"):
        assert_complete_planned_ledger(planned, pd.concat([planned, planned]))



def test_terminal_status_counts_accepts_complete_failure_inclusive_ledger():
    import pandas as pd

    ledger = pd.DataFrame({"status": ["SUCCESS"] * 9 + ["NON_EVALUABLE"]})
    counts = terminal_status_counts(ledger)
    assert counts["SUCCESS"] == 9
    assert counts["NON_EVALUABLE"] == 1
    assert counts["FAILED"] == counts["NOT_RUN"] == 0


def test_terminal_status_counts_rejects_running_or_unknown_states():
    import pandas as pd

    with pytest.raises(AssertionError, match="non-terminal"):
        terminal_status_counts(pd.DataFrame({"status": ["SUCCESS", "RUNNING"]}))
    with pytest.raises(AssertionError, match="non-terminal"):
        terminal_status_counts(pd.DataFrame({"status": ["SUCCESS", "MYSTERY"]}))


def test_failure_inclusive_coverage_bounds_include_non_success_denominator():
    lower, upper = failure_inclusive_coverage_bounds(
        covered_success=50, successful=99, planned=100
    )
    assert lower == pytest.approx(0.50)
    assert upper == pytest.approx(0.51)


def test_failure_inclusive_coverage_bounds_reject_invalid_counts():
    with pytest.raises(ValueError):
        failure_inclusive_coverage_bounds(covered_success=4, successful=3, planned=5)
    with pytest.raises(ValueError):
        failure_inclusive_coverage_bounds(covered_success=1, successful=6, planned=5)



def test_cross_gamma_pairing_accepts_different_validation_sizes_with_common_prefix():
    import pandas as pd

    parent = pd.DataFrame({"object_id": [1, 2], "redshift": [0.1, 0.2]})
    events = pd.DataFrame({"event_id": [1], "distance_mpc": [500.0]})
    short = pd.DataFrame({
        "validation_id": [1, 2], "redshift": [0.1, 0.2], "retained": [0, 1],
    })
    long = pd.DataFrame({
        "validation_id": [1, 2, 3], "redshift": [0.1, 0.2, 0.3],
        "retained": [1, 0, 1],
    })
    assert_cross_gamma_pairing(parent, parent, events, events, short, long)


def test_cross_gamma_pairing_rejects_mismatch_inside_common_validation_prefix():
    import pandas as pd

    parent = pd.DataFrame({"object_id": [1, 2], "redshift": [0.1, 0.2]})
    events = pd.DataFrame({"event_id": [1], "distance_mpc": [500.0]})
    short = pd.DataFrame({
        "validation_id": [1, 2], "redshift": [0.1, 0.2], "retained": [0, 1],
    })
    long = pd.DataFrame({
        "validation_id": [1, 9, 3], "redshift": [0.1, 0.2, 0.3],
        "retained": [1, 0, 1],
    })
    with pytest.raises(AssertionError, match="validation IDs/redshifts"):
        assert_cross_gamma_pairing(parent, parent, events, events, short, long)
