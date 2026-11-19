"""Regression tests for independent H0 pilot summary reconstruction."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_quadratic_logit_h0_calibration_pilot import (
    _assert_summary_matches_ledger,
    _assert_validation_samples_paired,
    _assert_nested_validation_frames,
    _completion_status,
)


def fixtures():
    ledger = pd.DataFrame([
        {"gamma":0.0,"n_validation":100,"status":"SUCCESS","median_bias":-2.0,
         "covered68":1,"covered90":1,"covered95":1,"width68":8.0,
         "min_catalogue_kish_ess":12.0,"max_normalized_weight":0.2},
        {"gamma":0.0,"n_validation":100,"status":"SUCCESS","median_bias":4.0,
         "covered68":0,"covered90":1,"covered95":1,"width68":10.0,
         "min_catalogue_kish_ess":10.0,"max_normalized_weight":0.3},
        {"gamma":0.0,"n_validation":100,"status":"NON_EVALUABLE","median_bias":float("nan"),
         "covered68":float("nan"),"covered90":float("nan"),"covered95":float("nan"),
         "width68":float("nan"),"min_catalogue_kish_ess":float("nan"),
         "max_normalized_weight":float("nan")},
    ])
    summary = pd.DataFrame([{
        "gamma":0.0,"n_validation":100,"planned":3,"success":2,"non_evaluable":1,
        "failed":0,"not_run":0,"bias":1.0,"rmse":10.0**0.5,"coverage68":0.5,
        "coverage90":1.0,"coverage95":1.0,"mean_width68":9.0,
        "median_min_kish_ess":11.0,"median_max_normalized_weight":0.25,
    }])
    return ledger, summary


def test_summary_is_recomputed_from_ledger():
    ledger, summary = fixtures()
    _assert_summary_matches_ledger(ledger, summary)


@pytest.mark.parametrize("column,bad_value", [
    ("bias",1.1),("rmse",2.0),("coverage68",1.0),("mean_width68",8.0),
    ("median_min_kish_ess",99.0),("median_max_normalized_weight",0.01),
    ("success",1),("non_evaluable",0),
])
def test_tampered_summary_is_rejected(column,bad_value):
    ledger, summary = fixtures()
    summary.loc[0,column] = bad_value
    with pytest.raises(AssertionError):
        _assert_summary_matches_ledger(ledger, summary)


def test_summary_rejects_missing_and_duplicate_cells():
    ledger, summary = fixtures()
    with pytest.raises(AssertionError,match="cell keys"):
        _assert_summary_matches_ledger(ledger,summary.iloc[0:0])
    with pytest.raises(AssertionError,match="duplicate"):
        _assert_summary_matches_ledger(ledger,pd.concat([summary,summary],ignore_index=True))



def test_validation_observations_must_match_across_curvature(tmp_path):
    import numpy as np
    from pathlib import Path

    dirs = {}
    rows = []
    for gamma in (-1.0, 0.0, 1.0):
        cell = f"g{gamma}"
        root = tmp_path / cell
        (root / "observable").mkdir(parents=True)
        pd.DataFrame({"validation_id":[1,2,3],"redshift":[.1,.2,.3],
                      "retained":[1,0,1]}).to_csv(
            root / "observable/validation_observed.csv", index=False)
        dirs[cell] = root
        rows.append({"cell_id":cell,"replicate":0,"gamma":gamma,"n_validation":100})
    ledger = pd.DataFrame(rows)
    assert _assert_validation_samples_paired(ledger, dirs) == 1

    # Retention is the outcome being varied by gamma and must not be
    # compared as part of the shared SRS design.
    changed = pd.read_csv(dirs["g1.0"] / "observable/validation_observed.csv")
    changed.loc[0,"retained"] = 0
    changed.to_csv(dirs["g1.0"] / "observable/validation_observed.csv",index=False)
    assert _assert_validation_samples_paired(ledger, dirs) == 1

    changed.loc[0,"validation_id"] = 999
    changed.to_csv(dirs["g1.0"] / "observable/validation_observed.csv",index=False)
    with pytest.raises(AssertionError):
        _assert_validation_samples_paired(ledger, dirs)


def _nested_frame_fixture(tmp_path, *, tamper_core=False, tamper_prefix=False):
    import numpy as np

    dirs = {}
    rows = []
    base_ids = [0, 2, 1, 4, 3]
    redshift = np.array([0.05, 0.10, 0.20, 0.30, 0.38])
    probability = np.array([0.8, 0.7, 0.6, 0.5, 0.4])
    retained = np.array([1, 0, 1, 1, 0])
    for n in (2, 3, 5):
        cell = f"n{n}"
        root = tmp_path / cell
        (root / "truth").mkdir(parents=True)
        (root / "observable").mkdir()
        ids = base_ids[:n]
        frame = pd.DataFrame({
            "frame_id": np.arange(5),
            "redshift": redshift,
            "inclusion_probability": probability,
            "retained": retained,
            "selected_for_validation": np.isin(np.arange(5), ids).astype(int),
        })
        if tamper_core and n == 3:
            frame.loc[0, "retained"] = 0
        if tamper_prefix and n == 3:
            ids = [0, 1, 2]
            frame["selected_for_validation"] = np.isin(np.arange(5), ids).astype(int)
        frame.to_csv(root / "truth/validation_frame_truth.csv", index=False)
        pd.DataFrame({
            "object_id": [10, 11, 12], "redshift": [0.11, 0.22, 0.33],
        }).to_csv(root / "observable/analysis_catalogue.csv", index=False)
        pd.DataFrame({
            "object_id": [10, 11, 12, 13, 14], "redshift": [0.02, 0.08, 0.16, 0.25, 0.39],
        }).to_csv(root / "truth/analysis_parent_truth.csv", index=False)
        pd.DataFrame({
            "event_id": [0, 1], "observed_distance": [120.0, 240.0],
        }).to_csv(root / "observable/gw_events_observed.csv", index=False)
        pd.DataFrame({
            "validation_id": ids,
            "redshift": redshift[ids],
            "retained": retained[ids],
        }).to_csv(root / "observable/validation_observed.csv", index=False)
        dirs[cell] = root
        rows.append({
            "cell_id": cell, "replicate": 0, "gamma": 1.0,
            "n_validation": n,
        })
    return pd.DataFrame(rows), dirs


def test_nested_validation_audit_allows_size_specific_selection_flags(tmp_path):
    ledger, dirs = _nested_frame_fixture(tmp_path)
    assert _assert_nested_validation_frames(ledger, dirs) == 1


def test_nested_validation_audit_rejects_changed_underlying_frame(tmp_path):
    ledger, dirs = _nested_frame_fixture(tmp_path, tamper_core=True)
    with pytest.raises(AssertionError):
        _assert_nested_validation_frames(ledger, dirs)


def test_nested_validation_audit_rejects_non_nested_sample_ids(tmp_path):
    ledger, dirs = _nested_frame_fixture(tmp_path, tamper_prefix=True)
    with pytest.raises(AssertionError, match="not nested prefixes"):
        _assert_nested_validation_frames(ledger, dirs)


def test_nested_validation_audit_rejects_event_changed_with_sample_size(tmp_path):
    ledger, dirs = _nested_frame_fixture(tmp_path)
    changed = pd.read_csv(dirs["n3"] / "observable/gw_events_observed.csv")
    changed.loc[0, "observed_distance"] = 121.0
    changed.to_csv(dirs["n3"] / "observable/gw_events_observed.csv", index=False)
    with pytest.raises(AssertionError):
        _assert_nested_validation_frames(ledger, dirs)


def test_nested_validation_audit_rejects_catalogue_changed_with_sample_size(tmp_path):
    ledger, dirs = _nested_frame_fixture(tmp_path)
    changed = pd.read_csv(dirs["n3"] / "observable/analysis_catalogue.csv")
    changed.loc[0, "redshift"] = 0.19
    changed.to_csv(dirs["n3"] / "observable/analysis_catalogue.csv", index=False)
    with pytest.raises(AssertionError):
        _assert_nested_validation_frames(ledger, dirs)



def test_summary_accepts_schema_without_not_run_column():
    ledger, summary = fixtures()
    _assert_summary_matches_ledger(ledger, summary.drop(columns=["not_run"]))



def test_interrupted_ledger_preserves_planned_denominator_and_is_incomplete():
    ledger = pd.DataFrame([
        {"cell_id":"a","status":"SUCCESS"},
        {"cell_id":"b","status":"RUNNING"},
        {"cell_id":"c","status":"NOT_RUN"},
    ])
    assert _completion_status(ledger, planned_count=3, summary_available=False) == "INCOMPLETE"


def test_complete_ledger_requires_summary_and_reports_pass():
    ledger = pd.DataFrame([
        {"cell_id":"a","status":"SUCCESS"},
        {"cell_id":"b","status":"NON_EVALUABLE"},
        {"cell_id":"c","status":"FAILED"},
    ])
    assert _completion_status(ledger, planned_count=3, summary_available=True) == "PASS"
    with pytest.raises(AssertionError, match="missing calibration_summary"):
        _completion_status(ledger, planned_count=3, summary_available=False)


def test_missing_planned_rows_are_incomplete_even_if_attempts_are_terminal():
    ledger = pd.DataFrame([
        {"cell_id":"a","status":"SUCCESS"},
        {"cell_id":"b","status":"FAILED"},
    ])
    assert _completion_status(ledger, planned_count=3, summary_available=True) == "INCOMPLETE"


def test_interpretation_complete_with_non_evaluable_cells_is_not_called_incomplete():
    from scripts.audit_quadratic_logit_h0_calibration_pilot import _interpretation
    message = _interpretation("PASS", {"SUCCESS": 595, "NON_EVALUABLE": 5})
    assert "Complete planned attempt ledger" in message
    assert "5 NON_EVALUABLE" in message
    assert "not joint SBC" in message
    assert "INCOMPLETE" not in message


def test_interpretation_incomplete_ledger_is_explicit():
    from scripts.audit_quadratic_logit_h0_calibration_pilot import _interpretation
    message = _interpretation("INCOMPLETE", {"SUCCESS": 10, "NOT_RUN": 2})
    assert message.startswith("INCOMPLETE:")
    assert "unattempted or running" in message
