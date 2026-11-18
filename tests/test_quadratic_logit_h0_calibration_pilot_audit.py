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
    _assert_nested_validation_samples,
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



def test_nested_validation_audit_ignores_size_dependent_frame_indicator(tmp_path):
    dirs = {}
    rows = []
    base_ids = [10, 11, 12, 13, 14]
    for nval, selected in [(2, [1, 1, 0, 0, 0]), (4, [1, 1, 1, 1, 0]), (5, [1, 1, 1, 1, 1])]:
        cell = f"n{nval}"
        root = tmp_path / cell
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        pd.DataFrame({"object_id":[1,2],"redshift":[.1,.2]}).to_csv(
            root / "observable/analysis_catalogue.csv",index=False)
        pd.DataFrame({"event_id":[1],"distance":[100.]}).to_csv(
            root / "observable/gw_events_observed.csv",index=False)
        pd.DataFrame({"frame_id":base_ids,"redshift":[.1,.2,.3,.4,.5],
                      "selected_for_validation":selected}).to_csv(
            root / "truth/validation_frame_truth.csv",index=False)
        pd.DataFrame({"validation_id":base_ids[:nval],
                      "redshift":[.1,.2,.3,.4,.5][:nval]}).to_csv(
            root / "observable/validation_observed.csv",index=False)
        dirs[cell] = root
        rows.append({"cell_id":cell,"replicate":0,"gamma":0.0,"n_validation":nval})
    ledger = pd.DataFrame(rows)
    assert _assert_nested_validation_samples(ledger, dirs) == 1

    bad_ids = pd.read_csv(dirs["n4"] / "observable/validation_observed.csv")
    bad_ids.loc[0,"validation_id"] = 999
    bad_ids.to_csv(dirs["n4"] / "observable/validation_observed.csv",index=False)
    with pytest.raises(AssertionError,match="not nested prefixes"):
        _assert_nested_validation_samples(ledger, dirs)



def test_legacy_summary_without_not_run_is_accepted_only_when_zero():
    ledger, summary = fixtures()
    legacy = summary.drop(columns=["not_run"])
    _assert_summary_matches_ledger(ledger, legacy)

    ledger.loc[2, "status"] = "NOT_RUN"
    summary.loc[0, "non_evaluable"] = 0
    with pytest.raises(AssertionError, match="missing required count column: not_run"):
        _assert_summary_matches_ledger(ledger, legacy)
