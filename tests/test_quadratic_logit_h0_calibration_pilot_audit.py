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
