"""Regression tests for independent H0 pilot summary reconstruction."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_quadratic_logit_h0_calibration_pilot import _assert_summary_matches_ledger


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
