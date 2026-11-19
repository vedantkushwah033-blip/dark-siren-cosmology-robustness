"""Focused tests for quadratic five-arm replicate aggregation."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_quadratic_five_arm_replicates import (
    ARMS, _aggregate, _paired_contrasts, _seed,
)


def test_child_seed_schedule_is_deterministic_and_distinct():
    first=[(_seed(2026111929,r,0),_seed(2026111929,r,1)) for r in range(20)]
    second=[(_seed(2026111929,r,0),_seed(2026111929,r,1)) for r in range(20)]
    assert first==second
    flattened=[x for pair in first for x in pair]
    assert len(set(flattened))==40


def test_quadratic_five_arm_aggregate_and_paired_contrasts():
    summaries=[]
    ledger=[]
    errors={
        "A_catalogue_only":-4.0,
        "B_parent_oracle":1.0,
        "C_known_probability_ipw":-1.0,
        "D_map_plugin_ipw":2.0,
        "E_uncertainty_integrated_ipw":1.5,
    }
    for arm in ARMS:
        e=errors[arm]
        summaries.append({
            "replicate":0,"gamma":0.0,"validation_size":100,"arm":arm,
            "status":"SUCCESS","h0_true":70.0,"median_error":e,
            "posterior_median":70.0+e,"q16":65.0,"q84":75.0,
        })
        ledger.append({"replicate":0,"gamma":0.0,"validation_size":100,
            "arm":arm,"status":"SUCCESS"})
        ledger.append({"replicate":1,"gamma":0.0,"validation_size":100,
            "arm":arm,"status":"FAILED"})
    summary=pd.DataFrame(summaries)
    lifecycle=pd.DataFrame(ledger)
    contrasts=_paired_contrasts(summary)
    assert len(contrasts)==8
    assert set(contrasts.reference_arm)=={"A_catalogue_only","C_known_probability_ipw"}
    aggregate=_aggregate(summary,lifecycle,replicates=2,max_abs_error=20.0)
    row=aggregate.loc[aggregate.arm=="A_catalogue_only"].iloc[0]
    assert row.successful_replicates==1
    assert row.failed_or_non_evaluable==1
    assert row.success_only_bias==-4.0
    assert row.success_only_rmse==4.0
    assert row.failure_inclusive_rmse_lower==np.sqrt(8.0)
    assert row.failure_inclusive_rmse_upper==np.sqrt(208.0)
    assert row.failure_inclusive_coverage_lower==0.5
    assert row.failure_inclusive_coverage_upper==1.0
    assert row.failure_inclusive_width_lower==5.0
    assert row.failure_inclusive_width_upper==25.0


def test_all_failed_cells_remain_in_aggregate_denominator():
    ledger=pd.DataFrame([
        {"replicate":r,"gamma":g,"validation_size":n,"arm":a,"status":"FAILED"}
        for r in (0,1) for g in (-1.0,0.0) for n in (100,300) for a in ARMS
    ])
    summary=pd.DataFrame(columns=[
        "replicate","gamma","validation_size","arm","status","h0_true",
        "median_error","q16","q84",
    ])
    aggregate=_aggregate(summary,ledger,replicates=2,max_abs_error=20.0)
    assert len(aggregate)==2*2*len(ARMS)
    assert (aggregate.successful_replicates==0).all()
    assert (aggregate.failed_or_non_evaluable==2).all()
    assert aggregate.success_only_bias.isna().all()
