from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

MODULE_PATH=Path(__file__).resolve().parents[1]/"scripts"/"analyze_full_parent_h0_cluster_uncertainty.py"
SPEC=importlib.util.spec_from_file_location("full_parent_cluster_uncertainty",MODULE_PATH)
MODULE=importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_cluster_analysis_uses_parent_populations_as_independent_units(tmp_path):
    root=tmp_path/"baseline"; root.mkdir()
    rows=[]
    # Strong within-parent similarity: the independent-unit count is 3,
    # not 12 event draws.
    parent_bias=[-2.0,1.0,5.0]
    for parent,bias in enumerate(parent_bias):
        for event in range(4):
            rows.append({"parent_replicate":parent,"event_replicate":event,
                "median_bias":bias+0.1*event,"covered68":int(event<3),
                "covered90":1,"covered95":1,"width68":20.0+parent})
    pd.DataFrame(rows).to_csv(root/"replicates.csv",index=False)
    out=tmp_path/"readout"
    report=MODULE.analyze(root,out,seed=11)
    cluster=report["uncertainty"]["parent_mean_median_bias"]
    expected=np.mean([np.mean([b+0.1*e for e in range(4)]) for b in parent_bias])
    assert cluster["n"]==3
    assert cluster["mean"]==pytest.approx(expected)
    assert report["parent_populations"]==3
    assert report["total_event_realizations"]==12
    assert (out/"parent_cluster_summary.csv").is_file()
    assert (out/"cluster_uncertainty.json").is_file()
    assert report["uncertainty"]["event_weighted_median_bias_conditional_on_archived_parents"]["n"]==12


def test_cluster_bootstrap_is_deterministic_and_finite():
    values=np.array([-1.0,0.5,2.0,3.0])
    a=MODULE.cluster_bootstrap_interval(values,np.random.default_rng(7),draws=2000)
    b=MODULE.cluster_bootstrap_interval(values,np.random.default_rng(7),draws=2000)
    assert a==b
    assert np.isfinite(a).all()


def test_duplicate_parent_event_keys_are_rejected(tmp_path):
    root=tmp_path/"baseline"; root.mkdir()
    row={"parent_replicate":0,"event_replicate":0,"median_bias":1.0,
         "covered68":1,"covered90":1,"covered95":1,"width68":20.0}
    pd.DataFrame([row,row]).to_csv(root/"replicates.csv",index=False)
    with pytest.raises(ValueError,match="duplicate"):
        MODULE.analyze(root,tmp_path/"out")
