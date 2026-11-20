# tests for prior-predictive SBC summaries
import json
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from summarize_prior_predictive_sbc import summarize

def test_rank_and_interval_coverage_summary(tmp_path):
    root=tmp_path/"artifact";root.mkdir();out=tmp_path/"out"
    (root/"config.json").write_text(json.dumps({"truth_mode":"prior_predictive"}))
    ledger=[];rows=[]
    for rep,rank,truth in [(0,.2,40.0),(1,.8,80.0),(2,.5,60.0)]:
        for arm in ("A_catalogue_only",):
            ledger.append({"replicate":rep,"gamma":0.0,"validation_size":100,"arm":arm,"status":"SUCCESS"})
            rows.append({"replicate":rep,"gamma":0.0,"validation_size":100,"arm":arm,
                "posterior_rank":rank,"h0_true":truth,"q025":20,"q05":25,"q10":30,
                "q16":35,"q25":40,"q75":70,"q84":75,"q90":80,"q95":85,"q975":90})
    pd.DataFrame(ledger).to_csv(root/"replicate_ledger.csv",index=False)
    pd.DataFrame(rows).to_csv(root/"replicate_posterior_summaries.csv",index=False)
    result=summarize(root,out,bins=3)
    assert result["status"]=="PASS"
    summary=pd.read_csv(out/"prior_predictive_sbc_summary.csv").iloc[0]
    assert summary.successful_replicates==3
    assert summary.mean_rank==pytest.approx(.5)
    assert summary.coverage_50==pytest.approx(2/3)
    hist=pd.read_csv(out/"prior_predictive_rank_histograms.csv")
    assert hist["count"].sum()==3

def test_sbc_summary_rejects_fixed_truth_mode(tmp_path):
    root=tmp_path/"artifact";root.mkdir()
    (root/"config.json").write_text(json.dumps({"truth_mode":"fixed"}))
    with pytest.raises(ValueError,match="prior_predictive"):
        summarize(root,tmp_path/"out")


def test_sbc_summary_marks_incomplete_lifecycle(tmp_path):
    root=tmp_path/"artifact";root.mkdir();out=tmp_path/"out"
    (root/"config.json").write_text(json.dumps({"truth_mode":"prior_predictive"}))
    pd.DataFrame([{"replicate":0,"gamma":0.0,"validation_size":100,"arm":"A_catalogue_only",
        "status":"FAILED","audit_status":"NOT_RUN"}]).to_csv(root/"replicate_ledger.csv",index=False)
    pd.DataFrame(columns=["replicate","gamma","validation_size","arm","posterior_rank","h0_true",
        "q025","q05","q10","q16","q25","q75","q84","q90","q95","q975"]).to_csv(
            root/"replicate_posterior_summaries.csv",index=False)
    result=summarize(root,out,bins=3)
    assert result["status"]=="INCOMPLETE"
    assert result["failed_or_non_successful_cells"]==1
