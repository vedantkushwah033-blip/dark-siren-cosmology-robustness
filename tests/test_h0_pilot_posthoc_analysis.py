from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "analyze_quadratic_logit_h0_pilot_uncertainty.py"
SPEC = importlib.util.spec_from_file_location("h0_pilot_posthoc", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_posthoc_reports_wilson_and_paired_contrasts(tmp_path):
    rows=[]
    for rep in range(4):
        for gamma in (0.0,1.0):
            for nval in (100,300):
                median=70.0 + (rep-1.5)*0.5 + gamma*1.0 - (nval==300)*0.2
                rows.append({
                    "replicate":rep,"gamma":gamma,"n_validation":nval,
                    "status":"SUCCESS","h0_true":70.0,"h0_median":median,
                    "h0_q16":median-5.0,"h0_q84":median+5.0,
                    "covered68":int(abs(median-70.0)<=5.0),
                    "covered90":1,"covered95":1,
                    "min_catalogue_kish_ess":100.0+rep,
                    "max_normalized_weight":0.03,
                })
    root=tmp_path/"artifact"; root.mkdir()
    pd.DataFrame(rows).to_csv(root/"replicate_ledger.csv",index=False)
    output=tmp_path/"posthoc"
    report=MODULE.analyze(root,output,seed=12)
    assert report["status"]=="PASS"
    cells=pd.read_csv(output/"cell_uncertainty_summary.csv")
    assert len(cells)==4
    assert {"coverage68_wilson95_low","coverage68_wilson95_high"} <= set(cells.columns)
    gamma=pd.read_csv(output/"paired_gamma_contrasts.csv")
    assert set(gamma.metric)=={"signed_median_bias","absolute_median_error","width68"}
    size=pd.read_csv(output/"paired_validation_size_contrasts.csv")
    assert len(size)==6


def test_wilson_interval_contains_observed_rate():
    low,high=MODULE.wilson(7,10)
    assert low < 0.7 < high
