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
                    "analysis_seed":1000+rep,"validation_seed":2000+rep,
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


def test_posthoc_rejects_replicate_id_seed_mismatch(tmp_path):
    rows=[]
    for gamma in (0.0, 1.0):
        rows.append({
            "replicate":0,"gamma":gamma,"n_validation":100,
            "analysis_seed":1000 + int(gamma),"validation_seed":2000,
            "status":"SUCCESS","h0_true":70.0,"h0_median":71.0 + gamma,
            "h0_q16":66.0,"h0_q84":76.0,"covered68":1,"covered90":1,
            "covered95":1,"min_catalogue_kish_ess":100.0,
            "max_normalized_weight":0.03,
        })
    root=tmp_path/"artifact"; root.mkdir()
    pd.DataFrame(rows).to_csv(root/"replicate_ledger.csv",index=False)
    import pytest
    with pytest.raises(ValueError,match="analysis_seed"):
        MODULE.analyze(root,tmp_path/"posthoc",seed=12)



def test_posthoc_reports_failure_inclusive_coverage_bounds(tmp_path):
    rows = []
    statuses = ["SUCCESS", "SUCCESS", "FAILED", "NON_EVALUABLE", "NOT_RUN"]
    for rep, status in enumerate(statuses):
        success = status == "SUCCESS"
        rows.append({
            "replicate": rep, "gamma": 0.0, "n_validation": 100,
            "analysis_seed": 3000 + rep, "validation_seed": 4000 + rep,
            "status": status, "h0_true": 70.0,
            "h0_median": 71.0 + rep if success else float("nan"),
            "h0_q16": 66.0 if success else float("nan"),
            "h0_q84": 76.0 if success else float("nan"),
            "covered68": 1 if rep == 0 else (0 if success else float("nan")),
            "covered90": 1 if success else float("nan"),
            "covered95": 1 if success else float("nan"),
            "min_catalogue_kish_ess": 80.0 + rep if success else float("nan"),
            "max_normalized_weight": 0.04 if success else float("nan"),
        })
    root = tmp_path / "artifact"
    root.mkdir()
    pd.DataFrame(rows).to_csv(root / "replicate_ledger.csv", index=False)
    output = tmp_path / "posthoc"
    MODULE.analyze(root, output, seed=13)
    cell = pd.read_csv(output / "cell_uncertainty_summary.csv").iloc[0]
    assert cell["coverage68"] == 0.5
    assert cell["coverage68_all_attempt_lower"] == 0.2
    assert cell["coverage68_all_attempt_upper"] == 0.8
    assert cell["non_success"] == 3


def test_posthoc_zero_success_keeps_rmse_interval_undefined(tmp_path):
    rows = []
    for rep in range(3):
        rows.append({
            "replicate": rep, "gamma": 0.0, "n_validation": 100,
            "analysis_seed": 5000 + rep, "validation_seed": 6000 + rep,
            "status": "FAILED", "h0_true": 70.0,
            "h0_median": float("nan"), "h0_q16": float("nan"),
            "h0_q84": float("nan"), "covered68": float("nan"),
            "covered90": float("nan"), "covered95": float("nan"),
            "min_catalogue_kish_ess": float("nan"),
            "max_normalized_weight": float("nan"),
        })
    root = tmp_path / "artifact"
    root.mkdir()
    pd.DataFrame(rows).to_csv(root / "replicate_ledger.csv", index=False)
    output = tmp_path / "posthoc"
    MODULE.analyze(root, output, seed=14)
    cell = pd.read_csv(output / "cell_uncertainty_summary.csv").iloc[0]
    assert pd.isna(cell["rmse_posterior_median"])
    assert pd.isna(cell["rmse_bootstrap95_low"])
    assert pd.isna(cell["rmse_bootstrap95_high"])
    assert cell["coverage68_all_attempt_lower"] == 0.0
    assert cell["coverage68_all_attempt_upper"] == 1.0
