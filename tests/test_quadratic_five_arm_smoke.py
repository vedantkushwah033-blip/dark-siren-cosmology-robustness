from pathlib import Path
from types import SimpleNamespace
import json
import pandas as pd

from scripts.run_quadratic_five_arm_smoke import run
from scripts.audit_quadratic_five_arm_smoke import audit

def test_quadratic_five_arm_smoke_roundtrip(tmp_path: Path):
    root=tmp_path/"quadratic_five_arm"
    args=SimpleNamespace(seed=8123,validation_seed=8124,n_parent=300,n_detected=4,
      validation_frame=120,validation_sizes=[50,100],catalogue_floor=0.10,
      curvature_gamma=0.5,curvature_beta=-3.9,h0_true=70.0,h0_points=41,
      scott_multiplier=0.7,quadrature_order=3,sigma_mpc=140.0,
      detection_threshold_mpc=250.0,output_dir=root)
    result=run(args)
    assert result["status"]=="PASS"
    summary=pd.read_csv(root/"results"/"posterior_summaries.csv")
    assert len(summary)==10
    assert set(summary.arm)=={
      "A_catalogue_only","B_parent_oracle","C_known_probability_ipw",
      "D_map_plugin_ipw","E_uncertainty_integrated_ipw"}
    cfg=json.loads((root/"config.json").read_text())
    assert cfg["catalogue_model"]=="quadratic_logit"
    assert cfg["curvature_gamma"]==0.5
    report=audit(root)
    assert report["status"]=="PASS"
    assert report["posterior_summary_rows_reconstructed"]==10
