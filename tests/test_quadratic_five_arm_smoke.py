from pathlib import Path
from types import SimpleNamespace
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hashlib
import json
import numpy as np
import pandas as pd
import pytest

from run_quadratic_five_arm_smoke import run, _make_h0_grid
from audit_quadratic_five_arm_smoke import audit

def test_quadratic_five_arm_smoke_roundtrip(tmp_path: Path):
    root=tmp_path/"quadratic_five_arm"
    args=SimpleNamespace(seed=8123,validation_seed=8124,n_parent=300,n_detected=4,
      validation_frame=120,validation_sizes=[50,100],catalogue_floor=0.10,
      curvature_gamma=0.5,curvature_beta=-3.9,h0_true=70.0,
      h0_prior_min=30.0,h0_prior_max=110.0,h0_points=41,h0_grid_type="geometric",
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
    assert cfg["h0_prior_support"]==[30.0,110.0]
    with np.load(root/"results"/"log_likelihoods.npz") as archive:
      assert archive["H0_grid"][0]==30.0 and archive["H0_grid"][-1]==110.0
      assert np.allclose(archive["H0_grid"], np.geomspace(30.0,110.0,41))
    assert cfg["h0_grid_type"]=="geometric"
    report=audit(root)
    assert report["status"]=="PASS"
    assert report["posterior_summary_rows_reconstructed"]==10



def test_h0_grid_builder_supports_both_measures_and_rejects_unknown_type():
    linear=_make_h0_grid(30.0,110.0,11,"linear")
    geometric=_make_h0_grid(30.0,110.0,11,"geometric")
    assert np.allclose(linear,np.linspace(30.0,110.0,11))
    assert np.allclose(geometric,np.geomspace(30.0,110.0,11))
    assert linear[0]==geometric[0]==30.0
    assert linear[-1]==geometric[-1]==110.0
    with pytest.raises(ValueError,match="unsupported H0 grid type"):
        _make_h0_grid(30.0,110.0,11,"logarithmic")


def test_auditor_rejects_tampered_inclusion_quadrature_weights(tmp_path: Path):
    root=tmp_path/"quadratic_five_arm"
    args=SimpleNamespace(seed=8125,validation_seed=8126,n_parent=300,n_detected=4,
      validation_frame=120,validation_sizes=[50,100],catalogue_floor=0.10,
      curvature_gamma=0.5,curvature_beta=-3.9,h0_true=70.0,h0_points=41,
      scott_multiplier=0.7,quadrature_order=3,sigma_mpc=140.0,
      detection_threshold_mpc=250.0,output_dir=root)
    run(args)
    archive_path=root/"results"/"log_likelihoods.npz"
    with np.load(archive_path) as archive:
        arrays={name:archive[name] for name in archive.files}
    weights=arrays["posterior_weights_n50"].copy()
    weights[0] *= 0.5
    weights /= weights.sum()
    arrays["posterior_weights_n50"]=weights
    np.savez_compressed(archive_path,**arrays)

    manifest_path=root/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    for relative,metadata in manifest["files"].items():
        file_path=root/relative
        metadata["bytes"]=file_path.stat().st_size
        metadata["sha256"]=hashlib.sha256(file_path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    with pytest.raises(AssertionError,match="independent inclusion quadrature weights mismatch"):
        audit(root)
