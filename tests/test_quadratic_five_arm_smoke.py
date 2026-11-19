from pathlib import Path
from types import SimpleNamespace
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hashlib
import json
import numpy as np
import pandas as pd
import pytest

from run_quadratic_five_arm_smoke import run
from audit_quadratic_five_arm_smoke import audit

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


def test_expanded_common_h0_prior_roundtrip(tmp_path: Path):
    root = tmp_path / "quadratic_five_arm_expanded_prior"
    args = SimpleNamespace(
        seed=8127, validation_seed=8128, n_parent=300, n_detected=4,
        validation_frame=120, validation_sizes=[50, 100], catalogue_floor=0.10,
        curvature_gamma=0.5, curvature_beta=-3.9, h0_true=70.0,
        h0_prior_min=30.0, h0_prior_max=110.0, h0_points=81,
        scott_multiplier=0.7, quadrature_order=3, sigma_mpc=140.0,
        detection_threshold_mpc=250.0, output_dir=root,
    )
    result = run(args)
    assert result["status"] == "PASS"
    config = json.loads((root / "config.json").read_text())
    assert config["h0_prior_support"] == [30.0, 110.0]
    with np.load(root / "results" / "log_likelihoods.npz") as archive:
        assert archive["H0_grid"][0] == 30.0
        assert archive["H0_grid"][-1] == 110.0
    report = audit(root)
    assert report["status"] == "PASS"
    assert report["posterior_summary_rows_reconstructed"] == 10
