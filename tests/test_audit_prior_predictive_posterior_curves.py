from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_prior_predictive_posterior_curves import ARMS, _posterior_cdf_summary, audit

def _write_manifest(root):
    files = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "manifest.json":
            continue
        payload = path.read_bytes()
        files[path.relative_to(root).as_posix()] = {
            "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
    (root / "manifest.json").write_text(json.dumps({
        "status": "PASS", "file_count": len(files), "files": files}, indent=2))

def _fixture(root):
    child = root / "replicates" / "rep000_gp0"
    (child / "results").mkdir(parents=True)
    (root / "config.json").write_text(json.dumps({
        "truth_mode": "prior_predictive", "replicates": 1, "gammas": [0.0],
        "validation_sizes": [100], "h0_points": 5, "h0_prior_support": [30.0, 110.0]}))
    (child / "audit_report.json").write_text(json.dumps({"status": "PASS"}))
    (child / "config.json").write_text(json.dumps({"curvature_gamma": 0.0, "h0_true": 50.0}))
    grid = np.linspace(30.0, 110.0, 5)
    arrays, ledger, summaries = {"H0_grid": grid}, [], []
    for arm in ARMS:
        arrays[f"{arm}_n100"] = np.zeros(5)
        posterior = _posterior_cdf_summary(grid, arrays[f"{arm}_n100"], 50.0)
        ledger.append({"replicate": 0, "gamma": 0.0, "validation_size": 100,
            "arm": arm, "status": "SUCCESS", "audit_status": "PASS"})
        summaries.append({"replicate": 0, "gamma": 0.0, "validation_size": 100,
            "arm": arm, "status": "SUCCESS", "audit_status": "PASS",
            "h0_true": 50.0, **posterior})
    np.savez_compressed(child / "results" / "log_likelihoods.npz", **arrays)
    pd.DataFrame(ledger).to_csv(root / "replicate_ledger.csv", index=False)
    pd.DataFrame(summaries).to_csv(root / "replicate_posterior_summaries.csv", index=False)
    _write_manifest(root)

def test_flat_likelihood_has_uniform_cdf_and_expected_quantiles():
    result = _posterior_cdf_summary(np.linspace(30.0, 110.0, 5), np.zeros(5), 50.0)
    assert result["posterior_rank"] == pytest.approx(.25)
    assert result["q025"] == pytest.approx(32.0)
    assert result["q975"] == pytest.approx(108.0)

def test_curve_auditor_reconstructs_rows_and_checks_manifest(tmp_path):
    _fixture(tmp_path)
    result = audit(tmp_path)
    assert result["status"] == "PASS"
    assert result["planned_likelihood_curves"] == 5
    assert result["independently_reconstructed_posterior_ranks"] == 5
    assert result["successful_child_audits_rechecked"] == 1
    assert result["maximum_absolute_rank_difference"] < 1e-14
    assert result["maximum_absolute_quantile_difference"] < 1e-12

def test_curve_auditor_rejects_tampered_rank_after_manifest_refresh(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "replicate_posterior_summaries.csv"
    frame = pd.read_csv(path)
    frame.loc[0, "posterior_rank"] = .99
    frame.to_csv(path, index=False)
    _write_manifest(tmp_path)
    with pytest.raises(AssertionError, match="posterior rank mismatch"):
        audit(tmp_path)
