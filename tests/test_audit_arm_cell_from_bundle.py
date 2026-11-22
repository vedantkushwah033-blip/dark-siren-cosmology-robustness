import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts import audit_arm_cell_from_bundle as auditor


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(root: Path):
    bundle = root / "bundle"
    cell = root / "cell"
    (bundle / "results").mkdir(parents=True)
    cell.mkdir()
    (bundle / "shared_input_manifest.json").write_text('{"status":"FROZEN"}')
    config = {"validation_sizes": [100], "curvature_gamma": 0.5}
    (bundle / "config.json").write_text(json.dumps(config))
    h0 = np.linspace(50.0, 90.0, 5)
    z = np.linspace(0.01, 0.40, 7)
    ll = np.array([-1.0, -0.5, 0.0, -0.5, -1.0])
    density = np.ones(7) / 0.39
    np.savez_compressed(
        bundle / "results/log_likelihoods.npz",
        H0_grid=h0, redshift_grid=z,
        A_catalogue_only_n100=ll,
        density_A_catalogue_only_n100=density,
    )
    pd.DataFrame([{
        "validation_size": 100, "arm": "A_catalogue_only", "status": "SUCCESS",
        "posterior_median": 70.0, "q025": 55.0, "q05": 57.0,
        "q16": 62.0, "q84": 78.0, "q95": 83.0, "q975": 85.0,
    }]).to_csv(bundle / "results/posterior_summaries.csv", index=False)
    pd.DataFrame([{
        "validation_size": 100, "arm": "A_catalogue_only", "kish_ess": 2.0,
    }]).to_csv(bundle / "results/weight_diagnostics.csv", index=False)
    pd.DataFrame([{
        "validation_size": 100, "arm": "A_catalogue_only", "node": -1,
        "realized_bandwidth_redshift": 0.03,
    }]).to_csv(bundle / "results/realized_bandwidths.csv", index=False)

    np.savez_compressed(
        cell / "cell_arrays.npz",
        H0_grid=h0, redshift_grid=z, marginal_log_likelihood=ll,
        component_density=density[None, :],
        component_log_likelihood=ll[None, :],
    )
    metadata = {
        "schema_version": 1, "status": "COMPUTED_NOT_AUDITED",
        "cell_key": {"replicate": 1, "gamma": 0.5,
                     "validation_size": 100, "arm": "A_catalogue_only"},
        "bundle_manifest_sha256": _hash(bundle / "shared_input_manifest.json"),
        "config_sha256": _hash(bundle / "config.json"),
        "posterior_summary": {"posterior_median": 70.0, "q025": 55.0,
                              "q05": 57.0, "q16": 62.0, "q84": 78.0,
                              "q95": 83.0, "q975": 85.0},
        "weight_diagnostics": {"kish_ess": 2.0},
        "realized_bandwidths": [0.03],
    }
    (cell / "cell_result.json").write_text(json.dumps(metadata))
    files = {}
    for name in ("cell_arrays.npz", "cell_result.json"):
        path = cell / name
        files[name] = {"bytes": path.stat().st_size, "sha256": _hash(path)}
    (cell / "manifest.json").write_text(json.dumps({
        "schema_version": 1, "status": "COMPUTED_NOT_AUDITED", "files": files,
    }))
    return bundle, cell


def test_cell_auditor_cross_checks_independent_bundle_reconstruction(tmp_path, monkeypatch):
    bundle, cell = _fixture(tmp_path)
    monkeypatch.setattr(auditor, "audit_bundle", lambda _path: {"status": "PASS"})
    report = auditor.audit_cell(bundle, cell)
    assert report["status"] == "PASS"
    assert report["audit_status"] == "PASS"
    assert report["cell_status"] == "COMPUTED_NOT_AUDITED"


def test_cell_auditor_rejects_tampered_cell_artifact(tmp_path, monkeypatch):
    bundle, cell = _fixture(tmp_path)
    monkeypatch.setattr(auditor, "audit_bundle", lambda _path: {"status": "PASS"})
    with (cell / "cell_arrays.npz").open("ab") as stream:
        stream.write(b"tamper")
    with pytest.raises(AssertionError, match="byte count mismatch|SHA-256 mismatch"):
        auditor.audit_cell(bundle, cell)
