from __future__ import annotations
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


ROOT=Path(__file__).resolve().parents[1]

def _load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/path)
    module=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module

def test_absolute_bandwidth_density_smoke_independent_reconstruction_and_corruption(tmp_path):
    runner=_load("absolute_density_runner","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density"
    args=SimpleNamespace(
        output_dir=output, seed=20261114, validation_seed=20261115,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3, grid_points=101,
        absolute_widths=(0.012,0.024),
    )
    result=runner.run(args)
    assert result["status"]=="PASS"
    raw=np.load(output/"results"/"density_reconstruction_inputs.npz")
    grid=raw["redshift_grid"]
    assert np.allclose(raw["base_measure_grid"],1.0+1.5*(grid-0.01)/0.39)
    assert np.array_equal(raw["merger_weight_grid"],np.ones_like(grid))
    report=auditor.audit(output)
    assert report["status"]=="PASS"
    assert report["metric_rows_reconstructed"]==len(pd.read_csv(output/"results"/"density_metrics.csv"))
    assert report["scott_width_rows_reconstructed"]==len(pd.read_csv(output/"results"/"scott_widths.csv"))

    metrics_path=output/"results"/"density_metrics.csv"
    metrics=pd.read_csv(metrics_path)
    metrics.loc[0,"iae"] += 0.01
    metrics.to_csv(metrics_path,index=False)
    with pytest.raises(AssertionError,match="manifest verification failed"):
        auditor.audit(output)


def test_absolute_bandwidth_auditor_rejects_non_nested_validation_ids(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_ids","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_ids","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_ids"
    args=SimpleNamespace(
        output_dir=output, seed=20261116, validation_seed=20261117,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3, grid_points=101,
        absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    path=output/"observable"/"validation_n80.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"validation_id"] = frame.loc[1,"validation_id"]
    frame.to_csv(path,index=False)

    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(path.relative_to(output))
    payload=path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")

    with pytest.raises(AssertionError,match="validation IDs are not the declared nested prefix"):
        auditor.audit(output)

def test_absolute_bandwidth_auditor_rejects_modified_scott_width_even_with_updated_manifest(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_scott","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_scott","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_scott"
    args=SimpleNamespace(
        output_dir=output, seed=20261118, validation_seed=20261119,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    path=output/"results"/"scott_widths.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"realized_bandwidth_redshift"] += 0.01
    frame.to_csv(path,index=False)
    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(path.relative_to(output))
    payload=path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    with pytest.raises(AssertionError,match="independently reconstructed realized Scott widths"):
        auditor.audit(output)


def _rehash(root, path):
    import hashlib
    import json
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    payload = path.read_bytes()
    rel = str(path.relative_to(root))
    manifest["files"][rel] = {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def test_posterior_mean_density_and_expected_node_loss_are_separate_estimands(tmp_path):
    runner = _load("absolute_density_runner_estimands", "run_absolute_bandwidth_density_smoke.py")
    auditor = _load("absolute_density_auditor_estimands", "audit_absolute_bandwidth_density_smoke.py")
    output = tmp_path / "absolute_density_estimands"
    args = SimpleNamespace(
        output_dir=output, seed=20261122, validation_seed=20261123,
        n_parent=500, validation_frame=120, validation_sizes=(40, 80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012, 0.024),
    )
    runner.run(args)
    report = auditor.audit(output)
    assert report["status"] == "PASS"
    metrics = pd.read_csv(output / "results" / "density_metrics.csv")
    for (_, _), group in metrics.groupby(["validation_size", "absolute_bandwidth_sd"]):
        posterior_mean = group.loc[group.arm == "E_posterior_mean_density"].iloc[0]
        expected_node_loss = group.loc[group.arm == "E_posterior_weighted_node_metric"].iloc[0]
        for metric in ("iae", "ise", "boundary_iae"):
            assert posterior_mean[metric] <= expected_node_loss[metric] + 1e-12
        assert np.isnan(posterior_mean["kish_ess"])


def test_auditor_rejects_corrupted_posterior_mean_metric_after_manifest_refresh(tmp_path):
    runner = _load("absolute_density_runner_metric_corrupt", "run_absolute_bandwidth_density_smoke.py")
    auditor = _load("absolute_density_auditor_metric_corrupt", "audit_absolute_bandwidth_density_smoke.py")
    output = tmp_path / "absolute_density_metric_corrupt"
    args = SimpleNamespace(
        output_dir=output, seed=20261124, validation_seed=20261125,
        n_parent=500, validation_frame=120, validation_sizes=(40, 80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012, 0.024),
    )
    runner.run(args)
    path = output / "results" / "density_metrics.csv"
    frame = pd.read_csv(path)
    mask = frame.arm == "E_posterior_mean_density"
    assert mask.any()
    frame.loc[mask.idxmax(), "iae"] += 0.01
    frame.to_csv(path, index=False)
    _rehash(output, path)
    with pytest.raises(AssertionError, match="metric iae mismatch"):
        auditor.audit(output)


def test_density_smoke_does_not_apply_base_measure_twice(tmp_path):
    runner = _load("absolute_density_runner_base_once", "run_absolute_bandwidth_density_smoke.py")
    auditor = _load("absolute_density_auditor_base_once", "audit_absolute_bandwidth_density_smoke.py")
    output = tmp_path / "absolute_density_base_once"
    args = SimpleNamespace(
        output_dir=output, seed=20261126, validation_seed=20261127,
        n_parent=500, validation_frame=120, validation_sizes=(40, 80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012, 0.024),
    )
    runner.run(args)
    raw = np.load(output / "results" / "density_reconstruction_inputs.npz")
    grid = raw["redshift_grid"]
    base = 1.0 + 1.5 * (grid - 0.01) / 0.39
    assert np.allclose(raw["base_measure_grid"], base)
    assert np.array_equal(raw["merger_weight_grid"], np.ones_like(grid))
    assert auditor.audit(output)["status"] == "PASS"
