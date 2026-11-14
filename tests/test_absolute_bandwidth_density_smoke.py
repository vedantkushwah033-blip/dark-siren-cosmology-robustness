from __future__ import annotations
import importlib.util
from pathlib import Path
from types import SimpleNamespace

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
