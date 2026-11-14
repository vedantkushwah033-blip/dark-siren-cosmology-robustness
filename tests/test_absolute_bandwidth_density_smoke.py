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
        catalogue_floor=0.10, quadrature_order=3, grid_points=101,
        absolute_widths=(0.012,0.024),
    )
    result=runner.run(args)
    assert result["status"]=="PASS"
    report=auditor.audit(output)
    assert report["status"]=="PASS"
    assert report["metric_rows_reconstructed"]==len(pd.read_csv(output/"results"/"density_metrics.csv"))

    metrics_path=output/"results"/"density_metrics.csv"
    metrics=pd.read_csv(metrics_path)
    metrics.loc[0,"iae"] += 0.01
    metrics.to_csv(metrics_path,index=False)
    with pytest.raises(AssertionError,match="manifest verification failed"):
        auditor.audit(output)
