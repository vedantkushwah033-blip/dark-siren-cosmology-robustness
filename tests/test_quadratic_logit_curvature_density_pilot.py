from __future__ import annotations
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
def _load(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/path)
 mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod); return mod

def test_paired_curvature_density_pilot_and_outer_audit(tmp_path):
 runner=_load("curvature_pilot_runner","run_quadratic_logit_curvature_density_pilot.py")
 auditor=_load("curvature_pilot_auditor","audit_quadratic_logit_curvature_density_pilot.py")
 out=tmp_path/"curvature_pilot"
 args=SimpleNamespace(output_dir=out,replicates=1,seed=2026111500,n_parent=500,
   validation_frame=1500,quadrature_order=3,grid_points=101,scott_multiplier=0.7,
   absolute_widths="0.012,0.024")
 result=runner.run(args)
 assert result["status"]=="PASS"
 report=auditor.audit(out)
 assert report["status"]=="PASS"
 assert report["children_audited"]==7
 metrics=pd.read_csv(out/"results"/"density_metrics.csv")
 assert set(metrics.gamma)=={-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0}
 assert set(metrics.validation_size)=={100,300,1000}
 assert metrics.groupby(["gamma","replicate"]).size().eq(metrics.groupby(["gamma","replicate"]).size().iloc[0]).all()

def test_outer_auditor_rejects_ledger_corruption(tmp_path):
 runner=_load("curvature_pilot_runner_ledger","run_quadratic_logit_curvature_density_pilot.py")
 auditor=_load("curvature_pilot_auditor_ledger","audit_quadratic_logit_curvature_density_pilot.py")
 out=tmp_path/"curvature_pilot_ledger"
 runner.run(SimpleNamespace(output_dir=out,replicates=1,seed=2026111510,n_parent=500,
   validation_frame=1500,quadrature_order=3,grid_points=101,scott_multiplier=0.7,
   absolute_widths="0.012,0.024"))
 ledger=pd.read_csv(out/"results"/"replicate_ledger.csv")
 ledger.loc[0,"gamma"]=99.0
 ledger.to_csv(out/"results"/"replicate_ledger.csv",index=False)
 import hashlib,json
 mp=out/"manifest.json"; manifest=json.loads(mp.read_text())
 p=out/"results"/"replicate_ledger.csv"; raw=p.read_bytes()
 manifest["files"]["results/replicate_ledger.csv"]={"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
 mp.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="ledger gamma-replicate keys mismatch"):
  auditor.audit(out)


def test_outer_auditor_rejects_aggregate_metric_corruption_with_refreshed_manifest(tmp_path):
 runner=_load("curvature_pilot_runner_aggregate","run_quadratic_logit_curvature_density_pilot.py")
 auditor=_load("curvature_pilot_auditor_aggregate","audit_quadratic_logit_curvature_density_pilot.py")
 out=tmp_path/"curvature_pilot_aggregate"
 runner.run(SimpleNamespace(output_dir=out,replicates=1,seed=2026111520,n_parent=500,
   validation_frame=1500,quadrature_order=3,grid_points=101,scott_multiplier=0.7,
   absolute_widths="0.012,0.024"))
 import hashlib,json
 path=out/"results"/"density_metrics.csv"
 frame=pd.read_csv(path)
 frame.loc[0,"iae"]+=0.01
 frame.to_csv(path,index=False)
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 raw=path.read_bytes()
 manifest["files"]["results/density_metrics.csv"]={"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="aggregate numerical mismatch: iae"):
  auditor.audit(out)
