from __future__ import annotations
import hashlib, importlib.util, json, sys
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/"scripts"
sys.path.insert(0,str(SCRIPTS))

def _load(name,path):
    spec=importlib.util.spec_from_file_location(name,SCRIPTS/path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

def _rehash(root,path):
    rel=str(path.relative_to(root))
    payload=path.read_bytes()
    manifest_path=root/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")

def _args(output):
    return SimpleNamespace(
      output_dir=output,master_seed=20261120,seed_stride=10,replicates=2,
      n_parent=500,validation_frame=120,validation_sizes=(40,80),
      catalogue_floor=0.10,scott_multiplier=0.7,quadrature_order=3,grid_points=101,
      absolute_widths=(0.016,0.024))

def test_two_replicate_density_pilot_and_independent_orchestration_audit(tmp_path):
    runner=_load("absolute_density_replicates_runner","run_absolute_bandwidth_density_replicates.py")
    auditor=_load("absolute_density_replicates_auditor","audit_absolute_bandwidth_density_replicates.py")
    out=tmp_path/"replicate_pilot"
    result=runner.run(_args(out))
    assert result["planned_replicates"]==2
    assert result["successful_replicates"]==2
    report=auditor.audit(out)
    assert report["status"]=="PASS"
    assert report["child_raw_audits_passed"]==2
    assert report["aggregate_metric_rows_reconstructed"]==len(pd.read_csv(out/"density_metrics_all_replicates.csv"))
    assert report["scott_width_rows_reconstructed"]==len(pd.read_csv(out/"scott_widths_all_replicates.csv"))

def test_two_replicate_quadratic_logit_density_pilot_and_cross_gamma_audit(tmp_path):
    runner=_load("quadratic_density_replicates_runner","run_absolute_bandwidth_density_replicates.py")
    cross_auditor=_load("quadratic_density_cross_auditor","audit_quadratic_logit_curvature_pilot.py")
    root=tmp_path/"quadratic_replicate_pilot"
    for gamma in (0.0,1.0):
        args=_args(root/f"gamma_{gamma:g}")
        args.gamma=gamma
        args.beta=-3.9
        result=runner.run(args)
        assert result["planned_replicates"]==2
        assert result["successful_replicates"]==2

    report=cross_auditor.audit(
        root,expected_gammas=(0.0,1.0),output_dir=root/"summary"
    )
    assert report["status"]=="PASS"
    assert report["paired_seed_schedule_verified"] is True
    assert report["paired_parent_and_validation_redshifts_verified"] is True
    assert report["replicates_per_gamma"]==2
    assert (root/"summary"/"curvature_density_iae_summary.csv").is_file()
    assert (root/"summary"/"curvature_paired_iae_contrasts.csv").is_file()


def test_orchestration_auditor_rejects_duplicate_aggregate_key_even_with_updated_manifest(tmp_path):
    runner=_load("absolute_density_replicates_runner_badkey","run_absolute_bandwidth_density_replicates.py")
    auditor=_load("absolute_density_replicates_auditor_badkey","audit_absolute_bandwidth_density_replicates.py")
    out=tmp_path/"replicate_badkey"
    runner.run(_args(out))
    path=out/"density_metrics_all_replicates.csv"
    frame=pd.read_csv(path)
    frame.loc[1]=frame.loc[0]
    frame.to_csv(path,index=False)
    _rehash(out,path)
    with pytest.raises(AssertionError,match="aggregate metric row count mismatch|duplicate aggregate metric key"):
        auditor.audit(out)


def test_failed_replicate_is_preserved_in_planned_denominator(tmp_path, monkeypatch):
    runner=_load("absolute_density_replicates_runner_failure","run_absolute_bandwidth_density_replicates.py")
    auditor=_load("absolute_density_replicates_auditor_failure","audit_absolute_bandwidth_density_replicates.py")
    original=runner.run_one
    calls={"count":0}
    def fail_second(args):
        calls["count"]+=1
        if calls["count"]==2:
            raise RuntimeError("injected replicate failure for lifecycle test")
        return original(args)
    monkeypatch.setattr(runner,"run_one",fail_second)
    out=tmp_path/"replicate_failure"
    result=runner.run(_args(out))
    assert result["successful_replicates"]==1
    assert result["failed_replicates"]==1
    ledger=pd.read_csv(out/"replicate_ledger.csv")
    assert ledger.status.tolist()==["SUCCESS","FAILED"]
    assert "injected replicate failure" in ledger.loc[1,"error"]
    assert (out/"replicates"/"replicate_0001").is_dir()
    report=auditor.audit(out)
    assert report["status"]=="PASS"
    assert report["planned_replicates"]==2
    assert report["failed_replicates"]==1

def test_outer_auditor_rejects_modified_aggregate_scott_width_even_with_updated_manifest(tmp_path):
    runner=_load("absolute_density_replicates_runner_scott","run_absolute_bandwidth_density_replicates.py")
    auditor=_load("absolute_density_replicates_auditor_scott","audit_absolute_bandwidth_density_replicates.py")
    out=tmp_path/"replicate_bad_scott"
    runner.run(_args(out))
    path=out/"scott_widths_all_replicates.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"realized_bandwidth_redshift"] += 0.01
    frame.to_csv(path,index=False)
    _rehash(out,path)
    with pytest.raises(AssertionError,match="aggregate Scott width values mismatch"):
        auditor.audit(out)
