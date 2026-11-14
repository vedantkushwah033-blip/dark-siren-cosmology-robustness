from __future__ import annotations
import hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/"scripts"
sys.path.insert(0,str(SCRIPTS))

def _load(name,path):
 spec=importlib.util.spec_from_file_location(name,SCRIPTS/path)
 module=importlib.util.module_from_spec(spec)
 sys.modules[name]=module
 assert spec.loader is not None
 spec.loader.exec_module(module)
 return module

def _args(output):
 from argparse import Namespace
 return Namespace(output_dir=output,master_seed=20261220,seed_stride=10,replicates=2,
  gammas=(-1.0,0.0,1.0),n_parent=350,validation_frame=120,
  validation_sizes=(40,80),beta=-3.9,quadrature_order=3,grid_points=101,
  scott_multiplier=0.7,absolute_widths=(0.012,0.024))

def test_paired_curvature_runner_and_independent_outer_audit(tmp_path):
 runner=_load("curvature_multi_runner","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_pilot"
 result=runner.run(_args(out))
 assert result["status"]=="PASS"
 assert result["planned_cells"]==6 and result["successful_cells"]==6
 outer_config=json.loads((out/"config.json").read_text())
 assert outer_config["arm_E_primary"]=="loss of posterior-mean density"
 assert outer_config["arm_E_secondary"]=="posterior expected node-specific density loss"
 ledger=pd.read_csv(out/"cell_ledger.csv")
 assert len(ledger)==6
 assert (ledger.audit_status=="PENDING_OUTER_AUDIT").all()
 for rid,group in ledger.groupby("replicate"):
  assert group.population_seed.nunique()==1
  assert group.validation_seed.nunique()==1
  assert group.population_seed.iloc[0]==20261220+int(rid)*10
  assert group.validation_seed.iloc[0]==20261221+int(rid)*10
  first_parent=None; first_frame=None; first_selected=None; first_validation=None
  for gamma in (-1.0,0.0,1.0):
   label=f"gamma_{gamma:+.2f}".replace("+","p").replace("-","m").replace(".","_")
   child=out/"cells"/label/f"replicate_{int(rid):04d}"
   parent=pd.read_csv(child/"truth"/"parent_truth.csv")
   frame=pd.read_csv(child/"truth"/"validation_frame_truth.csv")
   validation=pd.read_csv(child/"observable"/"validation_max.csv")
   if first_parent is None:
    first_parent=parent.redshift.to_numpy(float)
    first_frame=frame.redshift.to_numpy(float)
    first_selected=frame.selected_for_validation.to_numpy(int)
    first_validation=validation.validation_id.to_numpy(int)
   else:
    np.testing.assert_allclose(parent.redshift.to_numpy(float),first_parent,rtol=0,atol=1e-14)
    np.testing.assert_allclose(frame.redshift.to_numpy(float),first_frame,rtol=0,atol=1e-14)
    np.testing.assert_array_equal(frame.selected_for_validation.to_numpy(int),first_selected)
    np.testing.assert_array_equal(validation.validation_id.to_numpy(int),first_validation)
 metrics=pd.read_csv(out/"density_metrics_all_cells.csv")
 assert set(metrics.gamma)=={-1.0,0.0,1.0}
 assert set(metrics.replicate)=={0,1}
 report=auditor.audit(out)
 assert report["status"]=="PASS"
 assert report["planned_cells"]==6
 assert report["child_raw_audits_passed"]==6
 assert report["aggregate_metric_rows_reconstructed"]==len(metrics)
 assert report["paired_contrast_rows_reconstructed"]==576
 assert report["paired_contrast_summary_rows_reconstructed"]==288
 contrasts=pd.read_csv(out/"paired_contrasts_all_cells.csv")
 assert set(contrasts["contrast_type"])=={"gamma_vs_zero","arm_pair"}
 assert len(contrasts)==576
 secondary=contrasts[(contrasts.contrast_type=="gamma_vs_zero") & (contrasts.arm=="E_posterior_weighted_node_metric")]
 assert len(secondary)==48
 assert not ((contrasts.contrast_type=="arm_pair") & (contrasts.arm=="E_posterior_weighted_node_metric")).any()

def test_curvature_outer_auditor_rejects_aggregate_corruption_even_with_updated_manifest(tmp_path):
 runner=_load("curvature_multi_runner_corrupt","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_corrupt","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_corrupt"
 runner.run(_args(out))
 path=out/"density_metrics_all_cells.csv"
 frame=pd.read_csv(path)
 frame.loc[0,"iae"]+=0.01
 frame.to_csv(path,index=False)
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(path.relative_to(out))
 manifest["files"][rel]={"bytes":path.stat().st_size,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="aggregate values mismatch iae"):
  auditor.audit(out)

def test_curvature_runner_preserves_failed_cell_and_auditor_reports_incomplete(tmp_path, monkeypatch):
 runner=_load("curvature_multi_runner_failure","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_failure","audit_quadratic_logit_curvature_density_replicates.py")
 original=runner.run_one
 def fail_one(args):
  if args.gamma == 1.0 and args.seed == 20261220:
   raise RuntimeError("injected lifecycle failure for regression test")
  return original(args)
 monkeypatch.setattr(runner,"run_one",fail_one)
 out=tmp_path/"curvature_one_failed_cell"
 result=runner.run(_args(out))
 assert result["status"]=="COMPLETED_WITH_FAILURES"
 report=auditor.audit(out)
 assert report["status"]=="INCOMPLETE"
 assert report["planned_cells"]==6
 assert report["successful_cells"]==5
 assert report["failed_cells"]==1
 ledger=pd.read_csv(out/"cell_ledger.csv")
 failed=ledger[ledger.status=="FAILED"]
 assert len(failed)==1
 assert failed.iloc[0].gamma==1.0 and failed.iloc[0].replicate==0
 assert "injected lifecycle failure" in failed.iloc[0].error
 assert report["paired_contrast_rows_reconstructed"]==0
 assert report["paired_contrast_summary_rows_reconstructed"]==0
 assert pd.read_csv(out/"paired_contrasts_all_cells.csv").empty
 assert pd.read_csv(out/"paired_contrast_summary.csv").empty

def test_curvature_runner_preserves_auditable_artifact_when_every_cell_fails(tmp_path, monkeypatch):
 runner=_load("curvature_multi_runner_all_failed","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_all_failed","audit_quadratic_logit_curvature_density_replicates.py")
 def fail_all(args):
  raise RuntimeError("injected all-cell failure")
 monkeypatch.setattr(runner,"run_one",fail_all)
 out=tmp_path/"curvature_all_failed"
 result=runner.run(_args(out))
 assert result["status"]=="COMPLETED_WITH_FAILURES"
 report=auditor.audit(out)
 assert report["status"]=="INCOMPLETE"
 assert report["planned_cells"]==6
 assert report["successful_cells"]==0
 assert report["failed_cells"]==6
 assert report["paired_contrast_rows_reconstructed"]==0
 assert report["paired_contrast_summary_rows_reconstructed"]==0


def test_curvature_outer_auditor_rejects_unmanifested_files(tmp_path):
 runner=_load("curvature_multi_runner_unmanifested","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_unmanifested","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_unmanifested"
 runner.run(_args(out))
 (out/"unexpected_extra.txt").write_text("not covered by manifest\\n")
 with pytest.raises(AssertionError,match="outer manifest inventory mismatch"):
  auditor.audit(out)


def test_curvature_outer_auditor_rejects_false_child_audit_claim_with_updated_manifest(tmp_path):
 runner=_load("curvature_multi_runner_false_audit","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_false_audit","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_false_audit_claim"
 runner.run(_args(out))
 ledger_path=out/"cell_ledger.csv"
 ledger=pd.read_csv(ledger_path)
 ledger.loc[ledger.status=="SUCCESS","audit_status"]="PASS"
 ledger.to_csv(ledger_path,index=False)
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(ledger_path.relative_to(out))
 manifest["files"][rel]={"bytes":ledger_path.stat().st_size,"sha256":hashlib.sha256(ledger_path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="successful cells must await"):
  auditor.audit(out)


def test_curvature_outer_manifest_covers_nested_child_manifest(tmp_path):
 runner=_load("curvature_multi_runner_child_manifest","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_child_manifest","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_child_manifest"
 runner.run(_args(out))
 child_manifest=out/"cells"/"gamma_m1_00"/"replicate_0000"/"manifest.json"
 assert child_manifest.is_file()
 child_manifest.write_text("{\\n  \\"tampered\\": true\\n}\\n")
 with pytest.raises(AssertionError,match="outer manifest verification failed"):
  auditor.audit(out)


def test_curvature_outer_auditor_rejects_child_config_drift_with_updated_manifest(tmp_path):
 runner=_load("curvature_multi_runner_config_drift","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_config_drift","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_config_drift"
 runner.run(_args(out))
 config_path=out/"config.json"
 config=json.loads(config_path.read_text())
 config["beta"]=-3.8
 config_path.write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(config_path.relative_to(out))
 manifest["files"][rel]={"bytes":config_path.stat().st_size,"sha256":hashlib.sha256(config_path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="child validation-frame or beta mismatch"):
  auditor.audit(out)


def test_curvature_outer_auditor_rejects_prior_drift_with_updated_manifest(tmp_path):
 runner=_load("curvature_multi_runner_prior_drift","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_prior_drift","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_prior_drift"
 runner.run(_args(out))
 config_path=out/"config.json"
 config=json.loads(config_path.read_text())
 config["coefficient_prior"]["sd"]=[8.0,4.0]
 config_path.write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(config_path.relative_to(out))
 manifest["files"][rel]={"bytes":config_path.stat().st_size,"sha256":hashlib.sha256(config_path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="coefficient prior mismatch"):
  auditor.audit(out)


def test_curvature_runner_rejects_gamma_values_with_colliding_child_labels(tmp_path):
 runner=_load("curvature_multi_runner_label_collision","run_quadratic_logit_curvature_density_replicates.py")
 args=_args(tmp_path/"gamma_label_collision")
 args.gammas=(0.0,0.001)
 with pytest.raises(ValueError,match="collide under the two-decimal child-directory"):
  runner.run(args)
 assert list(args.output_dir.iterdir())==[]


def test_curvature_outer_auditor_rejects_child_diagnostic_count_drift_even_with_updated_manifest(tmp_path):
 runner=_load("curvature_multi_runner_count_drift","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_count_drift","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_count_drift"
 runner.run(_args(out))
 ledger_path=out/"cell_ledger.csv"
 ledger=pd.read_csv(ledger_path)
 ledger.loc[0,"weight_diagnostic_rows"]+=1
 ledger.to_csv(ledger_path,index=False)
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(ledger_path.relative_to(out))
 manifest["files"][rel]={"bytes":ledger_path.stat().st_size,"sha256":hashlib.sha256(ledger_path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="cell weight-diagnostic count differs from ledger"):
  auditor.audit(out)


def test_curvature_outer_auditor_rejects_entropy_ess_aggregate_corruption_after_manifest_refresh(tmp_path):
 runner=_load("curvature_multi_runner_entropy_corrupt","run_quadratic_logit_curvature_density_replicates.py")
 auditor=_load("curvature_multi_auditor_entropy_corrupt","audit_quadratic_logit_curvature_density_replicates.py")
 out=tmp_path/"curvature_entropy_corrupt"
 runner.run(_args(out))
 path=out/"weight_diagnostics_all_cells.csv"
 frame=pd.read_csv(path)
 assert "entropy_ess" in frame.columns
 frame.loc[0,"entropy_ess"]+=0.25
 frame.to_csv(path,index=False)
 manifest_path=out/"manifest.json"
 manifest=json.loads(manifest_path.read_text())
 rel=str(path.relative_to(out))
 manifest["files"][rel]={"bytes":path.stat().st_size,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
 manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\\n")
 with pytest.raises(AssertionError,match="aggregate values mismatch entropy_ess"):
  auditor.audit(out)
