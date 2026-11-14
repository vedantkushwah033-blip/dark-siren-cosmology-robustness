from __future__ import annotations
import hashlib, importlib.util, json, sys
from pathlib import Path
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
 ledger=pd.read_csv(out/"cell_ledger.csv")
 assert len(ledger)==6
 for rid,group in ledger.groupby("replicate"):
  assert group.population_seed.nunique()==1
  assert group.validation_seed.nunique()==1
  assert group.population_seed.iloc[0]==20261220+int(rid)*10
  assert group.validation_seed.iloc[0]==20261221+int(rid)*10
 metrics=pd.read_csv(out/"density_metrics_all_cells.csv")
 assert set(metrics.gamma)=={-1.0,0.0,1.0}
 assert set(metrics.replicate)=={0,1}
 report=auditor.audit(out)
 assert report["status"]=="PASS"
 assert report["planned_cells"]==6
 assert report["child_raw_audits_passed"]==6
 assert report["aggregate_metric_rows_reconstructed"]==len(metrics)
 assert report["paired_contrast_rows_reconstructed"]==528
 assert report["paired_contrast_summary_rows_reconstructed"]==264
 contrasts=pd.read_csv(out/"paired_contrasts_all_cells.csv")
 assert set(contrasts["contrast_type"])=={"gamma_vs_zero","arm_pair"}
 assert len(contrasts)==528

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
