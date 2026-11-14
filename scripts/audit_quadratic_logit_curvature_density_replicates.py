"""Independent outer auditor for paired quadratic-logit curvature density pilots."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd

def _hash(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for block in iter(lambda:f.read(1<<20),b""): h.update(block)
 return h.hexdigest()

def _load_child_auditor():
 path=Path(__file__).with_name("audit_quadratic_logit_curvature_density_smoke.py")
 spec=importlib.util.spec_from_file_location("quadratic_logit_child_raw_auditor",path)
 module=importlib.util.module_from_spec(spec)
 assert spec.loader is not None
 spec.loader.exec_module(module)
 return module.audit

def _compare_aggregate(root,filename,frames,keys,values):
 expected=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=keys+values)
 actual=pd.read_csv(root/filename)
 if len(actual)!=len(expected): raise AssertionError(f"{filename}: aggregate row count mismatch")
 if len(expected)==0: return 0
 if expected.duplicated(keys).any() or actual.duplicated(keys).any():
  raise AssertionError(f"{filename}: duplicate aggregate key")
 expected=expected.sort_values(keys).reset_index(drop=True)
 actual=actual.sort_values(keys).reset_index(drop=True)
 for key in keys:
  if actual[key].tolist()!=expected[key].tolist(): raise AssertionError(f"{filename}: aggregate key mismatch {key}")
 for col in values:
  if not np.allclose(actual[col].to_numpy(float),expected[col].to_numpy(float),rtol=2e-10,atol=2e-10,equal_nan=True):
   raise AssertionError(f"{filename}: aggregate values mismatch {col}")
 return len(expected)

def audit(root:Path):
 root=Path(root).resolve()
 manifest=json.loads((root/"manifest.json").read_text())
 checked=0
 for rel,meta in manifest["files"].items():
  path=root/rel
  if not path.is_file() or path.stat().st_size!=int(meta["bytes"]) or _hash(path)!=meta["sha256"]:
   raise AssertionError(f"outer manifest verification failed: {rel}")
  checked+=1
 if checked!=int(manifest["file_count"]): raise AssertionError("outer manifest file count mismatch")
 cfg=json.loads((root/"config.json").read_text())
 if cfg.get("experiment")!="paired_quadratic_logit_curvature_density_pilot": raise AssertionError("wrong experiment identifier")
 if cfg.get("status")!="DEVELOPMENT_PILOT_NOT_CONFIRMATORY_NOT_BANDWIDTH_SELECTION": raise AssertionError("claim status mismatch")
 if cfg.get("analysis_model")!="degree_one_logistic": raise AssertionError("analysis model mismatch")
 gammas=[float(x) for x in cfg["gammas"]]
 reps=int(cfg["replicates"]); stride=int(cfg["seed_stride"]); master=int(cfg["master_seed"])
 if reps<2 or len(gammas)<2 or len(set(gammas))!=len(gammas) or stride<2: raise AssertionError("invalid planned design")
 ledger=pd.read_csv(root/"cell_ledger.csv").fillna("")
 expected=[(r,g) for r in range(reps) for g in gammas]
 got=list(zip(ledger.replicate.astype(int),ledger.gamma.astype(float)))
 if got!=expected: raise AssertionError("cell ledger does not enumerate exact replicate x gamma design")
 if set(ledger.status)-{"SUCCESS","FAILED"}: raise AssertionError("invalid lifecycle status")
 if int(cfg["attempted_cells"])!=len(expected): raise AssertionError("attempted cell count mismatch")
 if int(cfg["successful_cells"])!=int((ledger.status=="SUCCESS").sum()) or int(cfg["failed_cells"])!=int((ledger.status=="FAILED").sum()):
  raise AssertionError("config totals disagree with cell ledger")
 child_audit=_load_child_auditor()
 metric_frames=[]; diagnostic_frames=[]; scott_frames=[]; child_reports=[]
 for row in ledger.to_dict("records"):
  rid=int(row["replicate"]); gamma=float(row["gamma"])
  seed=master+rid*stride
  if int(row["population_seed"])!=seed or int(row["validation_seed"])!=seed+1: raise AssertionError("paired seed schedule mismatch")
  label=f"gamma_{gamma:+.2f}".replace("+","p").replace("-","m").replace(".","_")
  child=root/"cells"/label/f"replicate_{rid:04d}"
  if row["status"]=="FAILED":
   if int(row["metric_rows"])!=0 or not row["error"] or not child.is_dir(): raise AssertionError("failed cell lacks preserved lifecycle record")
   continue
  child_cfg=json.loads((child/"config.json").read_text())
  if int(child_cfg["population_seed"])!=seed or int(child_cfg["validation_seed"])!=seed+1: raise AssertionError("child seed mismatch")
  if float(child_cfg["gamma"])!=gamma or int(child_cfg["n_parent"])!=int(cfg["n_parent"]): raise AssertionError("child DGP/config mismatch")
  if list(child_cfg["validation_sizes"])!=list(cfg["validation_sizes"]): raise AssertionError("child validation-size mismatch")
  if list(child_cfg["absolute_bandwidth_candidates"])!=list(cfg["absolute_widths"]): raise AssertionError("child width-grid mismatch")
  report=child_audit(child)
  child_reports.append({"replicate":rid,"gamma":gamma,**report})
  for fn,frames in [("density_metrics.csv",metric_frames),("weight_diagnostics.csv",diagnostic_frames),("scott_widths.csv",scott_frames)]:
   frame=pd.read_csv(child/"results"/fn)
   if fn=="density_metrics.csv" and len(frame)!=int(row["metric_rows"]): raise AssertionError("cell metric count differs from ledger")
   frame.insert(0,"replicate",rid); frame.insert(1,"gamma",gamma); frames.append(frame)
 metric_rows=_compare_aggregate(root,"density_metrics_all_cells.csv",metric_frames,
  ["replicate","gamma","validation_size","absolute_bandwidth_sd","arm","node"],
  ["iae","ise","boundary_iae","status","kish_ess"])
 diag_rows=_compare_aggregate(root,"weight_diagnostics_all_cells.csv",diagnostic_frames,
  ["replicate","gamma","validation_size","arm","node"],
  ["kish_ess","kish_ess_fraction","maximum_normalized_weight","top_1_percent_mass","top_5_percent_mass"])
 scott_rows=_compare_aggregate(root,"scott_widths_all_cells.csv",scott_frames,
  ["replicate","gamma","validation_size","arm","node"],["realized_bandwidth_redshift"])
 return {"status":"PASS","manifest_files_verified":checked,"planned_cells":len(expected),
  "successful_cells":len(child_reports),"failed_cells":int((ledger.status=="FAILED").sum()),
  "child_raw_audits_passed":len(child_reports),"aggregate_metric_rows_reconstructed":metric_rows,
  "aggregate_weight_diagnostic_rows_reconstructed":diag_rows,
  "aggregate_scott_width_rows_reconstructed":scott_rows,
  "scope":"Independent child raw-input audits plus exact outer lifecycle and aggregation reconstruction; synthetic density stage only."}

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_dir",type=Path)
 p.add_argument("--output",type=Path); a=p.parse_args()
 report=audit(a.artifact_dir)
 if a.output: a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\\n")
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=="__main__": main()
