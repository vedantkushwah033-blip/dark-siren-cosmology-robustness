"""Outer auditor for paired multi-gamma curvature density pilot."""
from __future__ import annotations
import argparse,hashlib,importlib.util,json
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
def _hash(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1<<20),b""): h.update(b)
 return h.hexdigest()
def _load():
 p=HERE/"audit_quadratic_logit_curvature_density_smoke.py"
 s=importlib.util.spec_from_file_location("curvature_child_auditor_outer",p)
 m=importlib.util.module_from_spec(s); assert s.loader is not None; s.loader.exec_module(m); return m
def audit(root):
 root=Path(root).resolve(); manifest=json.loads((root/"manifest.json").read_text())
 checked=0
 for rel,meta in manifest["files"].items():
  p=root/rel
  if not p.is_file() or p.stat().st_size!=int(meta["bytes"]) or _hash(p)!=meta["sha256"]:
   raise AssertionError(f"outer manifest verification failed: {rel}")
  checked+=1
 if checked!=int(manifest["file_count"]): raise AssertionError("outer manifest file count mismatch")
 cfg=json.loads((root/"config.json").read_text())
 if cfg.get("status")!="ENGINEERING_PILOT_NOT_CONFIRMATORY": raise AssertionError("unexpected experiment claim status")
 ledger=pd.read_csv(root/"results"/"replicate_ledger.csv")
 expected={(float(g),int(r)) for g in cfg["gammas"] for r in range(int(cfg["replicates_per_gamma"]))}
 actual={(float(x.gamma),int(x.replicate)) for x in ledger.itertuples()}
 if actual!=expected or ledger.duplicated(["gamma","replicate"]).any(): raise AssertionError("ledger gamma-replicate keys mismatch")
 child_auditor=_load(); success=0
 for row in ledger.itertuples():
  child=root/row.child_path
  if row.status=="SUCCESS":
   child_report=child_auditor.audit(child)
   if child_report.get("status")!="PASS": raise AssertionError(f"child audit failed: {row.child_path}")
   success+=1
  elif row.status not in ("FAILED","NON_EVALUABLE","NOT_RUN"):
   raise AssertionError(f"unknown lifecycle status: {row.status}")
 if success!=int((ledger.status=="SUCCESS").sum()): raise AssertionError("success ledger count mismatch")
 for filename,key in [("density_metrics.csv",["gamma","replicate","validation_size","absolute_bandwidth_sd","arm","node"]),
                      ("weight_diagnostics.csv",["gamma","replicate","validation_size","arm","node"]),
                      ("scott_widths.csv",["gamma","replicate","validation_size","arm","node"])]:
  frame=pd.read_csv(root/"results"/filename)
  if success and frame.empty: raise AssertionError(f"aggregate {filename} unexpectedly empty")
  if not frame.empty and frame.duplicated(key).any(): raise AssertionError(f"duplicate aggregate keys in {filename}")
  if not frame.empty:
   valid=set((float(r.gamma),int(r.replicate)) for r in ledger[ledger.status=="SUCCESS"].itertuples())
   got=set((float(r.gamma),int(r.replicate)) for r in frame.itertuples())
   if got!=valid: raise AssertionError(f"{filename} rows do not match successful ledger children")
 return {"status":"PASS" if success==len(ledger) else "INCOMPLETE",
  "outer_manifest_files_verified":checked,"children_audited":success,
  "children_expected":len(ledger),"failed_children":int((ledger.status=="FAILED").sum()),
  "aggregate_metric_rows":len(pd.read_csv(root/"results"/"density_metrics.csv")),
  "scope":"Outer manifest, lifecycle ledger, child raw-input audits, and aggregate key/coverage checks; no H0 inference."}
def main():
 p=argparse.ArgumentParser(); p.add_argument("artifact_dir",type=Path); a=p.parse_args()
 print(json.dumps(audit(a.artifact_dir),indent=2,sort_keys=True))
if __name__=="__main__":main()
