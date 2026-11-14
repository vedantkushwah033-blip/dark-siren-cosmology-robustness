"""Multi-replicate common-absolute-bandwidth density-only pilot.

This orchestration preserves one complete raw-input smoke artifact per
replicate. It is a development pilot, not a confirmatory bandwidth selector.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from types import SimpleNamespace
import pandas as pd

from run_absolute_bandwidth_density_smoke import run as run_one

def _hash(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def run(args):
    root=args.output_dir.resolve()
    root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()):
        raise ValueError("output directory must be empty")
    if args.replicates<2: raise ValueError("multi-replicate pilot requires at least two replicates")
    if args.seed_stride<2: raise ValueError("seed stride must be >=2 to keep streams disjoint")
    widths=tuple(float(x) for x in args.absolute_widths)
    if len(widths)<2: raise ValueError("at least two absolute widths are required")
    cfg={
      "schema_version":1,
      "experiment":"multi_replicate_absolute_physical_bandwidth_density_pilot",
      "status":"DEVELOPMENT_PILOT_NOT_CONFIRMATORY_NOT_BANDWIDTH_SELECTION",
      "master_seed":int(args.master_seed),"seed_stride":int(args.seed_stride),
      "replicates":int(args.replicates),"n_parent":int(args.n_parent),
      "validation_frame":int(args.validation_frame),
      "validation_sizes":list(args.validation_sizes),
      "catalogue_floor":float(args.catalogue_floor),
      "quadrature_order":int(args.quadrature_order),
      "grid_points":int(args.grid_points),
      "absolute_widths":list(widths),
      "seed_schedule":"population_seed=master_seed+replicate_id*seed_stride; validation_seed=population_seed+1",
      "limitations":[
        "density recovery only; no H0 likelihood, posterior, calibration, PIT or coverage",
        "stored inclusion posterior nodes are consumed from the production fitter and are not independently refitted",
        "no bandwidth selection; no holdout split; no estimator-superiority claim",
        "synthetic low-redshift Euclidean methodology only"
      ]
    }
    (root/"config.json").write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    (root/"replicates").mkdir()
    ledger=[]
    for rid in range(args.replicates):
        seed=args.master_seed+rid*args.seed_stride
        validation_seed=seed+1
        row={"replicate_id":rid,"population_seed":seed,"validation_seed":validation_seed,
             "status":"NOT_RUN","audit_status":"NOT_RUN","metric_rows":0,"error":""}
        dest=root/"replicates"/f"replicate_{rid:04d}"
        try:
            one=SimpleNamespace(
              output_dir=dest,seed=seed,validation_seed=validation_seed,
              n_parent=args.n_parent,validation_frame=args.validation_frame,
              validation_sizes=tuple(args.validation_sizes),
              catalogue_floor=args.catalogue_floor,
              quadrature_order=args.quadrature_order,grid_points=args.grid_points,
              absolute_widths=widths)
            result=run_one(one)
            row.update(status="SUCCESS",metric_rows=int(result["metric_rows"]),
                       child_manifest_files=int(result["manifest_files"]))
        except Exception as exc:
            row.update(status="FAILED",error=f"{type(exc).__name__}: {exc}")
        ledger.append(row)
        pd.DataFrame(ledger).to_csv(root/"replicate_ledger.csv",index=False)
    metrics=[]
    for row in ledger:
        if row["status"]=="SUCCESS":
            frame=pd.read_csv(root/"replicates"/f"replicate_{row['replicate_id']:04d}"/"results"/"density_metrics.csv")
            frame["replicate"]=int(row["replicate_id"])
            metrics.append(frame)
    combined=pd.concat(metrics,ignore_index=True) if metrics else pd.DataFrame(
       columns=["replicate","validation_size","absolute_bandwidth_sd","arm","node","iae","ise","boundary_iae","status","kish_ess"])
    combined.to_csv(root/"density_metrics_all_replicates.csv",index=False)
    cfg["attempted"]=len(ledger)
    cfg["successful"]=sum(x["status"]=="SUCCESS" for x in ledger)
    cfg["failed"]=sum(x["status"]=="FAILED" for x in ledger)
    (root/"config.json").write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    files=sorted(p for p in root.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS","file_count":len(files),"files":{
       str(p.relative_to(root)):{
         "bytes":p.stat().st_size,"sha256":_hash(p)} for p in files}}
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return {"status":"COMPLETED_WITH_FAILURES" if cfg["failed"] else "PASS",
      "planned_replicates":args.replicates,"successful_replicates":cfg["successful"],
      "failed_replicates":cfg["failed"],"metric_rows":len(combined),"manifest_files":len(files)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-seed",type=int,default=20261120)
    p.add_argument("--seed-stride",type=int,default=10)
    p.add_argument("--replicates",type=int,default=20)
    p.add_argument("--n-parent",type=int,default=4000)
    p.add_argument("--validation-frame",type=int,default=1500)
    p.add_argument("--validation-sizes",type=lambda x:tuple(int(v) for v in x.split(",")),default=(100,300,1000))
    p.add_argument("--catalogue-floor",type=float,default=0.10)
    p.add_argument("--quadrature-order",type=int,default=5)
    p.add_argument("--grid-points",type=int,default=501)
    p.add_argument("--absolute-widths",type=lambda x:tuple(float(v) for v in x.split(",")),
                   default=(0.012,0.016,0.020,0.024,0.028,0.032,0.040))
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    print(json.dumps(run(args),indent=2,sort_keys=True))
if __name__=="__main__": main()
