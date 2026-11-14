"""Paired multi-gamma curvature density pilot with child-level raw audits.

The default is an engineering pilot (3 replicates per gamma), not a
confirmatory experiment. Every child is independently audited before its
metrics enter the aggregate. Failed/non-evaluable children remain in ledger.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
GAMMAS=(-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0)
WIDTHS=(0.008,0.012,0.016,0.024,0.032,0.040,0.050,0.060,0.080)
SIZES=(100,300,1000)

def _load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec); assert spec.loader is not None
    spec.loader.exec_module(mod); return mod

def _hash(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def run(a):
    out=Path(a.output_dir).resolve()
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise ValueError("output directory must be empty")
    if a.replicates<1 or a.replicates>100: raise ValueError("replicates must be in [1,100]")
    if a.n_parent<100 or a.validation_frame<max(SIZES): raise ValueError("invalid population/validation sizes")
    if a.seed<0: raise ValueError("seed must be non-negative")
    gammas=GAMMAS
    widths=tuple(float(x) for x in a.absolute_widths.split(","))
    if len(widths)<2 or len(set(widths))!=len(widths) or any(x<=0 for x in widths):
        raise ValueError("at least two distinct positive bandwidths required")
    run_child=_load(HERE/"run_quadratic_logit_curvature_density_smoke.py","curvature_child_runner")
    audit_child=_load(HERE/"audit_quadratic_logit_curvature_density_smoke.py","curvature_child_auditor")
    ledger=[]; metric_parts=[]; diagnostic_parts=[]; scott_parts=[]
    for gamma_index,gamma in enumerate(gammas):
        for rep in range(a.replicates):
            # Common seed schedule across gamma gives paired latent parent and validation frames.
            pop_seed=int(a.seed+100000*rep+17)
            val_seed=int(a.seed+100000*rep+29)
            key=f"g{gamma:+.1f}_r{rep:03d}".replace("+","p").replace("-","m").replace(".","")
            child=out/"children"/key
            child.parent.mkdir(parents=True,exist_ok=True)
            status="NOT_RUN"; error=""
            try:
                from argparse import Namespace
                child_args=Namespace(
                    output_dir=child, seed=pop_seed, validation_seed=val_seed,
                    n_parent=a.n_parent, validation_frame=a.validation_frame,
                    validation_sizes=SIZES, gamma=gamma, beta=-3.9,
                    quadrature_order=a.quadrature_order, grid_points=a.grid_points,
                    scott_multiplier=a.scott_multiplier, absolute_widths=widths,
                )
                child_result=run_child.run(child_args)
                audit=audit_child.audit(child)
                if audit.get("status")!="PASS": raise AssertionError("child raw-input audit did not PASS")
                status="SUCCESS"
                for filename,target in [
                    ("density_metrics.csv",metric_parts),
                    ("weight_diagnostics.csv",diagnostic_parts),
                    ("scott_widths.csv",scott_parts),
                ]:
                    frame=pd.read_csv(child/"results"/filename)
                    frame.insert(0,"replicate",rep)
                    frame.insert(0,"gamma",gamma)
                    frame["population_seed"]=pop_seed
                    frame["validation_seed"]=val_seed
                    target.append(frame)
                audit_summary=json.dumps(audit,sort_keys=True)
            except Exception as exc:
                status="FAILED"; error=f"{type(exc).__name__}: {exc}"
                audit_summary=""
            ledger.append({
                "gamma":gamma,"replicate":rep,"population_seed":pop_seed,
                "validation_seed":val_seed,"status":status,"error":error,
                "child_path":str(child.relative_to(out)),"audit_summary":audit_summary,
            })
    ledger_df=pd.DataFrame(ledger)
    (out/"results").mkdir()
    ledger_df.to_csv(out/"results"/"replicate_ledger.csv",index=False)
    for filename,parts in [
        ("density_metrics.csv",metric_parts),
        ("weight_diagnostics.csv",diagnostic_parts),
        ("scott_widths.csv",scott_parts),
    ]:
        if parts:
            pd.concat(parts,ignore_index=True).to_csv(out/"results"/filename,index=False)
        else:
            pd.DataFrame().to_csv(out/"results"/filename,index=False)
    success=int((ledger_df.status=="SUCCESS").sum())
    failed=int((ledger_df.status=="FAILED").sum())
    config={
        "schema_version":1,"experiment":"quadratic_logit_curvature_paired_density_pilot",
        "status":"ENGINEERING_PILOT_NOT_CONFIRMATORY",
        "gammas":list(gammas),"replicates_per_gamma":a.replicates,
        "n_parent":a.n_parent,"validation_frame":a.validation_frame,
        "validation_sizes":list(SIZES),"absolute_bandwidth_candidates":list(widths),
        "quadrature_order":a.quadrature_order,"grid_points":a.grid_points,
        "scott_multiplier":a.scott_multiplier,"base_seed":a.seed,
        "seed_pairing":"same population and validation seeds across gamma within replicate",
        "arm_C":"true quadratic inclusion probabilities; oracle diagnostic only",
        "arm_D_E":"degree-one logistic fitted to independent validation outcomes",
        "arm_E_primary":"loss of posterior-mean density",
        "lifecycle":"SUCCESS/FAILED recorded per gamma-replicate; no failed rows silently dropped from ledger",
        "limitations":["pilot replicate count","no holdout bandwidth selection","density-only; no H0 inference",
          "synthetic low-redshift Euclidean model only","oracle arm C is diagnostic, not deployable"],
    }
    (out/"config.json").write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
    files=sorted(p for p in out.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS" if failed==0 else "INCOMPLETE","file_count":len(files),
      "files":{str(p.relative_to(out)):{"bytes":p.stat().st_size,"sha256":_hash(p)} for p in files}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return {"status":"PASS" if failed==0 else "INCOMPLETE","children":len(ledger),
      "success":success,"failed":failed,"metric_rows":sum(len(x) for x in metric_parts),
      "manifest_files":len(files),"output_dir":str(out)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--replicates",type=int,default=3)
    p.add_argument("--seed",type=int,default=2026111500)
    p.add_argument("--n-parent",type=int,default=4000)
    p.add_argument("--validation-frame",type=int,default=1500)
    p.add_argument("--quadrature-order",type=int,default=5)
    p.add_argument("--grid-points",type=int,default=401)
    p.add_argument("--scott-multiplier",type=float,default=0.7)
    p.add_argument("--absolute-widths",default="0.008,0.012,0.016,0.024,0.032,0.040,0.050,0.060,0.080")
    a=p.parse_args(); print(json.dumps(run(a),indent=2,sort_keys=True))
if __name__=="__main__": main()
