"""Paired fixed-truth H0 calibration pilot for quadratic-logit catalogue selection.

The design crosses curvature gamma with nested validation sizes while holding
analysis-parent/GW random streams fixed within replicate. Every child artifact
is independently audited. This is a pilot, not confirmatory calibration.
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd


def seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([master, replicate, stream]).generate_state(1, dtype=np.uint32)[0])


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=10)
    p.add_argument("--seed",type=int,default=2026111905)
    p.add_argument("--gammas",default="-1,0,1")
    p.add_argument("--validation-sizes",default="100,300,500")
    p.add_argument("--n-validation-frame",type=int,default=500)
    p.add_argument("--n-parent",type=int,default=1200)
    p.add_argument("--n-detected",type=int,default=4)
    p.add_argument("--h0-true",type=float,default=70.0)
    p.add_argument("--h0-points",type=int,default=141)
    p.add_argument("--prior-min",type=float,default=20.0)
    p.add_argument("--prior-max",type=float,default=300.0)
    p.add_argument("--orders",default="3,5,7,9")
    p.add_argument("--density-tolerance",type=float,default=0.001)
    p.add_argument("--quantile-tolerance",type=float,default=0.05)
    p.add_argument("--consecutive-passes",type=int,default=2)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    gammas=tuple(float(x) for x in a.gammas.split(","))
    sizes=tuple(int(x) for x in a.validation_sizes.split(","))
    orders=tuple(int(x) for x in a.orders.split(","))
    if a.replicates<2 or len(set(gammas))!=len(gammas) or not gammas:
        p.error("need >=2 replicates and distinct gamma values")
    if tuple(sorted(set(sizes)))!=sizes or sizes[-1]>a.n_validation_frame or sizes[0]<5:
        p.error("validation sizes must be increasing and within the validation frame")
    if not (a.prior_min<a.h0_true<a.prior_max) or a.n_detected<1 or a.h0_points<21:
        p.error("invalid H0 support, truth, event count, or grid")
    root=Path(__file__).resolve().parents[1]
    runner=root/"scripts/run_quadratic_logit_inclusion_h0_smoke.py"
    auditor=root/"scripts/audit_quadratic_logit_inclusion_h0_smoke.py"
    env=os.environ.copy(); env["PYTHONPATH"]=str(root/"src")+os.pathsep+env.get("PYTHONPATH","")
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): p.error("output directory must be empty")
    # Persist the complete plan before any expensive child run starts.
    planned_rows=[]
    for rep in range(a.replicates):
        analysis_seed=seed(a.seed,rep,0); validation_seed=seed(a.seed,rep,1)
        for gamma in gammas:
            for nval in sizes:
                cell=f"g{gamma:+g}_n{nval}_r{rep:03d}".replace("+","p").replace("-","m")
                planned_rows.append({"cell_id":cell,"replicate":rep,"gamma":gamma,
                    "n_validation":nval,"analysis_seed":analysis_seed,
                    "validation_seed":validation_seed,"h0_true":a.h0_true})
    pd.DataFrame(planned_rows).to_csv(out/"planned_cells.csv",index=False)
    rows=[]
    for rep in range(a.replicates):
        analysis_seed=seed(a.seed,rep,0); validation_seed=seed(a.seed,rep,1)
        for gamma in gammas:
            for nval in sizes:
                cell=f"g{gamma:+g}_n{nval}_r{rep:03d}".replace("+","p").replace("-","m")
                dest=out/"replicates"/cell
                row={"cell_id":cell,"replicate":rep,"gamma":gamma,"n_validation":nval,
                     "analysis_seed":analysis_seed,"validation_seed":validation_seed,
                     "h0_true":a.h0_true,"artifact_dir":str(dest),"status":"NOT_RUN",
                     "audit_status":"","error":"","quadrature_status":"",
                     "orders_evaluated":"","converged_at_order":None,
                     "n_catalogue":None,"min_catalogue_kish_ess":float("nan"),
                     "max_normalized_weight":float("nan"),"h0_median":float("nan"),
                     "h0_q025":float("nan"),"h0_q05":float("nan"),
                     "h0_q16":float("nan"),"h0_q84":float("nan"),
                     "h0_q95":float("nan"),"h0_q975":float("nan"),
                     "width68":float("nan"),"covered68":float("nan"),
                     "covered90":float("nan"),"covered95":float("nan"),
                     "median_bias":float("nan") }
                rows.append(row)
                dest.parent.mkdir(parents=True,exist_ok=True)
                cmd=[sys.executable,str(runner),"--seed",str(analysis_seed),
                     "--validation-seed",str(validation_seed),"--n-parent",str(a.n_parent),
                     "--n-detected",str(a.n_detected),"--h0-true",str(a.h0_true),
                     "--n-validation-frame",str(a.n_validation_frame),"--n-validation",str(nval),
                     "--degree","1","--catalogue-model","quadratic_logit",
                     "--curvature-gamma",str(gamma),"--curvature-beta","-3.9",
                     "--orders",a.orders,"--h0-points",str(a.h0_points),
                     "--prior-min",str(a.prior_min),"--prior-max",str(a.prior_max),
                     "--density-tolerance",str(a.density_tolerance),
                     "--quantile-tolerance",str(a.quantile_tolerance),
                     "--consecutive-passes",str(a.consecutive_passes),
                     "--output-dir",str(dest)]
                try:
                    run=subprocess.run(cmd,cwd=root,env=env,capture_output=True,text=True)
                    (dest.parent/(cell+".stdout.txt")).write_text(run.stdout)
                    (dest.parent/(cell+".stderr.txt")).write_text(run.stderr)
                    if run.returncode: raise RuntimeError("runner exit "+str(run.returncode)+": "+run.stderr[-1500:])
                    cfg=json.loads((dest/"config.json").read_text())
                    conv=json.loads((dest/"quadrature_convergence.json").read_text())
                    if cfg["analysis_generator_config"]["h0_true_km_s_mpc"]!=a.h0_true:
                        raise RuntimeError("child H0 truth differs from pilot truth")
                    audit=subprocess.run([sys.executable,str(auditor),str(dest),"--output",str(dest/"audit_report.json")],
                                         cwd=root,env=env,capture_output=True,text=True)
                    (dest.parent/(cell+".audit.stdout.txt")).write_text(audit.stdout)
                    (dest.parent/(cell+".audit.stderr.txt")).write_text(audit.stderr)
                    if audit.returncode: raise RuntimeError("child audit failed: "+audit.stderr[-1500:])
                    report=json.loads((dest/"audit_report.json").read_text())
                    row["audit_status"]=report.get("status","")
                    row["quadrature_status"]=conv["status"]
                    row["orders_evaluated"]=",".join(map(str,conv["orders_evaluated"]))
                    row["converged_at_order"]=conv.get("converged_at_order")
                    row["n_catalogue"]=len(pd.read_csv(dest/"observable/analysis_catalogue.csv"))
                    row["min_catalogue_kish_ess"]=min(x["minimum_catalogue_kish_ess"] for x in json.loads((dest/"prior_scale_sensitivity.json").read_text())["records"])
                    row["max_normalized_weight"]=max(x["maximum_catalogue_normalized_weight"] for x in json.loads((dest/"prior_scale_sensitivity.json").read_text())["records"])
                    if conv["status"]!="CONVERGED": row["status"]="NON_EVALUABLE"
                    else:
                        final=conv["posterior_by_order"][-1]
                        row.update({"status":"SUCCESS","h0_median":final["h0_median"],
                          "h0_q025":final["h0_q025"],"h0_q05":final["h0_q05"],
                          "h0_q16":final["h0_q16"],"h0_q84":final["h0_q84"],
                          "h0_q95":final["h0_q95"],"h0_q975":final["h0_q975"],
                          "width68":final["h0_q84"]-final["h0_q16"],
                          "covered68":int(final["h0_q16"]<=a.h0_true<=final["h0_q84"]),
                          "covered90":int(final["h0_q05"]<=a.h0_true<=final["h0_q95"]),
                          "covered95":int(final["h0_q025"]<=a.h0_true<=final["h0_q975"]),
                          "median_bias":final["h0_median"]-a.h0_true})
                except Exception as exc:
                    row["status"]="FAILED"; row["error"]=str(exc)
                pd.DataFrame(rows).to_csv(out/"replicate_ledger.csv",index=False)
    ledger=pd.DataFrame(rows)
    ledger.to_csv(out/"replicate_ledger.csv",index=False)
    summary=[]
    for (gamma,nval),g in ledger.groupby(["gamma","n_validation"]):
        ok=g[g.status=="SUCCESS"]
        summary.append({"gamma":gamma,"n_validation":nval,"planned":len(g),
          "success":len(ok),"non_evaluable":int((g.status=="NON_EVALUABLE").sum()),
          "failed":int((g.status=="FAILED").sum()),
          "not_run":int((g.status=="NOT_RUN").sum()),
          "bias":float(ok.median_bias.mean()) if len(ok) else None,
          "rmse":float(np.sqrt(np.mean((ok.h0_median-a.h0_true)**2))) if len(ok) else None,
          "coverage68":float(ok.covered68.mean()) if len(ok) else None,
          "coverage90":float(ok.covered90.mean()) if len(ok) else None,
          "coverage95":float(ok.covered95.mean()) if len(ok) else None,
          "mean_width68":float(ok.width68.mean()) if len(ok) else None,
          "median_min_kish_ess":float(ok.min_catalogue_kish_ess.median()) if len(ok) else None,
          "median_max_normalized_weight":float(ok.max_normalized_weight.median()) if len(ok) else None})
    pd.DataFrame(summary).to_csv(out/"calibration_summary.csv",index=False)
    status=ledger.status.value_counts().to_dict()
    report={"status":"PILOT_COMPLETE" if len(ledger)==a.replicates*len(gammas)*len(sizes) else "INCOMPLETE",
      "planned_cells":len(ledger),"status_counts":status,"interpretation":
      "Fixed-truth repeated-sampling pilot only; not joint SBC or confirmatory calibration."}
    (out/"pilot_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    manifest={str(x.relative_to(out)):hashlib.sha256(x.read_bytes()).hexdigest()
      for x in sorted(out.rglob("*")) if x.is_file() and x.name!="manifest.json"}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()
