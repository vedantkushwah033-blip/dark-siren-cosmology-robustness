#!/usr/bin/env python3
"""Inspect frozen paired-degree failures and quantify saved MAP-mode discrepancies."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize, root
from scipy.special import expit

def label(rep: int, gamma: float) -> str:
    return f"rep{rep:03d}_g{gamma:+g}".replace("+","p").replace("-","m")

def mode_diagnostics(child: Path, n: int) -> dict:
    config=json.loads((child/"config.json").read_text())
    validation=pd.read_csv(child/"observable"/f"validation_n{n}.csv")
    archive=np.load(child/"results"/"log_likelihoods.npz")
    z=validation["redshift"].to_numpy(float)
    y=validation["retained"].to_numpy(float)
    center,scale=0.205,0.195
    qcenter=float(config.get("fitted_quadratic_basis_center",0.0))
    x=(z-center)/scale
    design=np.column_stack([np.ones_like(x),x,x*x-qcenter])
    prior_mean=np.zeros(3); prior_sd=np.asarray([5.0,2.5,2.5]); precision=1/prior_sd**2
    def objective(theta):
        eta=design@theta
        return float(np.sum(np.logaddexp(0.0,eta)-y*eta)+0.5*np.sum(((theta-prior_mean)/prior_sd)**2))
    def gradient(theta):
        return design.T@(expit(design@theta)-y)+(theta-prior_mean)*precision
    def hessian(theta):
        p=expit(design@theta); w=p*(1-p)
        return design.T@(w[:,None]*design)+np.diag(precision)
    saved=np.asarray(archive[f"posterior_mode_n{n}"],float)
    saved_score=float(np.linalg.norm(gradient(saved),ord=np.inf))
    fit=minimize(objective,np.zeros(3),jac=gradient,hess=hessian,method="trust-exact",
                 options={"maxiter":2000,"gtol":1e-10})
    trust=np.asarray(fit.x,float)
    trust_score=float(np.linalg.norm(gradient(trust),ord=np.inf))
    solved=root(gradient,saved,jac=hessian,method="hybr",options={"xtol":1e-12})
    root_mode=np.asarray(solved.x,float)
    root_score=float(np.linalg.norm(gradient(root_mode),ord=np.inf))
    return {
      "validation_size":int(n),
      "saved_mode":saved.tolist(),
      "trust_exact_recomputed_mode":trust.tolist(),
      "root_refined_mode":root_mode.tolist(),
      "saved_score_inf":saved_score,
      "trust_exact_score_inf":trust_score,
      "root_score_inf":root_score,
      "saved_vs_trust_max_abs":float(np.max(np.abs(saved-trust))),
      "saved_vs_root_max_abs":float(np.max(np.abs(saved-root_mode))),
      "trust_vs_root_max_abs":float(np.max(np.abs(trust-root_mode))),
      "saved_objective":objective(saved),
      "root_objective":objective(root_mode),
      "saved_minus_root_objective":float(objective(saved)-objective(root_mode)),
      "trust_success":bool(fit.success),
      "trust_message":str(fit.message),
      "root_success":bool(solved.success),
      "root_message":str(solved.message),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("artifact_dir",type=Path)
    p.add_argument("--output",type=Path,default=Path("degree2_failure_diagnostics.json"))
    a=p.parse_args(); root_dir=a.artifact_dir; d2=root_dir/"degree2"
    ledger=pd.read_csv(d2/"replicate_ledger.csv")
    failed=ledger.loc[ledger.status!="SUCCESS"].copy()
    children=[]
    candidate_children=[]
    for cfg_path in d2.rglob("config.json"):
        try:
            cfg=json.loads(cfg_path.read_text())
            candidate_children.append((cfg_path.parent,cfg))
        except Exception:
            continue
    for (rep,gamma),g in failed.groupby(["replicate","gamma"],sort=True):
        cell=label(int(rep),float(gamma)); parent=d2/"replicates"
        analysis_seed=int(g.analysis_seed.iloc[0])
        matches=[(path,cfg) for path,cfg in candidate_children
                 if int(cfg.get("seed",-1))==analysis_seed
                 and float(cfg.get("curvature_gamma",float("nan")))==float(gamma)
                 and int(cfg.get("inclusion_model_degree",1))==2]
        child=matches[0][0] if len(matches)==1 else parent/cell
        def read(path):
            return path.read_text(errors="replace")[-12000:] if path.exists() else None
        mode_results=[]
        if child.exists() and (child/"results"/"log_likelihoods.npz").exists():
            for n in (100,300,500):
                try: mode_results.append(mode_diagnostics(child,n))
                except Exception as exc: mode_results.append({"validation_size":n,"diagnostic_error":f"{type(exc).__name__}: {exc}"})
        children.append({
          "replicate":int(rep),"gamma":float(gamma),"analysis_seed":analysis_seed,
          "resolved_child_path":str(child),"matching_config_count":len(matches),
          "status_values":sorted(g.status.astype(str).unique().tolist()),
          "audit_status_values":sorted(g.audit_status.astype(str).unique().tolist()),
          "ledger_error":str(g.error.iloc[0])[:5000],
          "runner_stderr_tail":read(parent/(cell+".stderr.txt")),
          "child_audit_stderr_tail":read(parent/(cell+".audit.stderr.txt")),
          "child_audit_stdout_tail":read(parent/(cell+".audit.stdout.txt")),
          "mode_diagnostics":mode_results,
        })
    counts=failed.groupby(["gamma","error"],dropna=False).size().reset_index(name="failed_cells")
    report={
      "status":"DIAGNOSTIC_ONLY" if len(failed) else "NO_FAILURES",
      "failed_cells":int(len(failed)),"planned_cells":int(len(ledger)),
      "affected_replicate_gamma_runs":len(children),
      "failed_cells_by_gamma":{str(float(k)):int(v) for k,v in failed.groupby("gamma").size().items()},
      "error_groups":[{"gamma":float(r.gamma),"failed_cells":int(r.failed_cells),"error":str(r.error)[:5000]} for r in counts.itertuples()],
      "children":children,
      "interpretation":"Frozen artifact diagnosis only. Numerical residuals and optimizer differences are diagnostic; no incomplete result is promoted to a scientific result."
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=="__main__": main()
