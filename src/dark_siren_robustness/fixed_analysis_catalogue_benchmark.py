"""Fixed-analysis-catalogue finite-validation density benchmark.

Independent validation frame; analysis catalogue is frozen across validation sizes.
Synthetic density-recovery only: no H0 inference or posterior-calibration claim.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import trapezoid
from scipy.special import expit
from dark_siren_robustness.estimated_inclusion import fit_inclusion_logistic
from dark_siren_robustness.finite_validation_density_recovery import (
    Z_MIN, Z_MAX, _draw_parent_redshifts, _curvature_center,
    _calibrate_intercept, _analytic_target, _safe_density, _weight_diagnostics,
)

ARMS=("catalogue_only","known_inclusion_ipw","parent_sample_kde")
DEGREES=(1,2)

def run_fixed_catalogue_benchmark(*, n_replicates=20, n_analysis=4000,
        n_validation_max=1000, validation_sizes=(100,300,1000),
        gamma_values=(-1.2,0.0,1.2), retention_values=(0.25,),
        degrees=(1,2), bandwidth_multipliers=(0.15,), grid_points=801,
        seed=20261109):
    if n_replicates<1 or n_analysis<1100 or n_validation_max<5 or grid_points<201:
        raise ValueError("replicates>=1, analysis size>=1100, validation max>=5, grid>=201 required")
    if tuple(sorted(set(validation_sizes)))!=tuple(validation_sizes) or not validation_sizes:
        raise ValueError("validation_sizes must be unique and increasing")
    if validation_sizes[-1]>n_validation_max or validation_sizes[0]<5:
        raise ValueError("validation sizes must lie in [5,n_validation_max]")
    if any(d not in (1,2) for d in degrees) or not degrees:
        raise ValueError("degrees must contain 1 and/or 2")
    if any(not 0<r<1 for r in retention_values) or not retention_values:
        raise ValueError("retention values must lie in (0,1)")
    if any(not np.isfinite(g) for g in gamma_values) or not gamma_values:
        raise ValueError("gamma values must be finite and nonempty")
    if any(not np.isfinite(h) or h<=0 for h in bandwidth_multipliers):
        raise ValueError("bandwidth multipliers must be positive and finite")
    grid=np.linspace(Z_MIN,Z_MAX,grid_points)
    merger=1+1.5*(grid-Z_MIN)/(Z_MAX-Z_MIN)
    truth=_analytic_target(grid); center=_curvature_center()
    unit_rows=[]; attempt_rows=[]
    for rep in range(n_replicates):
      # Independent streams for the analysis and validation frames.
      az=_draw_parent_redshifts(np.random.default_rng(seed+1000003*rep),n_analysis)
      vz=_draw_parent_redshifts(np.random.default_rng(seed+2000003*rep),n_validation_max)
      order=np.random.default_rng(seed+3000017*rep).permutation(n_validation_max)
      rank=np.empty(n_validation_max,dtype=int); rank[order]=np.arange(n_validation_max)
      au=np.random.default_rng(seed+4000037*rep).random(n_analysis)
      vu=np.random.default_rng(seed+5000081*rep).random(n_validation_max)
      for gamma in gamma_values:
       for retention in retention_values:
        alpha=_calibrate_intercept(float(gamma),float(retention),center)
        ax=(az-Z_MIN)/(Z_MAX-Z_MIN); vx=(vz-Z_MIN)/(Z_MAX-Z_MIN)
        api=expit(alpha-ax+gamma*(ax**2-center))
        vpi=expit(alpha-vx+gamma*(vx**2-center))
        ay=(au<api).astype(int); vy=(vu<vpi).astype(int)
        aidx=np.flatnonzero(ay==1); acat=az[aidx]; acatpi=api[aidx]
        for i in range(n_analysis):
         unit_rows.append(dict(replicate=rep,gamma=float(gamma),target_retention=float(retention),
          frame="ANALYSIS",object_id=i,redshift=float(az[i]),retention_uniform=float(au[i]),
          true_inclusion_probability=float(api[i]),retained=int(ay[i]),validation_rank=-1))
        for i in range(n_validation_max):
         unit_rows.append(dict(replicate=rep,gamma=float(gamma),target_retention=float(retention),
          frame="VALIDATION",object_id=i,redshift=float(vz[i]),retention_uniform=float(vu[i]),
          true_inclusion_probability=float(vpi[i]),retained=int(vy[i]),validation_rank=int(rank[i])))
        for nval in validation_sizes:
         vidx=order[:nval]; vz_n=vz[vidx]; vy_n=vy[vidx]
         fits={d:fit_inclusion_logistic(vz_n,vy_n,degree=d) for d in degrees}
         for h in bandwidth_multipliers:
          base={"catalogue_only":(acat,None),"known_inclusion_ipw":(acat,acatpi),
                "parent_sample_kde":(az,None)}
          for arm,(sample,prob) in base.items():
           status="SUCCESS"; msg=""; density=None
           try: density=_safe_density(grid,sample,prob,merger,h)
           except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
            status="NON_EVALUABLE"; msg=str(exc)
           metrics=_weight_diagnostics(prob,len(sample)) if len(sample) else {}
           attempt_rows.append(dict(replicate=rep,gamma=float(gamma),target_retention=float(retention),
            n_validation=int(nval),bandwidth_multiplier=float(h),degree=0,arm=arm,status=status,
            message=msg,n_analysis_parent=n_analysis,n_analysis_catalogue=len(acat),
            n_validation_max=n_validation_max,iae=float(trapezoid(abs(density-truth),grid)) if density is not None else np.nan,
            density_json=json.dumps(density.tolist(),separators=(",",":")) if density is not None else "",
            **metrics))
          for d,fit in fits.items():
           arm=f"estimated_degree_{d}_ipw"; status=fit.status; msg=fit.message; density=None; prob=None
           if status=="SUCCESS":
            try:
             prob=fit.predict(acat); density=_safe_density(grid,acat,prob,merger,h)
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
             status="NON_EVALUABLE"; msg=str(exc)
           metrics=_weight_diagnostics(prob,len(acat)) if density is not None else {}
           attempt_rows.append(dict(replicate=rep,gamma=float(gamma),target_retention=float(retention),
            n_validation=int(nval),bandwidth_multiplier=float(h),degree=int(d),arm=arm,status=status,
            message=msg,n_analysis_parent=n_analysis,n_analysis_catalogue=len(acat),
            n_validation_max=n_validation_max,iae=float(trapezoid(abs(density-truth),grid)) if density is not None else np.nan,
            density_json=json.dumps(density.tolist(),separators=(",",":")) if density is not None else "",
            coefficients_json=json.dumps(fit.coefficients.tolist(),separators=(",",":")) if fit.coefficients is not None else "",
            fit_status=fit.status,fit_center=fit.center,fit_scale=fit.scale,
            fit_gradient_inf_norm=fit.gradient_inf_norm,fit_objective=(
             float(np.sum(np.logaddexp(0, np.column_stack([((vz_n-fit.center)/fit.scale)**p for p in range(d+1)])@fit.coefficients)
              -vy_n*(np.column_stack([((vz_n-fit.center)/fit.scale)**p for p in range(d+1)])@fit.coefficients)))
             if fit.coefficients is not None else np.nan),**metrics))
    attempts=pd.DataFrame(attempt_rows); units=pd.DataFrame(unit_rows)
    expected=n_replicates*len(gamma_values)*len(retention_values)*len(validation_sizes)*len(bandwidth_multipliers)*(3+len(degrees))
    if len(attempts)!=expected: raise AssertionError("Internal attempt-factorial construction mismatch")
    status_counts={str(k):int(v) for k,v in attempts.status.value_counts().items()}
    condition_status=[]
    for key,group in attempts.groupby(["gamma","target_retention","n_validation","bandwidth_multiplier","arm"],sort=True):
        condition_status.append({"gamma":float(key[0]),"target_retention":float(key[1]),
         "n_validation":int(key[2]),"bandwidth_multiplier":float(key[3]),"arm":str(key[4]),
         "n_planned":int(len(group)),"n_success":int((group.status=="SUCCESS").sum()),
         "n_failed":int((group.status=="FAILED").sum()),
         "n_non_evaluable":int((group.status=="NON_EVALUABLE").sum()),
         "n_not_run":int((group.status=="NOT_RUN").sum()),
         "success_rate":float((group.status=="SUCCESS").mean()),
         "conditional_success_mean_iae":float(group.loc[group.status=="SUCCESS","iae"].mean())
            if (group.status=="SUCCESS").any() else None})
    paired_contrasts=[]
    group_keys=["gamma","target_retention","n_validation","bandwidth_multiplier"]
    for key,group in attempts.groupby(group_keys,sort=True):
        wide=group.pivot(index="replicate",columns="arm",values="iae")
        if "catalogue_only" not in wide: continue
        for degree in degrees:
            arm=f"estimated_degree_{degree}_ipw"
            if arm not in wide: continue
            pair=wide[["catalogue_only",arm]].dropna()
            diff=(pair[arm]-pair["catalogue_only"]).to_numpy(float)
            paired_contrasts.append({"gamma":float(key[0]),"target_retention":float(key[1]),
             "n_validation":int(key[2]),"bandwidth_multiplier":float(key[3]),"degree":int(degree),
             "n_pairs":int(len(diff)),"n_replicates_planned":int(n_replicates),
             "mean_iae_difference":float(diff.mean()) if len(diff) else None,
             "paired_mc_se":float(diff.std(ddof=1)/np.sqrt(len(diff))) if len(diff)>1 else None})
    summary={"artifact_schema_version":1,"design":"fixed_analysis_catalogue_independent_validation_frame",
     "status":"synthetic_density_recovery_no_H0_inference",
     "estimand":"analytic merger-weighted parent-host redshift density",
     "n_replicates":n_replicates,"n_analysis_parent":n_analysis,"n_validation_max":n_validation_max,
     "validation_sizes":list(validation_sizes),"gamma_values":list(gamma_values),
     "retention_values":list(retention_values),"degrees":list(degrees),
     "bandwidth_multipliers":list(bandwidth_multipliers),"seed_start":seed,
     "grid":[Z_MIN,Z_MAX,grid_points],"curvature_center":center,
     "analytic_target_density":truth.tolist(),"attempt_count":len(attempts),
     "status_counts":status_counts,
     "success_rate_by_arm":{str(k):float((v.status=="SUCCESS").mean()) for k,v in attempts.groupby("arm")},
     "condition_status_ledger":condition_status,"paired_ipw_minus_catalogue_only":paired_contrasts,
     "conditional_success_mean_iae":{str(k):float(v.loc[v.status=="SUCCESS","iae"].mean())
       for k,v in attempts.groupby("arm") if (v.status=="SUCCESS").any()},
     "primary_interpretation":"Validation-size contrasts condition on the identical analysis catalogue within replicate and scenario.",
     "limitations":["Synthetic density recovery only; no H0 inference or posterior calibration.",
      "Validation frame is independently sampled from the parent law with complete binary-label ascertainment.",
      "Oracle inclusion probabilities are diagnostic only.","Model-based inverse-information covariance is not a Bayesian posterior."]}
    return attempts,units,summary

def write_fixed_catalogue_artifacts(attempts,units,summary,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    files={"attempts":outdir/"attempt_metrics.csv","units":outdir/"frame_units.csv",
     "summary":outdir/"summary.json","config":outdir/"config.json"}
    attempts.to_csv(files["attempts"],index=False); units.to_csv(files["units"],index=False)
    config={k:v for k,v in summary.items() if k not in ("status_counts","success_rate_by_arm","conditional_success_mean_iae")}
    files["config"].write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
    files["summary"].write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    manifest={"artifact_schema_version":1,"design":summary["design"],"files":{
     k:{"path":p.name,"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
     for k,p in files.items()},"runtime":{"python":sys.version,"platform":platform.platform(),
     "numpy":np.__version__,"scipy":scipy.__version__,"pandas":pd.__version__,
     "github_sha":os.environ.get("GITHUB_SHA")}}
    (outdir/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=20);p.add_argument("--n-analysis",type=int,default=4000)
    p.add_argument("--n-validation-max",type=int,default=1000);p.add_argument("--seed",type=int,default=20261109)
    p.add_argument("--grid-points",type=int,default=801);p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args(); x,u,s=run_fixed_catalogue_benchmark(n_replicates=a.replicates,n_analysis=a.n_analysis,
      n_validation_max=a.n_validation_max,seed=a.seed,grid_points=a.grid_points)
    m=write_fixed_catalogue_artifacts(x,u,s,a.output_dir)
    print(json.dumps({"summary":s,"manifest":m},indent=2))
if __name__=="__main__": main()
