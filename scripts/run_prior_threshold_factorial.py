"""Factorial prior-support × detection-threshold sensitivity experiment.

Synthetic toy model only. Each detected observation is reused across all prior,
threshold, and likelihood-method configurations for its true H0 value.
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected
from dark_siren_robustness.threshold_misspecification_helpers import wilson_interval

HOST_REDSHIFTS=np.array([0.02,0.04])
HOST_WEIGHTS=np.array([0.7,0.3])
SIGMA=20.0
TRUE_MIN=0.0
TRUE_MAX=120.0
ASSUMED_MAXIMA=(120.0,150.0,180.0,210.0)
TRUTHS=(52.0,70.0,88.0)
PRIORS=((50.0,90.0),(40.0,100.0),(30.0,110.0))
METHODS=("uncorrected","corrected")

def draw_detected(truth,rng,max_attempts=100000):
    for _ in range(max_attempts):
        host=int(rng.choice(len(HOST_REDSHIFTS),p=HOST_WEIGHTS))
        distance=float(rng.normal(C_KM_S*HOST_REDSHIFTS[host]/truth,SIGMA))
        if TRUE_MIN < distance < TRUE_MAX: return host,distance
    raise RuntimeError("Detection rejection sampler exceeded attempt limit.")

def wilson(successes,n,z=1.959963984540054):
    if n<=0: raise ValueError("n must be positive")
    p=successes/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (0.0 if successes==0 else max(0.0,float(center-half)),
            1.0 if successes==n else min(1.0,float(center+half)))

def run(replicates=300,seed=20261028,grid_size=4001):
    if replicates<100: raise ValueError("Use at least 100 replicates per truth.")
    if grid_size<1001: raise ValueError("grid_size must be at least 1001.")
    rng=np.random.default_rng(seed); grid=np.linspace(25.0,125.0,grid_size)
    rows=[]
    for truth in TRUTHS:
        obs=[draw_detected(truth,rng) for _ in range(replicates)]
        for rep,(host,distance) in enumerate(obs):
            for pmin,pmax in PRIORS:
                for assumed_max in ASSUMED_MAXIMA:
                    uncorrected=posterior_h0(grid,HOST_REDSHIFTS,distance,SIGMA,
                        weights=HOST_WEIGHTS,prior_min=pmin,prior_max=pmax)
                    for method in METHODS:
                        if method=="uncorrected":
                            post=uncorrected
                        else:
                            post=posterior_h0_selected(grid,HOST_REDSHIFTS,distance,SIGMA,
                                assumed_max,weights=HOST_WEIGHTS,prior_min=pmin,prior_max=pmax,
                                detection_distance_min=TRUE_MIN)
                        sm=posterior_summary(grid,post)
                        rows.append({"true_H0":truth,"replicate":rep,"host_index":host,
                            "observed_distance_mpc":distance,"true_detection_max_mpc":TRUE_MAX,
                            "assumed_detection_max_mpc":assumed_max,"prior_min":pmin,"prior_max":pmax,
                            "method":method,"posterior_median":float(sm["median"]),
                            "lower_68":float(sm["lower_68"]),"upper_68":float(sm["upper_68"]),
                            "lower_95":float(sm["lower_95"]),"upper_95":float(sm["upper_95"])})
    expected=len(TRUTHS)*replicates*len(PRIORS)*len(ASSUMED_MAXIMA)*len(METHODS)
    if len(rows)!=expected: raise RuntimeError(f"Expected {expected} rows, got {len(rows)}")
    return rows

def summarize(rows):
    keys=("true_H0","prior_min","prior_max","assumed_detection_max_mpc","method")
    groups={}
    for row in rows: groups.setdefault(tuple(row[k] for k in keys),[]).append(row)
    out=[]
    for key,cell in sorted(groups.items()):
        truth,pmin,pmax,threshold,method=key
        errors=np.array([r["posterior_median"]-truth for r in cell])
        h68=np.array([r["lower_68"]<=truth<=r["upper_68"] for r in cell])
        h95=np.array([r["lower_95"]<=truth<=r["upper_95"] for r in cell])
        l68,u68=wilson(int(h68.sum()),len(cell));l95,u95=wilson(int(h95.sum()),len(cell))
        out.append({"true_H0":truth,"n":len(cell),"prior_min":pmin,"prior_max":pmax,
            "assumed_detection_max_mpc":threshold,"method":method,
            "bias":float(errors.mean()),"mae":float(np.mean(np.abs(errors))),
            "rmse":float(np.sqrt(np.mean(errors**2))),
            "coverage_68":float(h68.mean()),"coverage_68_wilson_low":l68,"coverage_68_wilson_high":u68,
            "coverage_95":float(h95.mean()),"coverage_95_wilson_low":l95,"coverage_95_wilson_high":u95,
            "mean_width_68":float(np.mean([r["upper_68"]-r["lower_68"] for r in cell])),
            "mean_width_95":float(np.mean([r["upper_95"]-r["lower_95"] for r in cell]))})
    if len(out)!=len(TRUTHS)*len(PRIORS)*len(ASSUMED_MAXIMA)*len(METHODS):
        raise RuntimeError("Summary grid is incomplete.")
    return out

def paired_contrasts(rows,bootstrap_draws=3000,seed=20261029):
    rng=np.random.default_rng(seed)
    lookup={(r["true_H0"],r["replicate"],r["prior_min"],r["prior_max"],r["assumed_detection_max_mpc"],r["method"]):r for r in rows}
    baseline_prior=(50.0,90.0);baseline_threshold=120.0
    metrics={"median":lambda r:r["posterior_median"],
      "absolute_error":lambda r:abs(r["posterior_median"]-r["true_H0"]),
      "width_68":lambda r:r["upper_68"]-r["lower_68"],
      "width_95":lambda r:r["upper_95"]-r["lower_95"],
      "coverage_68":lambda r:float(r["lower_68"]<=r["true_H0"]<=r["upper_68"]),
      "coverage_95":lambda r:float(r["lower_95"]<=r["true_H0"]<=r["upper_95"])}
    out=[]
    for truth in TRUTHS:
      for pmin,pmax in PRIORS:
       for threshold in ASSUMED_MAXIMA:
        if (pmin,pmax)==baseline_prior and threshold==baseline_threshold: continue
        for method in METHODS:
         a=[lookup[(truth,i,*baseline_prior,baseline_threshold,method)] for i in range(max(r["replicate"] for r in rows if r["true_H0"]==truth)+1)]
         b=[lookup[(truth,i,pmin,pmax,threshold,method)] for i in range(len(a))]
         for metric,fn in metrics.items():
          delta=np.array([fn(y)-fn(x) for x,y in zip(a,b)],float)
          boot=np.empty(bootstrap_draws)
          for j in range(bootstrap_draws):
           boot[j]=np.mean(delta[rng.integers(0,len(delta),len(delta))])
          out.append({"true_H0":truth,"method":method,"metric":metric,
            "baseline_prior":"[50,90]","comparison_prior":f"[{pmin:g},{pmax:g}]",
            "baseline_threshold_mpc":baseline_threshold,"comparison_threshold_mpc":threshold,
            "n_pairs":len(delta),"mean_difference_vs_baseline":float(delta.mean()),
            "bootstrap_95_low":float(np.quantile(boot,.025)),
            "bootstrap_95_high":float(np.quantile(boot,.975)),
            "bootstrap_draws":bootstrap_draws,"bootstrap_seed":seed})
    expected=len(TRUTHS)*((len(PRIORS)*len(ASSUMED_MAXIMA))-1)*len(METHODS)*len(metrics)
    if len(out)!=expected: raise RuntimeError(f"Expected {expected} contrasts, got {len(out)}")
    return out

def write_csv(path,rows):
    if not rows: raise RuntimeError(f"Refusing empty artifact {path}")
    path.parent.mkdir(parents=True,exist_ok=True); fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f); check=list(reader)
        if reader.fieldnames!=fields or len(check)!=len(rows): raise RuntimeError(f"CSV round-trip failed: {path}")

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261028);p.add_argument("--grid-size",type=int,default=4001)
    p.add_argument("--output-dir",default="results/prior-threshold-factorial")
    a=p.parse_args();rows=run(a.replicates,a.seed,a.grid_size);summaries=summarize(rows);contrasts=paired_contrasts(rows)
    out=Path(a.output_dir)
    for name,data in (("per_run.csv",rows),("summary.csv",summaries),("paired_contrasts.csv",contrasts)):
        write_csv(out/name,data)
    print(f"PASS: {len(rows)} per-run rows; {len(summaries)} summaries; {len(contrasts)} paired contrasts")
    for r in summaries: print(r)
if __name__=="__main__": main()
