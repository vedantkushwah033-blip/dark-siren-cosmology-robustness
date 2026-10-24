"""Paired fixed-truth single-host selection control; synthetic diagnostic only."""
from __future__ import annotations
import argparse,csv,json,platform
from pathlib import Path
import numpy as np, scipy
from dark_siren_robustness.hardened_inference import C_KM_S,posterior_h0,posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected
Z=0.03; SIGMA=20.; DMIN=0.; DMAX=180.; PLO=50.; PHI=90.
TRUTHS=(52.,55.,60.,65.,70.,75.,80.,85.,88.)
def draw(truth,rng,max_attempts=100000):
    if not np.isfinite(truth) or truth<=0: raise ValueError("truth must be finite and positive")
    mu=C_KM_S*Z/truth
    for _ in range(max_attempts):
        d=float(rng.normal(mu,SIGMA))
        if DMIN<d<DMAX:return d
    raise RuntimeError("detection rejection sampler exceeded limit")
def wilson(k,n,z=1.959963984540054):
    p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (0. if k==0 else max(0.,float(center-half)),1. if k==n else min(1.,float(center+half)))
def run(replicates=300,seed=20261030,grid_size=4001,bootstrap_draws=5000):
    if replicates<100:raise ValueError("Use at least 100 replicates per truth.")
    if grid_size<1001:raise ValueError("grid_size must be at least 1001.")
    if bootstrap_draws<1000:raise ValueError("Use at least 1000 bootstrap draws.")
    rng=np.random.default_rng(seed);grid=np.linspace(25.,125.,grid_size);rows=[]
    for truth in TRUTHS:
        for rep in range(replicates):
            d=draw(truth,rng)
            base=posterior_h0(grid,np.array([Z]),d,SIGMA,weights=np.array([1.]),prior_min=PLO,prior_max=PHI)
            corr=posterior_h0_selected(grid,np.array([Z]),d,SIGMA,DMAX,weights=np.array([1.]),prior_min=PLO,prior_max=PHI,detection_distance_min=DMIN)
            bs,cs=posterior_summary(grid,base),posterior_summary(grid,corr)
            row={"true_H0":truth,"replicate":rep,"observed_distance_mpc":d}
            for name,sm in (("uncorrected",bs),("corrected",cs)):
                for k in ("median","lower_68","upper_68","lower_95","upper_95"):row[f"{name}_{k}"]=float(sm[k])
                err=float(sm["median"]-truth);row[f"{name}_error"]=err;row[f"{name}_absolute_error"]=abs(err)
            row["paired_abs_error_change_corrected_minus_uncorrected"]=row["corrected_absolute_error"]-row["uncorrected_absolute_error"];rows.append(row)
    summaries=[];brng=np.random.default_rng(seed+1)
    for truth in TRUTHS:
        cell=[r for r in rows if r["true_H0"]==truth]
        for name in ("uncorrected","corrected"):
            e=np.array([r[f"{name}_error"] for r in cell])
            h68=np.array([r[f"{name}_lower_68"]<=truth<=r[f"{name}_upper_68"] for r in cell])
            h95=np.array([r[f"{name}_lower_95"]<=truth<=r[f"{name}_upper_95"] for r in cell])
            l68,u68=wilson(int(h68.sum()),len(cell));l95,u95=wilson(int(h95.sum()),len(cell))
            summaries.append({"true_H0":truth,"method":name,"n":len(cell),"bias":float(e.mean()),"mae":float(np.abs(e).mean()),"rmse":float(np.sqrt(np.mean(e**2))),
              "coverage_68":float(h68.mean()),"coverage_68_wilson_low":l68,"coverage_68_wilson_high":u68,
              "coverage_95":float(h95.mean()),"coverage_95_wilson_low":l95,"coverage_95_wilson_high":u95,
              "mean_width_68":float(np.mean([r[f"{name}_upper_68"]-r[f"{name}_lower_68"] for r in cell])),
              "mean_width_95":float(np.mean([r[f"{name}_upper_95"]-r[f"{name}_lower_95"] for r in cell]))})
        delta=np.array([r["paired_abs_error_change_corrected_minus_uncorrected"] for r in cell]);boot=np.empty(bootstrap_draws)
        for start in range(0,bootstrap_draws,256):
            stop=min(start+256,bootstrap_draws);idx=brng.integers(0,len(delta),size=(stop-start,len(delta)));boot[start:stop]=delta[idx].mean(axis=1)
        summaries.append({"true_H0":truth,"method":"paired_corrected_minus_uncorrected_absolute_error","n":len(cell),"mean_paired_change":float(delta.mean()),"bootstrap_95_low":float(np.quantile(boot,.025)),"bootstrap_95_high":float(np.quantile(boot,.975)),"bootstrap_draws":bootstrap_draws,"bootstrap_seed":seed+1})
    manifest={"experiment":"single_host_fixed_truth_selection_control","scope":"synthetic diagnostic; not an observational measurement","seed":seed,"replicates_per_truth":replicates,"grid_size":grid_size,"grid_bounds":[25.,125.],"prior":[PLO,PHI],"truths":list(TRUTHS),"host_redshift":Z,"sigma_mpc":SIGMA,"detection_interval_mpc":[DMIN,DMAX],"bootstrap_draws":bootstrap_draws,"bootstrap_seed":seed+1,"raw_rows":len(rows),"summary_rows":len(summaries),"python_version":platform.python_version(),"numpy_version":np.__version__,"scipy_version":scipy.__version__}
    return rows,summaries,manifest
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--replicates",type=int,default=300);p.add_argument("--seed",type=int,default=20261030);p.add_argument("--grid-size",type=int,default=4001);p.add_argument("--bootstrap-draws",type=int,default=5000);p.add_argument("--output-dir",default="results/single-host-fixed-truth-control");a=p.parse_args()
    rows,summaries,manifest=run(a.replicates,a.seed,a.grid_size,a.bootstrap_draws);out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for name,data in (("per_run.csv",rows),("summary.csv",summaries)):
        with (out/name).open("w",newline="",encoding="utf-8") as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in data for k in r)));w.writeheader();w.writerows(data)
    (out/"experiment_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(f"PASS: {len(rows)} per-run rows; {len(summaries)} summary rows.")
    for r in summaries:print(r)
if __name__=="__main__":main()
