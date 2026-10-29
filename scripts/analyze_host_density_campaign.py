"""Run paired host-density estimation, bandwidth, grid, and selection misspecification sweeps."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from dark_siren_robustness.host_density_benchmark import run_benchmark


def wilson(k, n, z=1.959964):
    p=k/n; d=1+z*z/n
    c=(p+z*z/(2*n))/d
    h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return float(c-h), float(c+h)


def summarise(frame, h0_true, n_boot=4000, seed=1):
    rng=np.random.default_rng(seed)
    out=[]
    piv=frame.pivot(index="replicate",columns="arm",values="posterior_mean")-h0_true
    for arm,g in frame.groupby("arm",sort=False):
        e=g.posterior_mean.to_numpy()-h0_true
        n=len(e); k=int(g.covered_68.sum()); lo,hi=wilson(k,n)
        row={"arm":arm,"n":n,"bias":float(e.mean()),"bias_se":float(e.std(ddof=1)/np.sqrt(n)),
             "rmse":float(np.sqrt(np.mean(e**2))),"coverage68":k/n,"cov_wilson95":[lo,hi],
             "width":float(g.ci68_width.mean()),"l1_density":float(g.l1_host_density_error.mean()),
             "max_weight_ratio_max":float(g.max_weight_ratio.max()),"ess_fraction_min":float(g.ess_fraction.min())}
        if arm!="oracle":
            d=(piv[arm].abs()-piv["oracle"].abs()).to_numpy()
            b=np.empty(n_boot)
            for start in range(0,n_boot,500):
                stop=min(start+500,n_boot)
                ids=rng.integers(0,n,size=(stop-start,n))
                b[start:stop]=d[ids].mean(axis=1)
            row["paired_abs_err_minus_oracle"]=float(d.mean())
            row["paired_ci95"]=[float(np.quantile(b,.025)),float(np.quantile(b,.975))]
        out.append(row)
    return out


# Fine-grid setting is a practical convergence check relative to the 161 x 801 default,
# not a claim of asymptotic convergence.
CONFIGS={
 "bw_scott":{"bandwidth_factor":None},
 "bw_0.05":{"bandwidth_factor":0.05},
 "bw_0.1":{"bandwidth_factor":0.1},
 "bw_0.2":{"bandwidth_factor":0.2},
 "bw_0.5":{"bandwidth_factor":0.5},
 "bw_0.8":{"bandwidth_factor":0.8},
 "bw_1.2":{"bandwidth_factor":1.2},
 "grid_coarse":{"bandwidth_factor":None,"h0_points":81,"z_points":401},
 "grid_fine":{"bandwidth_factor":None,"h0_points":201,"z_points":1001},
 "misspec_mid_0.18":{"bandwidth_factor":None,"supplied_midpoint":0.18},
 "misspec_mid_0.26":{"bandwidth_factor":None,"supplied_midpoint":0.26},
 "misspec_width_0.10":{"bandwidth_factor":None,"supplied_width":0.10},
 "misspec_width_0.025":{"bandwidth_factor":None,"supplied_width":0.025},
}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates",type=int,default=100)
    parser.add_argument("--seed",type=int,default=20261029)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--only",nargs="*")
    args=parser.parse_args()
    if args.replicates<2: parser.error("--replicates must be at least two")
    args.out.mkdir(parents=True,exist_ok=True)
    selected=[name for name in CONFIGS if not args.only or name in args.only]
    unknown=set(args.only or [])-set(CONFIGS)
    if unknown: parser.error(f"unknown --only configurations: {sorted(unknown)}")
    results={}
    for name in selected:
        frame,meta=run_benchmark(n_replicates=args.replicates,seed=args.seed,**CONFIGS[name])
        frame.drop(columns=["posterior_density_grid_json"],errors="ignore").to_csv(args.out/f"{name}_replicates.csv",index=False)
        results[name]={"config":CONFIGS[name],"arms":summarise(frame,meta["h0_true_km_s_mpc"])}
        print(name,[(r["arm"],round(r["bias"],2),round(r["rmse"],2),round(r["coverage68"],3)) for r in results[name]["arms"]],flush=True)
    summary_path=args.out/"campaign_summary.json"
    summary={"status":"synthetic_exploratory_host_density_campaign_not_real_H0_measurement",
             "replicates_per_configuration":args.replicates,"seed_start":args.seed,
             "seed_end":args.seed+args.replicates-1,"n_parent":6000,"n_detected":40,
             "configuration_names":selected,"results":results,
             "interpretation_limit":"Results apply only to the declared low-redshift fixed-H0 Gaussian toy."}
    summary_path.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    file_records={}
    for path in sorted(args.out.glob("*")):
        if path.is_file() and path.name!="campaign_manifest.json":
            file_records[path.name]={"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"bytes":path.stat().st_size}
    manifest={"schema_version":1,"seed_start":args.seed,"seed_end":args.seed+args.replicates-1,
              "replicates_per_configuration":args.replicates,"configuration_names":selected,"files":file_records}
    (args.out/"campaign_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
