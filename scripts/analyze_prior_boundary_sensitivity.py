"""Paired bootstrap contrasts for prior-boundary sensitivity results."""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np

BASELINE="baseline_50_90"
COMPARISONS=("wide_40_100","very_wide_30_120")
METRICS=("median_error","covers_68","covers_95","width_68","width_95",
         "mass_within_1_of_lower_prior","mass_within_1_of_upper_prior",
         "mass_within_2_of_lower_prior","mass_within_2_of_upper_prior",
         "mass_below_50","mass_above_90","mass_outside_baseline_50_90")

def ci_mean_difference(delta, rng, n_boot):
    n=len(delta)
    idx=rng.integers(0,n,size=(n_boot,n))
    means=delta[idx].mean(axis=1)
    return float(delta.mean()),float(np.quantile(means,.025)),float(np.quantile(means,.975))

def ci_rmse_difference(err_a,err_b,rng,n_boot):
    n=len(err_a)
    idx=rng.integers(0,n,size=(n_boot,n))
    rmse_a=np.sqrt(np.mean(err_a[idx]**2,axis=1))
    rmse_b=np.sqrt(np.mean(err_b[idx]**2,axis=1))
    delta=rmse_a-rmse_b
    return float(np.sqrt(np.mean(err_a**2))-np.sqrt(np.mean(err_b**2))),float(np.quantile(delta,.025)),float(np.quantile(delta,.975))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",default="results/prior-boundary-sensitivity/per_run.csv")
    p.add_argument("--output",default="results/prior-boundary-sensitivity/paired_contrasts.csv")
    p.add_argument("--bootstrap",type=int,default=5000)
    p.add_argument("--seed",type=int,default=20261025)
    a=p.parse_args()
    with Path(a.input).open(newline="",encoding="utf-8") as f: rows=list(csv.DictReader(f))
    keys=("true_H0","scenario","prior","replicate")
    by={}
    for row in rows:
        key=(float(row["true_H0"]),row["scenario"],row["prior"],int(row["replicate"]))
        by[key]=row
    rng=np.random.default_rng(a.seed)
    output=[]
    for h0 in (55.,60.,65.,70.,75.,80.,85.):
      for scenario in ("known_host","unknown_host"):
       for comp in COMPARISONS:
        base_rows=[]; comp_rows=[]
        for rep in range(300):
          kb=(h0,scenario,BASELINE,rep); kc=(h0,scenario,comp,rep)
          if kb not in by or kc not in by: raise ValueError(f"missing paired row: {kb} or {kc}")
          base_rows.append(by[kb]);comp_rows.append(by[kc])
        item={"true_H0":h0,"scenario":scenario,"baseline_prior":BASELINE,"comparison_prior":comp,"replicates":300}
        for metric in METRICS:
          # Compare alternative minus baseline, preserving replicate pairing.
          da=np.array([float(r[metric]) for r in comp_rows])
          db=np.array([float(r[metric]) for r in base_rows])
          mean,lo,hi=ci_mean_difference(da-db,rng,a.bootstrap)
          item[f"delta_{metric}_mean"]=mean
          item[f"delta_{metric}_ci95_low"]=lo
          item[f"delta_{metric}_ci95_high"]=hi
        rm,lo,hi=ci_rmse_difference(
          np.array([float(r["median_error"]) for r in comp_rows]),
          np.array([float(r["median_error"]) for r in base_rows]),rng,a.bootstrap)
        item["delta_rmse_median"]=rm
        item["delta_rmse_median_ci95_low"]=lo
        item["delta_rmse_median_ci95_high"]=hi
        output.append(item)
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    with Path(a.output).open("w",newline="",encoding="utf-8") as f:
      w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    assert len(output)==28
    print(f"Wrote {len(output)} paired prior contrasts with {a.bootstrap} bootstrap replicates each.")
if __name__=="__main__": main()
