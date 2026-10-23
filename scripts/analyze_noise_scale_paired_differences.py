"""Paired bootstrap comparisons for the noise-scale audit output."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np

METRICS = {
    "mean_signed_error": lambda r: float(r["median_error"]),
    "coverage_68": lambda r: float(r["covers_68"]),
    "coverage_95": lambda r: float(r["covers_95"]),
    "mean_width_68": lambda r: float(r["width_68"]),
    "mean_width_95": lambda r: float(r["width_95"]),
}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",default="results/noise-scale-consistency/per_run.csv")
    p.add_argument("--output",default="results/noise-scale-consistency/paired_bootstrap.csv")
    p.add_argument("--bootstrap",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261023)
    a=p.parse_args()
    if a.bootstrap < 1000: p.error("--bootstrap must be at least 1000")
    rows=list(csv.DictReader(open(a.input,newline="",encoding="utf-8")))
    by={}
    for r in rows: by.setdefault((float(r["true_H0"]),int(r["replicate"])),{})[r["mode"]]=r
    rng=np.random.default_rng(a.seed)
    out=[]
    for h0 in sorted({k[0] for k in by}):
        pairs=[v for (h,rep),v in sorted(by.items()) if h==h0]
        assert all("fixed_reported_sigma" in x and "sigma_depends_on_h0" in x for x in pairs)
        # Positive/negative deltas are H0-dependent mode minus fixed-sigma mode.
        for metric,fn in METRICS.items():
            fixed=np.array([fn(x["fixed_reported_sigma"]) for x in pairs])
            dependent=np.array([fn(x["sigma_depends_on_h0"]) for x in pairs])
            delta=dependent-fixed
            draws=np.empty(a.bootstrap)
            for start in range(0,a.bootstrap,250):
                stop=min(start+250,a.bootstrap)
                idx=rng.integers(0,len(delta),size=(stop-start,len(delta)))
                draws[start:stop]=delta[idx].mean(axis=1)
            lo,hi=np.quantile(draws,[0.025,0.975])
            out.append({"true_H0":h0,"metric":metric,"replicates":len(pairs),
                        "delta_dependent_minus_fixed":float(delta.mean()),
                        "bootstrap_95_low":float(lo),"bootstrap_95_high":float(hi),
                        "bootstrap_draws":a.bootstrap,"seed":a.seed})
        # RMSE is nonlinear; bootstrap the difference in RMSE from paired squared errors.
        e0=np.array([float(x["fixed_reported_sigma"]["median_error"]) for x in pairs])
        e1=np.array([float(x["sigma_depends_on_h0"]["median_error"]) for x in pairs])
        draws=np.empty(a.bootstrap)
        for start in range(0,a.bootstrap,250):
            stop=min(start+250,a.bootstrap)
            idx=rng.integers(0,len(e0),size=(stop-start,len(e0)))
            draws[start:stop]=np.sqrt(np.mean(e1[idx]**2,axis=1))-np.sqrt(np.mean(e0[idx]**2,axis=1))
        lo,hi=np.quantile(draws,[0.025,0.975])
        out.append({"true_H0":h0,"metric":"rmse_median",
                    "replicates":len(pairs),
                    "delta_dependent_minus_fixed":float(np.sqrt(np.mean(e1**2))-np.sqrt(np.mean(e0**2))),
                    "bootstrap_95_low":float(lo),"bootstrap_95_high":float(hi),
                    "bootstrap_draws":a.bootstrap,"seed":a.seed})
    path=Path(a.output); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"Wrote {len(out)} paired metric comparisons to {path}")
    for r in out:
        if r["metric"] in ("mean_signed_error","rmse_median","coverage_68"):
            print(r)

if __name__=="__main__": main()
