"""Cluster-aware uncertainty analysis for the full-parent fixed-count H0 baseline.

The independent sampling units are archived parent populations, not the many
event realizations nested within each parent. This script reports both
conditional event-level summaries and parent-cluster uncertainty.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t


def mean_interval(x: np.ndarray) -> dict:
    values=np.asarray(x,dtype=float)
    values=values[np.isfinite(values)]
    n=len(values)
    if not n:
        return {"n":0,"mean":float("nan"),"sd":float("nan"),
                "mcse":float("nan"),"t95_low":float("nan"),"t95_high":float("nan")}
    mean=float(values.mean())
    sd=float(values.std(ddof=1)) if n>1 else float("nan")
    se=sd/np.sqrt(n) if n>1 else float("nan")
    crit=float(t.ppf(.975,n-1)) if n>1 else float("nan")
    return {"n":n,"mean":mean,"sd":sd,"mcse":float(se),
            "t95_low":float(mean-crit*se) if n>1 else float("nan"),
            "t95_high":float(mean+crit*se) if n>1 else float("nan")}


def cluster_bootstrap_interval(values: np.ndarray, rng: np.random.Generator,
                               draws: int=20000) -> tuple[float,float]:
    x=np.asarray(values,dtype=float)
    x=x[np.isfinite(x)]
    if len(x)<2:
        return float("nan"),float("nan")
    indices=rng.integers(0,len(x),size=(draws,len(x)))
    means=x[indices].mean(axis=1)
    return tuple(map(float,np.quantile(means,[.025,.975])))


def analyze(root: Path, output: Path, seed: int=20261119) -> dict:
    root=root.resolve()
    frame=pd.read_csv(root/"replicates.csv")
    required={"parent_replicate","event_replicate","median_bias",
              "covered68","covered90","covered95","width68"}
    missing=required-set(frame.columns)
    if missing: raise ValueError(f"replicates file missing columns: {sorted(missing)}")
    if frame.duplicated(["parent_replicate","event_replicate"]).any():
        raise ValueError("duplicate parent/event replicate key")
    if frame.empty: raise ValueError("replicates file is empty")
    rng=np.random.default_rng(seed)
    metrics=("median_bias","covered68","covered90","covered95","width68")
    grouped=frame.groupby("parent_replicate",sort=True)
    parent=grouped.agg(
        event_draws=("event_replicate","size"),
        mean_median_bias=("median_bias","mean"),
        coverage68=("covered68","mean"),
        coverage90=("covered90","mean"),
        coverage95=("covered95","mean"),
        mean_width68=("width68","mean"),
    ).reset_index()
    output.mkdir(parents=True,exist_ok=True)
    parent.to_csv(output/"parent_cluster_summary.csv",index=False)
    n_parent=len(parent)
    results={}
    for col,label in [
        ("mean_median_bias","parent_mean_median_bias"),
        ("coverage68","parent_coverage68"),
        ("coverage90","parent_coverage90"),
        ("coverage95","parent_coverage95"),
        ("mean_width68","parent_mean_width68"),
    ]:
        x=parent[col].to_numpy(float)
        result=mean_interval(x)
        lo,hi=cluster_bootstrap_interval(x,rng)
        result["cluster_bootstrap95_low"]=lo
        result["cluster_bootstrap95_high"]=hi
        results[label]=result
    # This event-level standard error is conditional on the archived parent
    # populations and is not the unconditional population-level MCSE.
    event_bias=frame.median_bias.to_numpy(float)
    event_result=mean_interval(event_bias)
    results["event_weighted_median_bias_conditional_on_archived_parents"]=event_result
    output.mkdir(parents=True,exist_ok=True)
    report={
        "status":"PASS",
        "parent_populations":n_parent,
        "total_event_realizations":len(frame),
        "event_draws_per_parent":{str(k):int(v) for k,v in
            grouped.size().to_dict().items()},
        "independent_unit":"parent population",
        "uncertainty":results,
        "interpretation":"Cluster-aware uncertainty uses parent populations as independent units. Event-level uncertainty is conditional on the archived parent populations. This is a fixed-count synthetic diagnostic, not joint SBC or observational H0 inference."
    }
    (output/"cluster_uncertainty.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("baseline_output",type=Path)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261119)
    a=p.parse_args()
    print(json.dumps(analyze(a.baseline_output,a.output_dir,a.seed),indent=2))


if __name__=="__main__": main()
