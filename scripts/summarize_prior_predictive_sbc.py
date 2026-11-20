#!/usr/bin/env python3
"""Summarize posterior ranks and central-interval coverage for prior-predictive SBC."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chisquare, kstest

INTERVALS={
    "50":("q25","q75"),
    "68":("q16","q84"),
    "80":("q10","q90"),
    "90":("q05","q95"),
    "95":("q025","q975"),
}
ARMS=["A_catalogue_only","B_parent_oracle","C_known_probability_ipw","D_map_plugin_ipw","E_uncertainty_integrated_ipw"]

def _wilson(successes:int,n:int,z:float=1.959963984540054)->tuple[float,float]:
    if n<=0:return (float("nan"),float("nan"))
    p=successes/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half,center+half

def summarize(artifact_dir:Path,output_dir:Path,bins:int=10)->dict:
    root=artifact_dir.resolve(); out=output_dir.resolve()
    cfg=json.loads((root/"config.json").read_text())
    if cfg.get("truth_mode")!="prior_predictive":
        raise ValueError("SBC rank summary requires truth_mode=prior_predictive")
    if bins<2:raise ValueError("at least two rank bins are required")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    summaries=pd.read_csv(root/"replicate_posterior_summaries.csv")
    key=["replicate","gamma","validation_size","arm"]
    if ledger.duplicated(key).any() or summaries.duplicated(key).any():
        raise AssertionError("duplicate replicate/gamma/size/arm key")
    required_ledger={"replicate","gamma","validation_size","arm","status","audit_status"}
    if not required_ledger.issubset(ledger.columns):
        raise AssertionError("replicate ledger is missing lifecycle columns")
    successful=ledger[ledger.status=="SUCCESS"]
    lifecycle_complete=bool((ledger.status=="SUCCESS").all() and (ledger.audit_status=="PASS").all())
    if set(map(tuple,successful[key].itertuples(index=False,name=None))) != set(map(tuple,summaries[key].itertuples(index=False,name=None))):
        raise AssertionError("successful ledger cells and posterior-summary rows differ")
    if not summaries.empty:
        ranks=summaries.posterior_rank.to_numpy(float)
        if not np.isfinite(ranks).all() or np.any((ranks<0)|(ranks>1)):
            raise AssertionError("posterior ranks must be finite and lie in [0,1]")
    summary_rows=[]; histogram_rows=[]
    edges=np.linspace(0.0,1.0,bins+1)
    for (gamma,nval,arm),planned in ledger.groupby(["gamma","validation_size","arm"],sort=True):
        group=summaries[(summaries.gamma==gamma)&(summaries.validation_size==nval)&(summaries.arm==arm)]
        rank=group.posterior_rank.to_numpy(float)
        hist,_=np.histogram(rank,bins=edges)
        for i,count in enumerate(hist):
            histogram_rows.append({"gamma":gamma,"validation_size":nval,"arm":arm,
                "bin":i,"lower":edges[i],"upper":edges[i+1],"count":int(count),
                "planned_replicates":len(planned),"successful_replicates":len(group)})
        expected=len(rank)/bins
        chi2=float(chisquare(hist).statistic) if len(rank)>0 else float("nan")
        chi_p=float(chisquare(hist).pvalue) if len(rank)>0 and expected>=5 else float("nan")
        ks_p=float(kstest(rank,"uniform").pvalue) if len(rank)>0 else float("nan")
        row={"gamma":gamma,"validation_size":nval,"arm":arm,
            "planned_replicates":len(planned),"successful_replicates":len(group),
            "failed_replicates":len(planned)-len(group),"mean_rank":float(rank.mean()) if len(rank) else np.nan,
            "rank_sd":float(rank.std(ddof=1)) if len(rank)>1 else np.nan,
            "rank_ks_statistic":float(kstest(rank,"uniform").statistic) if len(rank) else np.nan,
            "rank_ks_pvalue_descriptive_only":ks_p,"rank_chi_square":chi2,
            "rank_chi_square_pvalue_if_expected_bin_count_ge_5":chi_p}
        for level,(lo,hi) in INTERVALS.items():
            if len(group):
                covered=((group[lo].to_numpy(float)<=group.h0_true.to_numpy(float))&
                         (group.h0_true.to_numpy(float)<=group[hi].to_numpy(float)))
                count=int(covered.sum()); lower,upper=_wilson(count,len(group))
                row[f"coverage_{level}"]=count/len(group)
                row[f"coverage_{level}_wilson95_low"]=lower
                row[f"coverage_{level}_wilson95_high"]=upper
            else:
                row[f"coverage_{level}"]=np.nan
                row[f"coverage_{level}_wilson95_low"]=np.nan
                row[f"coverage_{level}_wilson95_high"]=np.nan
        summary_rows.append(row)
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(out/"prior_predictive_sbc_summary.csv",index=False)
    pd.DataFrame(histogram_rows).to_csv(out/"prior_predictive_rank_histograms.csv",index=False)
    summary_keys=set(map(tuple,summaries[key].itertuples(index=False,name=None)))
    successful_keys=set(map(tuple,successful[key].itertuples(index=False,name=None)))
    lifecycle_complete = bool((ledger.status == "SUCCESS").all()
        and (ledger.audit_status == "PASS").all())
    complete=(lifecycle_complete and len(summary_keys)==len(successful_keys)
        and summary_keys==successful_keys and len(successful)>0)
    result={"status":"PASS" if complete else "INCOMPLETE","truth_mode":"prior_predictive","rank_bins":bins,
        "planned_cells":len(ledger),"successful_cells":len(successful),
        "failed_or_non_successful_cells":int(len(ledger)-len(successful)),
        "summary_rows":len(summary_rows),"rank_histogram_rows":len(histogram_rows),
        "interpretation":"Descriptive prior-predictive SBC summaries; p-values are not confirmatory for small replicate counts. INCOMPLETE means lifecycle or summary coverage is not fully audited."}
    (out/"summary_report.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("artifact_dir",type=Path);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--bins",type=int,default=10)
    a=p.parse_args();result=summarize(a.artifact_dir,a.output_dir,a.bins)
    print(json.dumps(result,indent=2,sort_keys=True))
    if result["status"]!="PASS":
        raise SystemExit(2)
if __name__=="__main__":main()
