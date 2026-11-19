"""Paired narrow-vs-wide H0 prior sensitivity audit for identical replicate seeds.

Checks identical synthetic inputs, validation samples, weight diagnostics, and
selected likelihood curves on the common [50,90] support before attributing
posterior-summary changes to prior support. This is a synthetic sensitivity
analysis, not observational inference.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t

ARMS=["A_catalogue_only","B_parent_oracle","C_known_probability_ipw",
      "D_map_plugin_ipw","E_uncertainty_integrated_ipw"]
RAW_FILES=["truth/parent_truth.csv","truth/gw_event_truth.csv",
           "truth/gw_generation_trace_truth.csv","truth/validation_frame_truth.csv",
           "observable/catalogue_observed.csv","observable/gw_events_observed.csv"]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def interval(x):
    x=np.asarray(x,float); n=len(x); m=float(x.mean())
    if n<2: return m,float("nan"),float("nan"),float("nan")
    se=float(x.std(ddof=1)/np.sqrt(n)); q=float(t.ppf(.975,n-1))
    return m,se,m-q*se,m+q*se

def audit_pair(narrow:Path,wide:Path,out:Path):
    narrow=narrow.resolve();wide=wide.resolve()
    a=json.loads((narrow/"config.json").read_text());b=json.loads((wide/"config.json").read_text())
    for key in ["master_seed","replicates","gammas","validation_sizes","n_parent","n_detected",
                "validation_frame","curvature_beta","h0_true","scott_multiplier",
                "quadrature_order","sigma_mpc","detection_threshold_mpc"]:
        if a[key]!=b[key]: raise AssertionError(f"paired configuration mismatch: {key}")
    if list(map(float,a.get("h0_prior_support",[50,90]))) != [50.0,90.0]:
        raise AssertionError("first artifact is not the [50,90] prior")
    if list(map(float,b.get("h0_prior_support",[50,90]))) != [30.0,110.0]:
        raise AssertionError("second artifact is not the [30,110] prior")
    if int(a["h0_points"])!=281 or int(b["h0_points"])!=561:
        raise AssertionError("unexpected grid sizes for common-support pairing")
    reps=int(a["replicates"]);gammas=list(map(float,a["gammas"]));sizes=list(map(int,a["validation_sizes"]))
    summaries_a=pd.read_csv(narrow/"replicate_posterior_summaries.csv")
    summaries_b=pd.read_csv(wide/"replicate_posterior_summaries.csv")
    weight_a=pd.read_csv(narrow/"replicate_weight_diagnostics.csv")
    weight_b=pd.read_csv(wide/"replicate_weight_diagnostics.csv")
    rows=[];max_ll=0.0;raw_checks=0;weight_checks=0
    for rep in range(reps):
        for gamma in gammas:
            suffix=f"g{gamma:+g}".replace("+","p").replace("-","m")
            ca=narrow/"replicates"/f"rep{rep:03d}_{suffix}"
            cb=wide/"replicates"/f"rep{rep:03d}_{suffix}"
            for rel in RAW_FILES:
                if sha(ca/rel)!=sha(cb/rel): raise AssertionError(f"raw input mismatch {rel} rep={rep} gamma={gamma}")
                raw_checks+=1
            for nval in sizes:
                rel=f"observable/validation_n{nval}.csv"
                if sha(ca/rel)!=sha(cb/rel): raise AssertionError(f"validation sample mismatch {rel} rep={rep} gamma={gamma}")
                raw_checks+=1
            for rel in ["results/weight_diagnostics.csv","results/realized_bandwidths.csv"]:
                if sha(ca/rel)!=sha(cb/rel): raise AssertionError(f"weight/bandwidth mismatch {rel} rep={rep} gamma={gamma}")
                weight_checks+=1
            xa=np.load(ca/"results"/"log_likelihoods.npz")
            xb=np.load(cb/"results"/"log_likelihoods.npz")
            ga=xa["H0_grid"];gb=xb["H0_grid"]
            common_idx=np.arange(140,421,1)
            if not np.allclose(gb[common_idx],ga,rtol=0,atol=1e-12):
                raise AssertionError("wide grid does not contain the original narrow grid")
            for nval in sizes:
                for arm in ARMS:
                    key=f"{arm}_n{nval}"
                    diff=float(np.max(np.abs(xa[key]-xb[key][common_idx])))
                    max_ll=max(max_ll,diff)
                    if diff>1e-9: raise AssertionError(f"common-support likelihood mismatch {key}: {diff}")
    keys=["posterior_median","q025","q05","q16","q84","q95","q975","median_error"]
    join=["replicate","gamma","validation_size","arm"]
    paired=summaries_a.merge(summaries_b,on=join,suffixes=("_narrow","_wide"),validate="one_to_one")
    for _,r in paired.iterrows():
        row={k:r[k] for k in join}
        for key in keys: row[key+"_difference_wide_minus_narrow"]=float(r[key+"_wide"]-r[key+"_narrow"])
        row["absolute_median_error_difference_wide_minus_narrow"]=abs(float(r.median_error_wide))-abs(float(r.median_error_narrow))
        truth=float(r.h0_true_narrow)
        row["coverage68_narrow"]=int(r.q16_narrow<=truth<=r.q84_narrow)
        row["coverage68_wide"]=int(r.q16_wide<=truth<=r.q84_wide)
        row["width68_difference_wide_minus_narrow"]=float((r.q84_wide-r.q16_wide)-(r.q84_narrow-r.q16_narrow))
        rows.append(row)
    diff=pd.DataFrame(rows)
    diff.to_csv(out.with_suffix(".csv"),index=False)
    summaries=[]
    for (gamma,nval,arm),g in diff.groupby(["gamma","validation_size","arm"],sort=True):
        for metric in ["median_error_difference_wide_minus_narrow",
                       "absolute_median_error_difference_wide_minus_narrow",
                       "width68_difference_wide_minus_narrow"]:
            mean,se,lo,hi=interval(g[metric].to_numpy(float))
            summaries.append({"gamma":float(gamma),"validation_size":int(nval),"arm":arm,
                "metric":metric,"n":len(g),"mean_difference":mean,"paired_se":se,"t95_low":lo,"t95_high":hi})
        cov=(g.coverage68_wide-g.coverage68_narrow).to_numpy(float)
        mean,se,lo,hi=interval(cov)
        summaries.append({"gamma":float(gamma),"validation_size":int(nval),"arm":arm,
            "metric":"coverage68_difference_wide_minus_narrow","n":len(g),
            "mean_difference":mean,"paired_se":se,"t95_low":lo,"t95_high":hi})
    pd.DataFrame(summaries).to_csv(out.with_name(out.stem+"_summary.csv"),index=False)
    result={"status":"PASS","paired_replicate_gamma_children":reps*len(gammas),
        "paired_posterior_rows":len(diff),"raw_input_file_checks":raw_checks,
        "weight_and_bandwidth_file_checks":weight_checks,
        "maximum_common_support_log_likelihood_difference":max_ll,
        "interpretation":"Paired prior sensitivity under identical synthetic data and likelihoods; no observational inference."}
    out.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("narrow_artifact",type=Path);p.add_argument("wide_artifact",type=Path)
    p.add_argument("--output-prefix",type=Path,default=Path("paired_prior_sensitivity"));x=p.parse_args()
    print(json.dumps(audit_pair(x.narrow_artifact,x.wide_artifact,x.output_prefix),indent=2))
if __name__=="__main__":main()
