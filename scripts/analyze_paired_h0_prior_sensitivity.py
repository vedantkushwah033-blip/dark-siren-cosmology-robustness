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
def verify_artifact(root):
    manifest=json.loads((root/"manifest.json").read_text())
    for rel,meta in manifest["files"].items():
        path=root/rel
        if not path.is_file() or path.stat().st_size!=int(meta["bytes"]) or sha(path)!=meta["sha256"]:
            raise AssertionError(f"artifact manifest mismatch: {root.name}/{rel}")
    if len(manifest["files"])!=int(manifest["file_count"]):
        raise AssertionError("artifact manifest count mismatch")
    cfg=json.loads((root/"config.json").read_text())
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    planned=pd.read_csv(root/"planned_cells.csv")
    summaries=pd.read_csv(root/"replicate_posterior_summaries.csv")
    expected=int(cfg["replicates"])*len(cfg["gammas"])*len(cfg["validation_sizes"])*len(ARMS)
    if len(ledger)!=expected or len(planned)!=expected or len(summaries)!=expected:
        raise AssertionError(f"artifact lifecycle denominator mismatch: {root}")
    if not (ledger.status=="SUCCESS").all() or not (ledger.audit_status=="PASS").all():
        raise AssertionError(f"artifact lifecycle is incomplete: {root}")
    for rep in range(int(cfg["replicates"])):
        for gamma in cfg["gammas"]:
            suffix=f"g{float(gamma):+g}".replace("+","p").replace("-","m")
            child=root/"replicates"/f"rep{rep:03d}_{suffix}"
            report=json.loads((child/"audit_report.json").read_text())
            if report.get("status")!="PASS":
                raise AssertionError(f"child audit not PASS: {child}")
    return len(manifest["files"]),len(ledger)
def validate_paired_summary_rows(a, b, join, expected):
    """Require a complete, unique, identical replicate/arm key set before pairing."""
    for label, frame in [("narrow", a), ("wide", b)]:
        if frame.duplicated(join).any():
            raise AssertionError(f"{label} posterior summaries contain duplicate pairing keys")
        if len(frame) != expected:
            raise AssertionError(f"{label} posterior summary row count {len(frame)} != {expected}")
    ka=set(map(tuple,a[join].itertuples(index=False,name=None)))
    kb=set(map(tuple,b[join].itertuples(index=False,name=None)))
    if ka != kb:
        raise AssertionError(
            f"posterior summary pairing keys differ: narrow_only={len(ka-kb)}, wide_only={len(kb-ka)}"
        )
    if "h0_true" in a.columns and "h0_true" in b.columns:
        aa=a[join+["h0_true"]].sort_values(join).reset_index(drop=True)
        bb=b[join+["h0_true"]].sort_values(join).reset_index(drop=True)
        if not np.array_equal(aa.h0_true.to_numpy(),bb.h0_true.to_numpy()):
            raise AssertionError("paired posterior summaries disagree on H0 truth")
    # Validate each source summary internally before computing paired contrasts.
    required = {"posterior_median", "q025", "q05", "q16", "q84", "q95", "q975"}
    for label, frame in [("narrow", a), ("wide", b)]:
        present = required.intersection(frame.columns)
        if present != required:
            raise AssertionError(f"{label} posterior summaries are missing quantile columns: {sorted(required-present)}")
        values = frame[sorted(required)].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise AssertionError(f"{label} posterior summaries contain non-finite quantiles")
        ordered = frame[["q025", "q05", "q16", "posterior_median", "q84", "q95", "q975"]].to_numpy(dtype=float)
        if np.any(np.diff(ordered, axis=1) < -1e-10):
            raise AssertionError(f"{label} posterior quantiles are not monotonically ordered")
        if "h0_true" in frame.columns and "median_error" in frame.columns:
            expected_error = frame["posterior_median"].to_numpy(dtype=float) - frame["h0_true"].to_numpy(dtype=float)
            observed_error = frame["median_error"].to_numpy(dtype=float)
            if not np.allclose(observed_error, expected_error, rtol=0, atol=1e-8):
                raise AssertionError(f"{label} median_error does not equal posterior_median minus h0_true")
    return a.merge(b,on=join,suffixes=("_narrow","_wide"),validate="one_to_one")


def interval(x):
    x=np.asarray(x,float); n=len(x); m=float(x.mean())
    if n<2: return m,float("nan"),float("nan"),float("nan")
    se=float(x.std(ddof=1)/np.sqrt(n)); q=float(t.ppf(.975,n-1))
    return m,se,m-q*se,m+q*se

def audit_pair(narrow:Path,wide:Path,out:Path):
    narrow=narrow.resolve();wide=wide.resolve()
    manifest_a,ledger_a=verify_artifact(narrow)
    manifest_b,ledger_b=verify_artifact(wide)
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
    expected_rows=reps*len(gammas)*len(sizes)*len(ARMS)
    paired=validate_paired_summary_rows(summaries_a,summaries_b,join,expected_rows)
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
    result={"status":"PASS","artifact_manifest_entries_verified":[manifest_a,manifest_b],
        "lifecycle_cells_verified":[ledger_a,ledger_b],
        "paired_replicate_gamma_children":reps*len(gammas),
        "paired_posterior_rows":len(diff),"raw_input_file_checks":raw_checks,
        "weight_and_bandwidth_file_checks":weight_checks,
        "maximum_common_support_log_likelihood_difference":max_ll,
        "interpretation":"Paired prior sensitivity under identical synthetic data and likelihoods; no observational inference."}
    out.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("narrow_artifact",type=Path);p.add_argument("wide_artifact",type=Path)
    p.add_argument("--output-prefix",type=Path,default=Path("paired_prior_sensitivity"));x=p.parse_args()
    result=audit_pair(x.narrow_artifact,x.wide_artifact,x.output_prefix)
    print(json.dumps(result,indent=2))
    summary=pd.read_csv(x.output_prefix.with_name(x.output_prefix.stem+"_summary.csv"))
    print("\nPAIRED PRIOR-SUPPORT READOUT — all validation sizes")
    selected=summary[summary.metric.isin([
        "median_error_difference_wide_minus_narrow",
        "absolute_median_error_difference_wide_minus_narrow",
        "width68_difference_wide_minus_narrow",
        "coverage68_difference_wide_minus_narrow",
    ])]
    print(selected.to_string(index=False,float_format=lambda value:f"{value:.6f}"))
if __name__=="__main__":main()
