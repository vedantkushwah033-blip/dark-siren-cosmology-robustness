"""Independent raw-seed reconstruction audit for the 20,000-draw full-parent baseline."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import logsumexp, ndtr, ndtri
from scipy.stats import norm, t

C_KM_S=299792.458

def _verify_manifest(root:Path):
    manifest=json.loads((root/"manifest.json").read_text())
    bad=[name for name,digest in manifest.items()
         if not (root/name).is_file() or hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest]
    if bad: raise AssertionError(f"manifest mismatch: {bad[:10]}")
    return len(manifest)

def _draw_detected_independent(rng,z,host_weight,h0,sigma,threshold,n_events):
    """Independent transcription of the conditional finite-parent detected-event sampler."""
    mu=C_KM_S*np.asarray(z,float)/float(h0)
    host=np.asarray(host_weight,float)
    pdet=ndtr((mu-threshold)/sigma)
    probs=host*pdet
    probs=probs/probs.sum()
    host_ids=rng.choice(len(mu),size=n_events,replace=True,p=probs)
    selected_mu=mu[host_ids]
    cdf_lower=ndtr((threshold-selected_mu)/sigma)
    uniforms=cdf_lower+rng.random(n_events)*(1.0-cdf_lower)
    distance=selected_mu+sigma*ndtri(uniforms)
    if np.any(distance<threshold): raise AssertionError("independent sampler generated an undetected event")
    return distance

def _posterior_quantiles_independent(z,host_weight,distances,h,sigma,threshold):
    z=np.asarray(z,float); w=np.asarray(host_weight,float); d=np.asarray(distances,float)
    logw=np.log(w); logw-=logsumexp(logw)
    means=C_KM_S*z[None,:]/h[:,None]
    log_eff=logsumexp(norm.logcdf((means-threshold)/sigma)+logw[None,:],axis=1)
    log_num=np.zeros(len(h))
    for value in d:
        log_num+=logsumexp(norm.logpdf(value,loc=means,scale=sigma)+logw[None,:],axis=1)
    logpost=log_num-len(d)*log_eff
    density=np.exp(logpost-np.max(logpost))
    density/=trapezoid(density,h)
    cdf=np.r_[0.0,cumulative_trapezoid(density,h)]
    cdf/=cdf[-1]
    probs=np.array([.025,.05,.16,.50,.84,.95,.975])
    quantiles=np.interp(probs,cdf,h)
    return dict(zip(("q025","q05","q16","median","q84","q95","q975"),map(float,quantiles)))

def _close(a,b,label,atol=2e-8):
    if not np.isclose(float(a),float(b),rtol=0,atol=atol,equal_nan=True):
        raise AssertionError(f"{label}: archived={a}, rebuilt={b}")

def _wilson(k,n,z=1.959963984540054):
    p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half,center+half

def audit(source_root:Path,baseline_root:Path,output:Path|None=None,
          expected_parents=100,draws_per_parent=200):
    source_root=source_root.resolve(); baseline_root=baseline_root.resolve()
    source_manifest_count=_verify_manifest(source_root)
    baseline_manifest_count=_verify_manifest(baseline_root)
    ledger=pd.read_csv(source_root/"replicate_ledger.csv")
    baseline=pd.read_csv(baseline_root/"replicates.csv")
    if len(ledger)!=600 or ledger.cell_id.duplicated().any():
        raise AssertionError("source 600-cell planned ledger is incomplete or duplicated")
    if len(baseline)!=expected_parents*draws_per_parent:
        raise AssertionError("20,000-draw baseline denominator mismatch")
    if baseline.duplicated(["parent_replicate","event_replicate"]).any():
        raise AssertionError("duplicate parent/event replicate keys")
    if baseline.groupby("parent_replicate").size().to_dict()!={i:draws_per_parent for i in range(expected_parents)}:
        raise AssertionError("per-parent event-draw counts are incomplete")
    design=json.loads((baseline_root/"summary.json").read_text())
    if int(design["parent_populations"])!=expected_parents or int(design["replicates_per_parent"])!=draws_per_parent:
        raise AssertionError("baseline design metadata differs from planned audit")
    h=np.linspace(20.0,300.0,561)
    rebuilt=[]
    max_qdiff=0.0
    for rep in range(expected_parents):
        cells=ledger[(ledger.replicate==rep)&(ledger.gamma==-1.0)&(ledger.n_validation==100)]
        if len(cells)!=1: raise AssertionError(f"source representative cell missing for parent {rep}")
        cell=cells.iloc[0]
        child=source_root/"replicates"/str(cell.cell_id)
        cfg=json.loads((child/"config.json").read_text())
        parent=pd.read_csv(child/"truth/analysis_parent_truth.csv")
        z=parent.redshift.to_numpy(float)
        zmin,zmax=map(float,cfg["analysis_generator_config"]["redshift_domain"])
        sigma=float(cfg["analysis_generator_config"]["sigma_mpc"])
        threshold=float(cfg["analysis_generator_config"]["detection_threshold_mpc"])
        truth=float(cell.h0_true)
        n_events=int(cfg["analysis_generator_config"]["n_detected_events"])
        if n_events!=4 or truth!=70.0: raise AssertionError("unexpected fixed-truth/event-count source design")
        base=1.0+1.5*(z-zmin)/(zmax-zmin)
        rng=np.random.default_rng(990000+rep)
        block=baseline[baseline.parent_replicate==rep].sort_values("event_replicate")
        if not np.array_equal(block.event_replicate.to_numpy(int),np.arange(draws_per_parent)):
            raise AssertionError(f"event replicate indices incomplete for parent {rep}")
        for row in block.itertuples(index=False):
            distances=_draw_detected_independent(rng,z,base,truth,sigma,threshold,n_events)
            q=_posterior_quantiles_independent(z,base,distances,h,sigma,threshold)
            expected={**q,"median_bias":q["median"]-truth,
              "covered68":int(q["q16"]<=truth<=q["q84"]),
              "covered90":int(q["q05"]<=truth<=q["q95"]),
              "covered95":int(q["q025"]<=truth<=q["q975"]),
              "width68":q["q84"]-q["q16"]}
            for key,value in expected.items():
                diff=abs(float(getattr(row,key))-float(value))
                max_qdiff=max(max_qdiff,diff)
                _close(getattr(row,key),value,f"{key} parent={rep} event={row.event_replicate}")
            rebuilt.append({"parent_replicate":rep,"event_replicate":int(row.event_replicate),
                            "median_bias":expected["median_bias"],
                            "covered68":expected["covered68"],"covered90":expected["covered90"],
                            "covered95":expected["covered95"],"width68":expected["width68"]})
    rebuilt_frame=pd.DataFrame(rebuilt)
    archived_summary=json.loads((baseline_root/"summary.json").read_text())
    bias=rebuilt_frame.median_bias.to_numpy(float)
    mean=float(bias.mean()); sd=float(bias.std(ddof=1)); mcse=sd/math.sqrt(len(bias))
    crit=float(t.ppf(.975,len(bias)-1))
    for key,value in {
      "mean_median_bias":mean,"median_bias_sd":sd,"median_bias_mcse":mcse,
      "median_bias_t95_low":mean-crit*mcse,"median_bias_t95_high":mean+crit*mcse,
      "mean_width68":float(rebuilt_frame.width68.mean())
    }.items(): _close(archived_summary[key],value,f"pooled summary {key}")
    for level in (68,90,95):
        k=int(rebuilt_frame[f"covered{level}"].sum()); p=k/len(rebuilt_frame)
        lo,hi=_wilson(k,len(rebuilt_frame))
        for key,value in ((f"covered{level}",p),(f"covered{level}_wilson95_low",lo),
                          (f"covered{level}_wilson95_high",hi)):
            _close(archived_summary[key],value,f"pooled summary {key}")
    parent=rebuilt_frame.groupby("parent_replicate").agg(
      n=("median_bias","size"),mean_median_bias=("median_bias","mean"),
      coverage68=("covered68","mean"),coverage90=("covered90","mean"),
      coverage95=("covered95","mean")).reset_index()
    archived_parent=pd.read_csv(baseline_root/"parent_summary.csv").sort_values("parent_replicate").reset_index(drop=True)
    pd.testing.assert_frame_equal(parent,archived_parent,check_exact=False,rtol=0,atol=2e-8)
    cluster=json.loads((baseline_root/"cluster_summary.json").read_text())
    for key in ("mean_median_bias","coverage68","coverage90","coverage95"):
        vals=parent[key].to_numpy(float); n=len(vals); m=float(vals.mean())
        se=float(vals.std(ddof=1)/math.sqrt(n)); c=float(t.ppf(.975,n-1))
        for metric,value in (("mean",m),("cluster_mcse",se),("t95_low",m-c*se),("t95_high",m+c*se)):
            _close(cluster[key][metric],value,f"parent cluster {key} {metric}")
    report={"status":"PASS","source_manifest_entries_verified":source_manifest_count,
      "baseline_manifest_entries_verified":baseline_manifest_count,
      "planned_source_cells":len(ledger),"parent_populations_regenerated":expected_parents,
      "detected_event_realizations_regenerated":expected_parents*draws_per_parent,
      "posterior_curves_reconstructed":expected_parents*draws_per_parent,
      "parent_cluster_summaries_recomputed":expected_parents,
      "max_abs_quantile_or_metric_difference":max_qdiff,
      "scope":"Independent raw-seed reconstruction of the repeated four-event full-parent fixed-truth baseline; not joint SBC."}
    if output:
        output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source_artifact",type=Path); p.add_argument("baseline_artifact",type=Path)
    p.add_argument("--expected-parents",type=int,default=100)
    p.add_argument("--draws-per-parent",type=int,default=200)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(audit(a.source_artifact,a.baseline_artifact,a.output,a.expected_parents,a.draws_per_parent),indent=2))
if __name__=="__main__": main()
