"""Repeated-sampling full-parent oracle baseline for the selected-H0 toy likelihood.

Uses archived parent populations but generates fresh detected-event samples
from the exact fixed-count selection process. This isolates the event-sample
baseline before catalogue incompleteness or estimated inclusion probabilities.
Raw parent inputs remain usable when a fitted-inclusion cell is NON_EVALUABLE.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import logsumexp, log_ndtr, ndtr, ndtri
from scipy.stats import t

C_KM_S = 299792.458


def draw_detected_distances(rng, redshift, host_probability, h0_true,
                            sigma_mpc, threshold_mpc, n_events):
    """Draw iid events exactly conditional on the one-sided distance cut."""
    mu = C_KM_S * np.asarray(redshift, float) / float(h0_true)
    p_detect = ndtr((mu - threshold_mpc) / sigma_mpc)
    detected_host_probability = np.asarray(host_probability, float) * p_detect
    detected_host_probability /= detected_host_probability.sum()
    ids = rng.choice(len(mu), size=n_events, replace=True, p=detected_host_probability)
    selected_mu = mu[ids]
    lower_cdf = ndtr((threshold_mpc - selected_mu) / sigma_mpc)
    u = lower_cdf + rng.random(n_events) * (1.0 - lower_cdf)
    distances = selected_mu + sigma_mpc * ndtri(u)
    if np.any(distances < threshold_mpc):
        raise ArithmeticError("conditional sampler produced an undetected event")
    return distances


def prepare_full_parent_likelihood(redshift, host_weights, h0_grid,
                                    sigma_mpc, threshold_mpc):
    """Precompute the exact finite-parent selected likelihood for one parent."""
    z = np.asarray(redshift, float)
    w = np.asarray(host_weights, float)
    h = np.asarray(h0_grid, float)
    if z.ndim != 1 or w.shape != z.shape or z.size < 2 or np.any(w <= 0):
        raise ValueError("invalid parent redshifts or positive host weights")
    logw = np.log(w) - logsumexp(np.log(w))
    means = C_KM_S * z[None, :] / h[:, None]
    log_beta = logsumexp(
        log_ndtr((means-threshold_mpc)/sigma_mpc) + logw[None, :], axis=1
    )
    return means, logw, log_beta, h


def selected_posterior_quantiles(distances, prepared, sigma_mpc, threshold_mpc):
    """Independent full-parent selected-event likelihood and equal-tailed quantiles."""
    means, logw, log_beta, h = prepared
    d = np.asarray(distances, float)
    if d.ndim != 1 or not d.size or np.any(d < threshold_mpc):
        raise ValueError("distances must all satisfy the detection threshold")
    log_numerator = np.zeros(h.size)
    log_norm = math.log(sigma_mpc*math.sqrt(2.0*math.pi))
    for distance in d:
        logpdf = -0.5*((distance-means)/sigma_mpc)**2-log_norm
        log_numerator += logsumexp(logpdf+logw[None,:],axis=1)
    log_likelihood = log_numerator-d.size*log_beta
    density=np.exp(log_likelihood-np.max(log_likelihood))
    density/=trapezoid(density,h)
    cdf=np.r_[0.0,np.cumsum(0.5*(density[1:]+density[:-1])*np.diff(h))]
    cdf/=cdf[-1]
    qs=np.interp([.025,.05,.16,.50,.84,.95,.975],cdf,h)
    return dict(zip(("q025","q05","q16","median","q84","q95","q975"),map(float,qs)))


def wilson(k, n, z=1.959963984540054):
    p=k/n
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half,center+half


def run(artifact_root, output_dir, replicates_per_parent=200, seed=990000, h0_points=281):
    root=Path(artifact_root).resolve()
    manifest=json.loads((root/"manifest.json").read_text())
    bad=[name for name,digest in manifest.items()
         if not (root/name).is_file()
         or hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest]
    if bad:
        raise AssertionError(f"input aggregate manifest mismatch: {bad[:10]}")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    allowed_statuses={"SUCCESS","NON_EVALUABLE"}
    if ledger.empty or not set(ledger.status).issubset(allowed_statuses):
        raise AssertionError(
            "input pilot must contain only SUCCESS/NON_EVALUABLE cells; "
            f"observed statuses={sorted(set(ledger.status))}"
        )
    if ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate cell IDs in input ledger")
    output=Path(output_dir).resolve()
    output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for replicate in sorted(ledger.replicate.unique()):
        selected=ledger[(ledger.replicate==replicate)&
                        (ledger.gamma==-1.0)&(ledger.n_validation==100)]
        if len(selected)!=1:
            raise AssertionError(f"expected one representative cell for replicate {replicate}")
        cell=selected.iloc[0]
        directory=root/"replicates"/cell.cell_id
        cfg=json.loads((directory/"config.json").read_text())
        parent=pd.read_csv(directory/"truth/analysis_parent_truth.csv")
        z=parent.redshift.to_numpy(float)
        z_min,z_max=cfg["analysis_generator_config"]["redshift_domain"]
        sigma=float(cfg["analysis_generator_config"]["sigma_mpc"])
        threshold=float(cfg["analysis_generator_config"]["detection_threshold_mpc"])
        h0_true=float(cell.h0_true)
        n_events=int(cfg["analysis_generator_config"]["n_detected_events"])
        base=1.0+1.5*(z-z_min)/(z_max-z_min)
        host_probability=base/base.sum()
        h0_min,h0_max=cfg["inference"]["uniform_h0_prior"]
        h0_grid=np.linspace(h0_min,h0_max,h0_points)
        prepared=prepare_full_parent_likelihood(z,base,h0_grid,sigma,threshold)
        rng=np.random.default_rng(seed+int(replicate))
        for sim in range(replicates_per_parent):
            distances=draw_detected_distances(
                rng,z,host_probability,h0_true,sigma,threshold,n_events
            )
            q=selected_posterior_quantiles(distances,prepared,sigma,threshold)
            rows.append({
                "parent_replicate":int(replicate),"event_replicate":int(sim),
                "n_events":n_events,"h0_true":h0_true,
                "median":q["median"],"median_bias":q["median"]-h0_true,
                "q025":q["q025"],"q05":q["q05"],"q16":q["q16"],
                "q84":q["q84"],"q95":q["q95"],"q975":q["q975"],
                "covered68":int(q["q16"]<=h0_true<=q["q84"]),
                "covered90":int(q["q05"]<=h0_true<=q["q95"]),
                "covered95":int(q["q025"]<=h0_true<=q["q975"]),
                "width68":q["q84"]-q["q16"],
            })
    frame=pd.DataFrame(rows)
    frame.to_csv(output/"replicates.csv",index=False)
    n=len(frame); bias=frame.median_bias.to_numpy(float)
    mean=float(bias.mean()); sd=float(bias.std(ddof=1)); se=sd/math.sqrt(n)
    critical=float(t.ppf(.975,n-1))
    summary={"status":"DIAGNOSTIC_COMPLETE",
        "parent_populations":int(frame.parent_replicate.nunique()),
        "replicates_per_parent":int(replicates_per_parent),"total_replicates":n,
        "h0_true":float(frame.h0_true.iloc[0]),"n_events":int(frame.n_events.iloc[0]),
        "h0_grid_points":int(h0_points),"mean_median_bias":mean,
        "median_bias_sd":sd,"median_bias_mcse":se,
        "median_bias_t95_low":mean-critical*se,"median_bias_t95_high":mean+critical*se,
        "mean_width68":float(frame.width68.mean())}
    for col in ("covered68","covered90","covered95"):
        k=int(frame[col].sum()); p=k/n; lo,hi=wilson(k,n)
        summary[col]=p; summary[col+"_mcse"]=math.sqrt(p*(1-p)/n)
        summary[col+"_wilson95_low"]=lo; summary[col+"_wilson95_high"]=hi
    frame.groupby("parent_replicate").agg(
        n=("median_bias","size"),mean_median_bias=("median_bias","mean"),
        coverage68=("covered68","mean"),coverage90=("covered90","mean"),
        coverage95=("covered95","mean")).reset_index().to_csv(
            output/"parent_summary.csv",index=False)
    (output/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
    output_manifest={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.rglob("*")) if p.is_file() and p.name!="manifest.json"}
    (output/"manifest.json").write_text(json.dumps(output_manifest,indent=2,sort_keys=True))
    print(json.dumps(summary,indent=2))
    return frame,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_root",type=Path)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--replicates-per-parent",type=int,default=200)
    parser.add_argument("--seed",type=int,default=990000)
    parser.add_argument("--h0-points",type=int,default=281)
    args=parser.parse_args()
    run(args.artifact_root,args.output_dir,args.replicates_per_parent,args.seed,args.h0_points)


if __name__=="__main__":
    main()
