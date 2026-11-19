"""Prior-predictive simulation-based calibration for the full-parent H0 likelihood.

Every replicate draws H0 from the same uniform prior used for inference, then
generates a fresh finite parent population and 30 detected events. The
full-parent selected-event posterior is evaluated without catalogue thinning.
This isolates the baseline likelihood before inclusion correction is added.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import log_ndtr, logsumexp, ndtr
from scipy.stats import kstest, chisquare

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299792.458

def seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate), int(stream)]).generate_state(1, dtype=np.uint32)[0])

def posterior_on_grid(z, host_weights, distances, h0_grid, sigma, threshold):
    z=np.asarray(z,float); w=np.asarray(host_weights,float); d=np.asarray(distances,float); h=np.asarray(h0_grid,float)
    if z.ndim!=1 or w.shape!=z.shape or np.any(w<=0) or not np.isfinite(z).all() or not np.isfinite(w).all():
        raise ValueError("invalid full-parent arrays")
    if d.ndim!=1 or len(d)<1 or np.any(d<threshold) or not np.isfinite(d).all():
        raise ValueError("invalid detected-event distances")
    logw=np.log(w)-logsumexp(np.log(w))
    mu=C_KM_S*z[None,:]/h[:,None]
    log_beta=logsumexp(log_ndtr((mu-threshold)/sigma)+logw[None,:],axis=1)
    log_num=np.zeros(len(h))
    log_norm=np.log(sigma*np.sqrt(2*np.pi))
    for distance in d:
        logpdf=-0.5*((distance-mu)/sigma)**2-log_norm
        log_num+=logsumexp(logpdf+logw[None,:],axis=1)
    logpost=log_num-len(d)*log_beta
    density=np.exp(logpost-np.max(logpost))
    area=float(trapezoid(density,h))
    if not np.isfinite(area) or area<=0: raise ArithmeticError("posterior normalization failed")
    density/=area
    cdf=np.r_[0.0,cumulative_trapezoid(density,h)]
    cdf/=cdf[-1]
    return density,cdf

def run(output_dir: Path, replicates=200, master_seed=2026111907, n_parent=1200,
        n_detected=30, h0_min=20.0, h0_max=300.0, h0_points=561):
    if replicates<20 or n_parent<100 or n_detected<1 or h0_points<101 or not h0_min<h0_max:
        raise ValueError("invalid SBC design")
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise ValueError("output directory must be empty")
    h=np.linspace(h0_min,h0_max,h0_points)
    rows=[]
    for rep in range(replicates):
        truth_rng=np.random.default_rng(seed(master_seed,rep,0))
        truth=float(truth_rng.uniform(h0_min,h0_max))
        sim_seed=seed(master_seed,rep,1)
        sim=generate_synthetic_population(seed=sim_seed,n_parent=n_parent,n_detected=n_detected,h0_true=truth)
        parent=sim["parent_truth"]; events=sim["gw_events_observed"]; cfg=sim["config"]
        z=parent.redshift.to_numpy(float)
        zmin,zmax=map(float,cfg["redshift_domain"])
        host_w=1.0+1.5*(z-zmin)/(zmax-zmin)
        distances=events.distance_mpc.to_numpy(float)
        density,cdf=posterior_on_grid(z,host_w,distances,h,float(cfg["sigma_mpc"]),float(cfg["detection_threshold_mpc"]))
        rank=float(np.interp(truth,h,cdf))
        q=np.interp([.025,.05,.10,.16,.25,.50,.75,.84,.90,.95,.975],cdf,h)
        rows.append({"replicate":rep,"truth_seed":seed(master_seed,rep,0),"simulation_seed":sim_seed,
          "h0_true":truth,"n_parent":len(z),"n_detected":len(distances),"n_catalogue":len(sim["catalogue_observed"]),
          "posterior_rank":rank,"h0_median":float(q[5]),"median_bias":float(q[5]-truth),
          "q025":float(q[0]),"q05":float(q[1]),"q10":float(q[2]),"q16":float(q[3]),
          "q25":float(q[4]),"q75":float(q[6]),"q84":float(q[7]),"q90":float(q[8]),
          "q95":float(q[9]),"q975":float(q[10]),
          "covered50":int(q[4]<=truth<=q[6]),"covered68":int(q[3]<=truth<=q[7]),
          "covered80":int(q[2]<=truth<=q[8]),"covered90":int(q[1]<=truth<=q[9]),
          "covered95":int(q[0]<=truth<=q[10]),"mean_posterior_rank":rank})
    frame=pd.DataFrame(rows); frame.to_csv(out/"sbc_replicates.csv",index=False)
    ranks=frame.posterior_rank.to_numpy(float)
    bins=np.linspace(0,1,11); counts,_=np.histogram(ranks,bins=bins)
    chi=chisquare(counts,np.full(10,replicates/10))
    ks=kstest(ranks,"uniform")
    summary={"status":"COMPLETE","replicates":replicates,"n_parent":n_parent,"n_detected":n_detected,
      "h0_prior":[h0_min,h0_max],"h0_grid_points":h0_points,"master_seed":master_seed,
      "mean_rank":float(ranks.mean()),"rank_sd":float(ranks.std(ddof=1)),
      "rank_histogram_10_bins":counts.tolist(),"rank_histogram_chi_square":float(chi.statistic),
      "rank_histogram_chi_square_pvalue":float(chi.pvalue),"rank_ks_statistic":float(ks.statistic),
      "rank_ks_pvalue":float(ks.pvalue),"mean_median_bias":float(frame.median_bias.mean()),
      "median_bias_mcse":float(frame.median_bias.std(ddof=1)/np.sqrt(replicates)),
      "rmse":float(np.sqrt(np.mean(frame.median_bias.to_numpy()**2)))}
    for level in (50,68,80,90,95):
        k=int(frame[f"covered{level}"].sum()); p=k/replicates
        summary[f"coverage{level}"]=p
        summary[f"coverage{level}_mcse"]=float(np.sqrt(p*(1-p)/replicates))
    (out/"sbc_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    design={"status":"COMPLETE","replicates":replicates,"master_seed":master_seed,
      "truth_seed_stream":0,"population_seed_stream":1,"n_parent":n_parent,"n_detected":n_detected,
      "h0_prior":[h0_min,h0_max],"h0_grid_points":h0_points,
      "generation":"generate_synthetic_population; full parent and fixed-count detected events",
      "inference":"full finite-parent discrete-host likelihood with one-sided detection-efficiency normalization and matching uniform H0 prior",
      "catalogue":"catalogue outcomes generated but excluded from full-parent baseline inference"}
    (out/"design.json").write_text(json.dumps(design,indent=2,sort_keys=True)+"\n")
    manifest={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
      for p in sorted(out.rglob("*")) if p.is_file() and p.name!="manifest.json"}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,indent=2))
    return frame,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--replicates",type=int,default=200)
    p.add_argument("--master-seed",type=int,default=2026111907)
    p.add_argument("--n-parent",type=int,default=1200)
    p.add_argument("--n-detected",type=int,default=30)
    p.add_argument("--h0-min",type=float,default=20.0)
    p.add_argument("--h0-max",type=float,default=300.0)
    p.add_argument("--h0-points",type=int,default=561)
    a=p.parse_args()
    run(a.output_dir,a.replicates,a.master_seed,a.n_parent,a.n_detected,a.h0_min,a.h0_max,a.h0_points)
if __name__=="__main__": main()
