"""Diagnose fixed-truth H0 bias with exact host density (no KDE or finite catalogue).

Prior-averaged calibration and conditional bias at fixed truth are different
questions. This is a controlled toy model, not real gravitational-wave inference.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import chi2
from dark_siren_robustness.generative_selection import latent_galaxy_intensity, C_KM_S
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood

SIGMA, THRESHOLD, N_DETECTED = 140.0, 250.0, 40


def exact_density(z):
    density=latent_galaxy_intensity(z)*(1.0+1.5*(z-z[0])/(z[-1]-z[0]))
    return density/trapezoid(density,z)


def make_sampler(zmax=.40,zmin=.01,n=40001):
    grid=np.linspace(zmin,zmax,n); density=exact_density(grid)
    cdf=np.r_[0.,np.cumsum(.5*(density[1:]+density[:-1])*np.diff(grid))]
    cdf/=cdf[-1]
    return lambda rng,size:np.interp(rng.random(size),cdf,grid)


def simulate(rng,h0,sample_z):
    distances=[]
    while len(distances)<N_DETECTED:
        z=sample_z(rng,4*N_DETECTED)
        d=rng.normal(C_KM_S*z/h0,SIGMA)
        distances.extend(d[d>=THRESHOLD].tolist())
    return np.asarray(distances[:N_DETECTED])


def run(h0_mode,reps,seed,z_points,h0_points,h0_lo=50.,h0_hi=90.,key="selected_log_likelihood"):
    if reps<2: raise ValueError("reps must be at least two")
    rng=np.random.default_rng(seed); sample_z=make_sampler()
    z=np.linspace(.01,.40,z_points); density=exact_density(z); grid=np.linspace(h0_lo,h0_hi,h0_points)
    errors=[]; median_errors=[]; cover68=[]; cover90=[]; ranks=[]
    for _ in range(reps):
        h0=rng.uniform(h0_lo,h0_hi) if h0_mode=="prior" else float(h0_mode)
        d=simulate(rng,h0,sample_z)
        likelihood=selection_aware_toy_likelihood(z,density,grid,d,sigma=SIGMA,detection_threshold=THRESHOLD)
        if key not in {"selected_log_likelihood","unconditional_log_likelihood"}: raise ValueError("invalid likelihood key")
        loglike=likelihood[key]
        s=posterior_summary(grid,loglike)
        p=np.exp(loglike-loglike.max())
        cdf=np.r_[0.,np.cumsum(.5*(p[1:]+p[:-1])*np.diff(grid))]
        cdf/=cdf[-1]
        ranks.append(float(np.interp(h0,grid,cdf)))
        errors.append(s["posterior_mean"]-h0); median_errors.append(s["posterior_median"]-h0)
        cover68.append(s["ci68_low"]<=h0<=s["ci68_high"])
        q05,q95=np.interp([.05,.95],cdf,grid); cover90.append(q05<=h0<=q95)
    e=np.asarray(errors); m=np.asarray(median_errors); r=np.asarray(ranks)
    counts,_=np.histogram(r,bins=10,range=(0,1))
    x2=float(np.sum((counts-reps/10.)**2/(reps/10.)))
    return {"h0":h0_mode,"reps":reps,"z_points":z_points,"h0_points":h0_points,
            "mean_bias":float(e.mean()),"mean_bias_se":float(e.std(ddof=1)/np.sqrt(reps)),
            "median_bias":float(m.mean()),"median_bias_se":float(m.std(ddof=1)/np.sqrt(reps)),
            "rmse":float(np.sqrt(np.mean(e**2))),"cov68":float(np.mean(cover68)),
            "cov90":float(np.mean(cover90)),"pit_chi2_9df":x2,"pit_p":float(chi2.sf(x2,9)),
            "pit_counts":counts.tolist(),"likelihood_component":key}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps",type=int,default=2000); parser.add_argument("--seed",type=int,default=20261029)
    parser.add_argument("--out",type=Path,required=True); parser.add_argument("--modes",nargs="+",default=["70","60","80","prior"])
    parser.add_argument("--z-points",type=int,default=801); parser.add_argument("--h0-points",type=int,default=161)
    args=parser.parse_args(); results=[]
    for mode in args.modes:
        result=run(mode if mode=="prior" else float(mode),args.reps,args.seed,args.z_points,args.h0_points)
        results.append(result); print(result,flush=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
