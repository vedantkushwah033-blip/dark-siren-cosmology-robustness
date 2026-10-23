"""Paired prior-boundary sensitivity audit for fixed-truth H0 coverage."""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
import numpy as np

C_KM_S=299792.458
H0_VALUES=(55.,60.,65.,70.,75.,80.,85.)
PRIORS={"baseline_50_90":(50.,90.),"wide_40_100":(40.,100.),"very_wide_30_120":(30.,120.)}

def cdf_from_density(grid, density):
    areas=0.5*(density[1:]+density[:-1])*np.diff(grid)
    cdf=np.concatenate(([0.],np.cumsum(areas)))
    if not np.isfinite(cdf[-1]) or cdf[-1]<=0: raise ValueError("invalid posterior normalization")
    return cdf/cdf[-1]

def interval(grid,cdf,level):
    tail=(1-level)/2
    return float(np.interp(tail,cdf,grid)),float(np.interp(1-tail,cdf,grid))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dir",required=True)
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261024)
    p.add_argument("--output-dir",default="results/prior-boundary-sensitivity")
    p.add_argument("--catalogue-size",type=int,default=500)
    a=p.parse_args()
    if a.replicates<100:p.error("--replicates must be at least 100")
    sys.path.insert(0,str(Path(a.original_dir).resolve()))
    from src.simulations_v2 import simulate_galaxy_catalogue,choose_host
    from src.bayesian_inference import posterior_h0
    grid=np.linspace(30.,120.,1801)
    raw=[]
    for ih,h0 in enumerate(H0_VALUES):
      for rep in range(a.replicates):
        seed=a.seed+ih*10_000_000+rep*100
        rng=np.random.default_rng(seed)
        cat=simulate_galaxy_catalogue(n_galaxies=a.catalogue_size,H0_true=h0,seed=seed+1)
        host=choose_host(cat,weighting="uniform",seed=seed+2)
        zhost=float(cat["redshift"][host])
        true_d=C_KM_S*zhost/h0
        sigma=float(rng.uniform(15.,35.))
        observed=float(rng.normal(true_d,sigma))
        scenarios={"known_host":(np.array([zhost]),np.array([1.])),
                   "unknown_host":(np.asarray(cat["redshift"],float),np.ones(len(cat["redshift"])))}
        for prior_name,(lo,hi) in PRIORS.items():
          for scenario,(zs,weights) in scenarios.items():
            post=np.asarray(posterior_h0(grid,zs,observed,sigma,weights=weights,prior_min=lo,prior_max=hi),float)
            cdf=cdf_from_density(grid,post)
            med=float(np.interp(.5,cdf,grid))
            # Quantify posterior probability concentrated near each prior edge.
            norm=float(np.trapezoid(post,grid))
            if not np.isfinite(norm) or norm<=0:
              raise ValueError(f"invalid posterior integral for {prior_name}/{scenario}/H0={h0}/rep={rep}")
            post_norm=post/norm
            edge_mass={}
            for band in (1.0,2.0):
              lower_mask=(grid>=lo)&(grid<=min(lo+band,hi))
              upper_mask=(grid<=hi)&(grid>=max(hi-band,lo))
              edge_mass[f"mass_within_{int(band)}_of_lower_prior"]=float(np.trapezoid(post_norm[lower_mask],grid[lower_mask]))
              edge_mass[f"mass_within_{int(band)}_of_upper_prior"]=float(np.trapezoid(post_norm[upper_mask],grid[upper_mask]))
            below_baseline=float(np.trapezoid(post_norm[grid<50.0],grid[grid<50.0]))
            above_baseline=float(np.trapezoid(post_norm[grid>90.0],grid[grid>90.0]))
            edge_mass["mass_below_50"]=below_baseline
            edge_mass["mass_above_90"]=above_baseline
            edge_mass["mass_outside_baseline_50_90"]=below_baseline+above_baseline
            row={"prior":prior_name,"prior_min":lo,"prior_max":hi,"scenario":scenario,
                 "true_H0":h0,"replicate":rep,"seed":seed,"host_redshift":zhost,
                 "observed_distance_mpc":observed,"sigma_mpc":sigma,
                 "posterior_median":med,"median_error":med-h0,**edge_mass}
            for lev in (68,95):
              left,right=interval(grid,cdf,lev/100)
              row[f"covers_{lev}"]=int(left<=h0<=right)
              row[f"width_{lev}"]=right-left
            raw.append(row)
    summaries=[]
    for prior in PRIORS:
      for h0 in H0_VALUES:
       for scenario in ("known_host","unknown_host"):
        rows=[r for r in raw if r["prior"]==prior and r["true_H0"]==h0 and r["scenario"]==scenario]
        err=np.array([r["median_error"] for r in rows],float)
        item={"prior":prior,"prior_min":PRIORS[prior][0],"prior_max":PRIORS[prior][1],
              "scenario":scenario,"true_H0":h0,"replicates":len(rows),
              "mean_signed_error":float(err.mean()),"rmse_median":float(np.sqrt(np.mean(err**2)))}
        for lev in (68,95):
          cov=float(np.mean([r[f"covers_{lev}"] for r in rows]))
          item[f"coverage_{lev}"]=cov
          item[f"coverage_{lev}_mc_se"]=float(np.sqrt(cov*(1-cov)/len(rows)))
          item[f"mean_width_{lev}"]=float(np.mean([r[f"width_{lev}"] for r in rows]))
        for band in (1,2):
          for edge in ("lower","upper"):
            key=f"mass_within_{band}_of_{edge}_prior"
            item[f"mean_{key}"]=float(np.mean([r[key] for r in rows]))
        for key in ("mass_below_50","mass_above_90","mass_outside_baseline_50_90"):
          item[f"mean_{key}"]=float(np.mean([r[key] for r in rows]))
        summaries.append(item)
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for name,rows in (("per_run.csv",raw),("summary.csv",summaries)):
      with (out/name).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    assert len(raw)==len(H0_VALUES)*a.replicates*len(PRIORS)*2,len(raw)
    assert len(summaries)==len(H0_VALUES)*len(PRIORS)*2,len(summaries)
    print(f"Validated {len(raw)} posterior fits and {len(summaries)} summary rows.")
    for row in summaries:
      if row["scenario"]=="unknown_host" and row["true_H0"] in (55.,60.,70.,85.):
        print(row)
if __name__=="__main__":main()
