"""Matched fixed-H0 coverage audit under the same noise design as SBC.

Unlike the legacy fractional-distance-noise audit, this design draws reported
sigma independently of H0 and host identity, exactly as in run_sbc_diagnostic.py.
This is a prior-predictive companion experiment with repeated fixed truths.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np

C_KM_S = 299792.458
H0_VALUES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
H0_MIN, H0_MAX = 50.0, 90.0

def posterior_cdf(grid: np.ndarray, posterior: np.ndarray) -> np.ndarray:
    if posterior.shape != grid.shape or np.any(~np.isfinite(posterior)) or np.any(posterior < 0):
        raise ValueError("posterior must be finite, nonnegative, and match grid")
    areas = 0.5 * (posterior[1:] + posterior[:-1]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(areas)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("posterior normalization failed")
    return cdf / cdf[-1]

def central_interval(grid: np.ndarray, cdf: np.ndarray, level: float) -> tuple[float,float]:
    if not 0 < level < 1:
        raise ValueError("level must be between 0 and 1")
    tail=(1-level)/2
    return float(np.interp(tail,cdf,grid)),float(np.interp(1-tail,cdf,grid))

def summarize(rows: list[dict], scenario: str, h0: float) -> dict:
    errors=np.asarray([r["median_error"] for r in rows],dtype=float)
    out={"scenario":scenario,"true_H0":h0,"replicates":len(rows),
         "mean_signed_error":float(errors.mean()),
         "rmse_median":float(np.sqrt(np.mean(errors**2)))}
    for level in (50,68,90,95):
        cov=float(np.mean([r[f"covers_{level}"] for r in rows]))
        out[f"coverage_{level}"]=cov
        out[f"coverage_{level}_mc_se"]=float(np.sqrt(cov*(1-cov)/len(rows)))
        out[f"mean_width_{level}"]=float(np.mean([r[f"width_{level}"] for r in rows]))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dir",required=True)
    p.add_argument("--replicates",type=int,default=500)
    p.add_argument("--seed",type=int,default=20261023)
    p.add_argument("--output-dir",default="results/matched-fixed-h0-coverage")
    p.add_argument("--catalogue-size",type=int,default=500)
    a=p.parse_args()
    if a.replicates<100: p.error("--replicates must be at least 100")
    original=Path(a.original_dir).resolve()
    sys.path.insert(0,str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host
    from src.bayesian_inference import posterior_h0
    grid=np.linspace(H0_MIN,H0_MAX,1601)
    raw=[]; summaries=[]
    for h_index,h0 in enumerate(H0_VALUES):
        scenario_rows={"known_host":[],"unknown_host":[]}
        for rep in range(a.replicates):
            seed=a.seed+h_index*10_000_000+rep*100
            rng=np.random.default_rng(seed)
            cat=simulate_galaxy_catalogue(n_galaxies=a.catalogue_size,H0_true=h0,seed=seed+1)
            host=choose_host(cat,weighting="uniform",seed=seed+2)
            sigma=float(rng.uniform(15.0,35.0))
            z_host=float(cat["redshift"][host])
            true_distance=C_KM_S*z_host/h0
            observed=float(rng.normal(true_distance,sigma))
            scenarios={"known_host":(np.asarray([z_host]),np.asarray([1.0])),
                       "unknown_host":(np.asarray(cat["redshift"],dtype=float),
                                       np.ones(len(cat["redshift"]),dtype=float))}
            for name,(zs,weights) in scenarios.items():
                post=np.asarray(posterior_h0(grid,zs,observed,sigma,weights=weights,
                                             prior_min=H0_MIN,prior_max=H0_MAX),dtype=float)
                cdf=posterior_cdf(grid,post)
                median=float(np.interp(.5,cdf,grid))
                row={"scenario":name,"true_H0":h0,"replicate":rep,"seed":seed,
                     "host_redshift":z_host,"observed_distance_mpc":observed,
                     "sigma_mpc":sigma,"posterior_median":median,
                     "median_error":median-h0}
                for level in (50,68,90,95):
                    lo,hi=central_interval(grid,cdf,level/100)
                    row[f"covers_{level}"]=int(lo<=h0<=hi)
                    row[f"width_{level}"]=hi-lo
                scenario_rows[name].append(row); raw.append(row)
        for name,rows in scenario_rows.items():
            summaries.append(summarize(rows,name,h0))
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for file,rows in (("per_run.csv",raw),("summary.csv",summaries)):
        with (out/file).open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(f"Validated design outputs: {len(summaries)} summaries; {len(raw)} posterior fits.")
    for row in summaries: print(row)
if __name__=="__main__": main()
