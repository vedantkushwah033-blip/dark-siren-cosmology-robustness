"""Paired sensitivity of fixed-truth coverage to the prior support.

Uses identical synthetic observations across prior intervals. This is a toy-model
diagnostic; it does not alter the data-generating process to improve coverage.
"""
from __future__ import annotations
import argparse, csv, importlib.util
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

BASE = Path(__file__).resolve().parent / "run_fixed_truth_coverage_map.py"
SPEC = importlib.util.spec_from_file_location("coverage_map_source", BASE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Unable to load shared fixed-truth generator.")
coverage_source = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coverage_source)

TRUTHS = (52.0, 70.0, 88.0)
PRIORS = ((50.0, 90.0), (40.0, 100.0), (30.0, 110.0))
METHODS = ("uncorrected", "corrected")

def wilson(successes, n, z=1.959963984540054):
    if n <= 0: raise ValueError("n must be positive")
    p=successes/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (0.0 if successes == 0 else max(0.0,float(center-half)),
            1.0 if successes == n else min(1.0,float(center+half)))

def run_audit(replicates=300, seed=20261024, grid_size=4001):
    if replicates < 100: raise ValueError("Use at least 100 replicates per truth.")
    if grid_size < 1001: raise ValueError("grid_size must be at least 1001.")
    rng=np.random.default_rng(seed)
    grid=np.linspace(25.0,125.0,grid_size)
    rows=[]; summaries=[]
    for truth in TRUTHS:
        observations=[coverage_source.draw_detected_observation(truth,rng) for _ in range(replicates)]
        for prior_min,prior_max in PRIORS:
            cell=[]
            for rep,(host,distance) in enumerate(observations):
                base=posterior_h0(grid,coverage_source.HOST_REDSHIFTS,distance,
                    coverage_source.SIGMA_MPC,weights=coverage_source.HOST_WEIGHTS,
                    prior_min=prior_min,prior_max=prior_max)
                corr=posterior_h0_selected(grid,coverage_source.HOST_REDSHIFTS,distance,
                    coverage_source.SIGMA_MPC,coverage_source.DETECTION_MAX,
                    weights=coverage_source.HOST_WEIGHTS,prior_min=prior_min,
                    prior_max=prior_max,detection_distance_min=coverage_source.DETECTION_MIN)
                row={"true_H0":truth,"replicate":rep,"host_index":host,
                     "observed_distance_mpc":distance,"prior_min":prior_min,"prior_max":prior_max}
                for method,post in (("uncorrected",base),("corrected",corr)):
                    sm=posterior_summary(grid,post)
                    for key in ("median","lower_68","upper_68","lower_95","upper_95"):
                        row[f"{method}_{key}"]=float(sm[key])
                rows.append(row); cell.append(row)
            for method in METHODS:
                err=np.array([r[f"{method}_median"]-truth for r in cell])
                hit68=np.array([r[f"{method}_lower_68"]<=truth<=r[f"{method}_upper_68"] for r in cell])
                hit95=np.array([r[f"{method}_lower_95"]<=truth<=r[f"{method}_upper_95"] for r in cell])
                lo68,hi68=wilson(int(hit68.sum()),replicates)
                lo95,hi95=wilson(int(hit95.sum()),replicates)
                summaries.append({"true_H0":truth,"method":method,"n":replicates,
                    "prior_min":prior_min,"prior_max":prior_max,"seed":seed,"grid_size":grid_size,
                    "bias":float(err.mean()),"rmse":float(np.sqrt(np.mean(err**2))),
                    "coverage_68":float(hit68.mean()),"coverage_68_wilson_low":lo68,
                    "coverage_68_wilson_high":hi68,"coverage_95":float(hit95.mean()),
                    "coverage_95_wilson_low":lo95,"coverage_95_wilson_high":hi95})
    return rows,summaries

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261024)
    p.add_argument("--grid-size",type=int,default=4001)
    p.add_argument("--output-dir",default="results/prior-support-sensitivity")
    a=p.parse_args(); rows,summaries=run_audit(a.replicates,a.seed,a.grid_size)
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for name,data in (("per_run.csv",rows),("summary.csv",summaries)):
        path=out/name
        fields=list(dict.fromkeys(k for r in data for k in r))
        with path.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
        with path.open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f); read=list(reader)
            if reader.fieldnames!=fields or len(read)!=len(data): raise RuntimeError("Artifact validation failed.")
    expected= len(TRUTHS)*len(PRIORS)*replicates
    if len(rows)!=expected or len(summaries)!=len(TRUTHS)*len(PRIORS)*len(METHODS):
        raise RuntimeError("Incomplete prior-support experiment.")
    print(f"PASS: {len(rows)} paired observation/prior rows and {len(summaries)} summaries.")
    for r in summaries: print(r)
if __name__=="__main__": main()
