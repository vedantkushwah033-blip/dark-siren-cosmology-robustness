"""Selection-threshold misspecification audit for the synthetic dark-siren model.

Data are generated under one known finite detection interval; inference is
repeated with the correct upper threshold and deliberately overestimated
thresholds. Synthetic diagnostic only, not an LVK detector model.
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

HOST_REDSHIFTS=np.array([0.02,0.04])
HOST_WEIGHTS=np.array([0.7,0.3])
SIGMA=20.0
TRUE_MIN=0.0
TRUE_MAX=120.0
INFERRED_MAXIMA=(120.0,150.0,180.0,210.0)
TRUTHS=(52.0,60.0,70.0,80.0,88.0)

def draw_detected(truth,rng,max_attempts=100000):
    for _ in range(max_attempts):
        host=int(rng.choice(len(HOST_REDSHIFTS),p=HOST_WEIGHTS))
        mean=C_KM_S*HOST_REDSHIFTS[host]/truth
        d=float(rng.normal(mean,SIGMA))
        if TRUE_MIN < d < TRUE_MAX:
            return host,d
    raise RuntimeError("Detection rejection sampler exceeded attempt limit.")

def run(replicates=300,seed=20261027,grid_size=4001):
    if replicates<100: raise ValueError("Use at least 100 replicates per truth.")
    if grid_size<1001: raise ValueError("Use grid_size >= 1001.")
    rng=np.random.default_rng(seed)
    grid=np.linspace(25.0,125.0,grid_size)
    rows=[]; summaries=[]
    for truth in TRUTHS:
        observations=[draw_detected(truth,rng) for _ in range(replicates)]
        for rep,(host,distance) in enumerate(observations):
            base=posterior_h0(grid,HOST_REDSHIFTS,distance,SIGMA,weights=HOST_WEIGHTS,prior_min=50,prior_max=90)
            bsum=posterior_summary(grid,base)
            for inferred_max in INFERRED_MAXIMA:
                corr=posterior_h0_selected(grid,HOST_REDSHIFTS,distance,SIGMA,inferred_max,
                    weights=HOST_WEIGHTS,prior_min=50,prior_max=90,detection_distance_min=TRUE_MIN)
                sm=posterior_summary(grid,corr)
                rows.append({"true_H0":truth,"replicate":rep,"host_index":host,
                    "observed_distance_mpc":distance,"true_detection_max_mpc":TRUE_MAX,
                    "inferred_detection_max_mpc":inferred_max,
                    "uncorrected_median":float(bsum["median"]),
                    "corrected_median":float(sm["median"]),
                    "corrected_lower_68":float(sm["lower_68"]),
                    "corrected_upper_68":float(sm["upper_68"]),
                    "corrected_lower_95":float(sm["lower_95"]),
                    "corrected_upper_95":float(sm["upper_95"])})
        for inferred_max in INFERRED_MAXIMA:
            cell=[r for r in rows if r["true_H0"]==truth and r["inferred_detection_max_mpc"]==inferred_max]
            errors=np.array([r["corrected_median"]-truth for r in cell])
            hit68=np.array([r["corrected_lower_68"]<=truth<=r["corrected_upper_68"] for r in cell])
            hit95=np.array([r["corrected_lower_95"]<=truth<=r["corrected_upper_95"] for r in cell])
            summaries.append({"true_H0":truth,"n":replicates,"seed":seed,
                "true_detection_max_mpc":TRUE_MAX,"inferred_detection_max_mpc":inferred_max,
                "bias":float(errors.mean()),"mae":float(np.mean(np.abs(errors))),
                "rmse":float(np.sqrt(np.mean(errors**2))),
                "coverage_68":float(hit68.mean()),"coverage_95":float(hit95.mean()),
                "mean_width_68":float(np.mean([r["corrected_upper_68"]-r["corrected_lower_68"] for r in cell])),
                "mean_width_95":float(np.mean([r["corrected_upper_95"]-r["corrected_lower_95"] for r in cell]))})
    return rows,summaries

def paired_threshold_contrasts(rows, seed=20261028, bootstrap_draws=5000):
    """Paired wrong-minus-correct threshold contrasts with percentile bootstrap CIs."""
    if bootstrap_draws < 1000:
        raise ValueError("Use at least 1000 bootstrap draws")
    rng=np.random.default_rng(seed)
    lookup={(r["true_H0"],r["replicate"],r["inferred_detection_max_mpc"]):r for r in rows}
    metrics={
        "posterior_median":lambda r:r["corrected_median"],
        "absolute_error":lambda r:abs(r["corrected_median"]-r["true_H0"]),
        "squared_error":lambda r:(r["corrected_median"]-r["true_H0"])**2,
        "width_68":lambda r:r["corrected_upper_68"]-r["corrected_lower_68"],
        "width_95":lambda r:r["corrected_upper_95"]-r["corrected_lower_95"],
        "coverage_68_indicator":lambda r:float(r["corrected_lower_68"]<=r["true_H0"]<=r["corrected_upper_68"]),
        "coverage_95_indicator":lambda r:float(r["corrected_lower_95"]<=r["true_H0"]<=r["corrected_upper_95"]),
    }
    out=[]
    for truth in TRUTHS:
        for wrong in INFERRED_MAXIMA:
            if wrong==TRUE_MAX: continue
            baseline=[lookup[(truth,rep,TRUE_MAX)] for rep in range(max(r["replicate"] for r in rows if r["true_H0"]==truth)+1)]
            comparison=[lookup[(truth,rep,wrong)] for rep in range(len(baseline))]
            for metric,fn in metrics.items():
                deltas=np.asarray([fn(b)-fn(a) for a,b in zip(baseline,comparison)],dtype=float)
                if not np.all(np.isfinite(deltas)): raise RuntimeError(f"Non-finite contrast: {metric}")
                indices=rng.integers(0,len(deltas),size=(bootstrap_draws,len(deltas)))
                means=deltas[indices].mean(axis=1)
                out.append({"true_H0":truth,"correct_threshold_mpc":TRUE_MAX,
                    "assumed_threshold_mpc":wrong,"metric":metric,"n_pairs":len(deltas),
                    "mean_wrong_minus_correct":float(deltas.mean()),
                    "median_paired_difference":float(np.median(deltas)),
                    "bootstrap_95_low":float(np.quantile(means,.025)),
                    "bootstrap_95_high":float(np.quantile(means,.975)),
                    "bootstrap_draws":bootstrap_draws,"bootstrap_seed":seed})
    expected=len(TRUTHS)*(len(INFERRED_MAXIMA)-1)*len(metrics)
    if len(out)!=expected: raise RuntimeError(f"Expected {expected} contrast rows, got {len(out)}")
    return out

def write_csv(path,rows):
    if not rows: raise RuntimeError(f"Refusing empty artifact {path}")
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    with path.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f);check=list(rd)
        if rd.fieldnames!=fields or len(check)!=len(rows): raise RuntimeError("CSV round-trip failed")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261027)
    p.add_argument("--grid-size",type=int,default=4001)
    p.add_argument("--output-dir",default="results/threshold-misspecification")
    a=p.parse_args(); rows,summaries=run(a.replicates,a.seed,a.grid_size)
    out=Path(a.output_dir);write_csv(out/"per_run.csv",rows);write_csv(out/"summary.csv",summaries)
    expected=len(TRUTHS)*a.replicates*len(INFERRED_MAXIMA)
    if len(rows)!=expected: raise RuntimeError(f"Expected {expected} rows, got {len(rows)}")
    print(f"PASS: {len(rows)} per-run rows; {len(summaries)} summaries; {len(contrasts)} paired threshold contrasts")
    for r in summaries: print(r)
if __name__=="__main__": main()
