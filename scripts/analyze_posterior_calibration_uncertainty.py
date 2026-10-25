"""Stratified-bootstrap uncertainty for posterior-calibration metrics.

Requires an independently audited posterior-calibration artifact. Resamples
replicates within each fixed truth/noise cell, preserving the design allocation.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np

TRUTHS=(55.0,65.0,75.0,85.0)
NOISES=(0.03,0.10,0.20)
SUPPORTS=("narrow_50_90","broad_20_220")
MEASURES=("uniform_H0","log_uniform_H0")
WEIGHTS=("correct","reversed")


def read_csv(path):
    with path.open(newline="",encoding="utf-8") as f: return list(csv.DictReader(f))


def stratified_bootstrap(values_by_cell, seed, resamples=4000, batch_size=100):
    """Return 95% percentile intervals for mean, RMSE, and a binary rate."""
    x=np.asarray(values_by_cell,dtype=float)
    if x.ndim!=2 or x.shape[0]!=len(TRUTHS)*len(NOISES) or x.shape[1]<2:
        raise ValueError("expected a 12-cell array with at least two replicates per cell")
    if not np.all(np.isfinite(x)): raise ValueError("non-finite bootstrap input")
    rng=np.random.default_rng(seed)
    mean_draws=np.empty(resamples); rmse_draws=np.empty(resamples)
    for start in range(0,resamples,batch_size):
        stop=min(start+batch_size,resamples)
        idx=rng.integers(0,x.shape[1],size=(stop-start,x.shape[0],x.shape[1]))
        sampled=x[np.arange(x.shape[0])[None,:,None],idx]
        mean_draws[start:stop]=sampled.mean(axis=(1,2))
        rmse_draws[start:stop]=np.sqrt(np.mean(sampled**2,axis=(1,2)))
    return {
        "mean_ci95_low":float(np.quantile(mean_draws,.025)),
        "mean_ci95_high":float(np.quantile(mean_draws,.975)),
        "rmse_ci95_low":float(np.quantile(rmse_draws,.025)),
        "rmse_ci95_high":float(np.quantile(rmse_draws,.975)),
    }


def analyze(root, seed=20261031, resamples=4000):
    meta=json.loads((root/"run_metadata.json").read_text(encoding="utf-8"))
    audit=json.loads((root/"audit_report.json").read_text(encoding="utf-8"))
    if meta.get("audit_status")!="PASS" or audit.get("audit")!="PASS":
        raise ValueError("posterior-calibration artifact must pass independent audit")
    for name in ("per_run.csv","summary.csv"):
        if meta.get("sha256",{}).get(name)!=hashlib.sha256((root/name).read_bytes()).hexdigest():
            raise ValueError(f"input hash mismatch: {name}")
    raw=read_csv(root/"per_run.csv")
    cells={}
    for row in raw:
        key=(float(row["true_H0"]),float(row["noise_fraction"]))
        cells.setdefault(key,[]).append(row)
    expected={(t,n) for t in TRUTHS for n in NOISES}
    if set(cells)!=expected or len(raw)!=len(expected)*int(meta["replicates_per_truth_noise_cell"]):
        raise ValueError("raw design does not match prespecified truth/noise cells")
    reps=int(meta["replicates_per_truth_noise_cell"])
    if any(len(v)!=reps for v in cells.values()): raise ValueError("unequal or incomplete replicate cells")
    results=[]
    arm_idx=0
    for support in SUPPORTS:
        for measure in MEASURES:
            for weight in WEIGHTS:
                prefix=f"{support}_{measure}_{weight}"
                errors=[]; cov68=[]; cov95=[]; widths68=[]; widths95=[]
                for truth in TRUTHS:
                    for noise in NOISES:
                        block=sorted(cells[(truth,noise)],key=lambda r:int(r["replicate"]))
                        errors.append(np.array([float(r[f"{prefix}_median"])-truth for r in block]))
                        cov68.append(np.array([int(r[f"{prefix}_covered68"]) for r in block],dtype=float))
                        cov95.append(np.array([int(r[f"{prefix}_covered95"]) for r in block],dtype=float))
                        widths68.append(np.array([float(r[f"{prefix}_upper68"])-float(r[f"{prefix}_lower68"]) for r in block]))
                        widths95.append(np.array([float(r[f"{prefix}_upper95"])-float(r[f"{prefix}_lower95"]) for r in block]))
                e=np.stack(errors); c68=np.stack(cov68); c95=np.stack(cov95)
                w68=np.stack(widths68); w95=np.stack(widths95)
                # Distinct deterministic seed streams for each arm/metric.
                ci_bias=stratified_bootstrap(e,seed+arm_idx*101+1,resamples)
                ci_rmse=stratified_bootstrap(e,seed+arm_idx*101+2,resamples)
                ci_cov68=stratified_bootstrap(c68,seed+arm_idx*101+3,resamples)
                ci_cov95=stratified_bootstrap(c95,seed+arm_idx*101+4,resamples)
                ci_w68=stratified_bootstrap(w68,seed+arm_idx*101+5,resamples)
                ci_w95=stratified_bootstrap(w95,seed+arm_idx*101+6,resamples)
                results.append({
                    "prior_support":support,"prior_measure":measure,"inference_weights":weight,
                    "replicates_per_cell":reps,"cells":len(expected),"bootstrap_resamples":resamples,
                    "seed":seed+arm_idx*101,
                    "bias":float(e.mean()),"bias_ci95_low":ci_bias["mean_ci95_low"],"bias_ci95_high":ci_bias["mean_ci95_high"],
                    "rmse":float(np.sqrt(np.mean(e**2))),"rmse_ci95_low":ci_rmse["rmse_ci95_low"],"rmse_ci95_high":ci_rmse["rmse_ci95_high"],
                    "coverage68":float(c68.mean()),"coverage68_ci95_low":ci_cov68["mean_ci95_low"],"coverage68_ci95_high":ci_cov68["mean_ci95_high"],
                    "coverage95":float(c95.mean()),"coverage95_ci95_low":ci_cov95["mean_ci95_low"],"coverage95_ci95_high":ci_cov95["mean_ci95_high"],
                    "mean_width68":float(w68.mean()),"mean_width68_ci95_low":ci_w68["mean_ci95_low"],"mean_width68_ci95_high":ci_w68["mean_ci95_high"],
                    "mean_width95":float(w95.mean()),"mean_width95_ci95_low":ci_w95["mean_ci95_low"],"mean_width95_ci95_high":ci_w95["mean_ci95_high"],
                })
                arm_idx+=1
    return results


def write_report(path, rows):
    lines=["# Stratified-bootstrap uncertainty for posterior calibration","",
        "Intervals are 95% percentile-bootstrap intervals. Replicates are resampled within each of the 12 fixed truth/noise cells, preserving the design allocation. These quantify Monte Carlo uncertainty for this simulation design; they are not uncertainty intervals for real-event H₀ inference.","",
        "| Support | Prior measure | Weights | Bias (95% CI) | RMSE (95% CI) | 68% coverage (95% CI) | 95% coverage (95% CI) |",
        "|---|---|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['prior_support']} | {r['prior_measure']} | {r['inference_weights']} | {r['bias']:.2f} [{r['bias_ci95_low']:.2f}, {r['bias_ci95_high']:.2f}] | {r['rmse']:.2f} [{r['rmse_ci95_low']:.2f}, {r['rmse_ci95_high']:.2f}] | {100*r['coverage68']:.1f}% [{100*r['coverage68_ci95_low']:.1f}%, {100*r['coverage68_ci95_high']:.1f}%] | {100*r['coverage95']:.1f}% [{100*r['coverage95_ci95_low']:.1f}%, {100*r['coverage95_ci95_high']:.1f}%] |")
    lines += ["","Bootstrap: 4,000 resamples per reported metric and inference arm, deterministic seed schedule beginning at 20261031. The intervals are conditional on the fixed truth/noise grid and the specified synthetic data-generating process.",""]
    path.write_text("\n".join(lines),encoding="utf-8")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("results_dir",nargs="?",default="results/posterior-calibration")
    p.add_argument("--seed",type=int,default=20261031)
    p.add_argument("--bootstrap-resamples",type=int,default=4000)
    a=p.parse_args()
    root=Path(a.results_dir)
    rows=analyze(root,a.seed,a.bootstrap_resamples)
    with (root/"uncertainty.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    write_report(root/"uncertainty_report.md",rows)
    print((root/"uncertainty_report.md").read_text(encoding="utf-8"))

if __name__=="__main__": main()
