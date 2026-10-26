"""Independently reconstruct posterior-level H0 calibration outputs."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm

C = 299792.458
TRUTHS = (55.0, 65.0, 75.0, 85.0)
NOISES = (0.03, 0.10, 0.20)
Z1, Z2 = 0.02, 0.04
SUPPORTS = {"narrow_50_90": (50.0, 90.0), "broad_20_220": (20.0, 220.0)}
MEASURES = ("uniform_H0", "log_uniform_H0")
WEIGHTS = {"correct": (0.7, 0.3), "reversed": (0.3, 0.7)}
GRID = np.linspace(10.0, 250.0, 2401)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconstruct(distance, sigma, weights, support, measure, truth):
    mask = (GRID >= support[0]) & (GRID <= support[1])
    h = GRID[mask]
    likelihood = weights[0] * norm.pdf(distance, loc=C*Z1/h, scale=sigma)
    likelihood += weights[1] * norm.pdf(distance, loc=C*Z2/h, scale=sigma)
    prior = np.ones_like(h) if measure == "uniform_H0" else 1.0/h
    density = likelihood * prior
    area = float(trapezoid(density, h))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("posterior normalization failed")
    density /= area
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, h)))
    cdf /= cdf[-1]
    quant = lambda q: float(np.interp(q, cdf, h))
    l68, u68, l95, u95 = quant(.16), quant(.84), quant(.025), quant(.975)
    return {"median":quant(.5),"lower68":l68,"upper68":u68,"lower95":l95,"upper95":u95,
            "covered68":int(l68 <= truth <= u68),"covered95":int(l95 <= truth <= u95)}


def audit(root):
    raw_path, sum_path, meta_path = root/"per_run.csv", root/"summary.csv", root/"run_metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("audit_status") != "PENDING":
        raise ValueError("metadata must be PENDING before independent audit")
    for name, path in (("per_run.csv",raw_path),("summary.csv",sum_path)):
        if meta.get("sha256",{}).get(name) != sha(path):
            raise ValueError(f"SHA-256 mismatch for {name}")
    # JSON booleans are Python ints, and int(100.5) silently truncates.
    # Reject both so malformed provenance cannot pass by coercion.
    reps_value = meta.get("replicates_per_truth_noise_cell")
    seed_value = meta.get("base_seed")
    if type(reps_value) is not int or reps_value < 100:
        raise ValueError("replicates_per_truth_noise_cell must be an integer >= 100")
    if type(seed_value) is not int or seed_value < 0:
        raise ValueError("base_seed must be a non-negative integer")
    reps, seed = reps_value, seed_value
    raw, summaries = read_csv(raw_path), read_csv(sum_path)
    if len(raw) != len(TRUTHS)*len(NOISES)*reps: raise ValueError("raw row count mismatch")
    if len(summaries) != len(TRUTHS)*len(NOISES)*len(SUPPORTS)*len(MEASURES)*len(WEIGHTS):
        raise ValueError("summary row count mismatch")
    cells = {}
    metric_names = ("median","lower68","upper68","lower95","upper95","covered68","covered95")
    for i, row in enumerate(raw):
        truth, noise, rep = float(row["true_H0"]), float(row["noise_fraction"]), int(row["replicate"])
        if truth not in TRUTHS or noise not in NOISES or not 0 <= rep < reps:
            raise ValueError(f"invalid design key at raw row {i}")
        if int(row["base_seed"]) != seed: raise ValueError(f"seed mismatch at raw row {i}")
        host = int(row["true_host_index"])
        zhost = Z1 if host == 0 else Z2 if host == 1 else None
        if zhost is None or not math.isclose(float(row["true_host_redshift"]),zhost,abs_tol=1e-12):
            raise ValueError(f"host identity mismatch at raw row {i}")
        sigma = noise*C*Z2/truth
        if not math.isclose(float(row["sigma_mpc"]),sigma,rel_tol=1e-12,abs_tol=1e-12):
            raise ValueError(f"sigma mismatch at raw row {i}")
        d = float(row["observed_distance_mpc"])
        if not math.isfinite(d) or d <= 0: raise ValueError(f"invalid distance at raw row {i}")
        key=(truth,noise); cell=cells.setdefault(key,{})
        if rep in cell: raise ValueError(f"duplicate replicate in cell {key}")
        cell[rep]=row
        for sname,support in SUPPORTS.items():
            for measure in MEASURES:
                for wname,weights in WEIGHTS.items():
                    prefix=f"{sname}_{measure}_{wname}"
                    expected=reconstruct(d,sigma,weights,support,measure,truth)
                    for metric in metric_names:
                        actual=float(row[f"{prefix}_{metric}"])
                        target=expected[metric]
                        if not math.isfinite(actual) or (metric.startswith("covered") and int(actual)!=target) or (not metric.startswith("covered") and not math.isclose(actual,target,rel_tol=1e-9,abs_tol=1e-8)):
                            raise ValueError(f"{prefix}_{metric} reconstruction mismatch at raw row {i}")
    expected_cells={(t,n) for t in TRUTHS for n in NOISES}
    if set(cells)!=expected_cells or any(set(c)!=set(range(reps)) for c in cells.values()):
        raise ValueError("incomplete replicate design")
    expected_summaries=[]
    for truth in TRUTHS:
        for noise in NOISES:
            block=[cells[(truth,noise)][r] for r in range(reps)]
            for sname,support in SUPPORTS.items():
                for measure in MEASURES:
                    for wname in WEIGHTS:
                        prefix=f"{sname}_{measure}_{wname}"
                        errors=np.array([float(r[f"{prefix}_median"])-truth for r in block])
                        expected_summaries.append({
                            "true_H0":truth,"noise_fraction":noise,"prior_support":sname,
                            "prior_min":support[0],"prior_max":support[1],"prior_measure":measure,
                            "inference_weights":wname,"replicates":reps,
                            "mean_median_error":float(errors.mean()),
                            "median_absolute_error":float(np.median(np.abs(errors))),
                            "rmse_median":float(np.sqrt(np.mean(errors**2))),
                            "mean_width68":float(np.mean([float(r[f"{prefix}_upper68"])-float(r[f"{prefix}_lower68"]) for r in block])),
                            "mean_width95":float(np.mean([float(r[f"{prefix}_upper95"])-float(r[f"{prefix}_lower95"]) for r in block])),
                            "coverage68":float(np.mean([int(r[f"{prefix}_covered68"]) for r in block])),
                            "coverage95":float(np.mean([int(r[f"{prefix}_covered95"]) for r in block])),
                        })
    if len(summaries)!=len(expected_summaries): raise ValueError("summary row count mismatch")
    numeric=("true_H0","noise_fraction","prior_min","prior_max","mean_median_error","median_absolute_error","rmse_median","mean_width68","mean_width95","coverage68","coverage95")
    for i,(a,e) in enumerate(zip(summaries,expected_summaries)):
        for k in ("prior_support","prior_measure","inference_weights"):
            if a[k]!=e[k]: raise ValueError(f"{k} mismatch at summary row {i}")
        for k in numeric:
            x,y=float(a[k]),float(e[k])
            if not math.isfinite(x) or not math.isclose(x,y,rel_tol=1e-9,abs_tol=1e-8):
                raise ValueError(f"{k} mismatch at summary row {i}")
        if int(a["replicates"])!=reps: raise ValueError(f"replicate count mismatch at summary row {i}")
    meta["audit_status"]="PASS"
    meta["independent_audit"]={"raw_records":len(raw),"summary_rows":len(summaries),
        "posterior_summaries_reconstructed_from_observations":True,
        "calibration_metrics_reconstructed":True,"hashes_verified":True}
    meta_path.write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    report={"audit":"PASS","raw_records":len(raw),"summary_rows":len(summaries),
        "posterior_summaries_reconstructed_from_observations":True,
        "calibration_metrics_reconstructed":True,"hashes_verified":True}
    (root/"audit_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("results_dir",nargs="?",default="results/posterior-calibration")
    a=p.parse_args()
    print(json.dumps(audit(Path(a.results_dir)),indent=2))

if __name__=="__main__": main()
