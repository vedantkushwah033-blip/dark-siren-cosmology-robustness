"""Posterior-level calibration of a synthetic two-host H0 mixture.

Same synthetic observation is evaluated under both catalogue-weight choices,
two prior supports, and uniform/log-uniform prior measures. This is not a
real gravitational-wave analysis.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm

C = 299792.458
TRUTHS = (55.0, 65.0, 75.0, 85.0)
NOISES = (0.03, 0.10, 0.20)
Z1, Z2 = 0.02, 0.04
TRUE_WEIGHTS = (0.7, 0.3)
INFERENCE_WEIGHTS = {"correct": (0.7, 0.3), "reversed": (0.3, 0.7)}
SUPPORTS = {"narrow_50_90": (50.0, 90.0), "broad_20_220": (20.0, 220.0)}
MEASURES = ("uniform_H0", "log_uniform_H0")
GRID = np.linspace(10.0, 250.0, 2401)


def posterior_summary(distance, sigma, weights, support, measure, truth):
    lo, hi = support
    mask = (GRID >= lo) & (GRID <= hi)
    h = GRID[mask]
    if h.size < 2:
        raise ValueError("Prior support contains fewer than two grid points")
    likelihood = (
        weights[0] * norm.pdf(distance, loc=C * Z1 / h, scale=sigma)
        + weights[1] * norm.pdf(distance, loc=C * Z2 / h, scale=sigma)
    )
    prior = np.ones_like(h) if measure == "uniform_H0" else 1.0 / h
    if measure not in MEASURES:
        raise ValueError(f"Unknown prior measure {measure}")
    density = likelihood * prior
    area = float(trapezoid(density, h))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Posterior normalization failed")
    density = density / area
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, h)))
    cdf = cdf / cdf[-1]
    def quantile(q):
        return float(np.interp(q, cdf, h))
    median, lower68, upper68 = quantile(.5), quantile(.16), quantile(.84)
    lower95, upper95 = quantile(.025), quantile(.975)
    return {
        "median": median, "lower68": lower68, "upper68": upper68,
        "lower95": lower95, "upper95": upper95,
        "covered68": int(lower68 <= truth <= upper68),
        "covered95": int(lower95 <= truth <= upper95),
    }


def run(replicates=200, seed=20261030):
    if replicates < 100:
        raise ValueError("Protocol requires at least 100 replicates per truth/noise cell")
    rng = np.random.default_rng(seed)
    raw = []
    summaries = []
    for truth in TRUTHS:
        for noise in NOISES:
            sigma = noise * C * Z2 / truth
            for rep in range(replicates):
                host = int(rng.choice(2, p=TRUE_WEIGHTS))
                zhost = Z1 if host == 0 else Z2
                distance_true = C * zhost / truth
                distance = float(rng.normal(distance_true, sigma))
                while distance <= 0:
                    distance = float(rng.normal(distance_true, sigma))
                row = {
                    "true_H0": truth, "noise_fraction": noise, "replicate": rep,
                    "base_seed": seed, "true_host_index": host, "true_host_redshift": zhost,
                    "observed_distance_mpc": distance, "sigma_mpc": sigma,
                }
                for support_name, support in SUPPORTS.items():
                    for measure in MEASURES:
                        for weight_name, weights in INFERENCE_WEIGHTS.items():
                            key = f"{support_name}_{measure}_{weight_name}"
                            row.update({f"{key}_{k}": v for k, v in posterior_summary(
                                distance, sigma, weights, support, measure, truth).items()})
                raw.append(row)
            block = [r for r in raw if r["true_H0"] == truth and r["noise_fraction"] == noise]
            for support_name, support in SUPPORTS.items():
                for measure in MEASURES:
                    for weight_name in INFERENCE_WEIGHTS:
                        key = f"{support_name}_{measure}_{weight_name}"
                        errors = np.array([float(r[f"{key}_median"]) - truth for r in block])
                        summaries.append({
                            "true_H0": truth, "noise_fraction": noise,
                            "prior_support": support_name, "prior_min": support[0], "prior_max": support[1],
                            "prior_measure": measure, "inference_weights": weight_name,
                            "replicates": replicates, "mean_median_error": float(errors.mean()),
                            "median_absolute_error": float(np.median(np.abs(errors))),
                            "rmse_median": float(np.sqrt(np.mean(errors ** 2))),
                            "mean_width68": float(np.mean([float(r[f"{key}_upper68"])-float(r[f"{key}_lower68"]) for r in block])),
                            "mean_width95": float(np.mean([float(r[f"{key}_upper95"])-float(r[f"{key}_lower95"]) for r in block])),
                            "coverage68": float(np.mean([int(r[f"{key}_covered68"]) for r in block])),
                            "coverage95": float(np.mean([int(r[f"{key}_covered95"]) for r in block])),
                        })
    return raw, summaries


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates", type=int, default=200)
    p.add_argument("--seed", type=int, default=20261030)
    p.add_argument("--output-dir", default="results/posterior-calibration")
    a = p.parse_args()
    raw, summaries = run(a.replicates, a.seed)
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "per_run.csv", raw); write_csv(out / "summary.csv", summaries)
    meta = {"experiment":"posterior-level H0 calibration in a synthetic two-host mixture",
            "audit_status":"PENDING","base_seed":a.seed,"replicates_per_truth_noise_cell":a.replicates,
            "raw_records":len(raw),"summary_rows":len(summaries),
            "sha256":{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ("per_run.csv","summary.csv")}}
    (out/"run_metadata.json").write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    print(f"Generated {len(raw)} synthetic observations and {len(summaries)} posterior-calibration cells.")

if __name__ == "__main__": main()
