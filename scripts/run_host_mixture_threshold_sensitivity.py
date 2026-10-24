"""Paired host-mixture and detection-threshold sensitivity factorial.

Synthetic diagnostic only; not a real GW detector or catalogue analysis.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import ndtr
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

TRUTHS = (60.0, 70.0, 80.0, 88.0)
HOST_PAIRS = ((0.02, 0.04), (0.02, 0.03), (0.02, 0.05))
WEIGHT_PAIRS = ((0.5, 0.5), (0.7, 0.3), (0.9, 0.1))
THRESHOLDS = (90.0, 120.0, 150.0)
SIGMA_MPC = 20.0
PRIOR = (50.0, 90.0)
GRID_BOUNDS = (25.0, 125.0)
METHODS = ("two_host_uncorrected", "two_host_corrected", "collapsed_uncorrected", "collapsed_corrected")


def cell_seed(seed, truth, redshifts, weights, threshold):
    key = json.dumps([int(seed), float(truth), list(redshifts), list(weights), float(threshold)], separators=(",", ":"))
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")


def expected_detected_host_fractions(truth, redshifts, weights, threshold, sigma=SIGMA_MPC):
    means = C_KM_S * np.asarray(redshifts, dtype=float) / float(truth)
    efficiency = ndtr((float(threshold) - means) / sigma) - ndtr((0.0 - means) / sigma)
    unnormalized = np.asarray(weights, dtype=float) * efficiency
    return unnormalized / unnormalized.sum()


def wilson(successes, n, z=1.959963984540054):
    p = successes / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0.0, float(center-half)), min(1.0, float(center+half))


def run(replicates=300, seed=20261031, grid_size=4001, bootstrap_draws=5000):
    if replicates < 1 or grid_size < 1001 or bootstrap_draws < 1000:
        raise ValueError("replicates must be positive, grid_size >= 1001, bootstrap_draws >= 1000")
    grid = np.linspace(*GRID_BOUNDS, grid_size)
    rows, summaries = [], []
    for truth in TRUTHS:
        for redshifts in HOST_PAIRS:
            for weights in WEIGHT_PAIRS:
                for threshold in THRESHOLDS:
                    rng = np.random.default_rng(cell_seed(seed, truth, redshifts, weights, threshold))
                    host_ids, distances = [], []
                    for _ in range(replicates):
                        for attempt in range(100000):
                            host = int(rng.choice(2, p=weights))
                            mean = C_KM_S * redshifts[host] / truth
                            distance = float(rng.normal(mean, SIGMA_MPC))
                            if 0.0 < distance < threshold:
                                host_ids.append(host); distances.append(distance); break
                        else:
                            raise RuntimeError("Detection rejection sampler exceeded limit")
                    z = np.asarray(redshifts, dtype=float)
                    w = np.asarray(weights, dtype=float)
                    zeff = float(np.dot(z, w))
                    cell_rows = []
                    for replicate, (host, distance) in enumerate(zip(host_ids, distances)):
                        posts = {
                            "two_host_uncorrected": posterior_h0(grid, z, distance, SIGMA_MPC, weights=w, prior_min=PRIOR[0], prior_max=PRIOR[1]),
                            "two_host_corrected": posterior_h0_selected(grid, z, distance, SIGMA_MPC, threshold, weights=w, prior_min=PRIOR[0], prior_max=PRIOR[1]),
                            "collapsed_uncorrected": posterior_h0(grid, np.array([zeff]), distance, SIGMA_MPC, weights=np.array([1.0]), prior_min=PRIOR[0], prior_max=PRIOR[1]),
                            "collapsed_corrected": posterior_h0_selected(grid, np.array([zeff]), distance, SIGMA_MPC, threshold, weights=np.array([1.0]), prior_min=PRIOR[0], prior_max=PRIOR[1]),
                        }
                        row = {"true_H0": truth, "host_z1": z[0], "host_z2": z[1], "weight1": w[0], "weight2": w[1],
                               "effective_redshift": zeff, "detection_max_mpc": threshold, "replicate": replicate,
                               "cell_seed": cell_seed(seed, truth, redshifts, weights, threshold),
                               "generated_host_index": host, "observed_distance_mpc": distance}
                        for method, posterior in posts.items():
                            s = posterior_summary(grid, posterior)
                            for key in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
                                row[f"{method}_{key}"] = float(s[key])
                            error = float(s["median"] - truth)
                            row[f"{method}_error"] = error
                            row[f"{method}_absolute_error"] = abs(error)
                        rows.append(row); cell_rows.append(row)
                    cell_id = {"true_H0": truth, "host_z1": z[0], "host_z2": z[1], "weight1": w[0], "weight2": w[1], "detection_max_mpc": threshold}
                    for method in METHODS:
                        errors = np.array([r[f"{method}_error"] for r in cell_rows])
                        hit68 = np.array([r[f"{method}_lower_68"] <= truth <= r[f"{method}_upper_68"] for r in cell_rows])
                        hit95 = np.array([r[f"{method}_lower_95"] <= truth <= r[f"{method}_upper_95"] for r in cell_rows])
                        lo68, hi68 = wilson(int(hit68.sum()), len(cell_rows)); lo95, hi95 = wilson(int(hit95.sum()), len(cell_rows))
                        summaries.append({**cell_id, "method": method, "n": len(cell_rows), "bias": float(errors.mean()),
                            "mae": float(np.abs(errors).mean()), "rmse": float(np.sqrt(np.mean(errors**2))),
                            "coverage_68": float(hit68.mean()), "coverage_68_wilson_low": lo68, "coverage_68_wilson_high": hi68,
                            "coverage_95": float(hit95.mean()), "coverage_95_wilson_low": lo95, "coverage_95_wilson_high": hi95,
                            "mean_width_68": float(np.mean([r[f"{method}_upper_68"]-r[f"{method}_lower_68"] for r in cell_rows])),
                            "mean_width_95": float(np.mean([r[f"{method}_upper_95"]-r[f"{method}_lower_95"] for r in cell_rows]))})
                    for correction in ("uncorrected", "corrected"):
                        collapsed, mixture = f"collapsed_{correction}", f"two_host_{correction}"
                        delta = np.array([r[f"{collapsed}_absolute_error"] - r[f"{mixture}_absolute_error"] for r in cell_rows])
                        brng = np.random.default_rng(cell_seed(seed + 1, truth, redshifts, weights, threshold) % (2**63-1))
                        boot = np.empty(bootstrap_draws)
                        for start in range(0, bootstrap_draws, 128):
                            stop = min(start+128, bootstrap_draws)
                            idx = brng.integers(0, len(delta), size=(stop-start, len(delta)))
                            boot[start:stop] = delta[idx].mean(axis=1)
                        summaries.append({**cell_id, "method": f"paired_collapsed_minus_two_host_absolute_error_{correction}",
                            "n": len(cell_rows), "mean_paired_change": float(delta.mean()),
                            "bootstrap_95_low": float(np.quantile(boot, .025)), "bootstrap_95_high": float(np.quantile(boot, .975)),
                            "bootstrap_draws": bootstrap_draws, "bootstrap_seed": cell_seed(seed + 1, truth, redshifts, weights, threshold)})
                    expected = expected_detected_host_fractions(truth, z, w, threshold)
                    observed = np.bincount(host_ids, minlength=2) / len(host_ids)
                    summaries.append({**cell_id, "method": "detected_host_fraction_diagnostic", "n": len(cell_rows),
                        "host1_observed_fraction": float(observed[0]), "host1_expected_fraction": float(expected[0]),
                        "host2_observed_fraction": float(observed[1]), "host2_expected_fraction": float(expected[1])})
    manifest = {"experiment":"host_mixture_threshold_sensitivity_factorial", "status":"synthetic diagnostic only",
        "seed":seed, "replicates_per_cell":replicates, "grid_size":grid_size, "grid_bounds":list(GRID_BOUNDS),
        "prior":list(PRIOR), "truths":list(TRUTHS), "host_redshift_pairs":[list(x) for x in HOST_PAIRS],
        "intrinsic_weight_pairs":[list(x) for x in WEIGHT_PAIRS], "upper_detection_thresholds_mpc":list(THRESHOLDS),
        "sigma_mpc":SIGMA_MPC, "lower_detection_threshold_mpc":0.0, "bootstrap_draws":bootstrap_draws,
        "cell_count":len(TRUTHS)*len(HOST_PAIRS)*len(WEIGHT_PAIRS)*len(THRESHOLDS),
        "raw_rows":len(rows), "summary_rows":len(summaries), "methods":list(METHODS),
        "python_version":platform.python_version(),"numpy_version":np.__version__,"scipy_version":scipy.__version__}
    return rows, summaries, manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=300); p.add_argument("--seed",type=int,default=20261031)
    p.add_argument("--grid-size",type=int,default=4001); p.add_argument("--bootstrap-draws",type=int,default=5000)
    p.add_argument("--output-dir",default="results/host-mixture-threshold-sensitivity")
    a=p.parse_args(); rows,summaries,manifest=run(a.replicates,a.seed,a.grid_size,a.bootstrap_draws)
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    for name,data in (("per_run.csv",rows),("summary.csv",summaries)):
        with (out/name).open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in data for k in row))); w.writeheader(); w.writerows(data)
    (out/"experiment_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(f"PASS: {len(rows)} raw rows, {len(summaries)} summaries, {manifest['cell_count']} cells")


if __name__=="__main__": main()
