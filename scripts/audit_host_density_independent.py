"""Independent re-derivation audit for the toy host-density benchmark.

Re-simulates selected replicate seeds, computes a reflected weighted Gaussian KDE
and selected H0 posterior with separate numerical code, and compares posterior
mean/median/interval endpoints using tolerances fixed before the audit.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.special import ndtr
from dark_siren_robustness.generative_selection import generate_synthetic_population

C = 299_792.458
TOL_MEAN, TOL_Q = 0.05, 0.10


def inclusion(z, floor=0.10, mid=0.22, width=0.05):
    return floor + (1-floor)/(1+np.exp((z-mid)/width))


def naive_kde(zg, sample, weights, factor, lo, hi):
    total = float(weights.sum())
    mu = float(np.sum(weights*sample)/total)
    sd = float(np.sqrt(np.sum(weights*(sample-mu)**2)/total))
    neff = total**2/float(np.sum(weights**2))
    h = sd*(neff**(-1/5) if factor is None else factor)
    if not np.isfinite(h) or h <= 0:
        raise ValueError("independent KDE bandwidth must be finite and positive")
    out = np.zeros_like(zg)
    for value, weight in zip(sample, weights):
        for center in (value, 2*lo-value, 2*hi-value):
            out += weight*np.exp(-0.5*((zg-center)/h)**2)
    area = float(simpson(out, x=zg))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("independent KDE normalization failed")
    return out/area


def posterior(zg, density, distances, h0g, sigma, threshold):
    density = density/simpson(density, x=zg)
    loglike = np.empty(h0g.size)
    for i, H in enumerate(h0g):
        mu = C*zg/H
        per_event = [
            np.log(simpson(density*np.exp(-0.5*((distance-mu)/sigma)**2) /
                           (sigma*np.sqrt(2*np.pi)), x=zg))
            for distance in distances
        ]
        beta = simpson(density*ndtr((mu-threshold)/sigma), x=zg)
        loglike[i] = sum(per_event)-len(distances)*np.log(beta)
    p = np.exp(loglike-loglike.max())
    p /= simpson(p, x=h0g)
    cdf = np.r_[0.0, np.cumsum(0.5*(p[1:]+p[:-1])*np.diff(h0g))]
    cdf /= cdf[-1]
    return (float(simpson(h0g*p, x=h0g)),
            *(float(np.interp(q, cdf, h0g)) for q in (0.50, 0.16, 0.84)))


def _verify_manifest(outdir):
    path = outdir/"manifest.json"
    if not path.exists():
        return {"status": "NOT_PRESENT"}
    manifest = json.loads(path.read_text(encoding="utf-8"))
    records = manifest.get("files", {})
    for name, record in records.items():
        file_path = outdir/name
        if not file_path.is_file():
            raise ValueError(f"manifest-listed file missing: {name}")
        actual = hashlib.sha256(file_path.read_bytes()).hexdigest()
        if actual != record.get("sha256"):
            raise ValueError(f"SHA-256 mismatch for {name}")
    return {"status": "PASS", "files_checked": len(records)}


def audit(outdir, n_reps=3, z_points=1601, h0_points=801):
    outdir = Path(outdir)
    meta = json.loads((outdir/"summary.json").read_text(encoding="utf-8"))
    manifest = _verify_manifest(outdir)
    rows = list(csv.DictReader((outdir/"replicate_metrics.csv").open(newline="", encoding="utf-8")))
    if not rows:
        raise ValueError("replicate_metrics.csv is empty")
    required = {"replicate", "seed", "arm", "posterior_mean", "posterior_median", "ci68_low", "ci68_high"}
    if not required.issubset(rows[0]):
        raise ValueError("replicate CSV is missing required posterior fields")
    factor = meta.get("bandwidth_factor")
    factor = None if factor is None or factor == "scott" else float(factor)
    midpoint = float(meta.get("supplied_inclusion_midpoint", 0.22))
    width = float(meta.get("supplied_inclusion_width", 0.05))
    zlo, zhi = meta["redshift_grid"][:2]
    h0lo, h0hi = meta["h0_grid"][:2]
    zg, h0g = np.linspace(zlo, zhi, z_points), np.linspace(h0lo, h0hi, h0_points)
    merger = 1+1.5*(zg-zg[0])/(zg[-1]-zg[0])
    row_map = {}
    for row in rows:
        key = (int(row["replicate"]), row["arm"])
        if key in row_map:
            raise ValueError(f"duplicate replicate-arm row: {key}")
        row_map[key] = row
    failures, checked = [], 0
    replicate_ids = sorted({int(r["replicate"]) for r in rows})[:n_reps]
    for rep in replicate_ids:
        sim = generate_synthetic_population(
            seed=int(meta["seed_start"])+rep, n_parent=int(meta["n_parent"]),
            n_detected=int(meta["n_detected"]), h0_true=float(meta["h0_true_km_s_mpc"]),
            sigma_mpc=float(meta["sigma_mpc"]),
            detection_threshold_mpc=float(meta["detection_threshold_mpc"]),
            z_min=float(zlo), z_max=float(zhi))
        parent = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        observed = sim["catalogue_observed"]["redshift"].to_numpy(dtype=float)
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        p_include = inclusion(observed, mid=midpoint, width=width)
        arms = {
            "oracle": (parent, np.ones(parent.size)),
            "catalogue_only": (observed, np.ones(observed.size)),
            "known_completeness": (observed, 1.0/p_include),
        }
        for arm, (sample, weights) in arms.items():
            density = naive_kde(zg, sample, weights, factor, zlo, zhi)*merger
            density /= simpson(density, x=zg)
            values = posterior(zg, density, distances, h0g,
                               float(meta["sigma_mpc"]), float(meta["detection_threshold_mpc"]))
            if (rep, arm) not in row_map:
                failures.append({"replicate": rep, "arm": arm, "field": "row", "error": "missing arm"})
                continue
            row = row_map[(rep, arm)]
            checked += 1
            for key, value, tolerance in zip(
                ("posterior_mean", "posterior_median", "ci68_low", "ci68_high"),
                values, (TOL_MEAN, TOL_Q, TOL_Q, TOL_Q)):
                recorded = float(row[key])
                if abs(recorded-value) > tolerance:
                    failures.append({"replicate": rep, "arm": arm, "field": key,
                                     "recorded": recorded, "independent": value, "tolerance": tolerance})
    return {"status": "PASS" if not failures else "FAIL", "rows_reconstructed": checked,
            "replicates_checked": replicate_ids, "failures": failures, "manifest": manifest,
            "tolerances": {"mean": TOL_MEAN, "quantile": TOL_Q}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--n-reps", type=int, default=3)
    parser.add_argument("--z-points", type=int, default=1601)
    parser.add_argument("--h0-points", type=int, default=801)
    args = parser.parse_args()
    result = audit(args.outdir, args.n_reps, args.z_points, args.h0_points)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
