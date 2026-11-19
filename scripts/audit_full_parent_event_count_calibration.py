"""Independent raw-seed reconstruction audit for the paired full-parent event-count study."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import logsumexp
from scipy.stats import norm, t

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299792.458


def _seed(master: int, replicate: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate)]).generate_state(1, dtype=np.uint32)[0])


def _posterior_quantiles_independent(z, host_weights, distances, h0, sigma, threshold):
    z = np.asarray(z, float)
    w = np.asarray(host_weights, float)
    d = np.asarray(distances, float)
    h = np.asarray(h0, float)
    if z.ndim != 1 or w.shape != z.shape or len(z) < 2:
        raise ValueError("invalid parent redshift/host-weight arrays")
    if not np.isfinite(z).all() or not np.isfinite(w).all() or (w <= 0).any():
        raise ValueError("parent data must be finite with positive weights")
    if d.ndim != 1 or not len(d) or not np.isfinite(d).all() or (d < threshold).any():
        raise ValueError("invalid detected event distances")
    logw = np.log(w)
    logw -= logsumexp(logw)
    mu = C_KM_S * z[None, :] / h[:, None]
    log_detection = logsumexp(
        norm.logcdf((mu-threshold)/sigma) + logw[None, :], axis=1
    )
    log_numerator = np.zeros(len(h), dtype=float)
    for distance in d:
        log_numerator += logsumexp(
            norm.logpdf(distance, loc=mu, scale=sigma) + logw[None, :], axis=1
        )
    log_post = log_numerator - len(d)*log_detection
    density = np.exp(log_post - np.max(log_post))
    area = float(trapezoid(density, h))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("invalid posterior normalization")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h)]
    cdf /= cdf[-1]
    return {name: float(np.interp(q, cdf, h)) for name, q in (
        ("q025", .025), ("q050", .05), ("q160", .16),
        ("q500", .50), ("q840", .84), ("q950", .95), ("q975", .975),
    )}


def _wilson(k: int, n: int, z: float = 1.959963984540054):
    if n <= 0:
        return float("nan"), float("nan")
    p = k/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half, center+half


def _assert_close(actual, expected, label, atol=2e-9):
    if not np.isclose(float(actual), float(expected), rtol=0, atol=atol, equal_nan=True):
        raise AssertionError(f"{label} mismatch: archived={actual}, reconstructed={expected}")


def audit(root: Path, expected_replicates: int | None = None, output: Path | None = None):
    root = root.resolve()
    manifest = json.loads((root/"manifest.json").read_text())
    bad = []
    for name, digest in manifest.items():
        path = root/name
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != digest:
            bad.append(name)
    if bad:
        raise AssertionError(f"input manifest mismatch: {bad[:10]}")

    design = json.loads((root/"design.json").read_text())
    raw = pd.read_csv(root/"replicates.csv")
    archived_summary = pd.read_csv(root/"summary.csv")
    archived_paired = pd.read_csv(root/"paired_contrast.csv")
    nrep = int(design["replicates"])
    if expected_replicates is not None and nrep != expected_replicates:
        raise AssertionError(f"expected {expected_replicates} replicates, found {nrep}")
    if design.get("status") != "COMPLETE" or design.get("event_counts") != [4, 30]:
        raise AssertionError("unexpected event-count design metadata")
    if len(raw) != 2*nrep or raw.duplicated(["replicate", "event_count"]).any():
        raise AssertionError("raw replicate/event-count factorial is incomplete or duplicated")
    if raw.groupby("event_count").size().to_dict() != {4: nrep, 30: nrep}:
        raise AssertionError("event-count arm denominators are incomplete")
    if not raw.groupby("replicate").event_count.nunique().eq(2).all():
        raise AssertionError("replicates do not contain both event-count arms")
    if set(raw.event_count.astype(int)) != {4, 30}:
        raise AssertionError("unexpected event-count values")

    prior_min, prior_max = map(float, design["uniform_h0_prior"])
    h0 = np.linspace(prior_min, prior_max, int(design["h0_grid_points"]))
    truth = float(design["h0_true"])
    n_parent = int(design["n_parent"])
    master = int(design["master_seed"])
    reconstructed = []
    for rep in range(nrep):
        block = raw[raw.replicate == rep].sort_values("event_count")
        expected_seed = _seed(master, rep)
        if block.simulation_seed.nunique() != 1 or int(block.simulation_seed.iloc[0]) != expected_seed:
            raise AssertionError(f"simulation seed mismatch for replicate {rep}")
        sim = generate_synthetic_population(
            seed=expected_seed, n_parent=n_parent, n_detected=30, h0_true=truth
        )
        parent = sim["parent_truth"]
        events = sim["gw_events_observed"]
        cfg = sim["config"]
        if len(parent) != n_parent or len(events) != 30:
            raise AssertionError(f"regenerated parent/event dimensions differ at replicate {rep}")
        z = parent.redshift.to_numpy(float)
        z_min, z_max = map(float, cfg["redshift_domain"])
        host_weights = 1.0 + 1.5*(z-z_min)/(z_max-z_min)
        distances = events.distance_mpc.to_numpy(float)
        if len(distances) != 30 or not np.isfinite(distances).all():
            raise AssertionError(f"invalid regenerated event vector at replicate {rep}")
        if not np.array_equal(distances[:4], distances[:4]):
            raise AssertionError("event prefix is not reproducible")
        for count in (4, 30):
            row = block[block.event_count == count].iloc[0]
            q = _posterior_quantiles_independent(
                z, host_weights, distances[:count], h0,
                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]),
            )
            _assert_close(row.event_prefix_distance_sum, distances[:count].sum(),
                          f"event distance prefix sum rep={rep} count={count}", atol=1e-8)
            for key, value in q.items():
                archived_key = "h0_median" if key == "q500" else key
                _assert_close(row[archived_key], value,
                              f"{archived_key} rep={rep} count={count}")
            bias = q["q500"]-truth
            _assert_close(row.median_bias, bias, f"median bias rep={rep} count={count}")
            for level, low_key, high_key in (
                (68, "q160", "q840"), (90, "q050", "q950"), (95, "q025", "q975")
            ):
                covered = int(q[low_key] <= truth <= q[high_key])
                _assert_close(row[f"covered{level}"], covered,
                              f"coverage{level} rep={rep} count={count}", atol=0)
            _assert_close(row.width68, q["q840"]-q["q160"], "width68", atol=1e-8)
            _assert_close(row.width90, q["q950"]-q["q050"], "width90", atol=1e-8)
            _assert_close(row.width95, q["q975"]-q["q025"], "width95", atol=1e-8)
            reconstructed.append({"replicate": rep, "event_count": count,
                                   "median_bias": bias, **q})
    rebuilt = pd.DataFrame(reconstructed)

    if len(archived_summary) != 2 or set(archived_summary.event_count.astype(int)) != {4, 30}:
        raise AssertionError("summary does not contain exactly two event-count arms")
    for count in (4, 30):
        g = raw[raw.event_count == count]
        s = archived_summary[archived_summary.event_count == count].iloc[0]
        bias = rebuilt[rebuilt.event_count == count].median_bias.to_numpy(float)
        _assert_close(s.planned, nrep, "summary planned")
        _assert_close(s.success, nrep, "summary success")
        _assert_close(s.mean_median_bias, bias.mean(), "summary mean bias")
        _assert_close(s.median_bias_mcse, bias.std(ddof=1)/math.sqrt(nrep), "summary bias MCSE")
        _assert_close(s.rmse, np.sqrt(np.mean(bias*bias)), "summary RMSE")
        for level in (68, 90, 95):
            k = int(g[f"covered{level}"].sum())
            lo, hi = _wilson(k, nrep)
            _assert_close(s[f"coverage{level}"], k/nrep, f"summary coverage{level}")
            _assert_close(s[f"coverage{level}_wilson95_low"], lo, f"Wilson low {level}")
            _assert_close(s[f"coverage{level}_wilson95_high"], hi, f"Wilson high {level}")
            col = f"width{level}"
            _assert_close(s[f"mean_{col}"], g[col].mean(), f"mean {col}")
            _assert_close(s[f"mcse_{col}"], g[col].std(ddof=1)/math.sqrt(nrep), f"MCSE {col}")

    wide = rebuilt.pivot(index="replicate", columns="event_count", values="median_bias")
    diff = (wide[30]-wide[4]).to_numpy(float)
    n = len(diff)
    mean = float(diff.mean())
    sd = float(diff.std(ddof=1))
    mcse = sd/math.sqrt(n)
    critical = float(t.ppf(.975, n-1))
    if len(archived_paired) != 1:
        raise AssertionError("paired contrast table must contain one record")
    p = archived_paired.iloc[0]
    _assert_close(p.paired_replicates, nrep, "paired denominator")
    _assert_close(p.mean_difference, mean, "paired mean difference")
    _assert_close(p.sd_difference, sd, "paired SD")
    _assert_close(p.mcse, mcse, "paired MCSE")
    _assert_close(p.ci95_low, mean-critical*mcse, "paired CI low")
    _assert_close(p.ci95_high, mean+critical*mcse, "paired CI high")

    report = {
        "status": "PASS",
        "manifest_entries_verified": len(manifest),
        "planned_replicates": nrep,
        "raw_rows_reconstructed": len(rebuilt),
        "parent_populations_regenerated": nrep,
        "detected_events_regenerated": 30*nrep,
        "posterior_curves_reconstructed": 2*nrep,
        "summary_arms_reconstructed": 2,
        "paired_contrasts_reconstructed": 1,
        "method": "Independent discrete-host likelihood implementation and deterministic regeneration from archived seeds.",
        "scope": "Synthetic fixed-truth full-parent event-count diagnostic; not joint SBC or observational H0 inference.",
    }
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--expected-replicates", type=int, default=500)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.artifact_dir, a.expected_replicates, a.output), indent=2))


if __name__ == "__main__":
    main()
