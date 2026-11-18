"""Independent oracle-arm audit of a fixed-truth H0 calibration artifact.

Computes selected-H0 posteriors using true inclusion probabilities on the
observed catalogue and the complete finite parent population. Both are
simulation-only diagnostic oracles, not practical estimators.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm, t

C_KM_S = 299792.458


def selected_posterior_quantiles(z, raw_weights, distances, h0, sigma, threshold):
    """Independent discrete-host selected-event likelihood and H0 quantiles."""
    z = np.asarray(z, dtype=float)
    w = np.asarray(raw_weights, dtype=float)
    distances = np.asarray(distances, dtype=float)
    h0 = np.asarray(h0, dtype=float)
    if z.ndim != 1 or z.size < 2 or w.shape != z.shape:
        raise ValueError("z and positive host weights must be matching 1D arrays")
    if not np.isfinite(z).all() or not np.isfinite(w).all() or (w <= 0).any():
        raise ValueError("host inputs must be finite and weights strictly positive")
    if distances.ndim != 1 or not distances.size or not np.isfinite(distances).all():
        raise ValueError("distances must be a non-empty finite vector")
    if (distances < threshold).any() or sigma <= 0 or threshold <= 0:
        raise ValueError("invalid detected-event distances or selection parameters")
    if h0.ndim != 1 or h0.size < 2 or (np.diff(h0) <= 0).any():
        raise ValueError("H0 grid must be strictly increasing")
    logw = np.log(w)
    logw -= logsumexp(logw)
    means = C_KM_S * z[None, :] / h0[:, None]
    log_beta = logsumexp(
        norm.logcdf((means - threshold) / sigma) + logw[None, :], axis=1
    )
    log_numerator = np.zeros(h0.size)
    for distance in distances:
        log_numerator += logsumexp(
            norm.logpdf(distance, loc=means, scale=sigma) + logw[None, :], axis=1
        )
    log_likelihood = log_numerator - distances.size * log_beta
    density = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(trapezoid(density, h0))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("invalid H0 posterior integral")
    density /= area
    cdf = np.r_[0.0, np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(h0))]
    cdf /= cdf[-1]
    return {f"q{int(1000*q):03d}": float(np.interp(q, cdf, h0))
            for q in (0.025, 0.05, 0.16, 0.50, 0.84, 0.95, 0.975)}


def _paired(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    mean = float(x.mean()) if n else float("nan")
    sd = float(x.std(ddof=1)) if n > 1 else float("nan")
    se = sd / math.sqrt(n) if n > 1 else float("nan")
    crit = float(t.ppf(0.975, n - 1)) if n > 1 else float("nan")
    return {"n": n, "mean": mean, "sd": sd, "mcse": se,
            "ci95_low": mean - crit*se if n > 1 else float("nan"),
            "ci95_high": mean + crit*se if n > 1 else float("nan")}


def make_h0_grid(prior_min: float, prior_max: float, configured_points: int, requested_points: int | None = None):
    points = int(configured_points if requested_points is None else requested_points)
    if points < 2 or not np.isfinite([prior_min, prior_max]).all() or prior_min >= prior_max:
        raise ValueError("invalid H0 grid specification")
    return np.linspace(float(prior_min), float(prior_max), points)


def audit(root: Path, h0_points: int | None = None):
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    bad = []
    for name, expected in manifest.items():
        path = root / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != expected:
            bad.append(name)
    if bad:
        raise AssertionError(f"aggregate manifest mismatch: {bad[:10]}")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    if ledger.empty or set(ledger.status) != {"SUCCESS"}:
        raise AssertionError("input pilot must contain only successful cells")

    rows = []
    for (rep, gamma), block in ledger.groupby(["replicate", "gamma"]):
        block = block.sort_values("n_validation")
        first = block.iloc[0]
        directory = root / "replicates" / first.cell_id
        cfg = json.loads((directory / "config.json").read_text())
        h0_min, h0_max = cfg["inference"]["uniform_h0_prior"]
        h0 = make_h0_grid(h0_min, h0_max, cfg["inference"]["h0_grid_points"], h0_points)
        z_min, z_max = cfg["analysis_generator_config"]["redshift_domain"]
        sigma = float(cfg["analysis_generator_config"]["sigma_mpc"])
        threshold = float(cfg["analysis_generator_config"]["detection_threshold_mpc"])
        base = lambda z: 1.0 + 1.5*(np.asarray(z)-z_min)/(z_max-z_min)

        parent = pd.read_csv(directory / "truth/analysis_parent_truth.csv")
        catalogue = pd.read_csv(directory / "observable/analysis_catalogue.csv")
        events = pd.read_csv(directory / "observable/gw_events_observed.csv")
        if parent.object_id.duplicated().any() or catalogue.object_id.duplicated().any():
            raise AssertionError("duplicate object IDs")
        indexed = parent.set_index("object_id")
        if not catalogue.object_id.isin(indexed.index).all():
            raise AssertionError("catalogue IDs absent from parent")
        matched = indexed.loc[catalogue.object_id]
        if not np.allclose(matched.redshift, catalogue.redshift, rtol=0, atol=1e-14):
            raise AssertionError("catalogue redshifts disagree with parent truth")
        pi = matched.inclusion_probability.to_numpy(float)
        if not np.isfinite(pi).all() or (pi <= 0).any() or (pi > 1).any():
            raise AssertionError("invalid true inclusion probabilities")
        distances = events.distance_mpc.to_numpy(float)
        known = selected_posterior_quantiles(
            catalogue.redshift, base(catalogue.redshift)/pi, distances, h0, sigma, threshold
        )
        full = selected_posterior_quantiles(
            parent.redshift, base(parent.redshift), distances, h0, sigma, threshold
        )
        for cell in block.itertuples(index=False):
            rows.append({
                "replicate": int(rep), "gamma": float(gamma),
                "h0_grid_points": int(h0.size),
                "n_validation": int(cell.n_validation), "cell_id": cell.cell_id,
                "estimated_median": float(cell.h0_median),
                "estimated_bias": float(cell.h0_median-70.0),
                "known_pi_median": known["q500"], "known_pi_bias": known["q500"]-70.0,
                "known_pi_coverage68": int(known["q160"] <= 70.0 <= known["q840"]),
                "full_parent_median": full["q500"], "full_parent_bias": full["q500"]-70.0,
                "full_parent_coverage68": int(full["q160"] <= 70.0 <= full["q840"]),
                "estimated_minus_known_pi": float(cell.h0_median-known["q500"]),
                "known_pi_minus_full_parent": float(known["q500"]-full["q500"]),
            })
        # Catalogue, parent and event data are invariant across validation sizes.
        for cell in block.iloc[1:].itertuples(index=False):
            d = root / "replicates" / cell.cell_id
            pd.testing.assert_frame_equal(
                pd.read_csv(d/"observable/analysis_catalogue.csv"), catalogue)
            pd.testing.assert_frame_equal(
                pd.read_csv(d/"observable/gw_events_observed.csv"), events)
            pd.testing.assert_frame_equal(
                pd.read_csv(d/"truth/analysis_parent_truth.csv"), parent)

    cells = pd.DataFrame(rows)
    summaries = []
    contrasts = []
    for (gamma, nval), group in cells.groupby(["gamma", "n_validation"]):
        row = {"gamma": gamma, "n_validation": nval, "replicates": len(group),
               "h0_grid_points": int(group.h0_grid_points.iloc[0])}
        for col in ("estimated_bias", "known_pi_bias", "full_parent_bias",
                    "estimated_minus_known_pi", "known_pi_minus_full_parent"):
            row[col] = float(group[col].mean())
            row[col+"_mcse"] = float(group[col].std(ddof=1)/math.sqrt(len(group)))
        row["estimated_coverage68"] = float(
            ledger[(ledger.gamma==gamma)&(ledger.n_validation==nval)].covered68.mean())
        row["known_pi_coverage68"] = float(group.known_pi_coverage68.mean())
        row["full_parent_coverage68"] = float(group.full_parent_coverage68.mean())
        summaries.append(row)
        for label, col in (
            ("estimated minus known-pi IPW median", "estimated_minus_known_pi"),
            ("known-pi IPW minus full-parent median", "known_pi_minus_full_parent"),
        ):
            contrasts.append({"gamma": gamma, "n_validation": nval,
                              "contrast": label, **_paired(group[col].to_numpy())})
    return cells, pd.DataFrame(summaries), pd.DataFrame(contrasts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--h0-points", type=int, default=1121,
                        help="H0 posterior grid size for the independent oracle reconstruction (default: 1121)")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cells, summary, contrasts = audit(args.artifact_root, h0_points=args.h0_points)
    cells.to_csv(args.output_dir/"oracle_cell_readout.csv", index=False)
    summary.to_csv(args.output_dir/"oracle_summary.csv", index=False)
    contrasts.to_csv(args.output_dir/"oracle_paired_contrasts.csv", index=False)
    (args.output_dir/"oracle_grid_config.json").write_text(json.dumps({
        "h0_grid_points": int(args.h0_points),
        "prior_bounds": [20.0, 300.0],
        "interpretation": "independent oracle posterior reconstruction on a refined grid; not an observational result"
    }, indent=2, sort_keys=True))
    print(summary.to_string(index=False))
    print(contrasts.to_string(index=False))


if __name__ == "__main__":
    main()
