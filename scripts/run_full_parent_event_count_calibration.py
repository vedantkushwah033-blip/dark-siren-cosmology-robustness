"""Paired full-parent fixed-truth event-count calibration diagnostic.

For each independently generated parent population, compare the first four
detected events with a nested sample of thirty detected events. The same
full-parent host mixture and selected-event likelihood are used in both arms.
This isolates finite-event effects from catalogue thinning and inclusion fitting.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_h0_known_inclusion_oracle_sensitivity import selected_posterior_quantiles
from dark_siren_robustness.generative_selection import generate_synthetic_population


def _seed(master: int, replicate: int) -> int:
    return int(np.random.SeedSequence([master, replicate]).generate_state(1, dtype=np.uint32)[0])


def _wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1.0 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return center-half, center+half


def _paired_summary(values: np.ndarray) -> dict[str, float | int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    mean = float(x.mean()) if n else float("nan")
    sd = float(x.std(ddof=1)) if n > 1 else float("nan")
    mcse = sd/math.sqrt(n) if n > 1 else float("nan")
    crit = float(t.ppf(0.975, n-1)) if n > 1 else float("nan")
    return {
        "paired_replicates": n, "mean_difference": mean, "sd_difference": sd,
        "mcse": mcse,
        "ci95_low": mean-crit*mcse if n > 1 else float("nan"),
        "ci95_high": mean+crit*mcse if n > 1 else float("nan"),
    }


def run(replicates: int, master_seed: int, n_parent: int, h0_true: float,
        prior_min: float, prior_max: float, h0_points: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if replicates < 2 or n_parent < 100 or not prior_min < h0_true < prior_max or h0_points < 21:
        raise ValueError("invalid replicate, parent, prior, truth, or grid configuration")
    h0 = np.linspace(prior_min, prior_max, h0_points)
    rows = []
    for rep in range(replicates):
        sim_seed = _seed(master_seed, rep)
        sim = generate_synthetic_population(
            seed=sim_seed, n_parent=n_parent, n_detected=30, h0_true=h0_true
        )
        parent = sim["parent_truth"]
        events = sim["gw_events_observed"]
        cfg = sim["config"]
        z = parent.redshift.to_numpy(float)
        z_min, z_max = cfg["redshift_domain"]
        base = 1.0 + 1.5*(z-z_min)/(z_max-z_min)
        distances = events.distance_mpc.to_numpy(float)
        if len(distances) != 30:
            raise AssertionError("generator did not return the planned 30 detected events")
        for event_count in (4, 30):
            q = selected_posterior_quantiles(
                z, base, distances[:event_count], h0,
                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]),
            )
            rows.append({
                "replicate": rep, "simulation_seed": sim_seed,
                "n_parent": n_parent, "event_count": event_count,
                "h0_true": h0_true, "h0_median": q["q500"],
                "median_bias": q["q500"]-h0_true,
                "q025": q["q025"], "q050": q["q050"], "q160": q["q160"],
                "q840": q["q840"], "q950": q["q950"], "q975": q["q975"],
                "width68": q["q840"]-q["q160"],
                "width90": q["q950"]-q["q050"],
                "width95": q["q975"]-q["q025"],
                "covered68": int(q["q160"] <= h0_true <= q["q840"]),
                "covered90": int(q["q050"] <= h0_true <= q["q950"]),
                "covered95": int(q["q025"] <= h0_true <= q["q975"]),
                "event_prefix_distance_sum": float(distances[:event_count].sum()),
            })
    raw = pd.DataFrame(rows)
    summaries = []
    for event_count, group in raw.groupby("event_count"):
        n = len(group)
        bias = group.median_bias.to_numpy(float)
        row = {
            "event_count": int(event_count), "planned": replicates, "success": n,
            "mean_median_bias": float(bias.mean()),
            "median_bias_mcse": float(bias.std(ddof=1)/math.sqrt(n)),
            "rmse": float(np.sqrt(np.mean(bias*bias))),
        }
        for level in (68, 90, 95):
            k = int(group[f"covered{level}"].sum())
            lo, hi = _wilson(k, n)
            row[f"coverage{level}"] = k/n
            row[f"coverage{level}_wilson95_low"] = lo
            row[f"coverage{level}_wilson95_high"] = hi
        for level in (68, 90, 95):
            col = f"width{level}"
            row[f"mean_{col}"] = float(group[col].mean())
            row[f"mcse_{col}"] = float(group[col].std(ddof=1)/math.sqrt(n))
        summaries.append(row)
    wide = raw.pivot(index="replicate", columns="event_count", values="median_bias")
    paired = pd.DataFrame([{
        "contrast": "30-event minus nested 4-event posterior-median bias",
        **_paired_summary((wide[30]-wide[4]).to_numpy()),
    }])
    return raw, pd.DataFrame(summaries), paired


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates", type=int, default=100)
    p.add_argument("--seed", type=int, default=2026111928)
    p.add_argument("--n-parent", type=int, default=1200)
    p.add_argument("--h0-true", type=float, default=70.0)
    p.add_argument("--prior-min", type=float, default=20.0)
    p.add_argument("--prior-max", type=float, default=300.0)
    p.add_argument("--h0-points", type=int, default=561)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        p.error("output directory must be empty")
    raw, summary, paired = run(
        a.replicates, a.seed, a.n_parent, a.h0_true,
        a.prior_min, a.prior_max, a.h0_points,
    )
    raw.to_csv(out/"replicates.csv", index=False, float_format="%.17g")
    summary.to_csv(out/"summary.csv", index=False, float_format="%.17g")
    paired.to_csv(out/"paired_contrast.csv", index=False, float_format="%.17g")
    design = {
        "replicates": a.replicates, "master_seed": a.seed, "n_parent": a.n_parent,
        "event_counts": [4, 30], "event_design": "nested prefix of the same 30-event realization",
        "h0_true": a.h0_true, "uniform_h0_prior": [a.prior_min, a.prior_max],
        "h0_grid_points": a.h0_points, "status": "COMPLETE",
        "scope": "synthetic low-redshift full-parent selected-likelihood fixed-truth diagnostic only",
        "not_claimed": ["SBC", "catalogue-thinning effect", "observational H0 measurement", "LVK/GWTC analysis"],
    }
    (out/"design.json").write_text(json.dumps(design, indent=2, sort_keys=True)+"\n")
    manifest = {
        str(path.relative_to(out)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(out.rglob("*")) if path.is_file() and path.name != "manifest.json"
    }
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n")
    print(summary.to_string(index=False))
    print(paired.to_string(index=False))


if __name__ == "__main__":
    main()
