"""Audit the effect of H0 prior truncation on saved five-arm likelihoods.

For each saved child likelihood curve, compute posterior summaries on the full
expanded support and on a restricted [50,90] support using the same likelihood.
This is a within-realization prior-support sensitivity diagnostic.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid

ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]


def posterior_on_support(
    grid: np.ndarray, log_likelihood: np.ndarray, lower: float, upper: float
) -> dict[str, float]:
    grid = np.asarray(grid, dtype=float)
    log_likelihood = np.asarray(log_likelihood, dtype=float)
    if grid.ndim != 1 or log_likelihood.shape != grid.shape or len(grid) < 3:
        raise ValueError("grid and log-likelihood must be matching one-dimensional arrays")
    if np.any(np.diff(grid) <= 0) or not np.all(np.isfinite(log_likelihood)):
        raise ValueError("grid must increase and log-likelihood must be finite")
    if lower < grid[0] or upper > grid[-1] or lower >= upper:
        raise ValueError("requested prior support is outside the saved grid")
    mask = (grid >= lower - 1e-12) & (grid <= upper + 1e-12)
    h = grid[mask]
    logl = log_likelihood[mask]
    if len(h) < 3 or not np.isclose(h[0], lower, atol=1e-10) or not np.isclose(h[-1], upper, atol=1e-10):
        raise ValueError("saved grid must contain both prior endpoints and at least three nodes")
    density = np.exp(logl - np.max(logl))
    area = float(trapezoid(density, h))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h)]
    cdf /= cdf[-1]
    return {
        "median": float(np.interp(0.50, cdf, h)),
        "q025": float(np.interp(0.025, cdf, h)),
        "q16": float(np.interp(0.16, cdf, h)),
        "q84": float(np.interp(0.84, cdf, h)),
        "q975": float(np.interp(0.975, cdf, h)),
        "support_mass": float(trapezoid(density, h)),
    }


def audit(root: Path, restricted_support=(50.0, 90.0)) -> dict:
    root = Path(root)
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    orchestration = json.loads(
        (root / "independent_orchestration_audit.json").read_text(encoding="utf-8")
    )
    if orchestration.get("status") != "PASS":
        raise AssertionError("independent orchestration audit must PASS first")
    if config.get("h0_prior_support") != [30.0, 110.0]:
        raise AssertionError("expected the frozen expanded [30,110] source support")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    if len(ledger) != 900 or set(ledger.status) != {"SUCCESS"} or set(ledger.audit_status) != {"PASS"}:
        raise AssertionError("source lifecycle ledger is incomplete or unaudited")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    if len(summaries) != len(ledger):
        raise AssertionError("posterior summaries do not match planned lifecycle cells")

    full_lower, full_upper = map(float, config["h0_prior_support"])
    restricted_lower, restricted_upper = map(float, restricted_support)
    rows = []
    for (replicate, gamma), group in ledger.groupby(["replicate", "gamma"], sort=True):
        child = root / "replicates" / (
            f"rep{int(replicate):03d}_g{float(gamma):+g}".replace("+", "p").replace("-", "m")
        )
        child_config = json.loads((child / "config.json").read_text(encoding="utf-8"))
        if child_config.get("h0_prior_support") != [full_lower, full_upper]:
            raise AssertionError(f"child prior support mismatch: replicate={replicate}, gamma={gamma}")
        with np.load(child / "results" / "log_likelihoods.npz") as archive:
            grid = np.asarray(archive["H0_grid"], dtype=float)
            if not np.isclose(grid[0], full_lower) or not np.isclose(grid[-1], full_upper):
                raise AssertionError("saved likelihood grid does not span declared support")
            for nval in config["validation_sizes"]:
                for arm in ARMS:
                    key = f"{arm}_n{int(nval)}"
                    if key not in archive:
                        raise AssertionError(f"missing saved likelihood curve: {key}")
                    logl = np.asarray(archive[key], dtype=float)
                    full = posterior_on_support(grid, logl, full_lower, full_upper)
                    restricted = posterior_on_support(grid, logl, restricted_lower, restricted_upper)
                    saved = summaries.loc[
                        (summaries.replicate == replicate)
                        & (summaries.gamma == gamma)
                        & (summaries.validation_size == int(nval))
                        & (summaries.arm == arm)
                    ]
                    if len(saved) != 1:
                        raise AssertionError("expected exactly one archived posterior summary")
                    saved_median = float(saved.iloc[0].posterior_median)
                    if not np.isclose(saved_median, full["median"], rtol=0, atol=2e-7):
                        raise AssertionError("reconstructed expanded-support median differs from source")
                    full_grid = (grid >= full_lower - 1e-12) & (grid <= full_upper + 1e-12)
                    below = grid[full_grid] <= restricted_lower + 1e-12
                    full_density = np.exp(logl - np.max(logl))
                    full_density /= trapezoid(full_density, grid)
                    mass_below = float(trapezoid(full_density[below], grid[below]))
                    rows.append({
                        "replicate": int(replicate),
                        "gamma": float(gamma),
                        "validation_size": int(nval),
                        "arm": arm,
                        "full_median": full["median"],
                        "restricted_median": restricted["median"],
                        "median_shift_full_minus_restricted": full["median"] - restricted["median"],
                        "full_q025": full["q025"],
                        "restricted_q025": restricted["q025"],
                        "q025_shift_full_minus_restricted": full["q025"] - restricted["q025"],
                        "full_mass_below_50": mass_below,
                    })
    frame = pd.DataFrame(rows)
    if len(frame) != 900:
        raise AssertionError(f"expected 900 support-sensitivity rows, found {len(frame)}")
    aggregate = frame.groupby(["gamma", "validation_size", "arm"], sort=True).agg(
        n=("replicate", "size"),
        mean_full_median=("full_median", "mean"),
        mean_restricted_median=("restricted_median", "mean"),
        mean_median_shift=("median_shift_full_minus_restricted", "mean"),
        mean_full_q025=("full_q025", "mean"),
        mean_restricted_q025=("restricted_q025", "mean"),
        mean_q025_shift=("q025_shift_full_minus_restricted", "mean"),
        mean_mass_below_50=("full_mass_below_50", "mean"),
        fraction_mass_below_50_gt_0p10=("full_mass_below_50", lambda x: float((x > 0.10).mean())),
    ).reset_index()
    frame.to_csv(root / "prior_support_sensitivity_cells.csv", index=False)
    aggregate.to_csv(root / "prior_support_sensitivity_summary.csv", index=False)
    result = {
        "status": "PASS",
        "source_prior_support": [full_lower, full_upper],
        "restricted_prior_support": [restricted_lower, restricted_upper],
        "planned_cells_reconstructed": len(frame),
        "aggregate_rows": len(aggregate),
        "interpretation": (
            "Within-realization prior-support sensitivity only. Same saved likelihoods are "
            "renormalized on expanded and restricted supports; not a new data realization "
            "or a calibration/estimator-superiority result."
        ),
    }
    (root / "prior_support_sensitivity_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--restricted-min", type=float, default=50.0)
    parser.add_argument("--restricted-max", type=float, default=90.0)
    args = parser.parse_args()
    audit(args.artifact_dir, (args.restricted_min, args.restricted_max))


if __name__ == "__main__":
    main()
