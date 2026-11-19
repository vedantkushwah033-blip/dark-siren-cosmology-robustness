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


def assert_cross_gamma_pairing(parent_reference, parent_candidate,
                               events_reference, events_candidate,
                               validation_reference, validation_candidate):
    """Check common parent/event and validation ID-redshift streams across gamma.

    Inclusion probabilities and validation retention labels are intentionally
    allowed to differ across curvature conditions.
    """
    parent_columns = ["object_id", "redshift"]
    validation_columns = ["validation_id", "redshift"]
    if not parent_reference[parent_columns].equals(parent_candidate[parent_columns]):
        raise AssertionError("parent IDs/redshifts differ across gamma")
    if not events_reference.equals(events_candidate):
        raise AssertionError("observed GW events differ across gamma")
    if not validation_reference[validation_columns].equals(
        validation_candidate[validation_columns]
    ):
        raise AssertionError("validation IDs/redshifts differ across gamma")


def assert_complete_planned_ledger(planned: pd.DataFrame, ledger: pd.DataFrame) -> None:
    """Require the result ledger to contain every preregistered cell exactly once.

    A ledger containing only successful rows is not sufficient evidence of a
    complete experiment: a missing replicate could otherwise disappear from
    every denominator and paired contrast.
    """
    keys = ["cell_id", "replicate", "gamma", "n_validation",
            "analysis_seed", "validation_seed", "h0_true"]
    for label, frame in (("planned", planned), ("result", ledger)):
        missing = set(keys) - set(frame.columns)
        if missing:
            raise AssertionError(f"{label} ledger missing design columns: {sorted(missing)}")
        if frame.duplicated(["cell_id"]).any() or frame.duplicated(
            ["replicate", "gamma", "n_validation"]
        ).any():
            raise AssertionError(f"{label} ledger contains duplicate planned-cell keys")
    left = planned[keys].sort_values(keys[:4]).reset_index(drop=True)
    right = ledger[keys].sort_values(keys[:4]).reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
    except AssertionError as exc:
        raise AssertionError(
            "result ledger does not exactly match the preregistered planned-cell "
            "inventory, seeds, and H0 truth"
        ) from exc


def terminal_status_counts(ledger: pd.DataFrame) -> dict[str, int]:
    """Count final ledger states and reject unknown or interrupted states."""
    allowed = {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"}
    observed = set(ledger.status.astype(str))
    unknown = observed - allowed
    if unknown:
        raise AssertionError(f"unknown or non-terminal ledger statuses: {sorted(unknown)}")
    return {status: int((ledger.status == status).sum())
            for status in sorted(allowed)}


def failure_inclusive_coverage_bounds(covered_success: int, successful: int,
                                      planned: int) -> tuple[float, float]:
    """Worst/best coverage bounds when every non-success is uncovered/covered."""
    if planned <= 0 or successful < 0 or successful > planned:
        raise ValueError("invalid planned/successful denominator")
    if covered_success < 0 or covered_success > successful:
        raise ValueError("covered successes must be between zero and successes")
    non_success = planned - successful
    return covered_success / planned, (covered_success + non_success) / planned


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
    planned_path = root / "planned_cells.csv"
    if not planned_path.is_file():
        raise AssertionError("preregistered planned_cells.csv is missing")
    planned = pd.read_csv(planned_path)
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    if ledger.empty:
        raise AssertionError("result ledger is empty")
    assert_complete_planned_ledger(planned, ledger)
    status_counts = terminal_status_counts(ledger)
    successful_ledger = ledger[ledger.status == "SUCCESS"].copy()
    if successful_ledger.empty:
        raise AssertionError("no successful cells are available for oracle reconstruction")
    # Each replicate/curvature pair must retain at least one successful cell so
    # the common parent, catalogue and detected-event realization can be audited.
    planned_pairs = planned[["replicate", "gamma"]].drop_duplicates()
    successful_pairs = successful_ledger[["replicate", "gamma"]].drop_duplicates()
    if len(planned_pairs.merge(successful_pairs, how="left", indicator=True).query("_merge == 'left_only'")):
        raise AssertionError("a replicate/curvature pair has no successful cell for oracle reconstruction")

    rows = []
    paired_inputs = {}
    for (rep, gamma), block in successful_ledger.groupby(["replicate", "gamma"]):
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
        validation = pd.read_csv(directory / "observable/validation_observed.csv")
        paired_inputs[(int(rep), float(gamma))] = (parent, events, validation)
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
                "estimated_bias": float(cell.h0_median-cell.h0_true),
                "known_pi_median": known["q500"], "known_pi_bias": known["q500"]-cell.h0_true,
                "known_pi_coverage68": int(known["q160"] <= cell.h0_true <= known["q840"]),
                "full_parent_median": full["q500"], "full_parent_bias": full["q500"]-cell.h0_true,
                "full_parent_coverage68": int(full["q160"] <= cell.h0_true <= full["q840"]),
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
            larger_validation = pd.read_csv(d/"observable/validation_observed.csv")
            prefix_columns = ["validation_id", "redshift", "retained"]
            if not validation[prefix_columns].equals(
                larger_validation.iloc[:len(validation)][prefix_columns].reset_index(drop=True)
            ):
                raise AssertionError("validation sample is not a nested prefix across sizes")

    gammas = sorted(float(x) for x in ledger.gamma.unique())
    for rep in sorted(int(x) for x in ledger.replicate.unique()):
        reference_parent, reference_events, reference_validation = paired_inputs[(rep, gammas[0])]
        for gamma in gammas[1:]:
            parent_candidate, events_candidate, validation_candidate = paired_inputs[(rep, gamma)]
            assert_cross_gamma_pairing(
                reference_parent, parent_candidate,
                reference_events, events_candidate,
                reference_validation, validation_candidate,
            )

    cells = pd.DataFrame(rows)
    summaries = []
    contrasts = []
    for (gamma, nval), group in cells.groupby(["gamma", "n_validation"]):
        planned_group = planned[(planned.gamma == gamma) & (planned.n_validation == nval)]
        ledger_group = ledger[(ledger.gamma == gamma) & (ledger.n_validation == nval)]
        successful_group = successful_ledger[
            (successful_ledger.gamma == gamma) & (successful_ledger.n_validation == nval)
        ]
        covered_success = int(successful_group.covered68.fillna(0).sum())
        coverage_lower, coverage_upper = failure_inclusive_coverage_bounds(
            covered_success, len(successful_group), len(planned_group)
        )
        row = {"gamma": gamma, "n_validation": nval, "planned": len(planned_group),
               "success": len(successful_group),
               "non_evaluable": int((ledger_group.status == "NON_EVALUABLE").sum()),
               "failed": int((ledger_group.status == "FAILED").sum()),
               "not_run": int((ledger_group.status == "NOT_RUN").sum()),
               "replicates": len(group),
               "h0_grid_points": int(group.h0_grid_points.iloc[0])}
        for col in ("estimated_bias", "known_pi_bias", "full_parent_bias",
                    "estimated_minus_known_pi", "known_pi_minus_full_parent"):
            row[col] = float(group[col].mean())
            row[col+"_mcse"] = float(group[col].std(ddof=1)/math.sqrt(len(group)))
        row["estimated_coverage68_success_only"] = float(successful_group.covered68.mean())
        row["estimated_coverage68_lower_all_attempts"] = coverage_lower
        row["estimated_coverage68_upper_all_attempts"] = coverage_upper
        row["known_pi_coverage68_success_only"] = float(group.known_pi_coverage68.mean())
        row["full_parent_coverage68_success_only"] = float(group.full_parent_coverage68.mean())
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
    ledger = pd.read_csv(args.artifact_root / "replicate_ledger.csv")
    (args.output_dir/"oracle_grid_config.json").write_text(json.dumps({
        "h0_grid_points": int(args.h0_points),
        "prior_bounds": [20.0, 300.0],
        "planned_cells": int(len(pd.read_csv(args.artifact_root / "planned_cells.csv"))),
        "status_counts": terminal_status_counts(ledger),
        "interpretation": "independent same-grid oracle posterior reconstruction; successful cells only; failure-inclusive estimated-coverage bounds reported separately; synthetic diagnostic, not observational inference"
    }, indent=2, sort_keys=True))
    print(summary.to_string(index=False))
    print(contrasts.to_string(index=False))


if __name__ == "__main__":
    main()
