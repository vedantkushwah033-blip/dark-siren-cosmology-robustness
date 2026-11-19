"""Independent full-parent / known-inclusion oracle audit for the 100-replicate H0 pilot.

The oracle arms are simulation-only diagnostics. The audit retains all 600
planned cells, including NON_EVALUABLE estimated-inclusion fits, while
reconstructing oracle posteriors from their archived parent/catalogue/event
inputs. It does not modify the inference runner.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from audit_h0_known_inclusion_oracle_sensitivity import (
    _paired, assert_complete_planned_ledger, assert_cross_gamma_pairing,
    make_h0_grid, selected_posterior_quantiles,
)


def coverage_bounds(success_covered: int, success: int, planned: int) -> dict:
    if planned < 1 or success < 0 or success > planned or success_covered < 0 or success_covered > success:
        raise ValueError("invalid coverage counts")
    unresolved = planned - success
    return {
        "success_only": success_covered / success if success else float("nan"),
        "all_attempt_lower": success_covered / planned,
        "all_attempt_upper": (success_covered + unresolved) / planned,
        "success": success,
        "planned": planned,
        "unresolved": unresolved,
    }


def wilson_interval(k: int, n: int, z: float = 1.959963984540054):
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return center-half, center+half


def assert_oracle_eligible_statuses(ledger: pd.DataFrame) -> None:
    allowed = {"SUCCESS", "NON_EVALUABLE"}
    observed = set(ledger.status.astype(str))
    if not observed.issubset(allowed):
        raise AssertionError(f"unsupported statuses for oracle audit: {ledger.status.value_counts().to_dict()}")


def _verify_manifest(root: Path) -> int:
    manifest = json.loads((root / "manifest.json").read_text())
    bad = []
    for name, expected in manifest.items():
        path = root / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != expected:
            bad.append(name)
    if bad:
        raise AssertionError(f"aggregate manifest mismatch: {bad[:20]}")
    return len(manifest)


def audit(root: Path, h0_points: int = 1121):
    root = root.resolve()
    manifest_count = _verify_manifest(root)
    planned = pd.read_csv(root / "planned_cells.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    assert_complete_planned_ledger(planned, ledger)
    assert_oracle_eligible_statuses(ledger)
    if len(planned) != 600:
        raise AssertionError(f"expected 600 planned cells, found {len(planned)}")

    cell_rows = []
    paired_inputs = {}
    checked_children = 0
    nested_blocks = 0
    validation_pair_blocks = 0
    for (rep, gamma), block in ledger.groupby(["replicate", "gamma"], sort=True):
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
        paired_inputs[(int(rep), float(gamma))] = (
            parent[["object_id", "redshift"]].copy(), events.copy(),
            validation[["validation_id", "redshift"]].copy(),
        )
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
            catalogue.redshift, base(catalogue.redshift)/pi,
            distances, h0, sigma, threshold
        )
        full = selected_posterior_quantiles(
            parent.redshift, base(parent.redshift),
            distances, h0, sigma, threshold
        )

        previous_ids = None
        reference_catalogue = catalogue.copy()
        reference_parent = parent.copy()
        reference_events = events.copy()
        for cell in block.itertuples(index=False):
            d = root / "replicates" / cell.cell_id
            report = json.loads((d / "audit_report.json").read_text())
            if report.get("status") != "PASS":
                raise AssertionError(f"child audit not PASS: {cell.cell_id}")
            child_cfg = json.loads((d / "config.json").read_text())
            if not np.isclose(child_cfg["inference"]["curvature_gamma"], float(gamma)):
                raise AssertionError(f"gamma mismatch: {cell.cell_id}")
            if int(child_cfg["validation_generator_config"]["n_validation"]) != int(cell.n_validation):
                raise AssertionError(f"validation size mismatch: {cell.cell_id}")
            if int(child_cfg["analysis_generator_config"]["seed"]) != int(cell.analysis_seed):
                raise AssertionError(f"analysis seed mismatch: {cell.cell_id}")
            if int(child_cfg["validation_generator_config"]["seed"]) != int(cell.validation_seed):
                raise AssertionError(f"validation seed mismatch: {cell.cell_id}")
            if not np.isclose(
                float(child_cfg["analysis_generator_config"]["h0_true_km_s_mpc"]),
                float(cell.h0_true), rtol=0, atol=0,
            ):
                raise AssertionError(f"H0 truth mismatch: {cell.cell_id}")
            candidate_catalogue = pd.read_csv(d / "observable/analysis_catalogue.csv")
            candidate_parent = pd.read_csv(d / "truth/analysis_parent_truth.csv")
            candidate_events = pd.read_csv(d / "observable/gw_events_observed.csv")
            pd.testing.assert_frame_equal(candidate_catalogue, reference_catalogue)
            pd.testing.assert_frame_equal(candidate_parent, reference_parent)
            pd.testing.assert_frame_equal(candidate_events, reference_events)
            observed = pd.read_csv(d / "observable/validation_observed.csv")
            ids = observed.validation_id.to_numpy(dtype=int)
            if len(ids) != int(cell.n_validation):
                raise AssertionError(f"validation sample count mismatch: {cell.cell_id}")
            if previous_ids is not None and not np.array_equal(previous_ids, ids[:len(previous_ids)]):
                raise AssertionError(f"validation samples not nested: {cell.cell_id}")
            previous_ids = ids
            checked_children += 1
            estimated_median = float(cell.h0_median) if cell.status == "SUCCESS" else float("nan")
            cell_rows.append({
                "replicate": int(rep), "gamma": float(gamma),
                "n_validation": int(cell.n_validation), "cell_id": cell.cell_id,
                "status": cell.status, "h0_grid_points_oracle": int(h0.size),
                "estimated_median": estimated_median,
                "estimated_bias": estimated_median-70.0 if np.isfinite(estimated_median) else float("nan"),
                "covered68": float(cell.covered68) if cell.status == "SUCCESS" else float("nan"),
                "covered90": float(cell.covered90) if cell.status == "SUCCESS" else float("nan"),
                "covered95": float(cell.covered95) if cell.status == "SUCCESS" else float("nan"),
                "known_pi_median": known["q500"], "known_pi_bias": known["q500"]-70.0,
                "known_pi_coverage68": int(known["q160"] <= 70.0 <= known["q840"]),
                "known_pi_coverage90": int(known["q050"] <= 70.0 <= known["q950"]),
                "known_pi_coverage95": int(known["q025"] <= 70.0 <= known["q975"]),
                "full_parent_median": full["q500"], "full_parent_bias": full["q500"]-70.0,
                "full_parent_coverage68": int(full["q160"] <= 70.0 <= full["q840"]),
                "full_parent_coverage90": int(full["q050"] <= 70.0 <= full["q950"]),
                "full_parent_coverage95": int(full["q025"] <= 70.0 <= full["q975"]),
                "estimated_minus_known_pi": estimated_median-known["q500"] if np.isfinite(estimated_median) else float("nan"),
                "known_pi_minus_full_parent": known["q500"]-full["q500"],
            })
        if len(block) > 1:
            nested_blocks += 1

    gammas = sorted(float(x) for x in ledger.gamma.unique())
    for rep in sorted(int(x) for x in ledger.replicate.unique()):
        reference = paired_inputs[(rep, gammas[0])]
        for gamma in gammas[1:]:
            candidate = paired_inputs[(rep, gamma)]
            assert_cross_gamma_pairing(reference[0], candidate[0],
                                       reference[1], candidate[1],
                                       reference[2], candidate[2])
            validation_pair_blocks += 1

    cells = pd.DataFrame(cell_rows)
    summaries = []
    contrasts = []
    for (gamma, nval), group in cells.groupby(["gamma", "n_validation"], sort=True):
        success = group[group.status == "SUCCESS"]
        planned_n = len(group)
        row = {
            "gamma": gamma, "n_validation": nval, "planned": planned_n,
            "success": len(success),
            "non_evaluable": int((group.status == "NON_EVALUABLE").sum()),
            "failed": int((group.status == "FAILED").sum()),
            "estimated_median_bias_mean": float(success.estimated_bias.mean()) if len(success) else float("nan"),
            "estimated_median_bias_mcse": float(success.estimated_bias.std(ddof=1)/math.sqrt(len(success))) if len(success)>1 else float("nan"),
            "estimated_median_bias_rmse": float(np.sqrt(np.mean(success.estimated_bias.to_numpy(float)**2))) if len(success) else float("nan"),
            "known_pi_median_bias_mean": float(group.known_pi_bias.mean()),
            "known_pi_median_bias_mcse": float(group.known_pi_bias.std(ddof=1)/math.sqrt(planned_n)),
            "full_parent_median_bias_mean": float(group.full_parent_bias.mean()),
            "full_parent_median_bias_mcse": float(group.full_parent_bias.std(ddof=1)/math.sqrt(planned_n)),
        }
        for arm in ("estimated", "known_pi", "full_parent"):
            for level in (68, 90, 95):
                col = f"covered{level}" if arm == "estimated" else f"{arm}_coverage{level}"
                if arm == "estimated":
                    k = int(success[col].sum())
                    lo, hi = wilson_interval(k, len(success))
                    bounds = coverage_bounds(k, len(success), planned_n)
                    row.update({
                        f"{arm}_coverage{level}_success_only": bounds["success_only"],
                        f"{arm}_coverage{level}_wilson95_low": lo,
                        f"{arm}_coverage{level}_wilson95_high": hi,
                        f"{arm}_coverage{level}_all_attempt_low": bounds["all_attempt_lower"],
                        f"{arm}_coverage{level}_all_attempt_high": bounds["all_attempt_upper"],
                    })
                else:
                    k = int(group[col].sum())
                    lo, hi = wilson_interval(k, planned_n)
                    row[f"{arm}_coverage{level}"] = k/planned_n
                    row[f"{arm}_coverage{level}_wilson95_low"] = lo
                    row[f"{arm}_coverage{level}_wilson95_high"] = hi
        for arm, col in (
            ("estimated_minus_known_pi", "estimated_minus_known_pi"),
            ("known_pi_minus_full_parent", "known_pi_minus_full_parent"),
        ):
            vals = group[col].to_numpy(float)
            contrasts.append({"gamma": gamma, "n_validation": nval,
                              "contrast": arm, **_paired(vals)})
        summaries.append(row)

    # Replicate-paired curvature contrasts at each validation size.
    for nval, group in cells.groupby("n_validation"):
        wide = group.pivot(index="replicate", columns="gamma", values="estimated_bias")
        if {-1.0, 0.0, 1.0}.issubset(wide.columns):
            for left, right in ((-1.0, 0.0), (0.0, 1.0), (-1.0, 1.0)):
                d = wide[right]-wide[left]
                contrasts.append({"gamma": float("nan"), "n_validation": int(nval),
                                  "contrast": f"estimated gamma {right:g} minus {left:g}",
                                  **_paired(d.to_numpy(float))})
    report = {
        "status": "PASS",
        "planned_cells": len(planned),
        "ledger_status_counts": ledger.status.value_counts().to_dict(),
        "non_success_cells": ledger.loc[ledger.status != "SUCCESS", [
            "cell_id", "replicate", "gamma", "n_validation", "status",
            "quadrature_status", "orders_evaluated", "converged_at_order", "error",
        ]].replace({np.nan: None}).to_dict(orient="records"),
        "child_audits_verified": checked_children,
        "nested_validation_blocks_checked": nested_blocks,
        "cross_gamma_validation_pairing_blocks_checked": validation_pair_blocks,
        "aggregate_manifest_entries_verified": manifest_count,
        "oracle_grid_points": h0_points,
        "paired_replicate_blocks": int(cells.replicate.nunique()),
        "replicate_gamma_blocks": int(cells[["replicate", "gamma"]].drop_duplicates().shape[0]),
        "interpretation": "Oracle decomposition is synthetic diagnostic evidence, not a practical estimator or observational H0 result."
    }
    return cells, pd.DataFrame(summaries), pd.DataFrame(contrasts), report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_root", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--h0-points", type=int, default=1121)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    cells, summary, contrasts, report = audit(a.artifact_root, a.h0_points)
    cells.to_csv(a.output_dir / "oracle_cell_readout.csv", index=False)
    summary.to_csv(a.output_dir / "oracle_summary.csv", index=False)
    contrasts.to_csv(a.output_dir / "oracle_paired_contrasts.csv", index=False)
    (a.output_dir / "oracle_audit_report.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(summary.to_string(index=False))
    print(contrasts.to_string(index=False))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
