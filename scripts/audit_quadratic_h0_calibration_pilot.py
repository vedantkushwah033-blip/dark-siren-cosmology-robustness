"""Independent integrity and reconstruction audit for the quadratic H0 pilot.

The auditor intentionally does not import the pilot runner. It verifies the
outer and per-cell manifests, planned-cell ledger, posterior curves, derived
posterior metrics, and key aggregate summaries from serialized outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ARMS = {"degree_1": "density_d1_o7", "degree_2": "density_d2_o7", "catalogue_only": "density_catalogue_only", "oracle_true_inclusion": "density_oracle_true_inclusion"}
QUANTILES = (0.025, 0.05, 0.16, 0.50, 0.84, 0.95, 0.975)


def _fail(message: str) -> None:
    raise ValueError(message)


def _verify_manifest(root: Path, manifest_path: Path) -> int:
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict) or not manifest:
        _fail(f"empty or malformed manifest: {manifest_path}")
    for relative, expected in manifest.items():
        path = root / relative
        if not path.is_file():
            _fail(f"manifest file missing: {path}")
        data = path.read_bytes()
        if len(data) != int(expected["bytes"]):
            _fail(f"byte-count mismatch: {path}")
        digest = hashlib.sha256(data).hexdigest()
        if digest != expected["sha256"]:
            _fail(f"SHA-256 mismatch: {path}")
    actual = {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.name != "manifest.json"
    }
    if actual != set(manifest):
        _fail(f"manifest inventory mismatch under {root}")
    return len(manifest)


def _posterior_summary(grid: np.ndarray, density: np.ndarray, truth: float) -> dict:
    if (grid.ndim != 1 or density.shape != grid.shape or len(grid) < 2
            or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(density))
            or np.any(np.diff(grid) <= 0) or np.any(density < 0)):
        _fail("invalid posterior curve")
    area = float(np.trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        _fail("posterior curve has invalid area")
    d = density / area
    increments = 0.5 * (d[:-1] + d[1:]) * np.diff(grid)
    cdf = np.r_[0.0, np.cumsum(increments)]
    cdf /= cdf[-1]
    if not grid[0] < truth < grid[-1]:
        _fail("truth outside posterior support")
    q = {p: float(np.interp(p, cdf, grid)) for p in QUANTILES}
    return {
        "area": area,
        "pit": float(np.interp(truth, grid, cdf)),
        "h0_q025": q[0.025], "h0_q05": q[0.05],
        "h0_q16": q[0.16], "h0_median": q[0.50],
        "h0_q84": q[0.84], "h0_q95": q[0.95],
        "h0_q975": q[0.975],
        "covered_68": int(q[0.16] <= truth <= q[0.84]),
        "covered_90": int(q[0.05] <= truth <= q[0.95]),
        "covered_95": int(q[0.025] <= truth <= q[0.975]),
        "width_68": q[0.84] - q[0.16],
        "width_90": q[0.95] - q[0.05],
        "width_95": q[0.975] - q[0.025],
    }


def _check_close(label: str, actual: float, expected: float,
                 *, atol: float = 2e-10, rtol: float = 2e-9) -> None:
    if not np.isclose(float(actual), float(expected), atol=atol, rtol=rtol):
        _fail(f"{label} mismatch: stored={actual}, reconstructed={expected}")


def audit(root: Path) -> dict:
    root = root.resolve()
    config_path = root / "config.json"
    outer_manifest_path = root / "manifest.json"
    metrics_path = root / "replicate_metrics.csv"
    curves_path = root / "posterior_curves.csv"
    ledger_path = root / "failure_ledger.csv"
    aggregate_path = root / "aggregate_summary.json"
    for path in (config_path, outer_manifest_path, metrics_path,
                 curves_path, ledger_path, aggregate_path):
        if not path.is_file():
            _fail(f"required artifact missing: {path}")

    config = json.loads(config_path.read_text())
    if config.get("experiment") != "paired_quadratic_curvature_h0_calibration_pilot":
        _fail("unexpected experiment identifier")
    if config.get("status") != "synthetic_pilot_not_sbc_or_observational_measurement":
        _fail("missing synthetic-scope status")
    if config.get("arms") != ARMS:
        _fail("arm-to-posterior-column mapping mismatch")
    reps = int(config["replicates"])
    gammas = tuple(float(x) for x in config["gamma_values"])
    sizes = tuple(int(x) for x in config["validation_sizes"])
    planned_cells = [(r, g, n) for r in range(reps) for g in gammas for n in sizes]
    expected_rows = len(planned_cells) * len(ARMS)

    manifest_entries = _verify_manifest(root, outer_manifest_path)
    metrics = pd.read_csv(metrics_path)
    curves = pd.read_csv(curves_path)
    ledger = pd.read_csv(ledger_path)
    if len(metrics) != expected_rows or int(config["planned_metric_rows"]) != expected_rows:
        _fail("planned metric-row denominator mismatch")
    if int(config["actual_metric_rows"]) != len(metrics):
        _fail("config actual metric-row count mismatch")
    key = ["replicate", "gamma", "n_validation", "arm"]
    if metrics.duplicated(key).any():
        _fail("duplicate planned arm row")
    expected_keys = {
        (r, g, n, arm) for r, g, n in planned_cells for arm in ARMS
    }
    actual_keys = set(
        (int(x.replicate), float(x.gamma), int(x.n_validation), str(x.arm))
        for x in metrics[key].itertuples(index=False)
    )
    if actual_keys != expected_keys:
        _fail("planned arm-key inventory mismatch")
    if set(metrics.status) - {"SUCCESS", "FAILED"}:
        _fail("unexpected status value")
    counts = metrics.status.value_counts().to_dict()
    if counts != config.get("status_counts"):
        _fail("config status counts disagree with replicate metrics")

    ledger_keys = set(
        (int(x.replicate), float(x.gamma), int(x.n_validation), str(x.arm), str(x.status))
        for x in ledger.itertuples(index=False)
    )
    metric_ledger_keys = set(
        (int(x.replicate), float(x.gamma), int(x.n_validation), str(x.arm), str(x.status))
        for x in metrics[["replicate", "gamma", "n_validation", "arm", "status"]].itertuples(index=False)
    )
    if ledger_keys != metric_ledger_keys or len(ledger) != len(metrics):
        _fail("failure ledger does not exactly reconstruct metrics status keys")

    success = metrics[metrics.status == "SUCCESS"]
    failed = metrics[metrics.status == "FAILED"]
    curve_keys = ["replicate", "gamma", "n_validation", "arm"]
    if len(curves):
        if curves.duplicated(curve_keys + ["H0"]).any():
            _fail("duplicate posterior grid point")
    successful_keys = set(
        (int(x.replicate), float(x.gamma), int(x.n_validation), str(x.arm))
        for x in success[curve_keys].itertuples(index=False)
    )
    grouped_curves = {
        (int(k[0]), float(k[1]), int(k[2]), str(k[3])): group.sort_values("H0")
        for k, group in curves.groupby(curve_keys, sort=False)
    }
    if set(grouped_curves) != successful_keys:
        _fail("posterior curve groups do not match successful metric rows")
    if len(failed) and any(
        (int(x.replicate), float(x.gamma), int(x.n_validation), str(x.arm)) in grouped_curves
        for x in failed.itertuples(index=False)
    ):
        _fail("failed arm unexpectedly has posterior curve")

    for row in success.itertuples(index=False):
        k = (int(row.replicate), float(row.gamma), int(row.n_validation), str(row.arm))
        group = grouped_curves[k]
        if len(group) != int(config["h0_grid"][2]):
            _fail(f"wrong H0 grid length for {k}")
        reconstructed = _posterior_summary(
            group.H0.to_numpy(float), group.density.to_numpy(float), float(config["h0_true"])
        )
        _check_close(f"{k} posterior area", reconstructed["area"], float(row.area_before_normalization), atol=5e-8)
        for field in ("pit", "h0_q025", "h0_q05", "h0_q16", "h0_median",
                      "h0_q84", "h0_q95", "h0_q975", "width_68", "width_90", "width_95"):
            _check_close(f"{k} {field}", reconstructed[field], float(getattr(row, field)))
        for field in ("covered_68", "covered_90", "covered_95"):
            if int(getattr(row, field)) != reconstructed[field]:
                _fail(f"{k} {field} mismatch")

    # Reconstruct aggregate summaries independently from the stored replicate rows.
    stored_aggregate = json.loads(aggregate_path.read_text())
    aggregate_keys = {(float(x["gamma"]), int(x["n_validation"]), str(x["arm"])) for x in stored_aggregate}
    expected_aggregate_keys = {(g, n, arm) for g in gammas for n in sizes for arm in ARMS}
    if aggregate_keys != expected_aggregate_keys or len(stored_aggregate) != len(expected_aggregate_keys):
        _fail("aggregate summary group inventory mismatch")
    by_key = {(float(x["gamma"]), int(x["n_validation"]), str(x["arm"])): x for x in stored_aggregate}
    for (gamma, n, arm), group in metrics.groupby(["gamma", "n_validation", "arm"]):
        item = by_key[(float(gamma), int(n), str(arm))]
        valid = group[group.status == "SUCCESS"]
        if int(item["n_planned"]) != len(group) or int(item["n_success"]) != len(valid):
            _fail("aggregate planned/success counts mismatch")
        if int(item["n_failed"]) != int((group.status == "FAILED").sum()):
            _fail("aggregate failed count mismatch")
        _check_close("aggregate success fraction", item["success_fraction"], len(valid) / len(group))
        if len(valid):
            errors = valid.h0_median.to_numpy(float) - valid.h0_true.to_numpy(float)
            _check_close("aggregate bias", item["mean_posterior_median_bias"], np.mean(errors))
            _check_close("aggregate RMSE", item["rmse_posterior_median"], np.sqrt(np.mean(errors**2)))
            _check_close("aggregate mean PIT", item["mean_pit"], valid.pit.mean())
            for level in (68, 90, 95):
                coverage = float(valid[f"covered_{level}"].mean())
                _check_close(f"aggregate coverage {level}", item[f"coverage_{level}"], coverage)

    # Every raw cell must contain observable/truth separation and a valid cell manifest.
    raw_cells = 0
    for r, g, n in planned_cells:
        cell = root / "raw" / f"gamma_{g:+.2f}" / f"replicate_{r:04d}" / f"n_validation_{n:04d}"
        for relative in (
            "summary.json", "manifest.json", "observable/catalogue.csv",
            "observable/validation.csv", "observable/events.csv",
            "truth/parent_truth.csv", "truth/validation_frame_truth.csv",
            "truth/event_hosts.csv", "h0_posterior_curves.csv",
        ):
            if not (cell / relative).is_file():
                _fail(f"raw cell missing {relative}: {cell}")
        _verify_manifest(cell, cell / "manifest.json")
        summary = json.loads((cell / "summary.json").read_text())
        if summary.get("status") != "synthetic_smoke_not_measurement":
            _fail(f"raw cell has unexpected status: {cell}")
        if float(summary.get("gamma")) != float(g) or int(summary.get("n_validation")) != int(n):
            _fail(f"raw cell configuration does not match planned key: {cell}")
        observable_catalogue = pd.read_csv(cell / "observable/catalogue.csv")
        observable_validation = pd.read_csv(cell / "observable/validation.csv")
        observable_events = pd.read_csv(cell / "observable/events.csv")
        parent_truth = pd.read_csv(cell / "truth/parent_truth.csv")
        validation_truth = pd.read_csv(cell / "truth/validation_frame_truth.csv")
        event_truth = pd.read_csv(cell / "truth/event_hosts.csv")
        if list(observable_catalogue.columns) != ["object_id", "redshift"]:
            _fail(f"catalogue observable/truth leakage or schema drift: {cell}")
        if list(observable_validation.columns) != ["validation_id", "redshift", "retained"]:
            _fail(f"validation observable/truth leakage or schema drift: {cell}")
        if list(observable_events.columns) != ["event_id", "distance_mpc", "sigma_mpc"]:
            _fail(f"event observable/truth leakage or schema drift: {cell}")
        if not {"inclusion_probability", "retained"}.issubset(parent_truth.columns):
            _fail(f"parent truth table missing scoring fields: {cell}")
        if not {"inclusion_probability", "selected_for_validation"}.issubset(validation_truth.columns):
            _fail(f"validation truth table missing audit fields: {cell}")
        if not {"host_object_id", "host_redshift"}.issubset(event_truth.columns):
            _fail(f"event-host truth table missing scoring fields: {cell}")
        expected_catalogue = parent_truth.loc[
            parent_truth.retained.astype(bool), ["object_id", "redshift"]
        ].reset_index(drop=True)
        observed_catalogue = observable_catalogue.reset_index(drop=True)
        if not np.array_equal(expected_catalogue.object_id.to_numpy(), observed_catalogue.object_id.to_numpy()):
            _fail(f"catalogue rows do not match retained parent truth: {cell}")
        if not np.allclose(expected_catalogue.redshift, observed_catalogue.redshift, rtol=0, atol=0):
            _fail(f"catalogue redshifts do not match parent truth: {cell}")
        selected_truth = validation_truth.loc[
            validation_truth.selected_for_validation.astype(bool),
            ["frame_id", "redshift", "retained"],
        ].rename(columns={"frame_id": "validation_id"}).sort_values("validation_id").reset_index(drop=True)
        observed_validation = observable_validation.sort_values("validation_id").reset_index(drop=True)
        if not np.array_equal(selected_truth.validation_id, observed_validation.validation_id):
            _fail(f"validation IDs do not match selected frame: {cell}")
        if not np.allclose(selected_truth[["redshift", "retained"]], observed_validation[["redshift", "retained"]], rtol=0, atol=0):
            _fail(f"validation observables do not match selected truth frame: {cell}")
        if len(observable_events) != len(event_truth) or not np.array_equal(
            observable_events.event_id.to_numpy(), event_truth.event_id.to_numpy()
        ):
            _fail(f"event observable/truth IDs do not align: {cell}")
        raw_curve = pd.read_csv(cell / "h0_posterior_curves.csv")
        for arm, column in ARMS.items():
            metric = success[
                (success.replicate == r) & (success.gamma == g)
                & (success.n_validation == n) & (success.arm == arm)
            ]
            if len(metric) != 1:
                _fail(f"missing unique successful metric for raw cell arm {arm}: {cell}")
            for field in ("h0_q05", "h0_median", "h0_q95"):
                _check_close(f"raw summary {arm} {field}", summary["results"][arm][field], float(metric.iloc[0][field]))
            if column not in raw_curve:
                _fail(f"raw posterior curve missing {column}: {cell}")
            outer_group = grouped_curves[(int(r), float(g), int(n), str(arm))]
            raw_ordered = raw_curve.sort_values("H0")
            if not np.array_equal(raw_ordered.H0.to_numpy(float), outer_group.H0.to_numpy(float)):
                _fail(f"raw and aggregate H0 grids differ for {arm}: {cell}")
            if not np.allclose(raw_ordered[column].to_numpy(float), outer_group.density.to_numpy(float), rtol=0, atol=1e-14):
                _fail(f"raw and aggregate posterior curves differ for {arm}: {cell}")
        raw_cells += 1

    return {
        "status": "PASS",
        "planned_cells": len(planned_cells),
        "planned_arm_rows": expected_rows,
        "successful_arm_rows": int(len(success)),
        "failed_arm_rows": int(len(failed)),
        "posterior_curve_groups_reconstructed": len(grouped_curves),
        "raw_cell_manifests_verified": raw_cells,
        "outer_manifest_entries_verified": manifest_entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
