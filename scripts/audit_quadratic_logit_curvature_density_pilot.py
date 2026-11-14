"""Independent outer auditor for the paired multi-gamma curvature density pilot."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_child_auditor():
    path = HERE / "audit_quadratic_logit_curvature_density_smoke.py"
    spec = importlib.util.spec_from_file_location("curvature_child_auditor_outer", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.audit


def _compare_aggregate(root: Path, filename: str, frames: list[pd.DataFrame], keys: list[str]) -> int:
    expected = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    actual = pd.read_csv(root / "results" / filename)
    if len(actual) != len(expected):
        raise AssertionError(f"{filename}: aggregate row count mismatch")
    if expected.empty:
        if not actual.empty:
            raise AssertionError(f"{filename}: unexpected aggregate rows")
        return 0
    if list(actual.columns) != list(expected.columns):
        raise AssertionError(f"{filename}: aggregate columns mismatch")
    if expected.duplicated(keys).any() or actual.duplicated(keys).any():
        raise AssertionError(f"{filename}: duplicate aggregate keys")
    expected = expected.sort_values(keys).reset_index(drop=True)
    actual = actual.sort_values(keys).reset_index(drop=True)
    for column in expected.columns:
        left, right = actual[column], expected[column]
        if pd.api.types.is_numeric_dtype(right):
            if not np.allclose(left.to_numpy(float), right.to_numpy(float),
                               rtol=2e-10, atol=2e-10, equal_nan=True):
                raise AssertionError(f"{filename}: aggregate numerical mismatch: {column}")
        elif left.fillna("").astype(str).tolist() != right.fillna("").astype(str).tolist():
            raise AssertionError(f"{filename}: aggregate value mismatch: {column}")
    return len(expected)


def audit(root: Path) -> dict[str, object]:
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    checked = 0
    for relative, metadata in manifest["files"].items():
        path = root / relative
        if (not path.is_file()
                or path.stat().st_size != int(metadata["bytes"])
                or _hash(path) != metadata["sha256"]):
            raise AssertionError(f"outer manifest verification failed: {relative}")
        checked += 1
    if checked != int(manifest["file_count"]):
        raise AssertionError("outer manifest file count mismatch")

    cfg = json.loads((root / "config.json").read_text())
    if cfg.get("experiment") != "quadratic_logit_curvature_paired_density_pilot":
        raise AssertionError("unexpected experiment identifier")
    if cfg.get("status") != "ENGINEERING_PILOT_NOT_CONFIRMATORY":
        raise AssertionError("unexpected experiment claim status")
    gammas = [float(value) for value in cfg["gammas"]]
    replicates = int(cfg["replicates_per_gamma"])
    if len(gammas) != 7 or len(set(gammas)) != len(gammas):
        raise AssertionError("expected seven unique gamma values")
    if replicates < 1:
        raise AssertionError("invalid replicate count")

    ledger = pd.read_csv(root / "results" / "replicate_ledger.csv").fillna("")
    expected_keys = {(gamma, rep) for gamma in gammas for rep in range(replicates)}
    actual_keys = {(float(row.gamma), int(row.replicate)) for row in ledger.itertuples()}
    if actual_keys != expected_keys or ledger.duplicated(["gamma", "replicate"]).any():
        raise AssertionError("ledger gamma-replicate keys mismatch")
    if set(ledger.status) - {"SUCCESS", "FAILED"}:
        raise AssertionError("unknown lifecycle status")

    child_audit = _load_child_auditor()
    frames_by_file = {
        "density_metrics.csv": [],
        "weight_diagnostics.csv": [],
        "scott_widths.csv": [],
    }
    audited = 0
    for row in ledger.itertuples():
        gamma, replicate = float(row.gamma), int(row.replicate)
        expected_pop_seed = int(cfg["base_seed"]) + 100000 * replicate + 17
        expected_val_seed = int(cfg["base_seed"]) + 100000 * replicate + 29
        if int(row.population_seed) != expected_pop_seed or int(row.validation_seed) != expected_val_seed:
            raise AssertionError("paired seed schedule mismatch")
        label = f"g{gamma:+.1f}_r{replicate:03d}".replace("+", "p").replace("-", "m").replace(".", "")
        expected_child = f"children/{label}"
        if row.child_path != expected_child:
            raise AssertionError("ledger child path mismatch")
        child = root / expected_child

        if row.status == "FAILED":
            if not str(row.error).strip():
                raise AssertionError("failed child missing error record")
            continue

        if not child.is_dir():
            raise AssertionError(f"successful child missing: {expected_child}")
        child_cfg = json.loads((child / "config.json").read_text())
        if (float(child_cfg["gamma"]) != gamma
                or int(child_cfg["population_seed"]) != expected_pop_seed
                or int(child_cfg["validation_seed"]) != expected_val_seed
                or int(child_cfg["n_parent"]) != int(cfg["n_parent"])
                or list(child_cfg["validation_sizes"]) != list(cfg["validation_sizes"])
                or list(child_cfg["absolute_bandwidth_candidates"]) != list(cfg["absolute_bandwidth_candidates"])):
            raise AssertionError(f"child configuration mismatch: {expected_child}")
        child_report = child_audit(child)
        if child_report.get("status") != "PASS":
            raise AssertionError(f"child raw-input audit did not pass: {expected_child}")
        audited += 1

        for filename, frames in frames_by_file.items():
            frame = pd.read_csv(child / "results" / filename)
            frame.insert(0, "replicate", replicate)
            frame.insert(0, "gamma", gamma)
            frame["population_seed"] = expected_pop_seed
            frame["validation_seed"] = expected_val_seed
            frames.append(frame)

    metric_rows = _compare_aggregate(
        root, "density_metrics.csv", frames_by_file["density_metrics.csv"],
        ["gamma", "replicate", "validation_size", "absolute_bandwidth_sd", "arm", "node"],
    )
    diagnostic_rows = _compare_aggregate(
        root, "weight_diagnostics.csv", frames_by_file["weight_diagnostics.csv"],
        ["gamma", "replicate", "validation_size", "arm", "node"],
    )
    scott_rows = _compare_aggregate(
        root, "scott_widths.csv", frames_by_file["scott_widths.csv"],
        ["gamma", "replicate", "validation_size", "arm", "node"],
    )
    failed = int((ledger.status == "FAILED").sum())
    return {
        "status": "PASS" if failed == 0 else "INCOMPLETE",
        "outer_manifest_files_verified": checked,
        "children_audited": audited,
        "children_expected": len(ledger),
        "failed_children": failed,
        "aggregate_metric_rows_reconstructed": metric_rows,
        "aggregate_weight_diagnostic_rows_reconstructed": diagnostic_rows,
        "aggregate_scott_width_rows_reconstructed": scott_rows,
        "scope": "Independent child raw-input reconstruction plus exact outer aggregate-value, key, seed, lifecycle, and manifest audit; synthetic density stage only.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.artifact_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
