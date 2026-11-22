"""Independent aggregate auditor for the global multi-bundle cell lifecycle."""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _canonical_hash(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _timestamp(value, label):
    if value is None:
        return None
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise AssertionError(f"{label} lacks UTC offset")
    return parsed


def _verify_config_provenance(config: dict, provenance: dict) -> None:
    expected = {
        "seed": config["seed"],
        "event_seed": config.get("event_seed"),
        "validation_seed": config["validation_seed"],
        "validation_sample_seed": config.get("validation_sample_seed"),
        "n_parent": config["n_parent"],
        "n_detected": config["n_detected"],
        "h0_true": config["h0_true"],
        "sigma_mpc": config["sigma_mpc"],
        "detection_threshold_mpc": config["detection_threshold_mpc"],
        "n_validation_frame": config["validation_frame"],
        "curvature_gamma": config["curvature_gamma"],
        "curvature_beta": config["curvature_beta"],
        "inclusion_model_degree": config["inclusion_model_degree"],
        "quadratic_basis_center": config["quadratic_calibration"]["x2_center"],
    }
    for field, value in expected.items():
        if provenance.get(field) != value:
            raise AssertionError(f"bundle config/shared manifest provenance mismatch: {field}")


def audit(output_dir: Path) -> dict:
    output = output_dir.resolve()
    checkpoint = json.loads((output / "cell_lifecycle_checkpoint.json").read_text())
    plan, ledger = checkpoint.get("planned_cells"), checkpoint.get("cell_ledger")
    if checkpoint.get("schema_version") != 1 or not plan or not isinstance(ledger, list):
        raise AssertionError("invalid global checkpoint schema or empty plan")
    if checkpoint.get("plan_sha256") != _canonical_hash(plan):
        raise AssertionError("immutable global plan hash mismatch")
    fields = ("replicate", "gamma", "validation_size", "arm")
    key = lambda row: tuple(row[f] for f in fields)
    plan_map = {key(row): row for row in plan}
    if len(plan_map) != len(plan) or len({key(row) for row in ledger}) != len(ledger) or set(plan_map) != {key(row) for row in ledger}:
        raise AssertionError("global ledger keys do not exactly match immutable plan")
    counts = {"FAILED": 0, "NON_EVALUABLE": 0, "NOT_RUN": 0, "SUCCESS": 0}
    for row in ledger:
        planned = plan_map[key(row)]
        for field, value in planned.items():
            if row.get(field) != value:
                raise AssertionError(f"ledger provenance mismatch for {field}")
        status, audit_status = row["status"], row["audit_status"]
        if status not in counts:
            raise AssertionError(f"unknown lifecycle status {status}")
        if status == "NON_EVALUABLE":
            raise AssertionError("NON_EVALUABLE is prohibited until a preregistered policy exists")
        if status == "SUCCESS" and (audit_status != "PASS" or row.get("reason_code") != ""):
            raise AssertionError("SUCCESS cell lacks audit PASS or has a reason")
        if status == "FAILED" and audit_status == "PASS":
            raise AssertionError("FAILED cell cannot claim audit PASS")
        if status == "NOT_RUN" and audit_status != "NOT_RUN":
            raise AssertionError("NOT_RUN cell has inconsistent audit state")
        started, completed = _timestamp(row.get("started_at"), "started_at"), _timestamp(row.get("completed_at"), "completed_at")
        if status in {"SUCCESS", "FAILED"} and (started is None or completed is None):
            raise AssertionError("terminal cell lacks timestamps")
        if started and completed and completed < started:
            raise AssertionError("cell completed before it started")

        root = Path(row["input_dir"]).resolve()
        config_path, manifest_path = root / "config.json", root / "shared_input_manifest.json"
        config = json.loads(config_path.read_text())
        shared = load_shared_input_manifest(root, expected_validation_sizes=config["validation_sizes"])
        _verify_config_provenance(config, shared["provenance"])
        if _sha(config_path) != row["config_sha256"] or _sha(manifest_path) != row["shared_input_manifest_sha256"]:
            raise AssertionError("cell config/shared-input hash differs from immutable plan")
        if float(config["curvature_gamma"]) != float(row["gamma"]):
            raise AssertionError("bundle curvature differs from planned cell")
        if int(row["validation_size"]) not in shared["validation_sizes"]:
            raise AssertionError("cell validation size is absent from frozen bundle")
        if int(config["seed"]) != int(row["analysis_seed"]) or float(config["h0_true"]) != float(row["h0_true"]):
            raise AssertionError("bundle seed/H0 provenance differs from plan")
        if status == "SUCCESS":
            report_path = output / row["audit_report_path"]
            report = json.loads(report_path.read_text())
            if report.get("status") != "PASS" or tuple(report.get(f) for f in fields) != key(row):
                raise AssertionError("cell audit report identity/status mismatch")
        counts[status] += 1

    if (output / "cell_ledger.csv").read_text() != pd.DataFrame(ledger).to_csv(index=False):
        raise AssertionError("ledger CSV differs from checkpoint values")
    summary_columns = [
        "replicate", "gamma", "validation_size", "arm", "status", "h0_true",
        "median_error", "posterior_median", "q025", "q05", "q16", "q84", "q95", "q975",
        "weight_diagnostics", "realized_bandwidths", "input_sha256",
    ]
    summaries = []
    for row in ledger:
        if row["status"] == "SUCCESS":
            summary_path = (output / row["audit_report_path"]).parent / "result" / "summary.json"
            summary = json.loads(summary_path.read_text())
            if tuple(summary.get(f) for f in fields) != key(row) or summary.get("status") != "SUCCESS":
                raise AssertionError("cell summary identity/status mismatch")
            summaries.append(summary)
    expected = pd.DataFrame(summaries, columns=summary_columns).to_csv(index=False)
    if (output / "cell_summaries.csv").read_text() != expected:
        raise AssertionError("cell summaries CSV differs from source summaries")

    summary = json.loads((output / "run_summary.json").read_text())
    expected_status = "PASS" if counts["SUCCESS"] == len(plan) else "INCOMPLETE"
    if (summary.get("status") != expected_status or summary.get("status_counts") != counts
        or summary.get("planned_cells") != len(plan)
        or summary.get("replicate_count") != len({r["replicate"] for r in plan})
        or summary.get("curvature_count") != len({r["gamma"] for r in plan})
        or summary.get("bundle_count") != len({(r["replicate"], r["gamma"]) for r in plan})):
        raise AssertionError("global run summary does not reconcile with plan/checkpoint")

    manifest = json.loads((output / "run_manifest.json").read_text())
    actual = {str(p.relative_to(output)) for p in output.rglob("*")
              if p.is_file() and p.name != "run_manifest.json" and not p.name.endswith(".lock")}
    if manifest.get("status") != expected_status or set(manifest.get("files", {})) != actual:
        raise AssertionError("global run manifest file set/status mismatch")
    for relative, meta in manifest["files"].items():
        path = output / relative
        if path.stat().st_size != int(meta["bytes"]) or _sha(path) != meta["sha256"]:
            raise AssertionError(f"global run manifest hash mismatch: {relative}")

    cell_auditor = Path(__file__).resolve().with_name("audit_arm_cell.py")
    repeated = 0
    for row in ledger:
        if row["status"] != "SUCCESS":
            continue
        cell_dir = (output / row["audit_report_path"]).parent / "result"
        cmd = [sys.executable, str(cell_auditor), "--input-dir", row["input_dir"], "--cell-dir", str(cell_dir)]
        check = subprocess.run(cmd, capture_output=True, text=True)
        if check.returncode:
            raise AssertionError("repeat cell audit failed: " + (check.stderr[-4000:] or check.stdout[-4000:]))
        if json.loads(check.stdout).get("status") != "PASS":
            raise AssertionError("repeat cell audit did not PASS")
        repeated += 1
    return {"status": "PASS", "planned_cells": len(plan), "status_counts": counts,
            "successful_cells_reaudited": repeated,
            "claim_limit": "Global lifecycle and per-cell independent audit; not calibration."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
