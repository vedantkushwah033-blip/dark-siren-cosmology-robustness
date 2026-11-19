from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "pool_h0_100rep_calibration_replicates.py"
SPEC = importlib.util.spec_from_file_location("pool_h0_replicates", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def write_replication(root: Path, label: str, seed_base: int, non_success=False):
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for rep in range(100):
        for gamma in (-1.0, 0.0, 1.0):
            for nval in (100, 500):
                status = "NON_EVALUABLE" if non_success and rep == 0 and gamma == -1.0 and nval == 100 else "SUCCESS"
                bias = 0.2 + 0.1*gamma + 0.001*rep + (0.02 if nval == 500 else 0.0)
                rows.append({
                    "cell_id": f"{label}_r{rep}_g{gamma}_n{nval}",
                    "replicate": rep, "gamma": gamma, "n_validation": nval,
                    "analysis_seed": seed_base + rep,
                    "validation_seed": seed_base + 10000 + rep,
                    "h0_true": 70.0, "status": status,
                    "h0_median": 70.0+bias if status == "SUCCESS" else np.nan,
                    "median_bias": bias if status == "SUCCESS" else np.nan,
                    "h0_q16": 65.0+bias if status == "SUCCESS" else np.nan,
                    "h0_q84": 75.0+bias if status == "SUCCESS" else np.nan,
                    "covered68": int(rep % 3 != 0) if status == "SUCCESS" else np.nan,
                    "covered90": 1 if status == "SUCCESS" else np.nan,
                    "covered95": 1 if status == "SUCCESS" else np.nan,
                    "min_catalogue_kish_ess": 30.0+rep if status == "SUCCESS" else np.nan,
                    "max_normalized_weight": 0.05 if status == "SUCCESS" else np.nan,
                })
    ledger = pd.DataFrame(rows)
    planned = ledger[MODULE.KEYS].copy()
    planned.to_csv(root/"planned_cells.csv", index=False)
    ledger.to_csv(root/"replicate_ledger.csv", index=False)
    audit = {
        "status": "PASS", "planned_cells": 600,
        "status_counts": {str(k): int(v) for k,v in ledger.status.value_counts().items()},
    }
    (root/"independent_orchestration_audit.json").write_text(json.dumps(audit))
    manifest = {
        name: hashlib.sha256((root/name).read_bytes()).hexdigest()
        for name in ("planned_cells.csv", "replicate_ledger.csv")
    }
    (root/"manifest.json").write_text(json.dumps(manifest))
    return ledger


def test_pools_two_independent_600_cell_replications_and_keeps_non_success(tmp_path):
    a = tmp_path/"run1"; b = tmp_path/"run2"
    write_replication(a, "r1", 100000, non_success=True)
    write_replication(b, "r2", 200000)
    report = MODULE.analyze(a, b, tmp_path/"out", expected_cells=600)
    assert report["status"] == "PASS"
    assert report["pooled_planned_cells"] == 1200
    assert report["pooled_success"] == 1199
    assert report["pooled_non_success"] == 1
    pooled = pd.read_csv(tmp_path/"out/pooled_cell_summary.csv")
    assert len(pooled) == 6
    cell = pooled[(pooled.gamma == -1) & (pooled.n_validation == 100)].iloc[0]
    assert cell.planned == 200 and cell.success == 199
    assert cell.coverage68_all_attempt_upper-cell.coverage68_all_attempt_lower == pytest.approx(1/200)
    per_rep = pd.read_csv(tmp_path/"out/replication_cell_summary.csv")
    assert len(per_rep) == 12
    contrasts = pd.read_csv(tmp_path/"out/pooled_paired_contrasts.csv")
    assert len(contrasts) == 15


def test_rejects_overlapping_seed_streams(tmp_path):
    a = tmp_path/"run1"; b = tmp_path/"run2"
    write_replication(a, "r1", 100000)
    write_replication(b, "r2", 100000)
    with pytest.raises(AssertionError, match="seed streams overlap"):
        MODULE.analyze(a, b, tmp_path/"out", expected_cells=600)


def test_rejects_plan_ledger_mismatch_after_manifest_refresh(tmp_path):
    a = tmp_path/"run1"; b = tmp_path/"run2"
    write_replication(a, "r1", 100000)
    write_replication(b, "r2", 200000)
    ledger = pd.read_csv(a/"replicate_ledger.csv").iloc[:-1]
    ledger.to_csv(a/"replicate_ledger.csv", index=False)
    manifest = {
        name: hashlib.sha256((a/name).read_bytes()).hexdigest()
        for name in ("planned_cells.csv", "replicate_ledger.csv")
    }
    (a/"manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(AssertionError):
        MODULE.analyze(a, b, tmp_path/"out", expected_cells=600)


def test_rejects_nonpassing_orchestration_audit(tmp_path):
    a = tmp_path/"run1"; b = tmp_path/"run2"
    write_replication(a, "r1", 100000)
    write_replication(b, "r2", 200000)
    (b/"independent_orchestration_audit.json").write_text(json.dumps({"status":"INCOMPLETE","planned_cells":600}))
    with pytest.raises(AssertionError, match="orchestration audit"):
        MODULE.analyze(a, b, tmp_path/"out", expected_cells=600)


def test_rejects_manifest_corruption(tmp_path):
    a = tmp_path/"run1"; b = tmp_path/"run2"
    write_replication(a, "r1", 100000)
    write_replication(b, "r2", 200000)
    with (a/"replicate_ledger.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(AssertionError, match="manifest mismatch"):
        MODULE.analyze(a, b, tmp_path/"out", expected_cells=600)
