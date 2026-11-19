from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "analyze_h0_100rep_calibration_uncertainty.py"
SPEC = importlib.util.spec_from_file_location("h0_failure_bounds", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def row(rep, gamma, nval, status="SUCCESS", *, seed_offset=0):
    ok = status == "SUCCESS"
    median = 70.0 + 0.4 * rep + 1.2 * gamma - 0.1 * (nval == 300) if ok else float("nan")
    return {
        "cell_id": f"r{rep}_g{gamma}_n{nval}",
        "replicate": rep, "gamma": gamma, "n_validation": nval,
        "analysis_seed": 1000 + rep + seed_offset,
        "validation_seed": 2000 + rep + seed_offset,
        "h0_true": 70.0, "status": status,
        "h0_median": median,
        "h0_q16": median - 5 if ok else float("nan"),
        "h0_q84": median + 5 if ok else float("nan"),
        "covered68": int(abs(median-70) <= 5) if ok else float("nan"),
        "covered90": 1 if ok else float("nan"),
        "covered95": 1 if ok else float("nan"),
        "min_catalogue_kish_ess": 50 + rep if ok else float("nan"),
        "max_normalized_weight": 0.04 if ok else float("nan"),
    }


def write_artifact(root: Path, rows):
    root.mkdir(parents=True, exist_ok=True)
    ledger = pd.DataFrame(rows)
    planned = ledger[["cell_id", "replicate", "gamma", "n_validation",
                      "analysis_seed", "validation_seed", "h0_true"]].copy()
    planned.to_csv(root / "planned_cells.csv", index=False)
    ledger.to_csv(root / "replicate_ledger.csv", index=False)
    manifest = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in ("planned_cells.csv", "replicate_ledger.csv")
    }
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True))
    return ledger


def test_complete_artifact_reports_success_metrics_and_paired_contrasts(tmp_path):
    rows = [row(rep, gamma, nval)
            for rep in range(4) for gamma in (0.0, 1.0) for nval in (100, 300)]
    root = tmp_path / "artifact"
    write_artifact(root, rows)
    report = MODULE.analyze(root, tmp_path / "out", seed=12)
    assert report["status"] == "PASS"
    assert report["planned_cells"] == 16
    cells = pd.read_csv(tmp_path / "out/cell_uncertainty_summary.csv")
    assert len(cells) == 4
    assert {"coverage68_wilson95_low", "coverage68_all_attempt_lower",
            "coverage68_all_attempt_upper"} <= set(cells.columns)
    gamma = pd.read_csv(tmp_path / "out/paired_gamma_contrasts.csv")
    assert set(gamma.metric) == {"signed_median_bias", "absolute_median_error", "width68"}
    size = pd.read_csv(tmp_path / "out/paired_validation_size_contrasts.csv")
    assert len(size) == 6


def test_failure_inclusive_bounds_preserve_all_planned_attempts(tmp_path):
    statuses = ["SUCCESS", "SUCCESS", "FAILED", "NON_EVALUABLE", "NOT_RUN"]
    rows = [row(rep, 0.0, 100, status) for rep, status in enumerate(statuses)]
    rows[1]["covered68"] = 0
    root = tmp_path / "artifact"
    write_artifact(root, rows)
    report = MODULE.analyze(root, tmp_path / "out", seed=13)
    cell = pd.read_csv(tmp_path / "out/cell_uncertainty_summary.csv").iloc[0]
    assert report["status"] == "INCOMPLETE"
    assert cell["coverage68"] == 0.5
    assert cell["coverage68_all_attempt_lower"] == 0.2
    assert cell["coverage68_all_attempt_upper"] == 0.8
    assert cell["non_success"] == 3
    assert cell["not_run"] == 1


def test_plan_ledger_mismatch_is_rejected(tmp_path):
    root = tmp_path / "artifact"
    write_artifact(root, [row(0, 0.0, 100), row(1, 0.0, 100)])
    ledger = pd.read_csv(root / "replicate_ledger.csv").iloc[:1]
    ledger.to_csv(root / "replicate_ledger.csv", index=False)
    manifest = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in ("planned_cells.csv", "replicate_ledger.csv")
    }
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(AssertionError, match="planned cell inventory"):
        MODULE.analyze(root, tmp_path / "out")


def test_corrupt_manifest_is_rejected(tmp_path):
    root = tmp_path / "artifact"
    write_artifact(root, [row(0, 0.0, 100)])
    (root / "replicate_ledger.csv").write_text("tampered")
    with pytest.raises(AssertionError, match="manifest mismatch"):
        MODULE.analyze(root, tmp_path / "out")


def test_pairing_seed_mismatch_is_rejected(tmp_path):
    rows = [row(0, 0.0, 100), row(0, 1.0, 100, seed_offset=1)]
    root = tmp_path / "artifact"
    write_artifact(root, rows)
    with pytest.raises(ValueError, match="analysis_seed"):
        MODULE.analyze(root, tmp_path / "out")


def test_zero_success_keeps_success_only_metrics_undefined(tmp_path):
    rows = [row(rep, 0.0, 100, "FAILED") for rep in range(3)]
    root = tmp_path / "artifact"
    write_artifact(root, rows)
    MODULE.analyze(root, tmp_path / "out")
    cell = pd.read_csv(tmp_path / "out/cell_uncertainty_summary.csv").iloc[0]
    assert pd.isna(cell["rmse_posterior_median"])
    assert pd.isna(cell["rmse_bootstrap95_low"])
    assert pd.isna(cell["rmse_bootstrap95_high"])
    assert cell["coverage68_all_attempt_lower"] == 0.0
    assert cell["coverage68_all_attempt_upper"] == 1.0
