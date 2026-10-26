import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_host_misspecification_audit.py"
SPEC = importlib.util.spec_from_file_location("host_misspecification_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_generated_host_support_matches_declared_scenarios():
    rows = module.generate_dataset(200, 17)
    assert len(rows) == 800
    for row in rows:
        assert row["observed_distance_mpc"] > 0
        assert row["sigma_mpc"] > 0
        if row["scenario"] in ("correct_weights", "reversed_weights"):
            assert row["true_host_redshift"] in (0.02, 0.04)
        else:
            assert row["true_host_redshift"] in (0.02, 0.03, 0.04)


def test_stress_test_metrics_are_finite_and_bounded():
    rows, summaries = module.run_stress_test(replicates=200, seed=18)
    assert len(rows) == 800
    assert len(summaries) == 4
    for row in summaries:
        for key in ("mean_median_error", "rmse_median", "coverage_68", "coverage_95",
                    "mean_width_68", "mean_width_95"):
            assert np.isfinite(row[key])
        assert row["rmse_median"] >= 0
        assert 0 <= row["coverage_68"] <= 1
        assert 0 <= row["coverage_95"] <= 1
        assert row["mean_width_68"] > 0
        assert row["mean_width_95"] > row["mean_width_68"]

def test_independent_auditor_reconstructs_seeded_small_run(tmp_path):
    import csv
    import importlib.util

    audit_path = Path(__file__).resolve().parents[1] / "scripts" / "audit_host_misspecification_audit.py"
    audit_spec = importlib.util.spec_from_file_location("independent_host_audit", audit_path)
    auditor = importlib.util.module_from_spec(audit_spec)
    audit_spec.loader.exec_module(auditor)

    rows, summaries = module.run_stress_test(replicates=200, seed=29)
    for filename, records in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (tmp_path / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    result = auditor.audit_directory(tmp_path, replicates=200, seed=29)
    assert result["status"] == "PASS"
    assert result["raw_rows"] == 800
    assert result["posterior_reconstructions"] == 800

    import pytest

    raw_path = tmp_path / "per_run.csv"
    tampered = list(rows)
    tampered[0] = {**tampered[0], "width_68": tampered[0]["width_68"] + 1.0}
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tampered[0]))
        writer.writeheader()
        writer.writerows(tampered)
    with pytest.raises(ValueError, match="Per-run derived-field mismatch"):
        auditor.audit_directory(tmp_path, replicates=200, seed=29)


    # The audit must reject missing raw rows even when all remaining values
    # are internally consistent, and it must reconstruct aggregate metrics
    # rather than trusting the summary file.
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows[1:])
    with pytest.raises(ValueError, match="Raw/summary row count"):
        auditor.audit_directory(tmp_path, replicates=200, seed=29)

    # Restore raw records and corrupt an aggregate coverage metric.
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = tmp_path / "summary.csv"
    corrupted_summaries = [dict(s) for s in summaries]
    corrupted_summaries[0]["coverage_68"] = min(
        1.0, float(corrupted_summaries[0]["coverage_68"]) + 0.05
    )
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(corrupted_summaries[0]))
        writer.writeheader()
        writer.writerows(corrupted_summaries)
    with pytest.raises(ValueError, match="Summary reconstruction mismatch"):
        auditor.audit_directory(tmp_path, replicates=200, seed=29)
