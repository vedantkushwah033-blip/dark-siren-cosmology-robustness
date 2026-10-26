"""Tests for multi-seed paired grid-sensitivity aggregation."""
import csv
import json

import pytest

from scripts.run_multiseed_grid_sensitivity import run_multiseed_sensitivity


def test_multiseed_runner_preserves_seed_level_audits_and_aggregates(tmp_path):
    results = run_multiseed_sensitivity(
        grid_points=[1201, 2401],
        seeds=[101, 202],
        n_replicates=100,
        output_dir=tmp_path,
    )
    assert [row["grid_points"] for row in results] == [1201, 2401]
    assert all(row["n_seeds"] == 2 for row in results)
    assert all(row["replicates_per_seed"] == 100 for row in results)
    assert all(row["total_replicates_across_seeds"] == 200 for row in results)
    assert all(row["all_seed_grid_audits_passed"] == "PASS" for row in results)
    finest = next(row for row in results if row["grid_points"] == 2401)
    assert finest["mean_paired_pit_delta_vs_finest_across_seeds"] == 0.0
    assert finest["sd_paired_pit_delta_across_seeds"] == 0.0
    assert finest["se_paired_pit_delta_across_seeds"] == 0.0

    with (tmp_path / "per_seed_paired_contrasts.csv").open(newline="", encoding="utf-8") as handle:
        per_seed = list(csv.DictReader(handle))
    assert len(per_seed) == 4
    assert {int(row["seed"]) for row in per_seed} == {101, 202}
    assert all(row["independent_audit"] == "PASS" for row in per_seed)
    for seed in (101, 202):
        assert (tmp_path / f"seed_{seed}" / "grid_sensitivity.csv").exists()
        for grid in (1201, 2401):
            assert (tmp_path / f"seed_{seed}" / f"grid_{grid}" / "audit_report.json").exists()

    metadata = json.loads((tmp_path / "multiseed_metadata.json").read_text(encoding="utf-8"))
    assert metadata["status"] == "PASS"
    assert metadata["seed_is_uncertainty_unit_for_aggregate_se"] is True
    assert metadata["seeds"] == [101, 202]


@pytest.mark.parametrize(
    ("grids", "seeds", "replicates", "message"),
    [
        ([1201], [1, 2], 100, "at least two grid resolutions"),
        ([1201, 2401], [1], 100, "At least two distinct seeds"),
        ([1201, 2401], [1, 1], 100, "Seeds must be unique"),
        ([1201, 2401], [1, 2], 99, "at least 100"),
        ([1201, 2401], [1, -2], 100, "non-negative integers"),
        ([1201, 2401], [1, True], 100, "non-negative integers"),
        ([1201, 1201], [1, 2], 100, "Grid resolutions must be unique"),
        ([1201, 2401.5], [1, 2], 100, "at least two points"),
    ],
)
def test_multiseed_runner_rejects_invalid_design(tmp_path, grids, seeds, replicates, message):
    with pytest.raises(ValueError, match=message):
        run_multiseed_sensitivity(grids, seeds, replicates, tmp_path)
    assert not (tmp_path / "multiseed_grid_sensitivity.csv").exists()


def test_failed_seed_run_does_not_leave_stale_success_marker(tmp_path, monkeypatch):
    import scripts.run_multiseed_grid_sensitivity as runner

    marker = tmp_path / "multiseed_metadata.json"
    marker.write_text(json.dumps({"status": "PASS"}), encoding="utf-8")

    def fail_before_audit(**kwargs):
        raise RuntimeError("synthetic failure for fail-closed regression test")

    monkeypatch.setattr(runner, "run_sensitivity", fail_before_audit)
    with pytest.raises(RuntimeError, match="synthetic failure"):
        runner.run_multiseed_sensitivity(
            grid_points=[1201, 2401],
            seeds=[101, 202],
            n_replicates=100,
            output_dir=tmp_path,
        )
    assert not marker.exists()
    assert not (tmp_path / "multiseed_grid_sensitivity.csv").exists()
