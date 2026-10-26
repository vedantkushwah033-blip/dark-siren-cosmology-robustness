"""Tests for the independently audited PIT grid-sensitivity experiment."""
import csv
import json

import pytest

import scripts.run_calibration_grid_sensitivity as grid_sensitivity
from scripts.run_calibration_grid_sensitivity import run_sensitivity


def test_grid_sensitivity_replays_and_audits_each_resolution(tmp_path):
    results = run_sensitivity(
        grid_points=[1201, 2401],
        n_replicates=100,
        seed=86420,
        output_dir=tmp_path,
    )
    assert [row["grid_points"] for row in results] == [1201, 2401]
    assert all(row["independent_audit"] == "PASS" for row in results)
    assert all(row["n_replicates"] == 100 for row in results)
    assert all(row["reference_grid_points"] == 2401 for row in results)
    finest = next(row for row in results if row["grid_points"] == 2401)
    assert finest["mean_pit_delta_vs_finest"] == 0.0
    assert finest["rms_pit_delta_vs_finest"] == 0.0
    assert finest["max_abs_pit_delta_vs_finest"] == 0.0
    assert results[0]["rms_pit_delta_vs_finest"] >= 0.0
    for row in results:
        assert row["mean_pit_delta_normal_95_low"] <= row["mean_pit_delta_vs_finest"]
        assert row["mean_pit_delta_vs_finest"] <= row["mean_pit_delta_normal_95_high"]
        assert row["mean_pit_delta_monte_carlo_se"] >= 0.0
    raw_by_grid = {}
    for points in (1201, 2401):
        case_dir = tmp_path / f"grid_{points}"
        raw_path = case_dir / "pit_replicates.csv"
        assert raw_path.exists()
        assert (case_dir / "pit_summary.json").exists()
        report = json.loads((case_dir / "audit_report.json").read_text(encoding="utf-8"))
        assert report["status"] == "PASS"
        assert report["n_replicates_audited"] == 100
        with raw_path.open(newline="", encoding="utf-8") as handle:
            raw_by_grid[points] = list(csv.DictReader(handle))

    # The experiment is paired: only posterior-grid discretization should
    # change, not the seeded synthetic truth/host/observation records.
    for field in ("truth_h0", "host_index", "host_redshift", "observed_distance_mpc", "rejection_attempts"):
        assert [row[field] for row in raw_by_grid[1201]] == [
            row[field] for row in raw_by_grid[2401]
        ]
    with (tmp_path / "grid_sensitivity.csv").open(newline="", encoding="utf-8") as handle:
        table = list(csv.DictReader(handle))
    assert len(table) == 2
    assert [int(row["grid_points"]) for row in table] == [1201, 2401]


@pytest.mark.parametrize(
    ("grid_points", "message"),
    [
        ([2401], "at least two grid resolutions"),
        ([1201, 1201], "must be unique"),
        ([1201, 1], "at least two points"),
    ],
)
def test_grid_sensitivity_rejects_invalid_grid_sets(tmp_path, grid_points, message):
    with pytest.raises(ValueError, match=message):
        run_sensitivity(grid_points, n_replicates=100, seed=1, output_dir=tmp_path)


@pytest.mark.parametrize("n_replicates", [0, 1, 2, 99, -5, 1.5, True])
def test_grid_sensitivity_rejects_invalid_replicate_counts(tmp_path, n_replicates):
    with pytest.raises(ValueError, match="n_replicates must be an integer of at least 100"):
        run_sensitivity([1201, 2401], n_replicates=n_replicates, seed=1, output_dir=tmp_path)


@pytest.mark.parametrize("seed", [1.5, "seed", True])
def test_grid_sensitivity_rejects_invalid_seeds(tmp_path, seed):
    with pytest.raises(ValueError, match="seed must be a non-negative integer"):
        run_sensitivity([1201, 2401], n_replicates=100, seed=seed, output_dir=tmp_path)


def test_grid_sensitivity_rejects_non_integer_grid_sizes(tmp_path):
    with pytest.raises(ValueError, match="Grid resolutions must be integers"):
        run_sensitivity([1201, 2401.5], n_replicates=100, seed=1, output_dir=tmp_path)


def test_grid_sensitivity_stops_when_independent_audit_fails(tmp_path, monkeypatch):
    stale_table = tmp_path / "grid_sensitivity.csv"
    stale_table.write_text("stale success\n", encoding="utf-8")
    stale_audit = tmp_path / "grid_1201" / "audit_report.json"
    stale_audit.parent.mkdir(parents=True)
    stale_audit.write_text('{"status": "PASS"}', encoding="utf-8")
    monkeypatch.setattr(grid_sensitivity, "audit", lambda *_args, **_kwargs: {"status": "FAIL"})
    with pytest.raises(RuntimeError, match="Independent audit failed for 1201-point grid"):
        run_sensitivity([1201, 2401], n_replicates=100, seed=123, output_dir=tmp_path)
    assert (tmp_path / "grid_1201" / "pit_replicates.csv").exists()
    assert not stale_table.exists()
    assert not stale_audit.exists()


@pytest.mark.parametrize("seed", [-1, -100])
def test_grid_sensitivity_rejects_negative_seeds(tmp_path, seed):
    with pytest.raises(ValueError, match="seed must be a non-negative integer"):
        run_sensitivity([1201, 2401], n_replicates=100, seed=seed, output_dir=tmp_path)
