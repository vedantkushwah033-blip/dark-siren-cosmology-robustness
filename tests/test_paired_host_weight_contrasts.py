import csv, hashlib, importlib.util, json
from pathlib import Path
import pytest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

RUN = load("paired_contrast_runner", "run_posterior_calibration.py")
AUDIT = load("paired_contrast_auditor", "audit_posterior_calibration.py")
PAIR = load("paired_host_weight_contrasts", "analyze_paired_host_weight_contrasts.py")
AUDIT_INDEP = load("independent_paired_contrast_audit", "audit_paired_host_weight_contrasts.py")

def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def make_audited_fixture(root):
    raw, summaries = RUN.run(replicates=100, seed=91772)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summaries)
    meta = {
        "audit_status": "PENDING",
        "base_seed": 91772,
        "replicates_per_truth_noise_cell": 100,
        "raw_records": len(raw),
        "summary_rows": len(summaries),
        "sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                   for name in ("per_run.csv", "summary.csv")},
    }
    (root / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    assert AUDIT.audit(root)["audit"] == "PASS"

def test_paired_host_weight_contrasts_complete_and_ordered(tmp_path):
    make_audited_fixture(tmp_path)
    raw, reps = PAIR.verify_artifact(tmp_path)
    # The vectorized bootstrap path must match scalar paired-metric evaluation.
    first_cell = [row for row in raw if float(row["true_H0"]) == RUN.TRUTHS[0]
                  and float(row["noise_fraction"]) == RUN.NOISES[0]]
    prefix = "broad_20_220_uniform_H0_"
    matrix_indices = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
    vectorized = PAIR._paired_metrics(first_cell, matrix_indices, prefix, RUN.TRUTHS[0])
    scalar = [
        PAIR._paired_metrics(first_cell, indices, prefix, RUN.TRUTHS[0])
        for indices in matrix_indices
    ]
    assert vectorized.shape == (2, len(PAIR.METRICS))
    np.testing.assert_allclose(vectorized, np.asarray(scalar))
    cells, aggregate = PAIR.analyze(raw, seed=331, resamples=100)
    cells_again, aggregate_again = PAIR.analyze(raw, seed=331, resamples=100)
    assert cells_again == cells
    assert aggregate_again == aggregate
    assert reps == 100
    assert len(cells) == 4 * 3 * 2 * 2
    assert len(aggregate) == 2 * 2 * len(PAIR.METRICS)
    assert {(r["true_H0"], r["noise_fraction"]) for r in cells} == {
        (t, n) for t in RUN.TRUTHS for n in RUN.NOISES
    }
    for row in cells + aggregate:
        for metric in PAIR.METRICS:
            if metric + "_ci95_low" in row:
                assert row[metric + "_ci95_low"] <= row[metric + "_ci95_high"]
        if "ci95_low" in row:
            assert row["ci95_low"] <= row["ci95_high"]
    broad_uniform = [r for r in aggregate if r["prior_support"] == "broad_20_220"
                     and r["prior_measure"] == "uniform_H0"]
    by_metric = {r["metric"]: r for r in broad_uniform}
    assert by_metric["bias_difference"]["ci95_low"] > 0
    assert by_metric["rmse_difference"]["ci95_low"] > 0
    assert by_metric["coverage68_difference"]["ci95_high"] < 0
    PAIR.write_csv(tmp_path / "cell.csv", cells)
    PAIR.write_csv(tmp_path / "aggregate.csv", aggregate)
    PAIR.write_report(tmp_path / "report.md", aggregate)
    assert len(list(csv.DictReader((tmp_path / "cell.csv").open()))) == 48
    assert len(list(csv.DictReader((tmp_path / "aggregate.csv").open()))) == 24
    assert "equal-weight average" in (tmp_path / "report.md").read_text().lower()
    PAIR.write_csv(tmp_path / "paired_host_weight_cell_contrasts.csv", cells)
    PAIR.write_csv(tmp_path / "paired_host_weight_contrasts.csv", aggregate)
    independent = AUDIT_INDEP.audit(tmp_path, seed=331, resamples=100)
    assert independent["audit"] == "PASS"
    assert independent["cell_contrast_rows_reconstructed"] == 48
    assert independent["aggregate_contrast_rows_reconstructed"] == 24

def test_paired_contrasts_refuse_non_audited_or_modified_artifact(tmp_path):
    make_audited_fixture(tmp_path)
    meta_path = tmp_path / "run_metadata.json"
    meta = json.loads(meta_path.read_text())
    meta["audit_status"] = "PENDING"
    meta_path.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="independently audited"):
        PAIR.verify_artifact(tmp_path)


def test_paired_contrasts_reject_incomplete_replicate_blocks(tmp_path):
    make_audited_fixture(tmp_path)
    raw, _ = PAIR.verify_artifact(tmp_path)
    # Drop one observation from only one truth/noise cell: the design is unbalanced.
    damaged = [row for row in raw if not (
        float(row["true_H0"]) == RUN.TRUTHS[0]
        and float(row["noise_fraction"]) == RUN.NOISES[0]
        and int(row["replicate"]) == 99
    )]
    with pytest.raises(ValueError, match="unequal replicate counts|missing, duplicate"):
        PAIR.analyze(damaged, seed=331, resamples=100)

    # A duplicated replicate ID must be rejected rather than silently reweighted.
    duplicated = list(raw) + [dict(raw[0])]
    with pytest.raises(ValueError, match="missing, duplicate"):
        PAIR.analyze(duplicated, seed=331, resamples=100)
