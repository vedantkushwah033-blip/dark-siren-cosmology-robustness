import csv, hashlib, importlib.util, json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

RUN = load("posterior_runner_uncertainty", "run_posterior_calibration.py")
AUDIT = load("posterior_auditor_uncertainty", "audit_posterior_calibration.py")
UNC = load("posterior_uncertainty", "analyze_posterior_calibration_uncertainty.py")

def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def make_audited_fixture(root):
    raw, summaries = RUN.run(replicates=100, seed=91771)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summaries)
    meta = {
        "audit_status": "PENDING",
        "base_seed": 91771,
        "replicates_per_truth_noise_cell": 100,
        "raw_records": len(raw),
        "summary_rows": len(summaries),
        "sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                   for name in ("per_run.csv", "summary.csv")},
    }
    (root / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    assert AUDIT.audit(root)["audit"] == "PASS"

def test_cell_uncertainty_has_all_truth_noise_and_inference_arms(tmp_path):
    make_audited_fixture(tmp_path)
    rows = UNC.analyze_cell_uncertainty(tmp_path, seed=7001, resamples=100)
    assert len(rows) == 4 * 3 * 2 * 2 * 2
    assert {(r["true_H0"], r["noise_fraction"]) for r in rows} == {
        (t, n) for t in RUN.TRUTHS for n in RUN.NOISES
    }
    for row in rows:
        assert row["replicates"] == 100
        for metric in ("bias", "rmse", "coverage68", "coverage95"):
            assert row[f"{metric}_ci95_low"] <= row[f"{metric}_estimate"] <= row[f"{metric}_ci95_high"]
        assert 0 <= row["coverage68_estimate"] <= 1
        assert 0 <= row["coverage95_estimate"] <= 1
    report = tmp_path / "cell_uncertainty_report.md"
    UNC.write_cell_report(report, rows)
    lines = report.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("# Truth/noise-specific")
    assert sum(line.startswith("| ") for line in lines) == len(rows) + 2

def test_cell_uncertainty_refuses_non_audited_artifact(tmp_path):
    make_audited_fixture(tmp_path)
    meta_path = tmp_path / "run_metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["audit_status"] = "PENDING"
    meta_path.write_text(json.dumps(meta), encoding="utf-8")
    with pytest.raises(ValueError, match="must pass independent audit"):
        UNC.analyze_cell_uncertainty(tmp_path, resamples=20)
