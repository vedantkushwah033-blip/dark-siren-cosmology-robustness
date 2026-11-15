import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUN = load("paired_host_misspec_runner", "scripts/run_paired_host_misspecification.py")
AUDIT = load("paired_host_misspec_auditor", "scripts/audit_paired_host_misspecification.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(root):
    raw, summary, contrasts = RUN.run_experiment(replicates=200, seed=91772, bootstrap_resamples=100)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summary)
    write_csv(root / "paired_contrasts.csv", contrasts)
    return raw, summary, contrasts


def test_paired_experiment_and_independent_reconstruction(tmp_path):
    raw, summary, contrasts = make_fixture(tmp_path)
    assert len(raw) == 800
    assert len(summary) == 4
    assert len(contrasts) == 12
    # Every replicate has the same observation and latent truth across all arms.
    for replicate in range(200):
        block = [r for r in raw if int(r["replicate"]) == replicate]
        assert [r["inference_arm"] for r in block] == list(RUN.ARMS)
        assert len({r["observed_distance_mpc"] for r in block}) == 1
        assert len({r["true_host_index"] for r in block}) == 1
    result = AUDIT.audit(tmp_path, replicates=200, seed=91772, resamples=100)
    assert result["audit"] == "PASS"
    assert result["posterior_reconstructions"] == 800
    assert result["paired_contrasts_reconstructed"] == 12
    for row in contrasts:
        assert float(row["ci95_low"]) <= float(row["ci95_high"])


def test_auditor_rejects_changed_posterior_or_contrast(tmp_path):
    raw, summary, contrasts = make_fixture(tmp_path)
    damaged = [dict(r) for r in raw]
    damaged[0]["posterior_median"] = str(float(damaged[0]["posterior_median"]) + 0.1)
    write_csv(tmp_path / "per_run.csv", damaged)
    with pytest.raises(ValueError, match="Posterior/derived-field mismatch"):
        AUDIT.audit(tmp_path, replicates=200, seed=91772, resamples=100)

    write_csv(tmp_path / "per_run.csv", raw)
    damaged_contrasts = [dict(r) for r in contrasts]
    damaged_contrasts[0]["estimate"] = str(float(damaged_contrasts[0]["estimate"]) + 0.1)
    write_csv(tmp_path / "paired_contrasts.csv", damaged_contrasts)
    with pytest.raises(ValueError, match="Paired-bootstrap reconstruction mismatch"):
        AUDIT.audit(tmp_path, replicates=200, seed=91772, resamples=100)


def test_pairing_and_minimum_design_are_enforced(tmp_path):
    with pytest.raises(ValueError, match="at least 200"):
        RUN.generate_observations(199, 3)
    with pytest.raises(ValueError, match="at least 100"):
        RUN.run_experiment(200, 3, 99)
