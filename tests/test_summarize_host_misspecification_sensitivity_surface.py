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


RUN = load("sensitivity_runner_for_summary", "scripts/run_host_misspecification_sensitivity_surface.py")
SUMMARIZER = load("sensitivity_surface_summarizer", "scripts/summarize_host_misspecification_sensitivity_surface.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture(root):
    raw, summary, contrasts = RUN.run_experiment(
        replicates=200, seed=7712, bootstrap_resamples=100,
        truths=(70.0,), noise_levels=(0.10,), hidden_fractions=(0.0, 0.30),
        hidden_redshifts=(0.03,),
    )
    for name, rows in (("per_run.csv", raw), ("summary.csv", summary), ("paired_contrasts.csv", contrasts)):
        write_csv(root / name, rows)


def test_marginal_summary_and_bootstrap_are_reproducible(tmp_path):
    fixture(tmp_path)
    out = tmp_path / "out"
    first = SUMMARIZER.summarize(tmp_path, out, resamples=100, seed=91)
    saved_arm = (out / "marginal_arm_summary.csv").read_text()
    saved_contrast = (out / "marginal_paired_contrasts.csv").read_text()
    second = SUMMARIZER.summarize(tmp_path, out, resamples=100, seed=91)
    assert first == second
    assert saved_arm == (out / "marginal_arm_summary.csv").read_text()
    assert saved_contrast == (out / "marginal_paired_contrasts.csv").read_text()
    arms = list(csv.DictReader((out / "marginal_arm_summary.csv").open()))
    contrasts = list(csv.DictReader((out / "marginal_paired_contrasts.csv").open()))
    assert len(arms) == 8
    assert len(contrasts) == 12
    zero_fraction = [r for r in contrasts if float(r["hidden_fraction"]) == 0
                     and r["contrast"] == "omitted_minus_modelled_hidden"]
    assert all(float(r["equal_cell_mean_difference"]) == pytest.approx(0, abs=1e-12)
               for r in zero_fraction)
    assert all(float(r["ci95_low"]) <= float(r["ci95_high"]) for r in contrasts)


def test_summarizer_rejects_tampered_cell_summary(tmp_path):
    fixture(tmp_path)
    rows = list(csv.DictReader((tmp_path / "summary.csv").open(newline="", encoding="utf-8")))
    rows[0]["rmse_error"] = str(float(rows[0]["rmse_error"]) + 0.1)
    write_csv(tmp_path / "summary.csv", rows)
    with pytest.raises(ValueError, match="Raw-to-summary mismatch"):
        SUMMARIZER.summarize(tmp_path, tmp_path / "out", resamples=100, seed=91)
