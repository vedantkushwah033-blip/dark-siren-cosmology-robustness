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


RUN = load("hidden_assumption_summary_runner", "scripts/run_hidden_host_assumption_sensitivity.py")
SUMMARY = load("hidden_assumption_summary", "scripts/summarize_hidden_host_assumption_sensitivity.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture(root):
    kwargs = dict(
        replicates=200, seed=11321, bootstrap_resamples=100,
        truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), hidden_zs=(0.03,),
        assumed_fractions=(0.15, 0.30), assumed_redshifts=(0.03,),
    )
    raw, summaries, contrasts = RUN.run_experiment(**kwargs)
    for name, rows in (("per_run.csv", raw), ("summary.csv", summaries), ("paired_contrasts.csv", contrasts)):
        write_csv(root / name, rows)


def test_marginal_summary_and_bootstrap_are_reproducible(tmp_path):
    fixture(tmp_path)
    out = tmp_path / "out"
    first = SUMMARY.summarize(tmp_path, out, resamples=100, seed=42)
    saved = {p.name: p.read_text() for p in out.glob("marginal_*.csv")}
    second = SUMMARY.summarize(tmp_path, out, resamples=100, seed=42)
    assert first == second
    assert saved == {p.name: p.read_text() for p in out.glob("marginal_*.csv")}
    marginal = list(csv.DictReader((out / "marginal_summary.csv").open(newline="", encoding="utf-8")))
    contrasts = list(csv.DictReader((out / "marginal_contrasts.csv").open(newline="", encoding="utf-8")))
    arms = RUN.arm_names((0.15, 0.30), (0.03,))
    assert len(marginal) == 2 * len(arms)
    assert len(contrasts) == 2 * (len(arms) - 1) * len(RUN.METRICS)
    zero = [r for r in contrasts if float(r["true_hidden_fraction"]) == 0
            and r["contrast"] == "base_omitted_minus_oracle"]
    assert all(float(r["equal_cell_mean_difference"]) == pytest.approx(0, abs=1e-12) for r in zero)
    assert all(float(r["ci95_low"]) <= float(r["ci95_high"]) for r in contrasts)


def test_summary_rejects_tampered_cell_summary(tmp_path):
    fixture(tmp_path)
    rows = list(csv.DictReader((tmp_path / "summary.csv").open(newline="", encoding="utf-8")))
    rows[0]["rmse_error"] = str(float(rows[0]["rmse_error"]) + 0.1)
    write_csv(tmp_path / "summary.csv", rows)
    with pytest.raises(ValueError, match="Raw-to-summary mismatch"):
        SUMMARY.summarize(tmp_path, tmp_path / "out", resamples=100, seed=42)
