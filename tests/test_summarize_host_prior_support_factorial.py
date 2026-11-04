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


RUN = load("host_prior_summary_runner", "scripts/run_host_prior_support_factorial.py")
SUMMARY = load("host_prior_summary", "scripts/summarize_host_prior_support_factorial.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture(root):
    raw, summaries, arm, prior = RUN.run_experiment(
        replicates=200, seed=9921, bootstrap_resamples=100,
        truths=(70.0,), noise_levels=(0.10,), hidden_fractions=(0.0, 0.30),
        hidden_redshifts=(0.03,),
    )
    for name, rows in (("per_run.csv", raw), ("summary.csv", summaries),
                       ("arm_contrasts.csv", arm), ("prior_contrasts.csv", prior)):
        write_csv(root / name, rows)


def test_marginal_summary_and_paired_bootstrap_are_reproducible(tmp_path):
    fixture(tmp_path)
    out = tmp_path / "out"
    first = SUMMARY.summarize(tmp_path, out, resamples=100, seed=123)
    saved = {p.name: p.read_text() for p in out.glob("marginal_*.csv")}
    second = SUMMARY.summarize(tmp_path, out, resamples=100, seed=123)
    assert first == second
    assert saved == {p.name: p.read_text() for p in out.glob("marginal_*.csv")}
    arms = list(csv.DictReader((out / "marginal_summary.csv").open(newline="", encoding="utf-8")))
    arm_contrasts = list(csv.DictReader((out / "marginal_arm_contrasts.csv").open(newline="", encoding="utf-8")))
    prior_contrasts = list(csv.DictReader((out / "marginal_prior_contrasts.csv").open(newline="", encoding="utf-8")))
    assert len(arms) == 2 * len(RUN.PRIORS) * len(RUN.ARMS)
    assert len(arm_contrasts) == 2 * len(RUN.PRIORS) * len(RUN.ARM_COMPARISONS) * len(RUN.METRICS)
    assert len(prior_contrasts) == 2 * len(RUN.ARMS) * len(RUN.PRIOR_COMPARISONS) * len(RUN.METRICS)
    zero_host = [r for r in arm_contrasts if float(r["hidden_fraction"]) == 0
                 and r["contrast"] == "omitted_minus_modelled_hidden"]
    assert all(float(r["equal_cell_mean_difference"]) == pytest.approx(0, abs=1e-12)
               for r in zero_host)
    # Prior sensitivity is identical for the base and hidden-modelled arms when the hidden fraction is zero.
    zero_prior = [r for r in prior_contrasts if float(r["hidden_fraction"]) == 0]
    by_key = {(r["contrast"], r["metric"], r["inference_arm"]): r for r in zero_prior}
    for contrast in RUN.PRIOR_COMPARISONS:
        label = contrast[2]
        for metric in RUN.METRICS:
            a = by_key[(label, metric, "base_weights")]
            b = by_key[(label, metric, "hidden_host_modelled")]
            assert float(a["equal_cell_mean_difference"]) == pytest.approx(
                float(b["equal_cell_mean_difference"]), abs=1e-12
            )


def test_summary_script_rejects_tampered_cell_summary(tmp_path):
    fixture(tmp_path)
    rows = list(csv.DictReader((tmp_path / "summary.csv").open(newline="", encoding="utf-8")))
    rows[0]["rmse_error"] = str(float(rows[0]["rmse_error"]) + 0.1)
    write_csv(tmp_path / "summary.csv", rows)
    with pytest.raises(ValueError, match="Raw-to-summary mismatch"):
        SUMMARY.summarize(tmp_path, tmp_path / "out", resamples=100, seed=123)
