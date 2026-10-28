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


RUN = load("host_prior_factorial_runner", "scripts/run_host_prior_support_factorial.py")
AUDIT = load("host_prior_factorial_auditor", "scripts/audit_host_prior_support_factorial.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_design(root):
    raw, summaries, arm_contrasts, prior_contrasts = RUN.run_experiment(
        replicates=200, seed=8812, bootstrap_resamples=100,
        truths=(70.0,), noise_levels=(0.10,), hidden_fractions=(0.0, 0.30),
        hidden_redshifts=(0.03,),
    )
    for name, rows in (("per_run.csv", raw), ("summary.csv", summaries),
                       ("arm_contrasts.csv", arm_contrasts), ("prior_contrasts.csv", prior_contrasts)):
        write_csv(root / name, rows)
    return raw, summaries, arm_contrasts, prior_contrasts


def test_paired_prior_factorial_and_independent_reconstruction(tmp_path):
    raw, summaries, arm_contrasts, prior_contrasts = make_design(tmp_path)
    assert len(raw) == 2 * 200 * len(RUN.PRIORS) * len(RUN.ARMS)
    assert len(summaries) == 2 * len(RUN.PRIORS) * len(RUN.ARMS)
    assert len(arm_contrasts) == 2 * len(RUN.PRIORS) * len(RUN.ARM_COMPARISONS) * len(RUN.METRICS)
    assert len(prior_contrasts) == 2 * len(RUN.ARMS) * len(RUN.PRIOR_COMPARISONS) * len(RUN.METRICS)
    block_size = len(RUN.PRIORS) * len(RUN.ARMS)
    for start in range(0, len(raw), block_size):
        block = raw[start:start + block_size]
        assert len({r["observed_distance_mpc"] for r in block}) == 1
        assert [r["inference_arm"] for r in block[:len(RUN.ARMS)]] == list(RUN.ARMS)
    result = AUDIT.audit(
        tmp_path, replicates=200, seed=8812, resamples=100,
        truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), hidden_zs=(0.03,),
    )
    assert result["audit"] == "PASS"
    assert result["posterior_reconstructions"] == len(raw)
    # At hidden fraction zero, the base and hidden-modelled arms coincide for every prior.
    zero = [r for r in raw if float(r["hidden_fraction"]) == 0]
    by_key = {(float(r["prior_min"]), float(r["prior_max"]), int(r["replicate"]), r["inference_arm"]): r
              for r in zero}
    for lo, hi in RUN.PRIORS:
        for rep in range(200):
            a, b = by_key[(lo, hi, rep, "base_weights")], by_key[(lo, hi, rep, "hidden_host_modelled")]
            for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                assert float(a[field]) == pytest.approx(float(b[field]), abs=1e-10)


def test_auditor_rejects_tampered_prior_contrast(tmp_path):
    make_design(tmp_path)
    altered = list(csv.DictReader((tmp_path / "prior_contrasts.csv").open(newline="", encoding="utf-8")))
    altered[0]["ci95_low"] = str(float(altered[0]["ci95_low"]) - 0.2)
    write_csv(tmp_path / "prior_contrasts.csv", altered)
    with pytest.raises(ValueError, match="prior contrast reconstruction mismatch"):
        AUDIT.audit(
            tmp_path, replicates=200, seed=8812, resamples=100,
            truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), hidden_zs=(0.03,),
        )


def test_prior_factorial_rejects_undersized_design():
    with pytest.raises(ValueError, match="at least 200"):
        RUN.run_experiment(
            replicates=199, seed=1, bootstrap_resamples=100,
            truths=(70.0,), noise_levels=(0.1,), hidden_fractions=(0.0,), hidden_redshifts=(0.03,),
        )
