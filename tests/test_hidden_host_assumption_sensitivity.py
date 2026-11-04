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


RUN = load("hidden_assumption_runner", "scripts/run_hidden_host_assumption_sensitivity.py")
AUDIT = load("hidden_assumption_auditor", "scripts/audit_hidden_host_assumption_sensitivity.py")


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
    return kwargs, raw, summaries, contrasts


def test_assumption_sensitivity_pairing_and_independent_audit(tmp_path):
    kwargs, raw, summaries, contrasts = fixture(tmp_path)
    arms = RUN.arm_names(kwargs["assumed_fractions"], kwargs["assumed_redshifts"])
    assert len(raw) == 2 * 200 * len(arms)
    assert len(summaries) == 2 * len(arms)
    assert len(contrasts) == 2 * (len(arms) - 1) * len(RUN.METRICS)
    for start in range(0, len(raw), len(arms)):
        block = raw[start:start + len(arms)]
        assert [r["inference_arm"] for r in block] == list(arms)
        assert len({r["observed_distance_mpc"] for r in block}) == 1
    audit_kwargs = {
        "replicates": kwargs["replicates"], "seed": kwargs["seed"], "resamples": kwargs["bootstrap_resamples"],
        "truths": kwargs["truths"], "noises": kwargs["noises"], "fractions": kwargs["fractions"],
        "true_zs": kwargs["hidden_zs"], "assumed_fractions": kwargs["assumed_fractions"],
        "assumed_zs": kwargs["assumed_redshifts"],
    }
    result = AUDIT.audit(tmp_path, **audit_kwargs)
    assert result["audit"] == "PASS"
    assert result["posterior_reconstructions"] == len(raw)
    # At true hidden fraction zero, the oracle arm must equal the omitted-host model exactly.
    zero = [r for r in raw if float(r["true_hidden_fraction"]) == 0]
    by_key = {(int(r["replicate"]), r["inference_arm"]): r for r in zero}
    for rep in range(200):
        a, b = by_key[(rep, "base_omitted")], by_key[(rep, "oracle_hidden")]
        for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
            assert float(a[field]) == pytest.approx(float(b[field]), abs=1e-10)
    # A deliberately wrong positive hidden fraction should not collapse to the zero-fraction control.
    wrong = by_key[(0, RUN.assumption_label(0.15, 0.03))]
    base = by_key[(0, "base_omitted")]
    assert abs(float(wrong["posterior_median"]) - float(base["posterior_median"])) > 1e-6


def test_auditor_rejects_tampered_assumption_contrast(tmp_path):
    kwargs, raw, summaries, contrasts = fixture(tmp_path)
    altered = [dict(r) for r in contrasts]
    altered[0]["estimate"] = str(float(altered[0]["estimate"]) + 0.1)
    write_csv(tmp_path / "paired_contrasts.csv", altered)
    with pytest.raises(ValueError, match="paired contrast reconstruction mismatch"):
        AUDIT.audit(
            tmp_path, replicates=200, seed=11321, resamples=100,
            truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), true_zs=(0.03,),
            assumed_fractions=(0.15, 0.30), assumed_zs=(0.03,),
        )


def test_assumption_factorial_rejects_invalid_assumed_fraction():
    with pytest.raises(ValueError, match="strictly between zero and one"):
        RUN.run_experiment(
            replicates=200, seed=1, bootstrap_resamples=100,
            truths=(70.0,), noises=(0.1,), fractions=(0.0,), hidden_zs=(0.03,),
            assumed_fractions=(0.0,), assumed_redshifts=(0.03,),
        )
