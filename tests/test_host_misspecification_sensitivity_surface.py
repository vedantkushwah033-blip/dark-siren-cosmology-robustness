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


RUN = load("host_sensitivity_runner", "scripts/run_host_misspecification_sensitivity_surface.py")
AUDIT = load("host_sensitivity_auditor", "scripts/audit_host_misspecification_sensitivity_surface.py")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_small_design(root):
    raw, summary, contrasts = RUN.run_experiment(
        replicates=200, seed=7712, bootstrap_resamples=100,
        truths=(70.0,), noise_levels=(0.10,), hidden_fractions=(0.0, 0.30),
        hidden_redshifts=(0.03,),
    )
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summary)
    write_csv(root / "paired_contrasts.csv", contrasts)
    return raw, summary, contrasts


def test_factorial_shape_pairing_negative_control_and_independent_audit(tmp_path):
    assert len(list(RUN.cell_specs())) == 108
    raw, summary, contrasts = make_small_design(tmp_path)
    assert len(raw) == 2 * 200 * 4
    assert len(summary) == 2 * 4
    assert len(contrasts) == 2 * 2 * len(RUN.METRICS)
    for block_start in range(0, len(raw), 4):
        block = raw[block_start:block_start + 4]
        assert [r["inference_arm"] for r in block] == list(RUN.ARMS)
        assert len({r["observed_distance_mpc"] for r in block}) == 1
        assert len({r["true_host_index"] for r in block}) == 1
    result = AUDIT.audit(
        tmp_path, replicates=200, seed=7712, resamples=100,
        truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), hidden_zs=(0.03,),
    )
    assert result["audit"] == "PASS"
    assert result["posterior_reconstructions"] == len(raw)
    # At hidden fraction zero, base and modelled-host posteriors are the same model.
    f0 = [r for r in raw if float(r["hidden_fraction"]) == 0.0]
    by_rep_arm = {(int(r["replicate"]), r["inference_arm"]): r for r in f0}
    for rep in range(200):
        a = by_rep_arm[(rep, "base_weights")]
        b = by_rep_arm[(rep, "hidden_host_modelled")]
        for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
            assert float(a[field]) == pytest.approx(float(b[field]), abs=1e-10)


def test_independent_auditor_rejects_tampered_paired_output(tmp_path):
    raw, summary, contrasts = make_small_design(tmp_path)
    altered = [dict(r) for r in contrasts]
    altered[0]["estimate"] = str(float(altered[0]["estimate"]) + 0.1)
    write_csv(tmp_path / "paired_contrasts.csv", altered)
    with pytest.raises(ValueError, match="paired contrast reconstruction mismatch"):
        AUDIT.audit(
            tmp_path, replicates=200, seed=7712, resamples=100,
            truths=(70.0,), noises=(0.10,), fractions=(0.0, 0.30), hidden_zs=(0.03,),
        )


def test_minimum_replicate_and_bootstrap_counts():
    with pytest.raises(ValueError, match="at least 200"):
        RUN.run_experiment(replicates=199, seed=1, bootstrap_resamples=100,
                           truths=(70.0,), noise_levels=(0.1,),
                           hidden_fractions=(0.0,), hidden_redshifts=(0.03,))
    with pytest.raises(ValueError, match="at least 100"):
        RUN.run_experiment(replicates=200, seed=1, bootstrap_resamples=99,
                           truths=(70.0,), noise_levels=(0.1,),
                           hidden_fractions=(0.0,), hidden_redshifts=(0.03,))
