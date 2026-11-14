from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _load_runner():
    path = ROOT / "scripts" / "run_quadratic_logit_curvature_density_smoke.py"
    spec = importlib.util.spec_from_file_location("quadratic_logit_density_smoke", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_auditor():
    path = ROOT / "scripts" / "audit_quadratic_logit_curvature_density_smoke.py"
    spec = importlib.util.spec_from_file_location("quadratic_logit_density_auditor", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _rehash(root: Path, relative_path: Path):
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    path = root / relative_path
    manifest["files"][str(relative_path)]["bytes"] = path.stat().st_size
    manifest["files"][str(relative_path)]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\\n")


def _args(output_dir: Path, gamma: float = 0.0):
    from argparse import Namespace

    return Namespace(
        output_dir=output_dir,
        seed=77101,
        validation_seed=77102,
        n_parent=500,
        validation_frame=180,
        validation_sizes=(60, 120, 180),
        gamma=gamma,
        beta=-3.9,
        quadrature_order=3,
        grid_points=101,
        scott_multiplier=0.7,
        absolute_widths=(0.012, 0.024),
    )


def test_quadratic_logit_density_smoke_writes_five_arm_artifacts(tmp_path: Path):
    runner = _load_runner()
    output = tmp_path / "gamma_zero"
    result = runner.run(_args(output, gamma=0.0))

    assert result["status"] == "PASS"
    assert result["gamma"] == 0.0
    assert result["weight_diagnostic_rows"] == 42
    assert (output / "manifest.json").is_file()
    assert (output / "results" / "density_reconstruction_inputs.npz").is_file()

    metrics = pd.read_csv(output / "results" / "density_metrics.csv")
    assert set(metrics["validation_size"]) == {60, 120, 180}
    assert set(metrics["absolute_bandwidth_sd"]) == {0.012, 0.024}
    assert set(metrics["status"]) == {"SUCCESS"}
    for arm in (
        "A_catalogue_only",
        "B_parent_oracle",
        "C_known_probability_ipw",
        "D_map_plugin_ipw",
        "E_posterior_mean_density",
        "E_posterior_weighted_node_metric",
    ):
        assert arm in set(metrics["arm"])

    diagnostics = pd.read_csv(output / "results" / "weight_diagnostics.csv")
    assert len(diagnostics) == 42
    assert set(diagnostics["validation_size"]) == {60, 120, 180}
    audit = _load_auditor().audit(output)
    assert audit["status"] == "PASS"
    assert audit["weight_diagnostic_rows_reconstructed"] == 42

    with np.load(output / "results" / "density_reconstruction_inputs.npz") as archive:
        probability = archive["catalogue_true_inclusion_probability"]
        base = archive["catalogue_base_weights"]
        known_weights = archive["weights_C_known_probability_ipw_n180"]
        np.testing.assert_allclose(known_weights, base / probability, rtol=1e-13, atol=1e-13)
        assert archive["posterior_nodes_n180"].shape[1] == 2


def test_quadratic_logit_density_smoke_records_curvature_and_separate_validation_seed(tmp_path: Path):
    runner = _load_runner()
    output = tmp_path / "gamma_curved"
    result = runner.run(_args(output, gamma=1.0))

    assert result["status"] == "PASS"
    config = __import__("json").loads((output / "config.json").read_text())
    assert config["gamma"] == 1.0
    assert config["population_seed"] != config["validation_seed"]
    assert config["analysis_model"] == "degree_one_logistic"
    assert config["arm_C"].startswith("true quadratic inclusion probabilities")
    assert config["arm_D_E"].startswith("degree-one logistic")
    observable = pd.read_csv(output / "observable" / "validation_n180.csv")
    assert list(observable.columns) == ["validation_id", "redshift", "retained"]
    assert "inclusion_probability" not in observable.columns


def test_quadratic_logit_auditor_rejects_corrupted_weight_diagnostic_even_with_updated_manifest(tmp_path: Path):
    runner = _load_runner()
    auditor = _load_auditor()
    output = tmp_path / "corrupted_gamma_zero"
    runner.run(_args(output, gamma=0.0))
    path = Path("results") / "weight_diagnostics.csv"
    diagnostics = pd.read_csv(output / path)
    diagnostics.loc[0, "kish_ess"] += 3.0
    diagnostics.to_csv(output / path, index=False)
    _rehash(output, path)
    import pytest
    with pytest.raises(AssertionError, match="weight diagnostic numerical mismatch: kish_ess"):
        auditor.audit(output)
