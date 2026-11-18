from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_inclusion_quadrature_ladder import audit


def test_full_support_quadrature_gate_records_validation_extrapolation(tmp_path: Path):
    rng = np.random.default_rng(70118)
    root = tmp_path / "gamma_cell"
    child = root / "replicates" / "replicate_0000"
    (child / "observable").mkdir(parents=True)

    z_validation = rng.uniform(0.035, 0.365, 80)
    x = (z_validation - 0.205) / 0.195
    p = 1.0 / (1.0 + np.exp(-(-0.3 - 1.0 * x)))
    retained = rng.binomial(1, p)
    pd.DataFrame({"redshift": z_validation, "retained": retained}).to_csv(
        child / "observable" / "validation_n80.csv", index=False
    )
    pd.DataFrame({"redshift": np.linspace(0.01, 0.40, 120)}).to_csv(
        child / "observable" / "catalogue_observed.csv", index=False
    )
    (child / "config.json").write_text("{}")
    (root / "config.json").write_text(json.dumps({
        "catalogue_selection_model": "quadratic_logit",
        "replicates": 1,
        "validation_sizes": [80],
    }))

    report = audit(
        root,
        orders=(5, 7, 9),
        full_support_points=101,
        mean_tolerance=10.0,
        covariance_tolerance=10.0,
        predictive_tolerance=10.0,
        log_evidence_tolerance=1e6,
        mixture_max_abs_tolerance=1.0,
        mixture_total_variation_tolerance=1.0,
        consecutive_passes=1,
    )

    assert report["status"] == "PASS"
    assert report["completed_replicate_validation_cases"] == 1
    assert report["full_support"] == [0.01, 0.40]
    assert len(report["rows"]) >= 2
    assert all(row["full_support_min_node_inclusion"] > 0 for row in report["rows"])
    assert all(row["full_support_max_node_inclusion"] < 1 for row in report["rows"])
    assert all(
        row["catalogue_fraction_outside_validation_range"] > 0
        for row in report["rows"]
    )
    assert all(
        row["posterior_mean_outside_validation_ipw_mass"] >= 0
        for row in report["rows"]
    )


def test_nonconverged_quadrature_ladder_fails_closed(tmp_path: Path):
    rng = np.random.default_rng(70119)
    root = tmp_path / "gamma_cell"
    child = root / "replicates" / "replicate_0000"
    (child / "observable").mkdir(parents=True)

    z_validation = rng.uniform(0.01, 0.40, 120)
    x = (z_validation - 0.205) / 0.195
    p = 1.0 / (1.0 + np.exp(-(-0.2 - 0.8 * x)))
    retained = rng.binomial(1, p)
    pd.DataFrame({"redshift": z_validation, "retained": retained}).to_csv(
        child / "observable" / "validation_n120.csv", index=False
    )
    pd.DataFrame({"redshift": np.linspace(0.01, 0.40, 100)}).to_csv(
        child / "observable" / "catalogue_observed.csv", index=False
    )
    (child / "config.json").write_text("{}")
    (root / "config.json").write_text(json.dumps({
        "catalogue_selection_model": "quadratic_logit",
        "replicates": 1,
        "validation_sizes": [120],
    }))

    report = audit(
        root,
        orders=(5, 7, 9),
        full_support_points=101,
        mean_tolerance=1e-14,
        covariance_tolerance=1e-14,
        predictive_tolerance=1e-14,
        log_evidence_tolerance=1e-14,
        mixture_max_abs_tolerance=1.0,
        mixture_total_variation_tolerance=1.0,
        consecutive_passes=2,
    )

    assert report["status"] == "NOT_PASSED"
    assert report["completed_replicate_validation_cases"] == 1
    assert report["rows"]
    assert report["rows"][-1]["ladder_status"] == "NON_CONVERGED"
    assert any("did not converge" in message for message in report["failures"])
