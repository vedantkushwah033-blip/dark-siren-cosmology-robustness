"""Regression tests for the cross-gamma pilot's complete-cell gate."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_quadratic_logit_curvature_pilot import _require_complete_replicates


def test_complete_gamma_cell_is_accepted():
    _require_complete_replicates(
        {"planned_replicates": 20, "successful_replicates": 20, "failed_replicates": 0},
        gamma=0.0,
    )


@pytest.mark.parametrize(
    "report",
    [
        {"planned_replicates": 20, "successful_replicates": 19, "failed_replicates": 1},
        {"planned_replicates": 20, "successful_replicates": 18, "failed_replicates": 2},
    ],
)
def test_incomplete_gamma_cell_is_rejected(report):
    with pytest.raises(AssertionError, match="incomplete gamma=1.0 cell"):
        _require_complete_replicates(report, gamma=1.0)



def test_paired_curvature_contrasts_use_matched_replicates_against_gamma_zero():
    import pandas as pd

    from audit_quadratic_logit_curvature_pilot import _paired_curvature_contrasts

    rows = []
    baseline = [1.0, 2.0, 3.0]
    stressed = [2.0, 4.0, 6.0]
    for gamma, values in ((0.0, baseline), (1.0, stressed)):
        for replicate, iae in enumerate(values):
            rows.append({
                "gamma": gamma,
                "replicate": replicate,
                "validation_size": 100,
                "absolute_bandwidth_sd": 0.016,
                "arm": "A_catalogue_only",
                "node": -1,
                "iae": iae,
            })
    result = _paired_curvature_contrasts(pd.DataFrame(rows), reference_gamma=0.0)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["gamma"] == 1.0
    assert row["reference_gamma"] == 0.0
    assert row["replicate_count"] == 3
    assert row["mean_paired_difference"] == 2.0
    assert row["contrast"] == "gamma_minus_gamma0_iae"


def test_paired_curvature_contrasts_reject_unmatched_metric_keys():
    import pandas as pd

    from audit_quadratic_logit_curvature_pilot import _paired_curvature_contrasts

    rows = [
        {"gamma": 0.0, "replicate": 0, "validation_size": 100,
         "absolute_bandwidth_sd": 0.016, "arm": "A_catalogue_only", "node": -1, "iae": 1.0},
        {"gamma": 0.0, "replicate": 1, "validation_size": 100,
         "absolute_bandwidth_sd": 0.016, "arm": "A_catalogue_only", "node": -1, "iae": 2.0},
        {"gamma": 1.0, "replicate": 0, "validation_size": 100,
         "absolute_bandwidth_sd": 0.016, "arm": "A_catalogue_only", "node": -1, "iae": 1.5},
    ]
    with pytest.raises(AssertionError, match="identical paired metric keys"):
        _paired_curvature_contrasts(pd.DataFrame(rows), reference_gamma=0.0)


def test_paired_contrasts_include_posterior_mean_density_against_map_plugin():
    import pandas as pd

    from audit_quadratic_logit_curvature_pilot import _paired_contrasts

    rows = []
    values = {
        "D_map_plugin_ipw": [0.30, 0.40, 0.50],
        "E_posterior_mean_density": [0.20, 0.35, 0.50],
        "C_known_probability_ipw": [0.10, 0.10, 0.10],
        "A_catalogue_only": [0.70, 0.70, 0.70],
        "B_parent_oracle": [0.08, 0.08, 0.08],
    }
    for arm, replicate_values in values.items():
        for replicate, iae in enumerate(replicate_values):
            rows.append({
                "gamma": 0.5,
                "replicate": replicate,
                "validation_size": 300,
                "absolute_bandwidth_sd": 0.016,
                "arm": arm,
                "node": -1,
                "iae": iae,
            })

    result = _paired_contrasts(pd.DataFrame(rows))
    row = result[
        (result.arm == "E_posterior_mean_density")
        & (result.reference_arm == "D_map_plugin_ipw")
    ].iloc[0]
    assert row.replicate_count == 3
    assert row.mean_paired_difference == pytest.approx(-0.05)
    assert row.ci95_lower < row.mean_paired_difference < row.ci95_upper
