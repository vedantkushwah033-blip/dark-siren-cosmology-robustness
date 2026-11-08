import json

import pandas as pd
import pytest

from dark_siren_robustness.fixed_analysis_catalogue_analysis import analyze_fixed_catalogue_attempts
from dark_siren_robustness.fixed_analysis_catalogue_analysis_audit import audit_confirmatory_analysis


def _fixture():
    rows = []
    for rep in range(5):
        for gamma in (-1.2, 0.0, 1.2):
            for nval in (100, 300, 1000):
                base = 1.0 + 0.01 * rep + 0.02 * gamma
                for arm, degree, delta in [
                    ("catalogue_only", 0, 0.0),
                    ("known_inclusion_ipw", 0, 0.02),
                    ("parent_sample_kde", 0, -0.03),
                    ("estimated_degree_1_ipw", 1,
                     -0.1 + 0.01 * rep if gamma == -1.2 and nval == 300 else
                     0.2 + 0.02 * rep if gamma == 1.2 and nval == 300 else -0.03),
                    ("estimated_degree_2_ipw", 2, 0.01),
                ]:
                    rows.append({
                        "replicate": rep, "gamma": gamma, "target_retention": 0.25,
                        "n_validation": nval, "bandwidth_multiplier": 0.15,
                        "degree": degree, "arm": arm, "status": "SUCCESS",
                        "iae": base + delta,
                    })
    return pd.DataFrame(rows)


def test_independent_analysis_auditor_reconstructs_all_contrasts():
    attempts = _fixture()
    reported = analyze_fixed_catalogue_attempts(
        attempts, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
    )
    audit = audit_confirmatory_analysis(
        attempts, reported, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
    )
    assert audit["status"] == "PASS"
    assert audit["attempt_count"] == 225
    assert audit["primary_complete_pairs"] == 5
    assert audit["secondary_contrasts_reconstructed"] == 18


def test_independent_analysis_auditor_rejects_tampered_primary_mean():
    attempts = _fixture()
    reported = analyze_fixed_catalogue_attempts(
        attempts, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
    )
    reported["primary_interaction"]["mean_interaction"] += 0.01
    with pytest.raises(ValueError, match="primary.mean_interaction"):
        audit_confirmatory_analysis(
            attempts, reported, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
        )


def test_independent_analysis_auditor_rejects_tampered_holm_adjustment():
    attempts = _fixture()
    reported = analyze_fixed_catalogue_attempts(
        attempts, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
    )
    reported["secondary_contrasts"][0]["holm_adjusted_p_value"] = 0.999
    with pytest.raises(ValueError, match="holm_adjusted_p_value"):
        audit_confirmatory_analysis(
            attempts, reported, planned_replicates=5, bootstrap_resamples=1000, bootstrap_seed=3
        )
