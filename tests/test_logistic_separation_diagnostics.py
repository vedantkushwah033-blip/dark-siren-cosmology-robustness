import numpy as np

from dark_siren_robustness.logistic_separation_diagnostics import (
    diagnose_logistic_separation,
)


def test_lp_diagnostic_detects_complete_separation():
    x = np.array([-2.0, -1.0, 1.0, 2.0])
    design = np.column_stack([np.ones(x.size), x])
    y = np.array([0, 0, 1, 1])
    result = diagnose_logistic_separation(design, y)
    assert result["separation_status"] == "complete_separation"
    assert result["complete_margin"] > result["tolerance"]
    assert result["complete_lp_status"] == 0


def test_lp_diagnostic_detects_quasi_complete_separation():
    x = np.array([-2.0, -1.0, 0.0, 0.0, 1.0, 2.0])
    design = np.column_stack([np.ones(x.size), x])
    # Opposite labels at x=0 force zero margins there; all other rows
    # can be weakly separated by the slope direction.
    y = np.array([0, 0, 0, 1, 1, 1])
    result = diagnose_logistic_separation(design, y)
    assert result["separation_status"] == "quasi_complete_separation"
    assert result["complete_margin"] <= result["tolerance"]
    assert result["quasi_score"] > result["tolerance"]


def test_lp_diagnostic_reports_overlap_for_interleaved_outcomes():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    design = np.column_stack([np.ones(x.size), x])
    y = np.array([0, 1, 0, 1, 0])
    result = diagnose_logistic_separation(design, y)
    assert result["separation_status"] == "overlap"
    assert result["complete_margin"] <= result["tolerance"]
    assert result["quasi_score"] <= result["tolerance"]


def test_rank_deficiency_is_reported_separately_from_overlap():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
    design = np.column_stack([np.ones(x.size), x, 2.0 * x])
    y = np.array([0, 1, 0, 1, 0, 1])
    result = diagnose_logistic_separation(design, y)
    assert result["separation_status"] == "overlap"
    assert result["design_rank"] == 2
    assert result["parameter_dimension"] == 3
    assert result["design_rank_deficient"] is True


def test_lp_diagnostic_is_invariant_to_nonzero_column_rescaling():
    x = np.array([-2.0, -1.0, 0.0, 0.0, 1.0, 2.0])
    design = np.column_stack([np.ones(x.size), x])
    y = np.array([0, 0, 0, 1, 1, 1])
    original = diagnose_logistic_separation(design, y)
    rescaled = diagnose_logistic_separation(
        design * np.array([17.0, 0.003]), y
    )
    assert original["separation_status"] == rescaled["separation_status"]


def test_lp_diagnostic_rejects_invalid_binary_design_inputs():
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    with np.testing.assert_raises(ValueError):
        diagnose_logistic_separation(design, [0, 0, 0, 0])
    with np.testing.assert_raises(ValueError):
        diagnose_logistic_separation(design, [0, 1, 0])
    with np.testing.assert_raises(ValueError):
        diagnose_logistic_separation(design, [0, 1, 2, 1])
