import numpy as np

from scripts.validate_weight_diagnostic_holdout import (
    _auc, _stratified_auc_bootstrap_ci, evaluate_frozen_rule,
)


def test_auc_matches_perfect_and_reversed_rankings():
    labels = np.array([0, 0, 1, 1])
    assert _auc(np.array([0.1, 0.2, 0.8, 0.9]), labels) == 1.0
    assert _auc(np.array([0.9, 0.8, 0.2, 0.1]), labels) == 0.0


def test_frozen_rule_learns_positive_direction_and_evaluates_holdout():
    train_x = np.arange(100, dtype=float)
    train_error = train_x + 0.01 * np.sin(train_x)
    test_x = np.arange(100, dtype=float) + 0.2
    test_error = test_x + 0.01 * np.cos(test_x)

    result = evaluate_frozen_rule(train_x, train_error, test_x, test_error)

    assert result["status"] == "PASS"
    assert result["risk_direction"] == "higher diagnostic"
    assert result["holdout_auroc"] > 0.99
    assert result["holdout_balanced_accuracy"] > 0.9
    assert sum(result["holdout_confusion_tp_fp_fn_tn"]) == 100


def test_frozen_rule_learns_negative_direction_without_reversing_holdout():
    train_x = np.arange(100, dtype=float)
    train_error = -train_x + 0.01 * np.sin(train_x)
    test_x = np.arange(100, dtype=float) + 0.2
    test_error = -test_x + 0.01 * np.cos(test_x)

    result = evaluate_frozen_rule(train_x, train_error, test_x, test_error)

    assert result["status"] == "PASS"
    assert result["risk_direction"] == "lower diagnostic"
    assert result["holdout_auroc"] > 0.99
    assert result["holdout_balanced_accuracy"] > 0.9


def test_frozen_rule_marks_weak_training_direction_non_evaluable():
    x = np.arange(100, dtype=float)
    error = np.tile([0.0, 1.0], 50)

    result = evaluate_frozen_rule(x, error, x[::-1], error)

    assert result["status"] == "NON_EVALUABLE"
    assert result["reason"] == "training direction too weak"


def test_frozen_rule_rejects_nonfinite_data():
    x = np.arange(100, dtype=float)
    error = x.copy()
    error[0] = np.nan

    try:
        evaluate_frozen_rule(x, error, x, x)
    except ValueError as exc:
        assert "non-finite" in str(exc)
    else:
        raise AssertionError("expected non-finite inputs to be rejected")



def test_stratified_auc_bootstrap_is_reproducible_and_bounded():
    labels = np.array([0] * 40 + [1] * 30)
    scores = np.concatenate([np.linspace(0.0, 0.8, 40), np.linspace(0.2, 1.0, 30)])
    first = _stratified_auc_bootstrap_ci(scores, labels, seed=17, replicates=500)
    second = _stratified_auc_bootstrap_ci(scores, labels, seed=17, replicates=500)
    assert first == second
    assert 0.0 <= first[0] <= first[1] <= 1.0
