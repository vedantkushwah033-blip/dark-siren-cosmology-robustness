import numpy as np

from dark_siren_robustness.estimated_inclusion_benchmark import (
    run_probability_validation_benchmark,
)


def test_probability_benchmark_preserves_full_paired_denominator():
    result = run_probability_validation_benchmark(
        n_parent=500,
        validation_sizes=(100, 200),
        curvature_values=(0.0, 1.2),
        degrees=(1, 2),
        replicates=2,
        seed=17,
    )
    assert len(result) == 2 * 2 * 2 * 2
    assert result["fit_status"].isin(
        ["SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"]
    ).all()
    assert result.groupby(
        ["gamma", "replicate", "n_validation"]
    )["model_degree"].nunique().eq(2).all()
    successful = result.loc[result["fit_status"] == "SUCCESS"]
    assert len(successful) > 0
    assert {"correctly_specified", "deliberately_misspecified"}.issubset(
        set(result["model_specification"])
    )
    assert np.all(np.isfinite(successful["heldout_brier_score"]))
    assert np.all(successful["heldout_brier_score"] >= 0)
    assert np.all(successful["kish_ess"] > 0)


def test_probability_benchmark_rejects_invalid_validation_sizes():
    import pytest
    with pytest.raises(ValueError):
        run_probability_validation_benchmark(
            n_parent=100, validation_sizes=(100,), replicates=1
        )
