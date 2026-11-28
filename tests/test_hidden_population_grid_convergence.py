import numpy as np

from scripts.audit_hidden_population_grid_convergence import run


def test_continuous_hidden_grid_convergence_audit_reports_paired_changes():
    result = run(
        n_events=12, seed=91, coarse=(7, 5, 5), fine=(9, 7, 7),
        quadrature_orders=(16, 24),
    )
    assert result["n_events"] == 12
    assert set(result["comparisons"]) == {"H0", "fraction", "mu_h", "sigma_h"}
    for parameter in result["comparisons"].values():
        for comparison in parameter.values():
            assert np.isfinite(comparison["absolute_max_quantile_change"])
