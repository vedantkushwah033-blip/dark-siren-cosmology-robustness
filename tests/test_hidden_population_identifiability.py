import numpy as np

from scripts.diagnose_hidden_population_identifiability import numerical_hessian, run
from scripts.run_continuous_hidden_population_recovery import generate_dataset


def test_local_information_hessian_is_symmetric_and_finite():
    u, _, _ = generate_dataset(20, np.random.default_rng(37))
    hessian = numerical_hessian(u, quadrature_order=16)
    assert hessian.shape == (4, 4)
    assert np.all(np.isfinite(hessian))
    np.testing.assert_allclose(hessian, hessian.T, rtol=0, atol=1e-10)


def test_local_information_reports_scaling_and_curvature_status():
    result = run(replicates=1, n_events=30, seed=38, quadrature_order=16)
    assert result["parameters"] == ["H0", "fraction", "mu_h", "sigma_h"]
    assert len(result["standardized_eigenvalues"]) == 4
    assert result["status"] in {"FINITE", "NON_POSITIVE_CURVATURE"}
