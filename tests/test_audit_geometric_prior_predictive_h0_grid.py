from pathlib import Path
import sys
import zipfile

import numpy as np
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_geometric_prior_predictive_h0_grid import (
    _evaluate,
    _find_root,
    _posterior,
)


def test_posterior_rank_and_quantiles_are_normalized_on_dh0_grid():
    grid = np.geomspace(30.0, 110.0, 2241)
    log_likelihood = -0.5 * ((grid - 70.0) / 8.0) ** 2
    density, quantiles, rank = _posterior(grid, log_likelihood, 70.0)
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-10)
    assert quantiles.shape == (9,)
    assert np.all(np.diff(quantiles) > 0)
    assert 0.45 < rank < 0.55


def test_selected_likelihood_recomputation_is_finite_and_shaped():
    z = np.linspace(0.01, 0.04, 101)
    densities = np.vstack([
        np.ones_like(z) / (z[-1] - z[0]),
        (z - z[0]) / np.trapezoid(z - z[0], z),
    ])
    h0 = np.geomspace(30.0, 110.0, 51)
    result = _evaluate(
        h0,
        densities,
        z,
        np.array([350.0, 500.0]),
        sigma=140.0,
        threshold=250.0,
    )
    assert result.shape == (2, h0.size)
    assert np.all(np.isfinite(result))


def test_source_root_detection_requires_expected_replicate_path(tmp_path):
    archive_path = tmp_path / "source.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "experiment/replicates/rep000_gm1/config.json",
            "{}",
        )
    with zipfile.ZipFile(archive_path) as archive:
        assert _find_root(archive) == "experiment/"
