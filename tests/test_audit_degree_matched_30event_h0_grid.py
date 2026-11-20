import numpy as np
import pytest
from scipy.special import logsumexp

from scripts.audit_degree_matched_30event_h0_grid import (
    _posterior,
    _selected_loglike,
)


def test_posterior_normalizes_and_returns_truth_rank():
    grid = np.linspace(30.0, 110.0, 281)
    loglike = -0.5 * ((grid - 70.0) / 5.0) ** 2
    density, quantiles, rank = _posterior(grid, loglike, 70.0)
    assert np.isclose(np.trapz(density, grid), 1.0)
    assert quantiles[2] == pytest.approx(70.0, abs=0.01)
    assert rank == pytest.approx(0.5, abs=0.01)


def test_selected_loglike_is_finite_and_shape_stable():
    z = np.linspace(0.01, 0.40, 201)
    density = np.exp(-0.5 * ((z - 0.2) / 0.06) ** 2)
    density /= np.trapz(density, z)
    grid = np.linspace(30.0, 110.0, 41)
    result = _selected_loglike(
        grid, density[None, :], z, np.array([300.0, 500.0]),
        sigma=140.0, threshold=250.0,
    )
    assert result.shape == (1, 41)
    assert np.isfinite(result).all()


def test_degree_two_e_nodes_can_be_marginalized():
    grid = np.linspace(30.0, 110.0, 41)
    z = np.linspace(0.01, 0.40, 201)
    base = np.exp(-0.5 * ((z - 0.2) / 0.06) ** 2)
    base /= np.trapz(base, z)
    densities = np.vstack([base, base * (1 + 0.1 * (z - 0.2))])
    densities /= np.trapz(densities, z, axis=1)[:, None]
    node_loglikes = _selected_loglike(
        grid, densities, z, np.array([400.0]), 140.0, 250.0
    )
    weights = np.array([0.25, 0.75])
    marginal = logsumexp(node_loglikes + np.log(weights[:, None]), axis=0)
    assert marginal.shape == (41,)
    assert np.isfinite(marginal).all()

from pathlib import Path
import zipfile

from scripts.audit_degree_matched_30event_h0_grid import _bundle_prefix


def test_bundle_prefix_resolves_downloaded_artifact_root_layout(tmp_path: Path):
    archive_path = tmp_path / "source.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("gamma_0/manifest.json", "{}")
    with zipfile.ZipFile(archive_path) as archive:
        assert _bundle_prefix(archive, "gamma_0") == "gamma_0/"


def test_bundle_prefix_accepts_nested_artifact_layout(tmp_path: Path):
    archive_path = tmp_path / "source_nested.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("degree_quadrature/gamma_p1/manifest.json", "{}")
    with zipfile.ZipFile(archive_path) as archive:
        assert _bundle_prefix(archive, "gamma_p1") == "degree_quadrature/gamma_p1/"


def test_bundle_prefix_fails_closed_when_manifest_is_missing(tmp_path: Path):
    archive_path = tmp_path / "missing.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("other/file.txt", "not a bundle")
    with zipfile.ZipFile(archive_path) as archive:
        with pytest.raises(FileNotFoundError):
            _bundle_prefix(archive, "gamma_m1")

from scripts.audit_degree_matched_30event_h0_grid import _verify_degree_pairing


def _write_paired_archive(path, *, mismatch=False):
    paths = (
        "observable/catalogue_observed.csv",
        "observable/gw_events_observed.csv",
        "observable/validation_n100.csv",
        "observable/validation_n500.csv",
        "truth/parent_truth.csv",
        "truth/gw_event_truth.csv",
        "truth/validation_frame_truth.csv",
    )
    with zipfile.ZipFile(path, "w") as archive:
        for degree in (1, 2):
            for rel in paths:
                payload = b"same-input"
                if mismatch and degree == 2 and rel == paths[0]:
                    payload = b"changed-input"
                archive.writestr(f"gamma_0/degree{degree}_order9/{rel}", payload)
            config = {
                "seed": 11, "validation_seed": 12, "n_parent": 1200,
                "n_detected": 30, "validation_frame": 500,
                "validation_sizes": [100, 500], "n_catalogue_observed": 600,
                "curvature_gamma": 0.0, "curvature_beta": -3.9,
                "quadratic_calibration": {"gamma": 0.0}, "h0_true": 70.0,
                "h0_prior_support": [30.0, 110.0], "h0_grid_points": 281,
                "scott_multiplier": 0.7, "sigma_mpc": 140.0,
                "detection_threshold_mpc": 250.0,
            }
            archive.writestr(
                f"gamma_0/degree{degree}_order9/config.json",
                __import__("json").dumps(config),
            )


def test_degree_pairing_checks_all_shared_inputs(tmp_path: Path):
    archive_path = tmp_path / "paired.zip"
    _write_paired_archive(archive_path)
    with zipfile.ZipFile(archive_path) as archive:
        assert _verify_degree_pairing(archive, "gamma_0/") == 7


def test_degree_pairing_rejects_changed_raw_input(tmp_path: Path):
    archive_path = tmp_path / "unpaired.zip"
    _write_paired_archive(archive_path, mismatch=True)
    with zipfile.ZipFile(archive_path) as archive:
        with pytest.raises(ValueError, match="degree-paired input mismatch"):
            _verify_degree_pairing(archive, "gamma_0/")


def test_selected_loglike_matches_direct_selected_integral():
    from scipy.integrate import trapezoid
    from scipy.special import log_ndtr
    from scipy.stats import norm

    z = np.linspace(0.01, 0.40, 401) ** 1.03
    density = np.exp(-0.5 * ((z - 0.20) / 0.055) ** 2)
    density /= trapezoid(density, z)
    distances = np.array([410.0, 520.0])
    sigma, threshold, h0 = 135.0, 245.0, 72.0
    actual = _selected_loglike(
        np.array([h0]), density[None, :], z, distances, sigma, threshold
    )[0, 0]

    means = 299792.458 * z / h0
    efficiency = trapezoid(norm.cdf((means - threshold) / sigma) * density, z)
    expected = 0.0
    for distance in distances:
        numerator = trapezoid(norm.pdf(distance, loc=means, scale=sigma) * density, z)
        expected += np.log(numerator) - np.log(efficiency)
    assert actual == pytest.approx(expected, rel=1e-11, abs=1e-11)


def test_posterior_normalization_respects_nonuniform_h0_measure():
    grid = np.geomspace(30.0, 110.0, 1121)
    loglike = -0.5 * ((grid - 72.0) / 4.0) ** 2
    density, quantiles, rank = _posterior(grid, loglike, 72.0)
    from scipy.integrate import trapezoid
    assert trapezoid(density, grid) == pytest.approx(1.0, abs=1e-12)
    assert quantiles[2] == pytest.approx(72.0, abs=0.01)
    assert 0.49 < rank < 0.51
