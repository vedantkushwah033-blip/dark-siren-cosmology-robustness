from __future__ import annotations

import numpy as np

from dark_siren_robustness.finite_parent_sampling import run_finite_parent_benchmark


def test_finite_parent_benchmark_is_reproducible_and_accounts_for_attempts():
    kwargs = dict(
        n_scenarios=1, n_validation_draws=2, n_parent=1200,
        n_detected=6, validation_sizes=(200,), seed=20261033,
        h0_points=41, z_points=201,
    )
    first = run_finite_parent_benchmark(**kwargs)
    second = run_finite_parent_benchmark(**kwargs)
    assert len(first) == 2
    assert first["status"].tolist() == second["status"].tolist()
    assert first["validation_indices"].tolist() == second["validation_indices"].tolist()
    for encoded in first["validation_indices"]:
        indices = np.fromstring(encoded, sep=",", dtype=int)
        assert indices.size == 200
        assert np.unique(indices).size == 200
        assert np.all((indices >= 0) & (indices < 1200))
    assert (first["status"] == "success").sum() + (first["status"] == "failed").sum() == 2
    for hash_column in ("parent_frame_sha256", "catalogue_sha256", "gw_distance_sha256"):
        assert first[hash_column].nunique(dropna=True) == 1
        assert first[hash_column].str.fullmatch(r"[0-9a-f]{64}").all()
    successful = first.loc[first["status"].eq("success")]
    if len(successful):
        assert successful["covered_68"].isin([True, False]).all()
        assert np.isfinite(successful["density_ise"]).all()
        assert np.isfinite(successful["weight_ess"]).all()
