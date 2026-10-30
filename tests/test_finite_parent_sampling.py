from __future__ import annotations

import numpy as np

from dark_siren_robustness.finite_parent_sampling import run_finite_parent_benchmark


def test_finite_parent_benchmark_is_reproducible_and_accounts_for_attempts():
    kwargs = dict(
        n_scenarios=1, n_validation_draws=2, n_parent=1200,
        n_detected=6, validation_sizes=(100, 200), seed=20261033,
        h0_points=41, z_points=201,
    )
    first = run_finite_parent_benchmark(**kwargs)
    second = run_finite_parent_benchmark(**kwargs)
    assert len(first) == 4
    assert first["status"].tolist() == second["status"].tolist()
    assert first["validation_indices"].tolist() == second["validation_indices"].tolist()
    for _, row in first.iterrows():
        indices = np.fromstring(row["validation_indices"], sep=",", dtype=int)
        assert indices.size == int(row["validation_size"])
        assert np.unique(indices).size == int(row["validation_size"])
        assert np.all((indices >= 0) & (indices < 1200))
    assert (first["status"] == "success").sum() + (first["status"] == "failed").sum() == 4
    for draw in range(2):
        small = first.loc[(first["draw"] == draw) & (first["validation_size"] == 100)].iloc[0]
        large = first.loc[(first["draw"] == draw) & (first["validation_size"] == 200)].iloc[0]
        assert large["validation_indices"].split(",")[:100] == small["validation_indices"].split(",")
    for hash_column in ("parent_frame_sha256", "catalogue_sha256", "gw_distance_sha256"):
        assert first[hash_column].nunique(dropna=True) == 1
        assert first[hash_column].str.fullmatch(r"[0-9a-f]{64}").all()
    successful = first.loc[first["status"].eq("success")]
    if len(successful):
        assert successful["covered_68"].isin([True, False]).all()
        assert np.isfinite(successful["density_ise"]).all()
        assert np.isfinite(successful["weight_ess"]).all()
        assert np.isfinite(successful["finite_parent_target_posterior_mean"]).all()
        expected_error = (
            successful["posterior_mean"]
            - successful["finite_parent_target_posterior_mean"]
        )
        assert np.allclose(successful["error_vs_finite_parent_target"], expected_error)
