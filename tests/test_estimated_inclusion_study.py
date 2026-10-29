from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot


def test_estimated_inclusion_pilot_is_reproducible_and_retains_failures():
    kwargs = dict(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100, 200),
        seed=20261031,
        h0_points=41,
        z_points=201,
    )
    first, summary = run_estimated_inclusion_pilot(**kwargs)
    second, _ = run_estimated_inclusion_pilot(**kwargs)

    assert len(first) == 4  # catalogue-only, oracle, and two estimated arms
    assert first[["replicate", "arm", "status"]].equals(
        second[["replicate", "arm", "status"]]
    )
    assert summary["status"] == "NON_CONFIRMATORY_SYNTHETIC_PILOT"
    assert {"catalogue_only", "oracle", "estimated_100", "estimated_200"} == set(first["arm"])
    assert all(first["status"].isin(["success", "failed"]))
    assert len(summary["contrasts"]) <= 3
