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
    assert first["status"].eq("success").any(), first[["arm", "error"]].to_dict("records")
    assert first[["replicate", "arm", "status"]].equals(
        second[["replicate", "arm", "status"]]
    )
    assert summary["status"] == "NON_CONFIRMATORY_SYNTHETIC_PILOT"
    assert {"catalogue_only", "oracle", "estimated_100", "estimated_200"} == set(first["arm"])
    assert all(first["status"].isin(["success", "failed"]))
    assert len(summary["contrasts"]) <= 3


def test_estimated_inclusion_pilot_rejects_non_integer_validation_sizes():
    import pytest

    for invalid_sizes in ((100.5,), (100.0,), (True,), ("100",)):
        with pytest.raises(ValueError, match="integer counts"):
            run_estimated_inclusion_pilot(
                n_replicates=1,
                n_parent=1000,
                n_detected=5,
                validation_sizes=invalid_sizes,
                h0_points=31,
                z_points=101,
            )


def test_inclusion_probability_mae_is_scored_on_full_parent_population(monkeypatch):
    import dark_siren_robustness.estimated_inclusion_study as study

    original = study.predict_inclusion_probability
    scored_sizes = []

    def tracking_predict(redshift, model):
        scored_sizes.append(len(redshift))
        return original(redshift, model)

    monkeypatch.setattr(study, "predict_inclusion_probability", tracking_predict)
    frame, _summary = study.run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100,),
        seed=20261043,
        h0_points=41,
        z_points=201,
    )
    estimated_successes = frame.loc[
        frame["arm"].eq("estimated_100") & frame["status"].eq("success")
    ]
    assert not estimated_successes.empty
    assert scored_sizes == [1000]
    assert estimated_successes["inclusion_probability_mae"].notna().all()


def test_quadratic_selection_stress_is_recorded_and_reproducible():
    import pytest

    kwargs = dict(
        n_replicates=2,
        n_parent=1200,
        n_detected=6,
        validation_sizes=(150,),
        seed=20261044,
        h0_points=41,
        z_points=201,
        true_quadratic=1.8,
    )
    first, summary = run_estimated_inclusion_pilot(**kwargs)
    second, _ = run_estimated_inclusion_pilot(**kwargs)
    assert summary["true_inclusion_model"]["family"] == "quadratic_logistic"
    assert summary["true_inclusion_model"]["quadratic"] == pytest.approx(1.8)
    assert first[["replicate", "arm", "status"]].equals(
        second[["replicate", "arm", "status"]]
    )
    assert set(first["arm"]) == {"catalogue_only", "oracle", "estimated_150"}


def test_quadratic_selection_coefficient_must_be_finite():
    import pytest

    with pytest.raises(ValueError, match="coefficients must be finite"):
        run_estimated_inclusion_pilot(
            n_replicates=1, n_parent=1000, n_detected=5,
            validation_sizes=(100,), true_quadratic=float("nan"),
            h0_points=31, z_points=101,
        )
