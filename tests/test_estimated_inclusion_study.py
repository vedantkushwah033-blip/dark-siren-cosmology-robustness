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

    assert len(first) == 6  # catalogue-only, oracle, and two model families at two sizes
    assert first["status"].eq("success").any(), first[["arm", "error"]].to_dict("records")
    assert first[["replicate", "arm", "status"]].equals(
        second[["replicate", "arm", "status"]]
    )
    assert summary["status"] == "NON_CONFIRMATORY_SYNTHETIC_PILOT"
    assert {
        "catalogue_only", "oracle",
        "estimated_linear_100", "estimated_quadratic_100",
        "estimated_linear_200", "estimated_quadratic_200",
    } == set(first["arm"])
    assert all(first["status"].isin(["success", "failed"]))
    assert len(summary["contrasts"]) <= 5
    assert len(summary["paired_model_contrasts"]) <= 2


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
        frame["arm"].str.startswith("estimated_") & frame["status"].eq("success")
    ]
    assert not estimated_successes.empty
    assert scored_sizes == [1000] * len(estimated_successes)
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
    assert set(first["arm"]) == {
        "catalogue_only", "oracle",
        "estimated_linear_150", "estimated_quadratic_150",
    }


def test_quadratic_selection_coefficient_must_be_finite():
    import pytest

    with pytest.raises(ValueError, match="coefficients must be finite"):
        run_estimated_inclusion_pilot(
            n_replicates=1, n_parent=1000, n_detected=5,
            validation_sizes=(100,), true_quadratic=float("nan"),
            h0_points=31, z_points=101,
        )



def test_positivity_and_validation_diagnostics_are_retained_for_each_arm():
    import pytest
    frame, _summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1200,
        n_detected=6,
        validation_sizes=(100, 300),
        seed=20261048,
        h0_points=41,
        z_points=201,
        true_quadratic=1.5,
    )
    diagnostics = [
        "true_probability_min", "true_probability_max", "true_probability_p01",
        "true_probability_p05", "true_probability_below_005",
        "true_probability_below_010", "true_max_inverse_probability",
        "true_retention_rate",
    ]
    assert set(diagnostics + ["audit_included_count"]).issubset(frame.columns)
    for column in diagnostics:
        assert frame[column].nunique(dropna=False) == 1
    assert frame["true_probability_min"].between(0, 1, inclusive="neither").all()
    assert frame["true_probability_max"].between(0, 1, inclusive="neither").all()
    assert frame["true_probability_below_005"].between(0, 1).all()
    assert frame["true_probability_below_010"].between(0, 1).all()
    assert (frame["true_probability_below_005"] <= frame["true_probability_below_010"]).all()
    assert frame["true_max_inverse_probability"].to_numpy() == pytest.approx(
        1.0 / frame["true_probability_min"].iloc[0]
    )
    for arm, size in (
        ("estimated_linear_100", 100), ("estimated_quadratic_100", 100),
        ("estimated_linear_300", 300), ("estimated_quadratic_300", 300),
    ):
        count = frame.loc[frame["arm"] == arm, "audit_included_count"].iloc[0]
        assert 0 <= count <= size
    assert frame.loc[frame["arm"] == "catalogue_only", "audit_included_count"].isna().all()


def test_pilot_records_preposterior_features_before_posterior_outputs():
    import pytest
    frame, _summary = run_estimated_inclusion_pilot(
        n_replicates=1, n_parent=1200, n_detected=6, validation_sizes=(150,),
        seed=20261051, h0_points=41, z_points=201,
    )
    estimated = frame.loc[frame["arm"].str.startswith("estimated_") & frame["status"].eq("success")]
    assert len(estimated) == 2, frame[["arm", "status", "error"]].to_dict("records")
    assert set(estimated["inclusion_model_degree"]) == {1, 2}
    required = {
        "diagnostic_schema", "n_validation", "validation_inclusion_fraction",
        "catalogue_outside_validation_range_fraction", "fitted_probability_p05",
        "inverse_weight_ess", "inverse_weight_ess_fraction",
        "inverse_weight_entropy_ess_fraction", "inverse_weight_max_share",
        "inverse_weight_top_1pct_mass", "inverse_weight_top_5pct_mass",
        "inverse_weight_max", "inverse_weight_p95", "inverse_weight_p99",
        "catalogue_near_validation_boundary_fraction",
    }
    assert required.issubset(estimated.columns)
    row = estimated.iloc[0]
    assert row["diagnostic_schema"] == "preposterior_weight_stability_v1"
    assert row["n_validation"] == 150
    assert 0 < row["inverse_weight_ess"] <= row["n_catalogue"]
    assert 0 < row["inverse_weight_ess_fraction"] <= 1
    assert row["inverse_weight_entropy_ess_fraction"] >= row["inverse_weight_ess_fraction"]
    assert 0 < row["inverse_weight_max_share"] <= row["inverse_weight_top_1pct_mass"]
    assert row["inverse_weight_top_1pct_mass"] <= row["inverse_weight_top_5pct_mass"] <= 1
    assert row["inverse_weight_ess"] == pytest.approx(
        row["inverse_weight_sum"] ** 2 / row["inverse_weight_sum_squares"]
    )
    assert not any("true_" in key or "posterior" in key or "coverage" in key for key in required)

def test_too_few_retained_catalogue_objects_are_recorded_without_unbound_features():
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=100,
        n_detected=1,
        validation_sizes=(20,),
        seed=20261056,
        true_intercept=-30.0,
        true_slope=0.0,
        h0_points=31,
        z_points=101,
    )
    assert set(frame["arm"]) == {
        "catalogue_only", "oracle",
        "estimated_linear_20", "estimated_quadratic_20",
    }
    assert frame["status"].eq("failed").all()
    assert frame["fit_status"].eq("not_attempted").all()
    assert frame["fit_boundary_limited"].eq(False).all()
    assert frame["error"].eq("too few retained catalogue objects").all()
    estimated = frame.loc[frame["arm"].str.startswith("estimated_")]
    assert estimated["inference_audit_payload_json"].notna().all()
    assert estimated["scoring_audit_payload_json"].notna().all()
    assert summary["status"] == "NON_CONFIRMATORY_SYNTHETIC_PILOT"

def test_pilot_rejects_numerically_saturated_true_inclusion_probabilities():
    import pytest

    with pytest.raises(ValueError, match="numerically saturated"):
        run_estimated_inclusion_pilot(
            n_replicates=1,
            n_parent=1000,
            n_detected=5,
            validation_sizes=(100,),
            seed=20261057,
            true_intercept=-1000.0,
            true_slope=0.0,
            h0_points=31,
            z_points=101,
        )


def test_pilot_reports_paired_secondary_outcomes_for_shared_replicates():
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=2,
        n_parent=1200,
        n_detected=6,
        validation_sizes=(150,),
        seed=20261067,
        h0_points=41,
        z_points=201,
        true_quadratic=1.2,
    )
    contrasts = summary["paired_secondary_contrasts"]
    assert contrasts
    keys = {
        (item["metric"], item["reference_arm"], item["comparison_arm"])
        for item in contrasts
    }
    expected = {
        (metric, "catalogue_only", arm)
        for arm in ("estimated_linear_150", "estimated_quadratic_150")
        for metric in ("coverage_68", "interval_width_68", "host_density_ise")
    }
    expected.add((
        "inclusion_probability_mae",
        "estimated_linear_150",
        "estimated_quadratic_150",
    ))
    assert keys == expected
    assert all(item["n_paired_replicates"] <= 2 for item in contrasts)
    assert all(item["difference_direction"] == "reference_minus_comparison" for item in contrasts)
    assert set(frame["arm"]) >= {
        "catalogue_only", "estimated_linear_150", "estimated_quadratic_150"
    }
