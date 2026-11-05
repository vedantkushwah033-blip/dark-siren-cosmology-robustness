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
    assert {"catalogue_only", "oracle", "estimated_linear_100", "estimated_linear_200", "estimated_quadratic_100", "estimated_quadratic_200"} == set(first["arm"])
    assert all(first["status"].isin(["success", "failed"]))
    assert len(summary["contrasts"]) <= 5


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
        frame["arm"].eq("estimated_linear_100") & frame["status"].eq("success")
    ]
    assert not estimated_successes.empty
    assert scored_sizes == [1000, 1000]
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
    assert set(first["arm"]) == {"catalogue_only", "oracle", "estimated_linear_150", "estimated_quadratic_150"}


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
    for arm, size in ((f"estimated_{kind}_{m}", m) for kind in ("linear", "quadratic") for m in (100, 300)):
        count = frame.loc[frame["arm"] == arm, "audit_included_count"].iloc[0]
        assert 0 <= count <= size
    assert frame.loc[frame["arm"] == "catalogue_only", "audit_included_count"].isna().all()


def test_pilot_records_preposterior_features_before_posterior_outputs():
    import pytest
    frame, _summary = run_estimated_inclusion_pilot(
        n_replicates=1, n_parent=1200, n_detected=6, validation_sizes=(150,),
        seed=20261051, h0_points=41, z_points=201,
    )
    estimated = frame.loc[frame["arm"].eq("estimated_linear_150") & frame["status"].eq("success")]
    assert len(estimated) == 1, frame[["arm", "status", "error"]].to_dict("records")
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
    assert set(frame["arm"]) == {"catalogue_only", "oracle", "estimated_linear_20", "estimated_quadratic_20"}
    assert frame["status"].eq("failed").all()
    assert frame["fit_status"].eq("not_attempted").all()
    assert frame["fit_boundary_limited"].eq(False).all()
    assert frame["error"].eq("too few retained catalogue objects").all()
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


def test_validation_bootstrap_propagates_into_fixed_toy_h0_posterior():
    import numpy as np
    from scipy.special import expit
    from scipy.integrate import trapezoid
    from dark_siren_robustness.estimated_inclusion_study import (
        run_validation_bootstrap_posterior_sensitivity,
    )

    rng = np.random.default_rng(20261073)
    z_grid = np.linspace(0.01, 0.40, 101)
    parent_z = rng.uniform(z_grid[0], z_grid[-1], 2500)
    x = (parent_z - (z_grid[0] + z_grid[-1]) / 2) / ((z_grid[-1] - z_grid[0]) / 2)
    true_p = expit(0.15 - 0.8 * x)
    retained = rng.random(parent_z.size) < true_p
    observed = parent_z[retained]
    validation_idx = rng.choice(parent_z.size, 350, replace=False)

    result = run_validation_bootstrap_posterior_sensitivity(
        redshift_grid=z_grid,
        observed_redshift=observed,
        validation_redshift=parent_z[validation_idx],
        validation_included=retained[validation_idx].astype(int),
        h0_grid=np.linspace(50, 90, 51),
        observed_gw_distances_mpc=np.array([300.0, 360.0, 420.0]),
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
        n_bootstrap=20,
        seed=20261074,
    )

    assert result["schema"] == "validation_bootstrap_posterior_sensitivity_v1"
    assert result["n_attempted"] == 20
    assert result["n_successful"] + result["n_failed"] == 20
    assert result["n_successful"] >= 2
    assert result["posterior_mean_draws"].shape == (result["n_successful"],)
    assert result["posterior_median_draws"].shape == (result["n_successful"],)
    assert result["inverse_weight_ess_draws"].shape == (result["n_successful"],)
    assert result["inverse_weight_maximum_draws"].shape == (result["n_successful"],)
    assert result["inverse_weight_p99_draws"].shape == (result["n_successful"],)
    assert np.all(result["inverse_weight_ess_draws"] > 0)
    assert result["posterior_ci68_low_draws"].shape == (result["n_successful"],)
    assert result["posterior_ci68_high_draws"].shape == (result["n_successful"],)
    assert np.all(
        result["posterior_ci68_low_draws"] < result["posterior_ci68_high_draws"]
    )
    assert result["posterior_mean_bootstrap_percentile_95"] is not None
    assert "not an H0" in result["interpretation"]
    assert np.allclose(
        trapezoid(result["density_draws"], z_grid, axis=1), 1.0, atol=1e-8
    )


def test_validation_bootstrap_posterior_rejects_events_below_detection_threshold():
    import numpy as np
    import pytest
    from dark_siren_robustness.estimated_inclusion_study import (
        run_validation_bootstrap_posterior_sensitivity,
    )

    with pytest.raises(ValueError, match="satisfy the detection threshold"):
        run_validation_bootstrap_posterior_sensitivity(
            redshift_grid=np.linspace(0.01, 0.40, 101),
            observed_redshift=np.linspace(0.02, 0.39, 30),
            validation_redshift=np.linspace(0.01, 0.40, 100),
            validation_included=np.tile([0, 1], 50),
            h0_grid=np.linspace(50, 90, 51),
            observed_gw_distances_mpc=np.array([200.0]),
            sigma_mpc=140.0,
            detection_threshold_mpc=250.0,
            n_bootstrap=5,
        )


def test_quadratic_inclusion_term_is_centered_under_parent_redshift_law():
    import json
    import numpy as np
    import pytest
    from scipy.special import expit

    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1200,
        n_detected=6,
        validation_sizes=(150,),
        seed=20261104,
        h0_points=41,
        z_points=201,
        true_quadratic=1.2,
    )
    model = summary["true_inclusion_model"]
    center = model["quadratic_center"]
    assert center == pytest.approx(0.28450056, abs=2e-7)
    assert center != pytest.approx(1.0 / 3.0, abs=1e-3)
    assert frame["true_quadratic_center"].to_numpy() == pytest.approx(
        np.full(len(frame), center)
    )

    successful = frame.loc[
        frame["status"].eq("success") & frame["arm"].str.startswith("estimated_")
    ]
    assert not successful.empty
    scoring = json.loads(successful.iloc[0]["scoring_audit_payload_json"])
    parent_z = np.asarray(scoring["parent_redshift"], dtype=float)
    stored_probability = np.asarray(scoring["true_inclusion_probability"], dtype=float)
    z_min, z_max = summary["redshift_grid"][:2]
    x = (parent_z - (z_min + z_max) / 2.0) / ((z_max - z_min) / 2.0)
    reconstructed = expit(
        model["intercept"] + model["slope"] * x
        + model["quadratic"] * (x**2 - center)
    )
    assert stored_probability == pytest.approx(reconstructed)


def test_linear_and_quadratic_arms_share_paired_validation_records():
    import json

    frame, _summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1200,
        n_detected=6,
        validation_sizes=(150,),
        seed=20261105,
        h0_points=41,
        z_points=201,
    )
    successful = frame.loc[frame["status"].eq("success")].set_index("arm")
    linear = successful.loc["estimated_linear_150"]
    quadratic = successful.loc["estimated_quadratic_150"]
    linear_payload = json.loads(linear["inference_audit_payload_json"])
    quadratic_payload = json.loads(quadratic["inference_audit_payload_json"])

    assert linear_payload["fitted_inclusion_model"]["model"] == "logistic_inclusion_v1"
    assert quadratic_payload["fitted_inclusion_model"]["model"] == "quadratic_logistic_inclusion_v1"
    assert linear_payload["validation_redshift"] == quadratic_payload["validation_redshift"]
    assert linear_payload["validation_included"] == quadratic_payload["validation_included"]
    assert len(linear_payload["validation_redshift"]) == 150



def test_estimated_inclusion_pilot_runs_v2_end_to_end_and_records_solver():
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1400,
        n_detected=6,
        validation_sizes=(350,),
        seed=20261105,
        h0_points=41,
        z_points=201,
        solver_version="v2",
    )
    assert summary["solver_version"] == "v2"
    assert frame["solver_version"].eq("v2").all()
    estimated = frame.loc[
        frame["arm"].eq("estimated_linear_350") & frame["status"].eq("success")
    ]
    assert len(estimated) == 1, frame[["arm", "status", "error"]].to_dict("records")
    import json
    payload = json.loads(estimated.iloc[0]["inference_audit_payload_json"])
    assert payload["solver_version"] == "v2"
    assert payload["fitted_inclusion_model"]["model"] == "logistic_inclusion_v2"
    assert payload["fitted_inclusion_model"]["solver_version"] == "v2"



def test_v2_pilot_artifact_passes_independent_auditor(tmp_path):
    import json
    import subprocess
    import sys

    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1, n_parent=1400, n_detected=6,
        validation_sizes=(350,), seed=20261107,
        h0_points=41, z_points=201, solver_version="v2",
    )
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    (tmp_path / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    result = subprocess.run(
        [sys.executable, "scripts/audit_estimated_inclusion_pilot.py", str(tmp_path)],
        check=False, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    audit = json.loads(result.stdout)
    assert audit["status"] == "PASS"



def test_v1_v2_pilot_comparison_uses_identical_synthetic_inputs():
    """Pair both solvers on the same parent, catalogue, validation and GW data."""
    import json
    import numpy as np

    common = dict(
        n_replicates=1,
        n_parent=1200,
        n_detected=8,
        validation_sizes=(100,),
        seed=20261108,
        h0_points=41,
        z_points=201,
    )
    v1, summary_v1 = run_estimated_inclusion_pilot(**common, solver_version="v1")
    v2, summary_v2 = run_estimated_inclusion_pilot(**common, solver_version="v2")

    assert summary_v1["solver_version"] == "v1"
    assert summary_v2["solver_version"] == "v2"
    assert set(v1["arm"]) == set(v2["arm"])
    left = v1.set_index("arm").sort_index()
    right = v2.set_index("arm").sort_index()

    invariant_columns = (
        "replicate_seed", "catalogue_seed", "audit_seed",
        "n_catalogue", "n_detected", "true_quadratic_center",
        "true_probability_min", "true_probability_max",
        "true_probability_p01", "true_probability_p05",
        "true_probability_below_005", "true_probability_below_010",
        "true_max_inverse_probability", "true_retention_rate",
        "audit_included_count",
    )
    for column in invariant_columns:
        assert np.array_equal(
            left[column].to_numpy(), right[column].to_numpy(), equal_nan=True
        ), column

    for arm in left.index:
        payload_v1 = json.loads(left.loc[arm, "inference_audit_payload_json"])
        payload_v2 = json.loads(right.loc[arm, "inference_audit_payload_json"])
        for key in (
            "redshift_grid", "catalogue_redshift", "density_inclusion_probability",
            "merger_weight_grid", "validation_redshift", "validation_included",
            "h0_grid", "observed_distances_mpc",
        ):
            assert payload_v1[key] == payload_v2[key], (arm, key)
        assert payload_v1["solver_version"] == "v1"
        assert payload_v2["solver_version"] == "v2"

        scoring_v1 = json.loads(left.loc[arm, "scoring_audit_payload_json"])
        scoring_v2 = json.loads(right.loc[arm, "scoring_audit_payload_json"])
        for key in ("parent_redshift", "true_inclusion_probability"):
            assert scoring_v1.get(key) == scoring_v2.get(key), (arm, key)

    # This is a paired-input regression test, not evidence that either solver
    # improves H0 calibration; downstream performance must be estimated over
    # a prespecified repeated-sampling campaign.
