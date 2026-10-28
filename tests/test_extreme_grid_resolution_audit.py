from scripts.audit_extreme_grid_resolution import run_audit


def test_extreme_grid_audit_detects_narrow_posterior_failure():
    report = run_audit()
    rows = {row["sigma_distance_mpc"]: row for row in report["rows"]}

    assert report["not_real_gw_data"] is True
    assert rows[2.0]["status"] == "WITHIN_TOLERANCE"
    assert rows[0.02]["status"] == "UNDERRESOLVED"
    assert rows[0.01]["status"] == "UNDERRESOLVED"
    assert rows[0.02]["absolute_sd_difference"] > 0.005
    assert rows[0.01]["maximum_absolute_quantile_difference"] > 0.01
    assert report["underresolved_case_count"] >= 3


def test_extreme_grid_audit_rejects_invalid_resolution_and_tolerances():
    import pytest

    from scripts.audit_extreme_grid_resolution import run_audit

    with pytest.raises(ValueError, match="fine_points"):
        run_audit(coarse_points=4001, fine_points=4001)
    with pytest.raises(ValueError, match="tolerances"):
        run_audit(quantile_tolerance=0.0)



def test_raw_coverage_audit_separates_prior_boundaries(tmp_path):
    import csv

    from scripts.audit_full_run_coverage_heterogeneity import audit

    arms = [
        "A_fixed_realized_sigma",
        "B_fractional_sigma_H0",
        "C_fixed_reference_sigma",
    ]
    rows = []
    for arm in arms:
        for truth in (50, 55, 90):
            for replicate in range(2):
                lower = 50.1 if truth in (50, 90) else 54.0
                upper = 89.9 if truth in (50, 90) else 56.0
                rows.append({
                    "true_h0": truth,
                    "fractional_error": 0.1,
                    "replicate": replicate,
                    "arm": arm,
                    "lower_68": lower,
                    "upper_68": upper,
                    "lower_95": lower,
                    "upper_95": upper,
                })
    source = tmp_path / "raw_replicates.csv"
    with source.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    report = audit(source)
    assert report["raw_record_count"] == 18
    for arm in arms:
        assert report["overall_by_arm"][arm]["coverage_68"]["estimate"] == 1 / 3
        assert report["interior_by_arm"][arm]["coverage_68"]["estimate"] == 1.0
        assert report["boundary_by_arm"][arm]["coverage_68"]["estimate"] == 0.0
