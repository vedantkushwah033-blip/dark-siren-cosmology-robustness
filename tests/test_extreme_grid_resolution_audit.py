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
