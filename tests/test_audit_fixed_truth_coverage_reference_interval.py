import csv
import math

import pytest

from scripts.audit_fixed_truth_coverage_reference_interval import audit


def write_fixture(path, rows):
    fields = [
        "true_h0", "fractional_error", "replicate", "arm", "host_redshift",
        "observed_distance_mpc", "generated_sigma_mpc",
        "lower_68", "upper_68", "lower_95", "upper_95",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fixture_rows():
    c = 299_792.458
    truth, redshift = 70.0, 0.02
    distance = c * redshift / truth
    sigma = 0.1 * distance
    base = {
        "true_h0": truth,
        "fractional_error": 0.1,
        "arm": "A_fixed_realized_sigma",
        "host_redshift": redshift,
        "generated_sigma_mpc": sigma,
        "lower_68": 60.0,
        "upper_68": 80.0,
        "lower_95": 50.0,
        "upper_95": 90.0,
    }
    first = {**base, "replicate": 0, "observed_distance_mpc": distance}
    second = {**base, "replicate": 1, "observed_distance_mpc": distance + 2.5 * sigma}
    return [first, second]


def test_reference_interval_audit_distinguishes_bayesian_and_frequentist_coverage(tmp_path):
    source = tmp_path / "replicates.csv"
    write_fixture(source, fixture_rows())

    report = audit(source)

    assert report["status"] == "PASS"
    assert report["arm_a_record_count"] == 2
    interior68 = report["groups"]["interior_truths_55_to_85"]["68"]
    assert interior68["bayesian_equal_tail"]["successes"] == 2
    assert interior68["inverted_gaussian_distance"]["successes"] == 1
    interior95 = report["groups"]["interior_truths_55_to_85"]["95"]
    assert interior95["bayesian_equal_tail"]["successes"] == 2
    assert interior95["inverted_gaussian_distance"]["successes"] == 1


def test_reference_interval_audit_rejects_duplicate_record_keys(tmp_path):
    source = tmp_path / "duplicates.csv"
    rows = fixture_rows()
    rows[1]["replicate"] = 0
    write_fixture(source, rows)

    with pytest.raises(ValueError, match="duplicate record key"):
        audit(source)


def test_reference_interval_audit_rejects_reversed_posterior_interval(tmp_path):
    source = tmp_path / "reversed.csv"
    rows = fixture_rows()
    rows[0]["lower_68"] = 81.0
    write_fixture(source, rows)

    with pytest.raises(ValueError, match="reversed 68% posterior interval"):
        audit(source)
