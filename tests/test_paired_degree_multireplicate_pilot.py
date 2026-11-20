from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_paired_degree_multireplicate_pilot import (
    _aggregate_paired, _paired_summary, _verify_child_pair, _write_manifest,
)


def _summaries():
    rows = []
    for rep, e1, e2 in [(0, 2.0, 1.0), (1, 4.0, 2.0), (2, 0.0, -1.0)]:
        for degree, error in [(1, e1), (2, e2)]:
            rows.append({
                "replicate": rep, "gamma": 1.0, "validation_size": 100,
                "arm": "E_uncertainty_integrated_ipw", "median_error": error,
                "q16": 60.0 + error, "q84": 80.0 + error, "h0_true": 70.0,
            })
    return pd.DataFrame(rows)


def test_degree_contrasts_are_paired_and_use_d2_minus_d1():
    data = _summaries()
    left = data[data.replicate.isin([0, 1, 2])].copy()
    # Construct the two input frames in the same shape as the replicate runner.
    d1 = left[left.index % 2 == 0].copy()
    d2 = left[left.index % 2 == 1].copy()
    paired = _paired_summary(d1, d2)
    assert len(paired) == 3
    assert paired.median_error_difference_d2_minus_d1.tolist() == [-1.0, -2.0, -1.0]
    aggregate = _aggregate_paired(paired)
    row = aggregate.iloc[0]
    assert row.paired_replicates == 3
    assert row.degree1_bias == 2.0
    assert row.degree2_bias == pytest.approx(2.0 / 3.0)
    assert row.median_error_difference_d2_minus_d1_mean == pytest.approx(-4.0 / 3.0)
    assert row.median_error_difference_d2_minus_d1_mcse > 0


def test_degree_pair_rejects_changed_observable_inputs(tmp_path: Path):
    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\n1,2\n")
        (root / "observable/validation_n100.csv").write_text("id,value\n1,2\n")
        np.savez(root / "results/log_likelihoods.npz", H0_grid=np.array([30., 31.]),
                 A_catalogue_only_n100=np.array([-1., -2.]),
                 B_parent_oracle_n100=np.array([-1., -2.]),
                 C_known_probability_ipw_n100=np.array([-1., -2.]))
    assert _verify_child_pair(left, right)["common_arm_likelihoods_within_tolerance"]
    (right / "observable/gw_events_observed.csv").write_text("id,value\n1,3\n")
    with pytest.raises(AssertionError):
        _verify_child_pair(left, right)


def test_degree_pair_rejects_common_arm_likelihood_mismatch(tmp_path: Path):
    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\n1,2\n")
        (root / "observable/validation_n100.csv").write_text("id,value\n1,2\n")
    for root, value in ((left, -1.0), (right, -1.01)):
        np.savez(root / "results/log_likelihoods.npz", H0_grid=np.array([30., 31.]),
                 A_catalogue_only_n100=np.array([value, -2.]),
                 B_parent_oracle_n100=np.array([-1., -2.]),
                 C_known_probability_ipw_n100=np.array([-1., -2.]))
    with pytest.raises(AssertionError, match="common arm"):
        _verify_child_pair(left, right)



def test_independent_auditor_verifies_manifest_and_rejects_tampering(tmp_path: Path):
    import hashlib
    import json
    from audit_paired_degree_multireplicate_pilot import _verify_manifest

    root = tmp_path / "artifact"
    root.mkdir()
    payload = root / "payload.txt"
    payload.write_text("integrity-check")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    manifest = {
        "status": "PASS", "file_count": 1,
        "files": {"payload.txt": {"bytes": payload.stat().st_size, "sha256": digest}},
    }
    (root / "manifest.json").write_text(json.dumps(manifest))
    assert _verify_manifest(root) == 1
    payload.write_text("tampered")
    with pytest.raises(AssertionError, match="hash/size mismatch"):
        _verify_manifest(root)


def test_independent_auditor_rejects_manifest_path_traversal(tmp_path: Path):
    import hashlib
    import json
    from audit_paired_degree_multireplicate_pilot import _verify_manifest

    root = tmp_path / "artifact"
    root.mkdir()
    outside = tmp_path / "outside.txt"
    outside.write_text("not inside artifact")
    manifest = {
        "status": "PASS", "file_count": 1,
        "files": {"../outside.txt": {
            "bytes": outside.stat().st_size,
            "sha256": hashlib.sha256(outside.read_bytes()).hexdigest(),
        }},
    }
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(AssertionError, match="escapes artifact root"):
        _verify_manifest(root)



def test_degree_one_prior_matches_degree_two_centered_basis_marginal():
    from audit_degree_prior_parameterization import audit_centered_prior_match

    result = audit_centered_prior_match()
    assert result["status"] == "PASS"
    assert result["marginal_prior_match"] is True
    assert result["degree1_prior_covariance_alpha_beta"] == [
        [25.0, 0.0], [0.0, 6.25]
    ]
    assert result["degree2_marginal_prior_covariance_alpha_beta"] == [
        [25.0, 0.0], [0.0, 6.25]
    ]
    assert result["degree2_alpha_beta_gamma_cross_covariance"] == [
        [0.0], [0.0]
    ]
    with pytest.raises(ValueError, match="finite and positive"):
        audit_centered_prior_match(prior_sd_curvature=0.0)



def test_independent_pair_auditor_uses_recorded_check_schema(tmp_path: Path):
    from audit_paired_degree_multireplicate_pilot import _verify_pair

    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\n1,2\n")
        (root / "observable/validation_n100.csv").write_text("id,value\n1,2\n")
        np.savez(
            root / "results/log_likelihoods.npz",
            H0_grid=np.array([30.0, 31.0]),
            A_catalogue_only_n100=np.array([-1.0, -2.0]),
            B_parent_oracle_n100=np.array([-1.0, -2.0]),
            C_known_probability_ipw_n100=np.array([-1.0, -2.0]),
        )
    result = _verify_pair(left, right)
    assert result == {
        "input_files_compared": 6,
        "common_arm_likelihoods_within_tolerance": True,
        "maximum_common_arm_log_likelihood_difference": 0.0,
    }



def test_independent_auditor_reconstructs_and_rejects_aggregate_tampering():
    from audit_paired_degree_multireplicate_pilot import _verify_aggregate_metrics

    data = _summaries()
    d1 = data[data.index % 2 == 0].copy()
    d2 = data[data.index % 2 == 1].copy()
    paired = _paired_summary(d1, d2)
    aggregate = _aggregate_paired(paired)

    assert _verify_aggregate_metrics(paired, aggregate) == 1

    tampered = aggregate.copy()
    tampered.loc[0, "degree2_rmse"] += 0.25
    with pytest.raises(AssertionError):
        _verify_aggregate_metrics(paired, tampered)


def test_independent_auditor_rejects_missing_aggregate_group():
    from audit_paired_degree_multireplicate_pilot import _verify_aggregate_metrics

    data = _summaries()
    d1 = data[data.index % 2 == 0].copy()
    d2 = data[data.index % 2 == 1].copy()
    paired = _paired_summary(d1, d2)
    aggregate = _aggregate_paired(paired)
    with pytest.raises(AssertionError):
        _verify_aggregate_metrics(paired, aggregate.iloc[0:0])


def test_independent_auditor_rejects_tampered_paired_interval():
    from audit_paired_degree_multireplicate_pilot import _verify_aggregate_metrics

    data = _summaries()
    d1 = data[data.index % 2 == 0].copy()
    d2 = data[data.index % 2 == 1].copy()
    paired = _paired_summary(d1, d2)
    aggregate = _aggregate_paired(paired)
    tampered = aggregate.copy()
    tampered.loc[0, "squared_error_difference_d2_minus_d1_ci95_high"] += 0.1
    with pytest.raises(AssertionError):
        _verify_aggregate_metrics(paired, tampered)


def test_paired_degree_checks_reject_asymmetric_validation_file_sets(tmp_path: Path):
    from audit_paired_degree_multireplicate_pilot import _verify_pair

    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\\n1,2\\n")
        (root / "observable/validation_n100.csv").write_text("id,value\\n1,2\\n")
        np.savez(
            root / "results/log_likelihoods.npz",
            H0_grid=np.array([30.0, 31.0]),
            A_catalogue_only_n100=np.array([-1.0, -2.0]),
            B_parent_oracle_n100=np.array([-1.0, -2.0]),
            C_known_probability_ipw_n100=np.array([-1.0, -2.0]),
        )

    (right / "observable/validation_n100.csv").unlink()
    with pytest.raises(AssertionError, match="validation-file sets differ"):
        _verify_pair(left, right)
    from run_paired_degree_multireplicate_pilot import _verify_child_pair
    with pytest.raises(AssertionError, match="validation-file sets differ"):
        _verify_child_pair(left, right)

    (right / "observable/validation_n100.csv").write_text("id,value\\n1,2\\n")
    (right / "observable/validation_n300.csv").write_text("id,value\\n1,2\\n")
    with pytest.raises(AssertionError, match="validation-file sets differ"):
        _verify_pair(left, right)
    with pytest.raises(AssertionError, match="validation-file sets differ"):
        _verify_child_pair(left, right)



def test_outer_manifest_preserves_incomplete_lifecycle_status(tmp_path: Path):
    import json

    payload = tmp_path / "partial.csv"
    payload.write_text("value\\n1\\n")
    manifest = _write_manifest(tmp_path, "INCOMPLETE")
    saved = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE"
    assert saved["status"] == "INCOMPLETE"
    assert saved["file_count"] == 1


def test_paired_weight_diagnostic_aggregates_are_independently_reconstructed():
    from run_paired_degree_multireplicate_pilot import (
        WEIGHT_DIAGNOSTIC_METRICS, _aggregate_weight_diagnostics,
    )
    from audit_paired_degree_multireplicate_pilot import (
        _verify_weight_diagnostic_outputs,
    )

    rows = []
    for degree in (1, 2):
        for replicate, offset in ((0, 0.0), (1, 2.0)):
            row = {
                "degree": degree, "replicate": replicate, "gamma": 1.0,
                "validation_size": 100, "arm": "D_map_plugin_ipw",
            }
            for index, metric in enumerate(WEIGHT_DIAGNOSTIC_METRICS):
                row[metric] = 10.0 * degree + offset + index
            rows.append(row)
    raw = pd.DataFrame(rows)
    aggregate = _aggregate_weight_diagnostics(raw)

    assert _verify_weight_diagnostic_outputs(raw, raw.copy(), aggregate) == (
        len(raw), len(aggregate)
    )

    tampered_aggregate = aggregate.copy()
    tampered_aggregate.loc[0, "mean"] += 0.5
    with pytest.raises(AssertionError):
        _verify_weight_diagnostic_outputs(raw, raw.copy(), tampered_aggregate)

    tampered_raw = raw.copy()
    tampered_raw.loc[0, WEIGHT_DIAGNOSTIC_METRICS[0]] += 0.5
    with pytest.raises(AssertionError):
        _verify_weight_diagnostic_outputs(raw, tampered_raw, aggregate)

