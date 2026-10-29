from __future__ import annotations

import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from reanalyze_bias_compensation_followup import analyze_rows


def rows():
    out = []
    truths = [60.0, 60.0, 60.0, 60.0]
    correct_errors = [1.0, -1.0, 2.0, -2.0]
    shallow_errors = [0.0, 0.0, 1.0, -1.0]
    for i, truth in enumerate(truths):
        out.append({"replicate": str(i), "arm": "correct_weighting", "h0_true": str(truth), "posterior_mean": str(truth + correct_errors[i])})
        out.append({"replicate": str(i), "arm": "shallow_slope", "h0_true": str(truth), "posterior_mean": str(truth + shallow_errors[i])})
    return out


def test_paired_bootstrap_is_deterministic_and_sign_is_documented():
    first = analyze_rows(rows(), true_h0=60.0, n_boot=1000, seed=17)
    second = analyze_rows(rows(), true_h0=60.0, n_boot=1000, seed=17)
    assert first["delta_rmse"] == pytest.approx(1.5811388300841898 - 0.7071067811865476)
    assert first["positive_delta_means"].startswith("shallow_slope")
    assert first["paired_bootstrap_95_interval"] == second["paired_bootstrap_95_interval"]
    assert first["n_paired_replicates"] == 4


def test_duplicate_and_incomplete_pairs_are_rejected():
    data = rows()
    with pytest.raises(ValueError, match="duplicate"):
        analyze_rows(data + [data[0]], true_h0=60.0, n_boot=1000)
    with pytest.raises(ValueError, match="incomplete"):
        analyze_rows(rows()[:-1], true_h0=60.0, n_boot=1000)


def test_mixed_truth_cells_are_rejected():
    data = rows()
    for row in data[-2:]:
        row["h0_true"] = "80"
        row["posterior_mean"] = str(80 + (float(row["posterior_mean"]) - 60))
    with pytest.raises(ValueError, match="true H0 varies"):
        analyze_rows(data, n_boot=1000)
