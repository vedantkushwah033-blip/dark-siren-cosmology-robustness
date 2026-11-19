from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hashlib
import json
import pandas as pd
import pytest
from run_full_parent_h0_sbc import run
from audit_full_parent_h0_sbc import audit

def test_small_full_parent_sbc_roundtrip(tmp_path: Path):
    root=tmp_path/"sbc"
    frame,summary=run(root,replicates=20,master_seed=91357,n_parent=300,
        n_detected=5,h0_min=30,h0_max=120,h0_points=101)
    assert len(frame)==20
    assert frame.posterior_rank.between(0,1).all()
    assert frame.n_detected.eq(5).all()
    assert 0 <= summary["mean_rank"] <= 1
    report=audit(root,expected_replicates=20)
    assert report["status"]=="PASS"
    assert report["posterior_rank_and_quantile_rows_reconstructed"]==20


def _refresh_source_manifest(root: Path):
    manifest = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in ("sbc_replicates.csv", "sbc_summary.json", "design.json")
    }
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True))


def test_sbc_auditor_rejects_tampered_coverage_mcse_with_refreshed_hash(tmp_path: Path):
    root = tmp_path / "sbc"
    run(root, replicates=20, master_seed=91358, n_parent=300,
        n_detected=5, h0_min=30, h0_max=120, h0_points=101)
    summary_path = root / "sbc_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["coverage68_mcse"] += 0.01
    summary_path.write_text(json.dumps(summary, sort_keys=True))
    _refresh_source_manifest(root)
    with pytest.raises(AssertionError, match="coverage68 MCSE mismatch"):
        audit(root, expected_replicates=20)


def test_sbc_auditor_rejects_missing_replicate_even_if_row_count_is_preserved(tmp_path: Path):
    root = tmp_path / "sbc"
    run(root, replicates=20, master_seed=91359, n_parent=300,
        n_detected=5, h0_min=30, h0_max=120, h0_points=101)
    raw_path = root / "sbc_replicates.csv"
    frame = pd.read_csv(raw_path)
    frame.loc[frame.index[-1], "replicate"] = 0
    frame.to_csv(raw_path, index=False)
    _refresh_source_manifest(root)
    with pytest.raises(AssertionError, match="replicate inventory"):
        audit(root, expected_replicates=20)
