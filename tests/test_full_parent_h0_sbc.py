from pathlib import Path
import json
import pandas as pd
from scripts.run_full_parent_h0_sbc import run
from scripts.audit_full_parent_h0_sbc import audit

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
