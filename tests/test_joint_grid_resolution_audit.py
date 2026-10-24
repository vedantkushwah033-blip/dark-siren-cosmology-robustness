from scripts.run_joint_grid_resolution_audit import run, summarize

def test_grid_resolution_audit_pairs_same_observations():
    rows=run(replicates=1,grids=(1001,2001))
    assert len(rows)==27*1*2*4
    groups={}
    for r in rows:
        key=(r["true_H0"],r["host_z1"],r["host_z2"],r["weight1"],r["weight2"],r["replicate"],r["arm"])
        groups.setdefault(key,{})[r["grid_size"]]=r
    assert len(groups)==27*4
    assert all(set(v)=={1001,2001} for v in groups.values())
    diffs=summarize(rows,grids=(1001,2001))
    assert len(diffs)==len(rows)
    assert all(r["median_abs_diff_vs_reference"]>=0 for r in diffs)
