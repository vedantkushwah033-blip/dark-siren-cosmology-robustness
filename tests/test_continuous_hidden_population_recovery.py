import numpy as np

from scripts.run_continuous_hidden_population_recovery import run, summarize, paired_contrasts


def test_continuous_hidden_population_recovery_smoke_is_paired_and_finite():
    rows = run(
        replicates=2, sample_sizes=(12,), seed=911,
        h0_points=7, fraction_points=5, mu_points=5, quadrature_order=16,
    )
    assert len(rows) == 2*4
    by_rep = {}
    for row in rows:
        by_rep.setdefault(row["replicate"], set()).add(row["method"])
        assert np.isfinite(row["median"])
        assert row["lower_95"] <= row["lower_68"] <= row["median"]
        assert row["median"] <= row["upper_68"] <= row["upper_95"]
    assert all(methods == {
        "oracle", "fixed_wrong_f0.10_mu0.060",
        "plug_in_joint_medians", "fully_integrated_f_mu",
    } for methods in by_rep.values())
    summary = summarize(rows)
    assert len(summary) == 4
    assert all(row["replicates"] == 2 for row in summary)
    contrasts = paired_contrasts(rows, bootstrap_resamples=100, seed=92)
    assert len(contrasts) == 3*6
    assert all(np.isfinite(row["estimate"]) for row in contrasts)

def test_continuous_hidden_population_can_estimate_hidden_scale():
    rows = run(
        replicates=1, sample_sizes=(10,), seed=912,
        h0_points=7, fraction_points=5, mu_points=5, quadrature_order=16,
        estimate_sigma_h=True, sigma_h_points=3,
    )
    assert len(rows) == 4
    assert all(row["sigma_h_estimated"] for row in rows)
    plugin = next(row for row in rows if row["method"] == "plug_in_joint_medians")
    assert 0.006 <= plugin["plugin_sigma_h"] <= 0.020
