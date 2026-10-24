import numpy as np
from scripts.run_joint_distance_model_host_mixture import run

def test_joint_experiment_smoke_is_paired_and_finite():
    rows=run(replicates=2,seed=20261025,grid_size=1001,bootstrap_draws=1000)
    assert len(rows)==3*3*3*2*4
    keys={}
    for r in rows:
        key=(r["true_H0"],r["host_z1"],r["host_z2"],r["weight1"],r["weight2"],r["replicate"])
        keys.setdefault(key,{})[r["inference_arm"]]=r
        assert np.isfinite(r["median"])
        assert 40 <= r["median"] <= 100
        assert r["lower_95"] <= r["lower_68"] <= r["median"] <= r["upper_68"] <= r["upper_95"]
    assert len(keys)==27*2
    for arms in keys.values():
        assert set(arms)=={"lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed"}
        assert len({arms[a]["observed_distance_mpc"] for a in arms})==1
        assert len({arms[a]["sigma_mpc"] for a in arms})==1

def test_distance_models_diverge_at_high_redshift():
    from dark_siren_robustness.cosmological_inference import distance_grid_for_model
    h0=[70.0]
    lcdm=distance_grid_for_model(h0,[0.4],"flat_lcdm")[0,0]
    lowz=distance_grid_for_model(h0,[0.4],"low_z")[0,0]
    assert abs(lcdm-lowz)>1.0
