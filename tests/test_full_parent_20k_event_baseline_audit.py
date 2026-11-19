from pathlib import Path
import numpy as np
from scripts.audit_full_parent_20k_event_baseline import _draw_detected_independent, _posterior_quantiles_independent

def test_independent_detected_sampler_and_posterior():
    z=np.linspace(.02,.35,300)
    w=1+1.5*(z-.01)/.39
    rng=np.random.default_rng(91)
    d=_draw_detected_independent(rng,z,w,70.0,140.0,250.0,4)
    assert d.shape==(4,) and np.all(d>=250.0)
    q=_posterior_quantiles_independent(z,w,d,np.linspace(20,300,101),140.0,250.0)
    assert 20<=q["q025"]<=q["q05"]<=q["q16"]<=q["median"]
    assert q["median"]<=q["q84"]<=q["q95"]<=q["q975"]<=300
