from __future__ import annotations
import numpy as np
from scripts.audit_host_density_campaign import (
    draw_parent_redshifts, inclusion, latent_intensity, reflected_weighted_kde,
    selected_posterior, simulate,
)


def test_independent_parent_sampler_and_full_simulation_are_seed_deterministic():
    z1=draw_parent_redshifts(100,np.random.default_rng(12345))
    z2=draw_parent_redshifts(100,np.random.default_rng(12345))
    np.testing.assert_array_equal(z1,z2)
    p1,o1,d1=simulate(20261029)
    p2,o2,d2=simulate(20261029)
    np.testing.assert_array_equal(p1,p2)
    np.testing.assert_array_equal(o1,o2)
    np.testing.assert_array_equal(d1,d2)
    assert len(p1)==6000 and len(d1)==40 and np.all(d1>=250.0)


def test_inclusion_curve_is_valid_and_monotonically_decreasing():
    z=np.linspace(0.01,0.40,401)
    p=inclusion(z)
    assert np.all((p>0)&(p<=1))
    assert np.all(np.diff(p)<0)
    assert np.all(latent_intensity(z)>0)


def test_independent_reflected_kde_and_selected_posterior_are_finite():
    z=np.linspace(0.01,0.40,401)
    sample=np.linspace(0.02,0.39,120)
    p=inclusion(sample)
    merger=1+1.5*(z-z[0])/(z[-1]-z[0])
    raw,w=reflected_weighted_kde(z,sample,p,0.2)
    density=raw*merger
    density/=np.trapezoid(density,x=z)
    h0=np.linspace(50,90,81)
    posterior=selected_posterior(z,density,np.array([300.,360.,425.]),h0)
    assert w.shape==sample.shape and np.all(np.isfinite(raw)) and np.all(raw>0)
    assert np.isclose(posterior["posterior_integral"],1.0,atol=1e-10)
    assert 50<posterior["posterior_mean"]<90
    assert 0<posterior["ci68_width"]<40
