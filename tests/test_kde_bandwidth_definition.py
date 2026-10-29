import numpy as np
from scipy.integrate import trapezoid
from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.generative_selection import latent_galaxy_intensity, generate_synthetic_population

Z=np.linspace(.01,.40,801)


def test_bandwidth_is_relative_to_original_sample_not_reflected_sample():
    truth=latent_galaxy_intensity(Z)
    truth=truth/trapezoid(truth,Z)
    sample=generate_synthetic_population(seed=3,n_parent=6000,n_detected=3)["parent_truth"]["redshift"].to_numpy()
    e_small=trapezoid(np.abs(estimate_host_density(Z,sample,bandwidth_factor=.1)-truth),Z)
    e_large=trapezoid(np.abs(estimate_host_density(Z,sample,bandwidth_factor=1.2)-truth),Z)
    uniform=trapezoid(np.abs(1/(Z[-1]-Z[0])-truth),Z)
    assert e_small<.06 and e_large>.9*uniform and e_small<e_large


def test_matches_direct_reflected_kernel_sum():
    rng=np.random.default_rng(0); sample=rng.uniform(.02,.39,50); weights=rng.uniform(1,3,50)
    p=1/weights
    h=.2*np.sqrt(np.sum(weights*(sample-np.sum(weights*sample)/weights.sum())**2)/weights.sum())
    direct=np.zeros_like(Z)
    for centers in (sample,2*Z[0]-sample,2*Z[-1]-sample):
        direct+=(weights[None,:]*np.exp(-.5*((Z[:,None]-centers[None,:])/h)**2)).sum(1)
    direct/=trapezoid(direct,Z)
    actual=estimate_host_density(Z,sample,inclusion_probability=p,bandwidth_factor=.2)
    assert np.allclose(actual,direct,rtol=1e-10,atol=1e-12)


def test_density_is_normalized_and_deterministic():
    sample=generate_synthetic_population(seed=3,n_parent=6000,n_detected=3)["parent_truth"]["redshift"].to_numpy()
    first=estimate_host_density(Z,sample); second=estimate_host_density(Z,sample)
    assert np.array_equal(first,second) and abs(trapezoid(first,Z)-1)<1e-10
