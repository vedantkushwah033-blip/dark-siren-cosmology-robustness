"""Deterministic identities for the pinned original dark-siren posterior.

These checks do not validate astrophysical realism. They verify that the
implemented weighted host mixture obeys simple algebraic identities.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-dir",default="original-project")
    a=p.parse_args()
    root=Path(a.original_dir).resolve()
    sys.path.insert(0,str(root))
    from src.bayesian_inference import posterior_h0

    grid=np.linspace(30.0,120.0,9001)
    z1,z2=0.025,0.061
    distance,sigma=245.0,28.0
    def post(zs,weights=None,missing_weight=0.0,missing_z=None,missing_weights=None):
        return np.asarray(posterior_h0(
            grid,np.asarray(zs,float),distance,sigma,weights=weights,
            missing_weight=missing_weight,
            missing_redshifts=None if missing_z is None else np.asarray(missing_z,float),
            missing_weights=None if missing_weights is None else np.asarray(missing_weights,float),
            prior_min=30.0,prior_max=120.0),float)
    def assert_close(a,b,label):
        if not np.allclose(a,b,rtol=1e-10,atol=1e-12):
            raise AssertionError(f"{label}: max absolute difference={np.max(np.abs(a-b)):.3g}")
    def assert_normalized(p,label):
        area=float(np.trapezoid(p,grid))
        if not np.isclose(area,1.0,rtol=0,atol=2e-10):
            raise AssertionError(f"{label}: integral={area}")

    # One candidate: explicit unit weight and default unit weight agree.
    one_default=post([z1])
    one_explicit=post([z1],[1.0])
    assert_close(one_default,one_explicit,"one-host/default-weight identity")
    assert_normalized(one_default,"one-host posterior")

    # A two-component catalogue must equal the same two components represented
    # as observed and missing prior mass. This checks the mixture mass algebra.
    direct=post([z1,z2],[0.7,0.3])
    split=post([z1],[1.0],missing_weight=0.3,missing_z=[z2],missing_weights=[1.0])
    assert_close(direct,split,"observed-plus-missing mixture identity")
    assert_normalized(split,"split-mixture posterior")

    # Multiplying all weights by a positive constant must not change the result.
    scaled=post([z1,z2],[7.0,3.0])
    assert_close(direct,scaled,"weight scale invariance")

    # Reordering candidate galaxies and their corresponding weights is invariant.
    reordered=post([z2,z1],[0.3,0.7])
    assert_close(direct,reordered,"candidate-order invariance")

    print("PASS: 4 posterior identities and normalization checks passed.")
if __name__=="__main__":
    main()
