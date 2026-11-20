# 30-event H0-grid refinement: independent readout

**Date:** 2026-11-21  
**Status:** independent numerical reconstruction completed locally; GitHub Actions validation of the PR implementation is still a separate pending check.  
**Scope:** exact H0 integration conditional on the archived order-9 KDE density arrays and coefficient-posterior nodes. This is not an end-to-end re-generation, calibration study, SBC, estimator ranking, or observational inference.

## Frozen source and independent checks

- Source experiment: PR #141, workflow `44806809662`.
- Source artifact: `degree-matched-30event-quadrature-sensitivity`, artifact ID `14669447416`.
- Downloaded ZIP SHA-256: `d0a8817b0716b2e6168ad4d69ad27b939d1359033fe49e6567b75c75077b77d0`.
- Design: three curvature bundles; degrees 1 and 2; order-9 KDEs; validation sizes 100 and 500; 30 detected events; five inference arms.
- I independently rechecked all 342 manifest-listed payload files (114 per curvature bundle) against declared byte counts and SHA-256 hashes.
- I independently confirmed the seven raw observable/truth inputs are byte-identical between degree-one and degree-two order-9 children in each curvature bundle (21/21 comparisons), and the paired invariant configuration fields agree.
- I independently reconstructed the archived 281-point selected log likelihoods for all 60 posterior curves (3 curvature × 2 degrees × 2 validation sizes × 5 arms). Maximum absolute discrepancy was **1.14 × 10⁻¹³**, far below the predeclared 10⁻⁸ reconstruction tolerance.

## H0-grid refinement method

For every curve, I recomputed all 30 event numerator terms and the selection-efficiency denominator from the archived redshift-density arrays at:

- 281 uniformly spaced H0 points on [30, 110];
- 561 geometrically spaced points on [30, 110];
- 1,121 geometrically spaced points on [30, 110].

For arm E, I evaluated the complete multi-event likelihood at each saved coefficient-posterior node and then marginalized the joint likelihood using the archived shared posterior quadrature weights. Posterior normalization and CDF integration used **dH0**, consistent with the uniform-H0 prior; geometric spacing was only the evaluation mesh.

The 281→561 and 561→1,121 comparisons were computed for all five arms, not only the integrated-IPW arm.

## Results across all 60 posterior curves

| Metric (maximum absolute change) | 281 → 561 | 561 → 1,121 |
|---|---:|---:|
| Posterior rank at H0 truth | 0.00005547 | 0.00005381 |
| Posterior median (km s⁻¹ Mpc⁻¹) | 0.0006487 | 0.0001689 |
| Largest shift among 2.5%, 16%, 50%, 84%, 97.5% quantiles (km s⁻¹ Mpc⁻¹) | 0.008271 | 0.004041 |
| Posterior-density L1 distance | 0.0004918 | 0.0001933 |

The largest 281→561 quantile shift occurred for degree 1, gamma=+1, n=500, MAP-IPW arm D. The largest 561→1,121 quantile shift occurred for degree 2, gamma=+1, n=500, MAP-IPW arm D. The largest density L1 differences occurred for degree 2, gamma=0, n=500, arm D in both comparisons.

## Predeclared gate assessment

The refined 561→1,121 comparison passes all limits in the PR #142 implementation:

- rank shift ≤ 0.002: **PASS** (maximum 0.00005381);
- median shift ≤ 0.01 km s⁻¹ Mpc⁻¹: **PASS** (maximum 0.0001689);
- maximum five-quantile shift ≤ 0.05 km s⁻¹ Mpc⁻¹: **PASS** (maximum 0.004041);
- posterior-density L1 ≤ 0.001: **PASS** (maximum 0.0001933).

The legacy 281→561 comparison also passes the same limits: maxima were 0.00005547, 0.0006487 km s⁻¹ Mpc⁻¹, 0.008271 km s⁻¹ Mpc⁻¹, and 0.0004918 respectively.

**Interpretation:** for these frozen KDEs, event data, coefficient nodes, and five arms, the 281-point uniform grid is already close to both finer grids under the predeclared numerical criteria. The refined grid comparison is reassuring for H0 integration error conditional on these archived inputs. It does not establish KDE convergence, bandwidth robustness, repeated-sampling calibration, or estimator superiority.

## Next scientific gate

1. Keep the H0-grid result separate from KDE-generation and bandwidth sensitivity.
2. Preserve the severe node-wise IPW concentration warnings from PR #141; grid convergence cannot repair unstable weights.
3. Proceed to independently seeded prior-predictive truths only after the PR #142 workflow and its independent artifact audit complete.
4. Retain the full five-arm comparison and the SUCCESS / NON_EVALUABLE / FAILED / NOT_RUN lifecycle ledger.
5. Do not promote degree two or claim calibration from this fixed-truth numerical audit.
