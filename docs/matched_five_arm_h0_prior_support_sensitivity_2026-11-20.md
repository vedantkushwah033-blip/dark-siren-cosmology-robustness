# Matched five-arm H0 prior-support sensitivity: 100 paired replicates

**Date:** 2026-11-20  
**Synthetic design:** low-redshift Euclidean dark-siren benchmark; fixed truth H0=70 km s^-1 Mpc^-1; 1,200 parent galaxies; 30 detected events; curvature gamma in {-1,0,+1}; validation sizes 100/300/500; five inference arms.  
**Comparison:** identical master seed 2026111929 and identical replicate-level parent/event/validation inputs under prior supports [50,90] and [30,110].  
**Status:** paired artifact audit PASS; this is prior-support sensitivity, not observational inference or a new estimator comparison.

## Immutable source artifacts and audit

- Narrow-prior source: workflow [44460222957](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44460222957), artifact 14504568415.
- Wide-prior paired source: workflow [44484212311](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44484212311), artifact 14520551336.
- Paired audit with full printed summary: workflow [44601642485](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44601642485), artifact 14567087743.

The paired audit verified:
- Both source artifact manifests and all lifecycle denominators.
- 300 paired replicate-curvature child directories and 4,500 one-to-one posterior summary rows.
- 2,700 SHA-256 comparisons of raw parent/event/validation inputs.
- 600 SHA-256 comparisons of weight-diagnostic and realized-bandwidth files.
- Common-support likelihood curves for every arm and validation size.
- Maximum absolute log-likelihood difference on common [50,90] support: **2.2737367544323206e-13**.

Each child had its own independent audit report marked PASS. The child auditor reconstructs selected likelihoods and posterior quantiles from the archived inputs/likelihood curves before accepting the posterior summary. The paired audit then checks that the same likelihood is being conditioned on the two different prior supports.

## Paired effects at validation size 500

All differences below are **wide prior minus narrow prior**. Median-error differences and interval-width differences are in km s^-1 Mpc^-1. Intervals are paired t-based 95% intervals across 100 identical replicate seeds within each curvature/arm cell.

| Curvature | Arm | Median-error change | Paired 95% interval | 68% width change | Paired 95% interval |
|---:|---|---:|---:|---:|---:|
| -1 | A catalogue-only | -2.2806 | [-2.7626, -1.7986] | +2.9580 | [2.5111, 3.4050] |
| -1 | B parent oracle | +0.0013 | [-0.0003, 0.0030] | +0.0030 | [-0.0006, 0.0066] |
| -1 | C known-probability IPW | +0.0097 | [-0.0087, 0.0280] | +0.0204 | [-0.0105, 0.0514] |
| -1 | D MAP estimated IPW | -0.0013 | [-0.0063, 0.0038] | +0.0127 | [0.0024, 0.0230] |
| -1 | E integrated estimated IPW | -0.0010 | [-0.0071, 0.0051] | +0.0149 | [0.0027, 0.0271] |
| 0 | A catalogue-only | -2.4773 | [-3.0260, -1.9285] | +3.3390 | [2.7668, 3.9111] |
| 0 | B parent oracle | +0.0013 | [-0.0003, 0.0030] | +0.0030 | [-0.0006, 0.0066] |
| 0 | C known-probability IPW | +0.0041 | [-0.0037, 0.0119] | +0.0098 | [-0.0052, 0.0248] |
| 0 | D MAP estimated IPW | +0.0049 | [-0.0059, 0.0156] | +0.0144 | [-0.0056, 0.0343] |
| 0 | E integrated estimated IPW | +0.0055 | [-0.0064, 0.0174] | +0.0161 | [-0.0056, 0.0378] |
| +1 | A catalogue-only | -2.6613 | [-3.3384, -1.9842] | +3.5946 | [2.9063, 4.2828] |
| +1 | B parent oracle | +0.0013 | [-0.0003, 0.0030] | +0.0030 | [-0.0006, 0.0066] |
| +1 | C known-probability IPW | +0.0043 | [-0.0001, 0.0087] | +0.0095 | [0.0003, 0.0187] |
| +1 | D MAP estimated IPW | +0.0291 | [-0.0014, 0.0596] | +0.0499 | [-0.0001, 0.1000] |
| +1 | E integrated estimated IPW | +0.0285 | [-0.0026, 0.0596] | +0.0479 | [-0.0012, 0.0971] |

The paired 68% coverage difference (wide minus narrow) was exactly zero in every cell except gamma=0, n=100, arm E, where it was -0.01 (paired 95% t interval [-0.0298, 0.0098]). Thus this prior change did not materially alter the binary 68% coverage outcome in this fixed-truth experiment.

A/B/C effects are identical across validation sizes because those arms do not depend on the fitted validation inclusion model. D/E effects vary slightly with n, as expected from the validation-dependent fitted weights.

## Interpretation

1. **Catalogue-only A is strongly prior-boundary sensitive.** Expanding support from [50,90] to [30,110] shifts its posterior median downward by about 2.28–2.66 km s^-1 Mpc^-1 and increases its 68% interval width by about 2.96–3.59. All paired 95% intervals exclude zero for these changes.
2. **The corrected arms are comparatively stable to this prior expansion.** For B and C, median shifts are close to zero. D/E have small changes relative to A; the positive-curvature D/E shifts are around +0.03 with paired intervals that narrowly include zero.
3. **The wide prior does not rescue catalogue-only inference.** It exposes more of A's lower-H0 tail and makes its negative median error more pronounced. This is consistent with the separate expanded-support boundary audit.
4. **The prior expansion does not explain the curvature-dependent estimated-IPW pattern.** D/E remain sensitive to inclusion-model curvature and validation uncertainty in the five-arm experiments. The matched prior-support comparison changes their summaries only modestly.
5. **No calibration claim follows.** This is one fixed H0 truth, one matched 100-replicate seed stream, and a synthetic model. The t intervals quantify paired Monte Carlo variation for this design; they are not general posterior calibration intervals.

## Methodological note

The two H0 grids are nested with equal spacing 0.142857 km s^-1 Mpc^-1 on the common [50,90] region: 281 nodes for [50,90] and 561 nodes for [30,110]. Therefore the likelihood comparison is a point-for-point comparison on identical H0 nodes. The prior support changes the normalized posterior, not the common-support likelihood itself. The numerical discrepancy of 2.27e-13 is floating-point level.

## Next gate

Prior-support sensitivity is now quantified. The remaining high-value questions are:
- full posterior calibration under the wider prior (prior-predictive SBC and fixed-truth coverage, separately);
- whether pre-posterior weight-tail/ESS diagnostics predict H0 posterior error and calibration failure;
- whether degree-two inclusion modelling actually reduces curvature-induced bias without unacceptable variance;
- full geometric-grid density-L1 convergence and holdout KDE bandwidth sensitivity.

Do not claim that a wider prior improves the estimator. It corrects truncation of the posterior support; it does not correct a misspecified host/catalogue model.
