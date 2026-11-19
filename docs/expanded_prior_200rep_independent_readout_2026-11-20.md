# Expanded-prior five-arm H0 experiment: independent 200-replicate readout

**Status:** independent artifact and numerical reconstruction PASS. Exploratory synthetic fixed-truth sensitivity; not observational inference, prior-predictive SBC, or a universal estimator ranking.

## Frozen design and provenance

Two independently seeded streams were run with identical design settings:

- Master seeds: `2026111940` and `2026111929`; 100 parent/event replicates per stream.
- Parent population: 1,200 galaxies; 30 detected GW events per realization.
- Truth: (H_0=70) km s(^{-1}) Mpc(^{-1}).
- Curvature: (gamma=-1,0,+1); (eta=-3.9).
- Validation frame: 500; nested validation sizes 100, 300, 500.
- Five arms: catalogue-only A, full-parent oracle B, known-π IPW C, degree-one MAP IPW D, degree-one uncertainty-integrated IPW E.
- Common uniform prior support: [30,110]; 561-point H0 grid.
- Fixed Scott multiplier 0.7; coefficient quadrature order 5.

Each stream contains 4,500 planned lifecycle cells (100 replicates × 3 curvatures × 3 validation sizes × 5 arms) and 300 child replicate/curvature runs. The independent orchestration auditor was rerun on both extracted source archives. Both returned PASS; every planned cell was SUCCESS and every child run was independently audited. The audit reconstructed aggregate metrics and paired contrasts, checked cross-curvature shared inputs and nested validation samples, and verified the aggregate manifests.

The follow-up audit independently reconstructed all **9,000 saved posterior medians** from the child likelihood arrays and all **9,000 same-likelihood prior-support comparisons** (full [30,110] versus restricted [50,90]). Maximum saved-versus-reconstructed median discrepancy was required to be at most (2\times10^{-7}); all passed.

- Audit workflow: `44602415948` (SUCCESS).
- Audit artifact: `14568202518`; SHA-256 `e3e1c1359d43ed7805f3473bca184571a2ddd5a2da79e589d477a788322bccd9`.
- Two source runs: `44483923518` and `44484212311`.

## Posterior-median performance at validation size 500

All entries use 200 independent parent/event replicates per curvature. Bias is mean posterior-median error; RMSE is the root mean squared posterior-median error. Coverage is empirical equal-tail interval coverage.

| Curvature | Arm | Bias | RMSE | 68% coverage | 90% coverage | 95% coverage |
|---:|---|---:|---:|---:|---:|---:|
| -1 | Catalogue-only A | -17.436 | 17.911 | 0.0% | 1.0% | 1.0% |
| -1 | Parent oracle B | +0.714 | 3.398 | 71.5% | 93.5% | 97.5% |
| -1 | Known-π IPW C | -0.003 | 4.647 | 52.0% | 81.5% | 88.0% |
| -1 | MAP IPW D | -3.000 | 5.081 | 52.0% | 71.0% | 78.5% |
| -1 | Integrated IPW E | -2.750 | 4.910 | 54.0% | 75.0% | 82.5% |
| 0 | Catalogue-only A | -16.992 | 17.611 | 0.0% | 1.0% | 2.5% |
| 0 | Parent oracle B | +0.714 | 3.398 | 71.5% | 93.5% | 97.5% |
| 0 | Known-π IPW C | +0.360 | 4.086 | 62.0% | 87.0% | 93.5% |
| 0 | MAP IPW D | +0.238 | 4.109 | 62.5% | 85.0% | 91.5% |
| 0 | Integrated IPW E | +0.405 | 4.081 | 66.5% | 87.5% | 94.5% |
| +1 | Catalogue-only A | -16.069 | 16.951 | 1.0% | 4.5% | 8.0% |
| +1 | Parent oracle B | +0.714 | 3.398 | 71.5% | 93.5% | 97.5% |
| +1 | Known-π IPW C | +0.601 | 3.983 | 63.5% | 86.5% | 93.5% |
| +1 | MAP IPW D | +3.328 | 5.343 | 48.5% | 78.5% | 86.5% |
| +1 | Integrated IPW E | +3.320 | 5.311 | 51.0% | 83.5% | 86.5% |

The oracle B is a simulation diagnostic, not an implementable estimator. Its relatively good performance does not establish a real-data benchmark.

## Paired squared-error contrasts at validation size 500

Each contrast is the paired mean difference in squared posterior-median error: first arm minus reference arm. Positive values mean the first arm has higher squared error. Intervals are paired 95% Student-t intervals over 200 replicates.

| Curvature | Contrast | Mean difference | MCSE | 95% interval |
|---:|---|---:|---:|---:|
| -1 | E − C | +2.516 | 1.889 | [-1.209, +6.241] |
| -1 | D − C | +4.222 | 2.012 | [+0.254, +8.190] |
| -1 | E − D | -1.706 | 0.328 | [-2.352, -1.060] |
| 0 | E − C | -0.040 | 0.724 | [-1.467, +1.388] |
| 0 | D − C | +0.191 | 0.731 | [-1.250, +1.632] |
| 0 | E − D | -0.231 | 0.341 | [-0.902, +0.441] |
| +1 | E − C | +12.343 | 1.941 | [+8.515, +16.171] |
| +1 | D − C | +12.680 | 1.988 | [+8.760, +16.599] |
| +1 | E − D | -0.337 | 0.448 | [-1.221, +0.547] |

Interpretation:

- At (gamma=-1), integrating coefficient uncertainty E improves squared error relative to MAP D, but E is not demonstrably better than known-π C.
- At (gamma=0), E, D, and C are statistically indistinguishable by these paired squared-error contrasts.
- At (gamma=+1), both estimated degree-one arms D/E are substantially worse than known-π C. E does not significantly improve on D.
- Therefore posterior integration does not repair the omitted quadratic term. It can reduce MAP plug-in error in one curvature regime without delivering robust gains against known-π IPW.

## Validation sample-size sensitivity

Paired change in squared error from (n_{val}=100) to (n_{val}=500), computed within the same parent/event replicate. Negative values mean lower squared error at (n=500).

| Curvature | Arm | Mean MSE change | MCSE | 95% interval |
|---:|---|---:|---:|---:|
| -1 | D MAP | -9.941 | 1.680 | [-13.253, -6.629] |
| -1 | E integrated | -3.564 | 1.239 | [-6.007, -1.120] |
| 0 | D MAP | -6.486 | 1.555 | [-9.553, -3.420] |
| 0 | E integrated | -3.416 | 1.204 | [-5.790, -1.042] |
| +1 | D MAP | -1.230 | 1.958 | [-5.090, +2.631] |
| +1 | E integrated | +0.739 | 1.348 | [-1.919, +3.396] |

Larger validation samples improve MSE for the fitted degree-one methods in the negative and zero-curvature settings. At positive curvature, neither method has a statistically clear MSE improvement from 100 to 500 validation objects. This is consistent with a misspecified inclusion model approaching a pseudo-true approximation rather than recovering the omitted curvature; the table itself is not a proof of an asymptotic limit.

## Same-likelihood prior-support sensitivity

For every saved likelihood curve, the posterior was normalized once on [30,110] and again on [50,90]. No data, likelihood, or event realization was changed. Values below are at validation size 500.

| Curvature | Arm | Mean posterior mass below 50 | Fraction with mass >10% | Median shift (wide − narrow) | q2.5% shift (wide − narrow) |
|---:|---|---:|---:|---:|---:|
| -1 | A catalogue-only | 0.325 | 80.5% | -2.110 | -7.754 |
| 0 | A catalogue-only | 0.308 | 76.5% | -2.307 | -8.492 |
| +1 | A catalogue-only | 0.278 | 70.5% | -2.474 | -8.966 |

For B-E, mean mass below 50 is at most 0.000486, and absolute mean median shifts are at most 0.0203. Thus the support effect is strongly arm-specific. Widening the support exposes a lower tail in A, but A's median bias remains approximately -16 to -17 km s(^{-1}) Mpc(^{-1}); wider support does not cure catalogue-selection bias.

## Scientific conclusions and limits

1. **Degree-one misspecification is consequential.** The estimated-IPW arms retain signed median bias near -2.75 at (gamma=-1) and +3.32 at (gamma=+1), even with 500 validation objects.
2. **Coefficient-posterior integration is not a general fix.** It improves on MAP D at negative curvature, is similar at zero curvature, and remains far worse than known-π C at positive curvature.
3. **Validation size is not a substitute for model adequacy.** Increasing (n_{val}) helps in some regimes but does not reliably remove the positive-curvature penalty.
4. **Catalogue-only A has a distinct prior-boundary sensitivity.** The [50,90] support discards a substantial lower tail in A, while the IPW/oracle arms are nearly invariant. A remains biased under [30,110].
5. **Coverage is not calibrated across all arms.** Several 68% coverage rates are materially below nominal; the fixed-truth experiment is not SBC.

Limitations: synthetic low-redshift Euclidean setup only; fixed truth (H_0=70); 30 events; 561-point H0 grid; fixed KDE bandwidth multiplier 0.7; quadrature order 5. The geometric-grid and grid-refinement audits remain separate numerical gates. Oracle weights are diagnostic only. No LVK/GWTC/public event data, real galaxy catalogue, observational H0 estimate, or universal estimator-superiority claim is involved.
