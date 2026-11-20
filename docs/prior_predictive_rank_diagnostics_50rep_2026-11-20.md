# Prior-predictive rank diagnostics: 50-truth five-arm pilot

**Date:** 2026-11-20  
**Scope:** Synthetic low-redshift Euclidean dark-siren benchmark only.

## Frozen source and exact-grid audit

- Source experiment: 50-replicate prior-predictive five-arm run, seed `2026112022`, source artifact ID `14613021797`.
- Source artifact SHA-256: `1d8667e2f7391283237e27b0c6db26cac6e9dfd56d344649b6931ec047745c0c`.
- Exact-grid audit workflow: [44735380998](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44735380998).
- Exact-grid audit artifact ID: `14629769931`; SHA-256: `2304c2e653727c21728cd560e8b6f63b661b4348efa69f1204039218dc8578d7`.
- 2,250 unique posterior curves: 50 truths × 3 curvatures × 3 validation sizes × 5 arms.
- Maximum selected-likelihood reconstruction error: (4.55\times10^{-13}).

The refined 1,121 → 2,241 geometric-grid comparison passed: maximum rank shift (6.87\times10^{-5}), maximum nine-quantile shift 0.00113 km/s/Mpc, and maximum posterior-density L1 distance 0.0001245. The legacy 561 uniform → 1,121 geometric comparison failed the strict density-L1 tolerance: maximum 0.001918 against a threshold of 0.001. Rank and quantile changes in that legacy comparison were small, but the density-L1 warning is retained rather than relaxing the tolerance.

The exact-grid audit conditions on archived KDE densities; it does not independently regenerate KDEs from raw galaxy weights.

## Refined-grid prior-predictive ranks at validation size 500

The following ranks were recomputed on the 2,241-point geometric grid. Each curvature stratum has 50 independent prior-predictive truths; all five arms are paired on the same truth and data.

| Curvature | Arm | Mean rank | MCSE | 95% t interval |
|---:|---|---:|---:|---:|
| -1 | Catalogue-only A | 0.9966 | 0.0018 | [0.9930, 1.0002] |
| -1 | Parent oracle B | 0.4811 | 0.0401 | [0.4005, 0.5618] |
| -1 | Known-π IPW C | 0.5245 | 0.0501 | [0.4238, 0.6252] |
| -1 | Degree-one MAP IPW D | 0.7096 | 0.0436 | [0.6219, 0.7973] |
| -1 | Degree-one integrated IPW E | 0.6998 | 0.0423 | [0.6149, 0.7847] |
| 0 | Catalogue-only A | 0.9934 | 0.0033 | [0.9869, 1.0000] |
| 0 | Parent oracle B | 0.4811 | 0.0401 | [0.4005, 0.5618] |
| 0 | Known-π IPW C | 0.4825 | 0.0461 | [0.3898, 0.5752] |
| 0 | Degree-one MAP IPW D | 0.5067 | 0.0485 | [0.4093, 0.6041] |
| 0 | Degree-one integrated IPW E | 0.5006 | 0.0464 | [0.4073, 0.5939] |
| +1 | Catalogue-only A | 0.9828 | 0.0060 | [0.9708, 0.9949] |
| +1 | Parent oracle B | 0.4811 | 0.0401 | [0.4005, 0.5618] |
| +1 | Known-π IPW C | 0.4655 | 0.0442 | [0.3766, 0.5544] |
| +1 | Degree-one MAP IPW D | 0.3026 | 0.0407 | [0.2208, 0.3844] |
| +1 | Degree-one integrated IPW E | 0.3060 | 0.0400 | [0.2257, 0.3863] |

### Paired rank displacement: estimated IPW minus known-π IPW

| Curvature | E − C mean rank | MCSE | 95% paired t interval |
|---:|---:|---:|---:|
| -1 | +0.1753 | 0.0205 | [+0.1341, +0.2165] |
| 0 | +0.0181 | 0.0138 | [-0.0098, +0.0459] |
| +1 | -0.1595 | 0.0185 | [-0.1966, -0.1224] |

The MAP-minus-known-π differences are +0.1851 (95% interval [+0.1394,+0.2308]) at gamma -1, +0.0242 ([-0.0055,+0.0539]) at gamma 0, and -0.1629 ([-0.1968,-0.1290]) at gamma +1.

## Interpretation

1. Catalogue-only ranks are concentrated near one, indicating severe prior-predictive rank distortion when the observed catalogue is incomplete and no inclusion correction is applied in this synthetic design.
2. Known-π IPW and the full-parent oracle have mean ranks near 0.5 in these strata, with broad intervals that include 0.5.
3. Both estimated degree-one IPW arms have substantial curvature-dependent rank displacement: positive relative to known-π at gamma -1, little displacement at gamma 0, and negative displacement at gamma +1.
4. Coefficient integration (E) does not remove this pattern. The omitted quadratic term remains a model-misspecification issue; integrating over coefficients within the wrong degree-one family does not integrate over the omitted curvature.
5. This is evidence of a specific synthetic sensitivity, not universal estimator superiority or a real-catalogue result.

## Statistical cautions

- The independent units for these intervals are the 50 prior-predictive truths within a curvature stratum. Validation sizes 100/300/500 reuse the same truth and paired data; they are not additional independent replicates.
- Mean rank and its t interval are not a complete rank-uniformity/SBC test. A rank histogram and additional independent truths are needed for stronger calibration claims.
- The parent-oracle arm is a simulation diagnostic, not an available real-world correction.
- No observational H0, LVK/GWTC, or real galaxy-survey validation claim is made.
