# Refined-grid 20-truth five-arm prior-predictive pilot — descriptive readout

**Date:** 2026-11-20  
**Scope:** Synthetic low-redshift Euclidean dark-siren benchmark only. This is an engineering/statistical pilot, not a confirmatory calibration study or observational inference.

## Frozen run and independent audit

- Workflow: [44731915717](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44731915717)
- Source artifact ID: `14628789845`; artifact name `geometric-h0-prior-predictive-20rep`.
- Seed: `2026112024`; 20 independently drawn prior-predictive H0 truths.
- Design: 1,200 parent galaxies; 30 detected events; H0 prior Uniform(30,110); 1,121-point geometric evaluation grid; gamma {-1,0,+1}; validation sizes {100,300,500}; five arms A–E; quadrature order 5.
- 900/900 planned cells SUCCESS; 900/900 cell audits PASS; 60/60 child audits successful.
- Independent aggregate audit: 1,148 manifest entries verified; 45 aggregate metric rows and 1,440 paired contrast rows reconstructed; 60 cross-curvature paired validation blocks; status PASS.
- Exact geometric-grid audit: [workflow 44746786352](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44746786352), artifact `14633702303`. It reconstructed all 900 selected likelihood curves with maximum absolute error 9.09e-13; 1,121→2,241 geometric-grid maximum rank shift 0.00004596, maximum nine-quantile shift 0.0009942 km/s/Mpc, and maximum posterior-density L1 distance 0.00007995.
- The grid audit holds archived KDE densities fixed. It is a numerical H0-integration audit, not an independent end-to-end regeneration of the density estimates.

## Validation-size-500 descriptive results

Bias is the mean posterior-median error (km/s/Mpc); RMSE is across the 20 paired prior-predictive truths. Mean rank is the mean posterior CDF evaluated at the true H0.

| gamma | Arm | Median bias | RMSE | Mean rank |
|---:|---|---:|---:|---:|
| -1 | A catalogue-only | -15.977 | 17.845 | 0.976 |
| -1 | B parent oracle | +1.912 | 6.353 | 0.439 |
| -1 | C known-π IPW | +1.642 | 7.355 | 0.466 |
| -1 | D degree-one MAP IPW | -1.562 | 7.100 | 0.652 |
| -1 | E degree-one integrated IPW | -1.334 | 7.028 | 0.639 |
| 0 | A catalogue-only | -15.607 | 17.593 | 0.969 |
| 0 | B parent oracle | +1.912 | 6.353 | 0.439 |
| 0 | C known-π IPW | +1.978 | 7.314 | 0.441 |
| 0 | D degree-one MAP IPW | +1.662 | 7.169 | 0.461 |
| 0 | E degree-one integrated IPW | +1.783 | 7.267 | 0.456 |
| +1 | A catalogue-only | -14.569 | 16.639 | 0.957 |
| +1 | B parent oracle | +1.912 | 6.353 | 0.439 |
| +1 | C known-π IPW | +2.137 | 6.638 | 0.431 |
| +1 | D degree-one MAP IPW | +4.580 | 8.053 | 0.298 |
| +1 | E degree-one integrated IPW | +4.567 | 8.078 | 0.301 |

### Paired E − C rank displacement at validation size 500

| gamma | Mean rank difference | MCSE | Paired 95% t interval |
|---:|---:|---:|---:|
| -1 | +0.1729 | 0.0297 | [+0.1107, +0.2351] |
| 0 | +0.0145 | 0.0171 | [-0.0212, +0.0503] |
| +1 | -0.1298 | 0.0230 | [-0.1779, -0.0818] |

This 20-truth refined-grid run reproduces the direction of the curvature-dependent E−C rank displacement seen in the independent 50-truth pilot: positive at gamma=-1, near zero at gamma=0, and negative at gamma=+1. The 20-truth estimates are less precise and should be treated as a replication check, not a stand-alone calibration result.

### Paired squared-error comparisons

For E−C, mean paired squared-median-error differences at gamma {-1,0,+1} were {-4.70, -0.69, +21.19} (km/s/Mpc)^2. Their paired 95% t intervals were [-25.23,15.82], [-12.81,11.44], and [-1.24,43.62], respectively. These intervals are wide; this pilot does not establish an RMSE advantage or disadvantage for E over C.

For E−D, squared-error differences were {-1.02,+1.42,+0.42} (km/s/Mpc)^2, with all three paired 95% intervals including zero. Coefficient integration is not materially separating E from D in this 20-truth run.

## Interpretation and limits

1. The catalogue-only arm has strongly shifted ranks and substantial negative median bias in this synthetic incomplete-catalogue design.
2. The parent-oracle and known-inclusion-probability arms have mean ranks closer to 0.5 in this small pilot, but the result is not a formal rank-uniformity or calibration verdict.
3. The fitted degree-one IPW arms show curvature-dependent rank displacement relative to known-π IPW. Integrating degree-one coefficients does not remove the omitted-quadratic misspecification pattern.
4. Increasing validation size from 100 to 500 does not change A/B/C in this design; those arms do not estimate inclusion probabilities from the validation sample. D/E vary with validation size, but not monotonically in every metric.
5. Each gamma/size/arm stratum contains only 20 independent truths. Validation sizes reuse the same truths and events, so they are paired repeated measurements, not additional independent replicates.
6. The oracle is a simulation diagnostic unavailable in ordinary applications. No claim of universal estimator superiority, real-catalogue validity, LVK/GWTC analysis, or observational H0 measurement is made.

## Decision

Keep this as a refined-grid replication/engineering pilot. The matching rank-displacement direction across the independent 20- and 50-truth runs justifies scaling the prior-predictive study to a larger number of independent truths, with the 1,121-point geometric H0 grid and the 2,241-point convergence audit retained. Do not promote the current 20- or 50-truth summaries to confirmatory calibration claims.
