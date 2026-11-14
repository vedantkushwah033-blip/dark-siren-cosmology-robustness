# Expanded 11-width density pilot: independently audited readout

Date: 2026-11-14

## Status and scope

The expanded common-physical-bandwidth pilot (workflow `43617052462`, source head `d94fdd62e514c81c20444196e166d2ebbb341a40`) completed successfully for all three catalogue-inclusion floors.

This is a density-recovery sensitivity pilot only. It is not a bandwidth-selection exercise, a confirmatory estimator comparison, an H0 likelihood or posterior analysis, a calibration study, or an analysis of real galaxy catalogues or LVK/GWTC data.

## Design

All three floor cells use the same paired replicate seed schedule and numerical design:

- 20 planned replicates; 20 successful in every floor cell.
- Parent population size 4,000.
- Independent validation frame size 1,500; nested validation sizes 100, 300, and 1,000.
- Master seed 20261215; stride 10; validation seed is population seed + 1.
- Degree-one logistic analysis model; quadrature order 5.
- 401-point redshift grid.
- Common absolute physical KDE widths: 0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080.
- Catalogue-inclusion floors: 0.00, 0.10, 0.25.

The floor-zero DGP is exactly logistic in redshift and is correctly specified by the degree-one logistic analysis model. A positive floor adds an offset to the logistic curve, so the same degree-one logistic analysis model is functionally misspecified. This distinction is central to interpreting the results.

## Independent artifact audit

All three floor jobs completed with PASS:

| Floor | Planned | Successful | Failed | Child raw audits | Manifest files | Aggregate metric rows |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 20 | 20 | 0 | 20 | 244 | 20,460 |
| 0.10 | 20 | 20 | 0 | 20 | 244 | 20,460 |
| 0.25 | 20 | 20 | 0 | 20 | 244 | 20,460 |

Each audit independently reconstructed 1,740 realized Scott-width rows. The outer auditor verified the seed schedule, replicate lifecycle ledger, manifest hashes, every successful child, and aggregate metric reconstruction. The audits explicitly do not authorize bandwidth selection or H0 inference.

## Paired contrast at a fixed physical width

The following is a descriptive comparison at the fixed grid point h = 0.016. This is a comparison slice, not a selected or optimized bandwidth. Contrasts are paired by replicate seed. Intervals are percentile bootstrap intervals over the 20 paired replicate differences (20,000 resamples); with only 20 pilot replicates, they are exploratory rather than confirmatory.

Positive values mean the estimated arm has larger IAE (worse density recovery) than known-probability IPW (Arm C).

| Floor | Validation n | MAP plug-in D − C | Posterior-mean density E − C |
|---:|---:|---:|---:|
| 0.00 | 100 | +0.0462 [-0.0114, +0.1087] | +0.0728 [+0.0101, +0.1401] |
| 0.00 | 300 | +0.0174 [-0.0175, +0.0539] | +0.0279 [-0.0088, +0.0661] |
| 0.00 | 1,000 | +0.0033 [-0.0197, +0.0257] | +0.0064 [-0.0166, +0.0290] |
| 0.10 | 100 | +0.1659 [+0.0874, +0.2500] | +0.1978 [+0.1167, +0.2839] |
| 0.10 | 300 | +0.1299 [+0.0789, +0.1843] | +0.1419 [+0.0899, +0.1973] |
| 0.10 | 1,000 | +0.1096 [+0.0843, +0.1364] | +0.1135 [+0.0880, +0.1405] |
| 0.25 | 100 | +0.1189 [+0.0592, +0.1827] | +0.1429 [+0.0797, +0.2104] |
| 0.25 | 300 | +0.0781 [+0.0471, +0.1117] | +0.0864 [+0.0548, +0.1206] |
| 0.25 | 1,000 | +0.0927 [+0.0737, +0.1127] | +0.0954 [+0.0763, +0.1156] |

### Interpretation

1. **Correct specification:** At floor zero, the MAP plug-in contrast moves toward zero as validation size increases. At n = 1,000, both D − C and E − C intervals include zero. This is consistent with the fitted logistic model learning the true inclusion curve as validation information increases.
2. **Positive-floor misspecification:** At floors 0.10 and 0.25, the estimated arms remain worse than known-probability IPW at all three validation sizes. Increasing validation size reduces some of the gap but does not remove it.
3. **Uncertainty integration is not a cure for functional-form error:** Arm E's posterior-mean density has a larger IAE contrast than Arm D in these fixed-width summaries. Integrating coefficient uncertainty does not repair a misspecified probability family.
4. **Known-probability IPW remains close to the parent oracle:** At h = 0.016, mean IAE for Arm B is 0.18378. Arm C means are 0.20391, 0.18653, and 0.18065 for floors 0.00, 0.10, and 0.25 respectively. The oracle comparison is a simulation diagnostic, not an estimator available in ordinary observational analysis.

## Weight concentration is a separate diagnostic dimension

At n = 1,000, mean weight diagnostics across the 20 replicates are:

| Floor | Arm | Kish ESS | ESS fraction | Max normalized weight | Top 1% mass | Top 5% mass |
|---:|---|---:|---:|---:|---:|---:|
| 0.00 | C known | 541.6 | 0.2469 | 0.0125 | 0.1446 | 0.3350 |
| 0.00 | D MAP | 536.2 | 0.2445 | 0.0128 | 0.1466 | 0.3380 |
| 0.00 | E aggregate | 300.5 | 0.1370 | 0.0204 | 0.2125 | 0.4427 |
| 0.10 | C known | 1,126.1 | 0.4750 | 0.0029 | 0.0641 | 0.2337 |
| 0.10 | D MAP | 813.0 | 0.3429 | 0.0057 | 0.1086 | 0.2980 |
| 0.10 | E aggregate | 535.6 | 0.2258 | 0.0083 | 0.1510 | 0.3810 |
| 0.25 | C known | 1,726.4 | 0.6546 | 0.0013 | 0.0353 | 0.1551 |
| 0.25 | D MAP | 1,429.9 | 0.5421 | 0.0025 | 0.0608 | 0.2125 |
| 0.25 | E aggregate | 1,103.8 | 0.4185 | 0.0035 | 0.0827 | 0.2663 |

For Arm E, the aggregate ESS values are the minimum across coefficient nodes, while maximum normalized weight and tail masses are the maxima across nodes. They are conservative aggregate diagnostics, not averages over nodes.

A useful counterexample to a weight-instability-only explanation is visible here: the floor-0.25 MAP weights are less concentrated than the floor-zero MAP weights by all listed concentration diagnostics, yet the floor-0.25 MAP density estimate is substantially worse relative to Arm C. The probability-model discrepancy and weight concentration must therefore be monitored as separate dimensions.

## Bandwidth-grid caution

The mean-IAE minima on the 11-point grid are not stable across arms, floors, or validation sizes. For example, the catalogue-only minimum is at the upper boundary h = 0.080 for floors 0.00 and 0.10, but at the lower boundary h = 0.012 for floor 0.25. Other arms have minima around h = 0.012–0.024. This is not evidence for a universal optimal bandwidth.

No bandwidth is selected from this pilot. A training/holdout design and a frozen estimand are still required before any confirmatory bandwidth-selection claim.

## Reconciliation with the earlier five-width pilot

The earlier five-width floor-0.10 pilot (workflow `43619594762`, code head `17466326a4be59ca7698c6c1f34529801e6e5507`) and the expanded 11-width floor-0.10 pilot share the same master seed (20261215), stride (10), parent size, validation frame, nested validation sizes, quadrature order, grid resolution, and catalogue floor. They are therefore a re-evaluation of the same paired replicate design with a wider bandwidth grid, not independent replications.

I compared the raw aggregate files at the five overlapping widths (0.012, 0.016, 0.024, 0.032, 0.040):

- 9,300 metric rows and all replicate/validation-size/arm/node keys matched.
- Maximum absolute differences were (1.72\times10^{-15}) for IAE, (1.07\times10^{-14}) for ISE, (8.05\times10^{-16}) for boundary IAE, and (1.37\times10^{-12}) for Kish ESS.
- The 1,740 realized Scott-width rows matched in key and differed by at most (1.01\times10^{-16}).
- For replicate 0, object IDs, validation IDs, and all retention/selection indicators matched exactly; redshift and inclusion-probability values differed only at floating-point roundoff scale (at most (1.11\times10^{-16}) and (4.44\times10^{-16}), respectively).

Thus the expanded run is numerically consistent with the earlier five-width run on their shared grid. Do not count the two runs as 40 independent replicates or pool them as independent evidence.

## Next gates

1. Preserve the three-floor expanded-grid artifacts and their independent PASS reports as the current density-only pilot record.
2. Add the explicit truth-probability API and independent weight-diagnostic audit improvements on the curvature-fix PR before running the quadratic-logit family.
3. Freeze the seven-gamma density estimand, paired seed schedule, failure ledger, and any train/holdout separation before running that experiment.
4. Keep H0 propagation blocked until the density-stage audit and downstream selection-conditioned likelihood requirements are separately satisfied.
