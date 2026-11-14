# Five-width density pilot: aggregate diagnostics audit and readout

Date: 2026-11-15  
Source branch: `research/absolute-width-distribution-pilot-20261114`  
Source commit: `9894307b78c569adcc137c27b530afbdaf2247bb`  
Workflow: `43628956224`

## Status and scope

The five-width, 20-replicate density-only pilot completed successfully for all three catalogue-inclusion floors. The updated outer auditor independently reconstructed the aggregate weight-diagnostic table in addition to the existing raw-input child audits, density metric aggregation, Scott-width aggregation, seed schedule, lifecycle ledger, and manifest checks.

This is an exploratory synthetic density-recovery experiment. It is not a bandwidth-selection exercise, a confirmatory estimator comparison, an H0 likelihood or posterior analysis, a calibration study, or an analysis of real galaxy catalogues or LVK/GWTC data.

## Design

All three floor cells used the same paired replicate seed schedule and numerical design:

- 20 planned replicates; 20 successful in each floor cell.
- Parent population size 4,000.
- Independent validation frame size 1,500; nested validation sizes 100, 300, and 1,000.
- Master seed 20261215; stride 10; validation seed is population seed + 1.
- Degree-one logistic analysis model; quadrature order 5.
- 401-point redshift grid.
- Five common absolute physical KDE widths: 0.012, 0.016, 0.024, 0.032, and 0.040.
- Catalogue-inclusion floors: 0.00, 0.10, and 0.25.
- Five estimator arms, including catalogue-only, parent oracle, known-probability IPW, MAP plug-in IPW, and posterior-mean density.

The floor-zero DGP is exactly logistic in redshift and is correctly specified by the degree-one logistic analysis model. Positive floors alter the inclusion curve and create functional-form misspecification for the degree-one logistic analysis model. This distinction matters when interpreting the results.

## Independent audit

All three workflow jobs completed successfully. Each separate artifact audit reported PASS:

| Catalogue floor | Planned | Successful | Failed | Child raw audits | Manifest files verified | Metric rows reconstructed | Scott-width rows reconstructed | Weight-diagnostic rows reconstructed |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 20 | 20 | 0 | 20 | 245 | 9,300 | 1,740 | 300 |
| 0.10 | 20 | 20 | 0 | 20 | 245 | 9,300 | 1,740 | 300 |
| 0.25 | 20 | 20 | 0 | 20 | 245 | 9,300 | 1,740 | 300 |

The aggregate auditor now verifies the diagnostic keys (replicate, validation size, arm, and posterior node) and independently checks Kish ESS, ESS fraction, maximum normalized weight, top-1% mass, and top-5% mass. The child raw-input auditor independently checks the per-replicate diagnostics. The adversarial regression test deliberately alters an aggregate Kish ESS value and refreshes the manifest; it is designed to ensure that a self-consistent hash cannot conceal a numerical mismatch.

The three-floor pilot workflow at source commit `9894307` completed successfully. The full repository test matrix subsequently completed successfully on Python 3.10, 3.11, and 3.12 at the same source commit. Thus both the three-floor artifact audit and the aggregate-diagnostics regression suite are PASS at `9894307`. Later quadratic-curvature smoke-runner additions are separate changes and have their own CI validation pending.

## Fixed-width paired comparison

The table below reports mean IAE at the common grid point (h=0.016), validation size (n=1{,}000), across 20 paired replicates. This is a descriptive slice at a prespecified grid point, not a selected or optimized bandwidth.

| Floor | Catalogue-only A | Parent oracle B | Known-probability IPW C | MAP plug-in IPW D | Posterior-mean density E |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.53049 | 0.18378 | 0.20391 | 0.20720 | 0.21032 |
| 0.10 | 0.37809 | 0.18378 | 0.18653 | 0.29618 | 0.30006 |
| 0.25 | 0.21349 | 0.18378 | 0.18065 | 0.27331 | 0.27605 |

Paired mean differences (first arm minus second arm) and Monte Carlo standard errors:

| Floor | Contrast | Mean difference | MCSE |
|---:|---|---:|---:|
| 0.00 | D − A | -0.32330 | 0.01145 |
| 0.00 | D − C | +0.00329 | 0.01189 |
| 0.00 | E − D | +0.00312 | 0.00035 |
| 0.00 | C − B | +0.02013 | 0.00654 |
| 0.10 | D − A | -0.08192 | 0.01470 |
| 0.10 | D − C | +0.10964 | 0.01359 |
| 0.10 | E − D | +0.00388 | 0.00012 |
| 0.10 | C − B | +0.00275 | 0.00401 |
| 0.25 | D − A | +0.05983 | 0.00907 |
| 0.25 | D − C | +0.09266 | 0.01012 |
| 0.25 | E − D | +0.00274 | 0.00008 |
| 0.25 | C − B | -0.00313 | 0.00256 |

Positive differences mean the first estimator has larger IAE and therefore worse density recovery.

### Interpretation

1. **Correctly specified floor-zero case:** MAP estimated IPW substantially improves over catalogue-only density recovery and approaches known-probability IPW at this width and validation size. The D−C difference is small relative to its Monte Carlo standard error.
2. **Positive-floor cases:** The MAP estimated-IPW arm remains materially worse than known-probability IPW. At floor 0.10 it still improves over catalogue-only, but at floor 0.25 it is worse than catalogue-only at this fixed width.
3. **Uncertainty integration does not remove model discrepancy:** Posterior-mean density E has slightly larger IAE than MAP plug-in D in all three cells at this fixed width. Integrating coefficient uncertainty does not correct a misspecified inclusion-probability family.
4. **Known-probability IPW tracks the parent oracle in these cells:** C is close to B, but B is an oracle simulation diagnostic and is not available in ordinary observational analyses.
5. **The floor comparison changes both the selection shape and the density difficulty:** The numerical results should be read as a controlled synthetic stress comparison, not as an isolated causal effect of one parameter without considering the full DGP.

## Weight concentration is not the same as model accuracy

At (n=1{,}000), mean diagnostics across 20 replicates are:

| Floor | Arm | Kish ESS | ESS fraction | Max normalized weight | Top 1% mass | Top 5% mass |
|---:|---|---:|---:|---:|---:|---:|
| 0.00 | C known | 541.6 | 0.2469 | 0.0125 | 0.1446 | 0.3350 |
| 0.00 | D MAP | 536.2 | 0.2445 | 0.0128 | 0.1466 | 0.3380 |
| 0.10 | C known | 1,126.1 | 0.4750 | 0.0029 | 0.0641 | 0.2337 |
| 0.10 | D MAP | 813.0 | 0.3429 | 0.0057 | 0.1086 | 0.2980 |
| 0.25 | C known | 1,726.4 | 0.6546 | 0.0013 | 0.0353 | 0.1551 |
| 0.25 | D MAP | 1,429.9 | 0.5421 | 0.0025 | 0.0608 | 0.2125 |

The floor-0.25 estimated-IPW arm has less concentrated weights than the floor-zero arm by these summaries, yet its density result is worse relative to the catalogue-only arm. This is direct motivation to keep functional-form misspecification and weight instability as separate diagnostic axes. A high ESS alone cannot certify a valid inclusion model.

## Bandwidth-grid limitation

The five tested widths are not sufficient to choose a bandwidth. For the catalogue-only arm at (n=1{,}000), the lowest mean IAE in the tested grid occurs at the lower boundary (0.012) for floor 0.25; the floor-zero and floor-0.10 curves are comparatively flat across the tested range. An optimum at a grid boundary is not bracketed. Even where another arm has an interior minimum, selecting a common width from these same exploratory replicates would risk tuning to the evaluation sample.

No bandwidth is selected here. The next density stage must widen the grid where needed, freeze a training/holdout protocol, select a common physical bandwidth using training replicates only, and evaluate it once on held-out replicates.

## Limitations and next gates

- Only 20 paired replicates per floor; uncertainty is exploratory.
- The five-width grid is not guaranteed to bracket the optimum for every arm.
- Stored inclusion posterior nodes are consumed from the production fitter; this experiment does not independently refit the inclusion model from raw validation data in the outer auditor.
- The experiment does not yet establish finite-validation posterior calibration or repeated-sampling coverage.
- No H0 likelihood, posterior, PIT, coverage, or selection-normalization analysis is included.
- No observational data, real galaxy catalogue, LVK/GWTC public data, or measured H0 result is used.

Next gates, in order:

1. The full Python 3.10/3.11/3.12 matrix is PASS at source commit `9894307`; rerun it for the subsequent curvature-smoke commits.
2. Finish the expanded-width re-evaluation and verify its independent audit, including the new aggregate weight diagnostics.
3. Exercise the new quadratic-logit curvature smoke runner and independent raw-input auditor; then scale to the seven-gamma paired-replicate benchmark.
4. Independently validate the inclusion fitter and posterior quadrature from raw validation records.
5. Freeze bandwidth selection on training replicates and pass held-out density validation.
6. Only then begin a separately audited H0 propagation study with the shared-`H_0` joint selected likelihood and explicit separation of catalogue inclusion from GW detection selection.

**Conclusion:** The new aggregate diagnostic audit is now exercised on all three 20-replicate floor cells and reports PASS. The results reinforce the central synthetic research question: estimated inclusion correction can work well under correct specification, but functional-form misspecification can eliminate or reverse its benefit even when weight concentration is not extreme. This remains an exploratory density-only result, not a general estimator claim.
