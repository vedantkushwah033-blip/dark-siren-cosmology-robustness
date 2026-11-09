# Exploratory finite-validation IPW host-density recovery pilot

> **Methodology status update (2026-11-09): SUPERSEDED PENDING RERUN.** The KDE implementation used for the numerical results below computed the Scott covariance scale from the reflected/augmented sample. That does not match the documented physical bandwidth definition based on the original unreflected weighted sample. A correction is under review in draft PR #39. Until the corrected implementation and independent auditor pass CI and this experiment is rerun from raw seeds, the numerical IAE values, paired contrasts, and the candidate positive-curvature failure regime below must not be treated as validated evidence. The design description remains a record of the earlier exploratory protocol; the results are retained for provenance only.


**Date:** 2026-11-08  
**Status:** isolated local pilot; not a repository benchmark, not CI-verified, not confirmatory  
**Scope:** synthetic parent-host density recovery only; no H0 posterior was computed.

## Objective

Test whether finite-validation fitted inclusion probabilities improve recovery of the parent host-redshift density, and identify scenarios where estimated IPW fails to improve on the unweighted catalogue. This is the first downstream density-stage smoke after the probability-only benchmark. It is not the full estimated-inclusion H0 experiment.

## Design

- Parent size: 4,000 galaxies.
- Replicates: 20 paired parent frames, seeds 70,000–70,019.
- True inclusion curvature: gamma = -1.2, 0, +1.2.
- Target parent-mean retention: 0.25; the logistic intercept is solved under the declared parent intensity.
- Validation sizes: 100, 300, 1,000.
- Validation design: simple random sample from the full parent frame, selected independently of retention; validation objects are excluded from the analysis catalogue.
- Fitted models: unpenalized logistic degree 1 and degree 2, trained only on validation labels.
- DGP: the same two-overdensity parent intensity and curvature-centered logistic inclusion family as the finite-validation probability pilot.
- Host target: normalized analytic parent intensity multiplied by the declared merger weight.
- Bandwidth: fixed Scott multiplier 0.15 across every arm, carried forward from the separate training/held-out density-recovery pilot. It was not chosen using H0 coverage.
- Arms: analytic parent-sample KDE (oracle diagnostic), catalogue-only KDE, known-inclusion IPW KDE, estimated degree-1 IPW KDE, and estimated degree-2 IPW KDE.
- Score: integrated absolute error (IAE) against the analytic parent-host density.

The local script used an isolated mirror of the branch's KDE and logistic fitting algorithms. It implemented the BFGS fit plus the same damped-Newton score-polish logic, but did not execute the repository benchmark runner, compute the fitter covariance, write a repository manifest, or run an independent raw-artifact auditor. All 360 estimated model fits in this pilot returned SUCCESS. This is local exploratory evidence only.

## Mean density-recovery error

Mean IAE over the 20 paired replicates:

| True curvature | Validation n | Catalogue-only | Known IPW | Estimated IPW degree 1 | Estimated IPW degree 2 | Oracle parent |
|---:|---:|---:|---:|---:|---:|---:|
| -1.2 | 100 | 0.3817 | 0.0944 | 0.1564 | 0.2539 | 0.0449 |
| -1.2 | 300 | 0.3811 | 0.0984 | 0.1182 | 0.1704 | 0.0449 |
| -1.2 | 1,000 | 0.3768 | 0.1097 | 0.1263 | 0.1427 | 0.0449 |
| 0.0 | 100 | 0.1810 | 0.0854 | 0.1505 | 0.2209 | 0.0449 |
| 0.0 | 300 | 0.1813 | 0.0881 | 0.1073 | 0.1463 | 0.0449 |
| 0.0 | 1,000 | 0.1777 | 0.0969 | 0.1076 | 0.1165 | 0.0449 |
| +1.2 | 100 | 0.1050 | 0.0788 | 0.1775 | 0.1934 | 0.0449 |
| +1.2 | 300 | 0.1058 | 0.0810 | 0.1207 | 0.1244 | 0.0449 |
| +1.2 | 1,000 | 0.1114 | 0.0878 | 0.1129 | 0.1031 | 0.0449 |

## Paired estimated-IPW minus catalogue-only contrasts

Negative values favor estimated IPW; positive values mean estimated IPW has larger density error than catalogue-only. Percentile-bootstrap intervals resample the 20 paired replicate differences.

| Curvature | Validation n | Degree 1 difference (95% interval) | Degree 2 difference (95% interval) |
|---:|---:|---:|---:|
| -1.2 | 100 | -0.2253 [-0.2427, -0.2069] | -0.1277 [-0.1915, -0.0495] |
| -1.2 | 300 | -0.2628 [-0.2773, -0.2493] | -0.2107 [-0.2524, -0.1644] |
| -1.2 | 1,000 | -0.2505 [-0.2718, -0.2289] | -0.2341 [-0.2661, -0.2003] |
| 0.0 | 100 | -0.0306 [-0.0676, 0.0154] | +0.0399 [-0.0187, 0.1053] |
| 0.0 | 300 | -0.0740 [-0.0893, -0.0582] | -0.0350 [-0.0641, -0.0036] |
| 0.0 | 1,000 | -0.0701 [-0.0879, -0.0517] | -0.0612 [-0.0802, -0.0425] |
| +1.2 | 100 | +0.0725 [0.0182, 0.1463] | +0.0884 [0.0290, 0.1571] |
| +1.2 | 300 | +0.0150 [0.0018, 0.0305] | +0.0186 [0.0040, 0.0329] |
| +1.2 | 1,000 | +0.0014 [-0.0062, 0.0096] | -0.0083 [-0.0204, 0.0042] |

## Interpretation

1. **Direction of curvature matters.** At gamma=-1.2, both estimated models improve density recovery over catalogue-only at all validation sizes, although known-inclusion IPW remains better.
2. **Correctly specified linear model.** At gamma=0, degree 1 is the correctly specified model. Its improvement over catalogue-only is uncertain at n=100 but clear in this pilot at n=300 and n=1,000.
3. **Finite-sample complexity penalty.** Degree 2 is noisier at small validation sizes, including when the true curvature is zero. At n=100 and gamma=0, the degree-2 point estimate is worse than catalogue-only and its paired interval is wide.
4. **Candidate failure condition at positive curvature.** At gamma=+1.2, both estimated models are worse than catalogue-only at n=100 and n=300. At n=1,000 the differences are near zero. In contrast, known-inclusion IPW improves over catalogue-only at all three validation sizes.
5. **Not simply extreme-weight instability.** In the gamma=+1.2 cells, mean estimated-IPW Kish ESS was approximately 910–930 at n=100 and the maximum normalized weight was around 0.0013–0.0016. The density-recovery failure therefore does not coincide with obviously catastrophic weight concentration in this pilot; probability-model error and finite validation are plausible contributors.
6. **No posterior claim follows.** These are host-density errors, not H0 posterior calibration or coverage. The asymmetry is a candidate finding to test with a larger independently audited campaign.

## Limitations

- Twenty replicates per cell are still exploratory.
- Only three curvature settings and one retention target were included.
- The local fit mirror did not compute covariance or Bayesian coefficient posterior draws.
- The bandwidth was frozen at a candidate multiplier from a separate local density-recovery split, but the complete bandwidth protocol has not yet been run through the repository's independent artifact pipeline.
- The experiment has no event-level GW data, no GW detection selection factor, and no H0 likelihood.

## Next gates

1. Implement this density-stage experiment as a repository runner with complete raw unit records, fit-attempt ledger, manifest hashes, and an independent auditor.
2. Expand to at least 100 replicates per cell and include more curvature/retention settings.
3. Reassess the positive-curvature failure under held-out seeds and an alternative parent law.
4. Compare probability RMSE/ECE and IPW ESS diagnostics with downstream density IAE on separate simulation seeds.
5. Only after the density stage is independently audited should the project implement the modular coefficient-posterior integration and selection-conditioned H0 likelihood described in the propagation design.

**Decision:** exploratory evidence supports a candidate finite-validation failure regime at positive curvature with n_validation <= 300, but the result is not yet confirmatory and does not establish an H0 effect.
