# Independent review of the 2026-11-08 CI artifacts

**Reviewed:** GitHub Actions run [#42455270698](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42455270698)  
**Run head:** `23cc0333324133712706f3e07df5d8e179412260`  
**Workflow conclusion:** success  
**Scope:** artifact inspection and independent interpretation; no claim of real-data validation.

## 1. Exact-head CI status

The run completed successfully on Python 3.10, 3.11 and 3.12. All three matrix jobs completed successfully. The Python 3.12 job executed the configured long-running numerical experiments; the equivalent expensive steps were intentionally skipped on Python 3.10 and 3.11. The matrix jobs ran the relevant unit tests and smoke checks.

The Python 3.12 artifact was downloaded and inspected. It was approximately 72.0 MB and contained the finite-validation probability pilot, host-density grid sensitivity, full distance-uncertainty experiment, selected-likelihood PIT calibration, catalogue-thinning audit, and supporting grid diagnostics.

This exact-head CI pass does **not** cover later commits adding the explicit Scott-multiplier KDE API. Those changes require a new CI run before being described as CI-verified.

## 2. Finite-validation probability pilot

The independent raw-record auditor reported:

- **PASS**
- 450 experiment/model groups reconstructed.
- 1,800,000 parent-frame unit records checked.
- 450 complete parent partitions checked.
- Numerical reconstruction tolerance: (10^{-10}).

This confirms artifact consistency for the configured five-replicate probability-stage pilot. It does not establish estimator superiority, H0 calibration, or coverage.

## 3. Host-density H0 grid sensitivity

The independent grid auditor reported **PASS** for paired grids of 161, 321, 641, 1,281 and 3,201 H0 points, with 15 replicate/arm pairs and a predeclared 0.01 km s^-1 Mpc^-1 tolerance.

The maximum absolute differences versus the 3,201-point reference were:

| H0 grid points | Median | 68% lower | 68% upper | 68% width |
|---:|---:|---:|---:|---:|
| 161 | 0.000275 | 0.003181 | 0.003324 | 0.006320 |
| 321 | 0.000074 | 0.000812 | 0.000889 | 0.001592 |
| 641 | 0.000017 | 0.000200 | 0.000204 | 0.000404 |
| 1,281 | 0.000004 | 0.000049 | 0.000046 | 0.000083 |
| 3,201 | 0 | 0 | 0 | 0 |

All are in km s^-1 Mpc^-1. This establishes numerical H0-grid stability **conditional on the supplied host-density estimates**. It does not validate KDE density recovery, posterior calibration, or the scientific likelihood.

## 4. Full distance-uncertainty sensitivity experiment

The full experiment produced 32,400 raw records: 9 fixed truths × 4 fractional-noise settings × 300 replicates × 3 likelihood arms. The independent auditor reconstructed all 32,400 posteriors and 108 summary cells, verified all 10,800 paired A/B comparisons, and reported **PASS**.

A separate fixed-truth interval audit compared Arm A's Bayesian equal-tail intervals with an inverted Gaussian-distance reference interval. For interior truths 55 through 85 (8,400 replicates):

| Interval | Bayesian equal-tail coverage | Inverted-distance reference coverage |
|---|---:|---:|
| 68% | 0.7329 (95% Wilson: 0.7233–0.7422) | 0.6762 (0.6661–0.6861) |
| 95% | 0.9712 (95% Wilson: 0.9674–0.9746) | 0.9512 (0.9464–0.9556) |

At the central truth H0=70 (1,200 replicates), Bayesian equal-tail coverage was 0.7150 at 68% and 0.9817 at 95%; the inverted-distance reference was 0.6633 and 0.9550, respectively.

Interpretation limits:
- Arm A uses the realized noise scale generated from the true H0, so it is an oracle-style conditional benchmark, not an operational noise model.
- Fixed-truth Bayesian credible-interval coverage is not generally required to equal its nominal level at every truth; prior effects and truth location matter.
- The boundary truths H0=50 and 90 lie on the prior boundaries and must be reported separately. Pooling them with interior truths produces misleading aggregate coverage.
- These results concern a one-host distance toy model. They do not validate the dark-siren catalogue-inclusion method.

The paired A-minus-B bootstrap found small coverage differences across the seven interior truths; for fractional error 0.05, the mean 68% coverage difference was -0.00333 with a stratified-bootstrap 95% interval [-0.00905, 0.00238]. This contrast is specific to the defined toy arms and is not evidence about IPW.

## 5. Other artifact checks

### Selected-likelihood PIT calibration

The 1,000-replicate matched-model PIT artifact passed its independent raw-PIT reconstruction. PIT mean was 0.4966 (Monte Carlo SE 0.0089), and the descriptive KS statistic was 0.0249 (p=0.559). Decile counts ranged from 77 to 114 out of 1,000.

A non-significant KS test is not proof of exact calibration. This is a separate selected-likelihood toy diagnostic, not a finite-validation inclusion experiment.

### Catalogue thinning

The 40,000-object synthetic distance-dependent thinning artifact passed its arithmetic/raw-record audit. It retained 15,073 objects versus 15,260.3 expected. One distance stratum had a standardized residual of -2.69. The artifact explicitly limits itself to conditional Bernoulli-thinning arithmetic; it does not validate a real survey selection law, galaxy-redshift intensity, or H0 inference.

### Extreme numerical grid stress

The deterministic distance-uncertainty stress report flagged four very small-noise cases (sigma 0.2, 0.05, 0.02 and 0.01 Mpc) as **UNDERRESOLVED** on the 1,601-point grid. The maximum posterior-quantile difference against a 400,001-point calculation reached 0.01734 km s^-1 Mpc^-1. This is correctly retained as a warning, not hidden by the successful ordinary-grid checks.

## 6. Highest-priority unresolved issue: host-density bandwidth

The existing host-density benchmark still calls `estimate_host_density(..., bandwidth_factor=0.8)`. That argument is an absolute SciPy covariance factor, not 0.8 times Scott's rule. The local forensic pilot showed substantial parent-density recovery and H0-posterior sensitivity to this setting.

Consequently:
1. The current factor-0.8 KDE oracle is not yet an adequate calibration reference.
2. H0-grid convergence cannot repair or detect a biased host-density estimate.
3. No finite-validation H0 effect should be interpreted until the bandwidth protocol is fixed and validated independently of H0 coverage.
4. Bandwidth selection should use training simulation seeds and a declared density-recovery criterion, then be frozen and evaluated on held-out seeds and parent laws.

## 7. API improvement added after this CI run

The KDE function now has a separate `scott_multiplier` argument, designed to multiply the weighted Scott factor computed from the original unreflected sample. The legacy `bandwidth_factor` remains an absolute SciPy covariance factor. The two explicit options are mutually exclusive. Tests were added for equivalence at multiplier 1, the weighted effective-sample-size calculation, mutual exclusion, and invalid multipliers.

These new source and test commits are **not covered by run #42455270698**. Do not mark them validated until a fresh workflow or equivalent targeted test run completes.

## 8. Project decision

- **GO:** independently audited finite-validation probability-stage software pilot.
- **GO:** independently audited numerical H0-grid stability conditional on current host densities.
- **CONDITIONAL:** distance-noise toy sensitivity and matched-model PIT diagnostics, within their stated inferential scope.
- **NO-GO:** interpreting the current factor-0.8 host-density H0 benchmark as a calibrated reference.
- **NO-GO:** estimated-inclusion H0 calibration/coverage claims; estimated-inclusion H0 propagation remains unimplemented.
- **NEXT:** validate the new bandwidth API, design training/held-out density-recovery bandwidth selection, then re-run oracle and known-inclusion H0 baselines before connecting finite-validation estimates to H0.

No default-branch changes or PR merges were made.


## 9. Additional independent reconstruction performed during review

After downloading the Python 3.12 artifact, I independently re-read the 1,800,000-row raw probability-stage CSV and reconstructed all 450 experiment/model groups directly from the unit records. This calculation used a separate analysis script and did not import the repository benchmark or auditor.

Checks:
- Every group had exactly 4,000 parent-frame records.
- Validation count matched its declared size (100, 300 or 1,000); all remaining objects were partitioned into analysis-catalogue or non-retained holdout roles.
- All 450 fitted-model groups were marked SUCCESS.
- Held-out Brier score, log loss, probability RMSE, 10-bin ECE, maximum calibration-bin gap, and the five estimated-IPW weight diagnostics were independently recomputed.
- The independently recomputed metrics agreed with the stored summary to the following maximum absolute differences:

| Metric | Maximum absolute difference |
|---|---:|
| Held-out Brier score | 1.11e-16 |
| Held-out log loss | 1.50e-15 |
| Held-out probability RMSE | 1.11e-16 |
| 10-bin ECE | 1.39e-16 |
| Maximum 10-bin gap | 1.67e-16 |
| Kish ESS | 1.79e-12 |
| Entropy ESS | 6.59e-12 |
| Maximum normalized weight | 2.18e-15 |
| Top-1% weight mass | 1.45e-14 |
| Top-5% weight mass | 1.31e-14 |

The raw partition totals were 210,000 validation records, 452,824 analysis-catalogue records, and 1,137,176 non-retained holdout records. The manifest SHA-256 hashes for both raw and summary CSV files were also independently recomputed and matched.

This materially strengthens confidence in the probability-stage artifact's numerical integrity. It does not change the pilot's scientific status: five replicates per cell remain exploratory, and the probability-stage results still do not establish downstream H0 calibration or estimator superiority.


## 10. Descriptive structure in the five-replicate probability pilot

The independently reconstructed summary also permits a descriptive check of finite-sample model-complexity effects. Across all retention targets and five replicates per cell, held-out RMSE by curvature and validation size was:

| True curvature | n=100: linear / quadratic | n=300: linear / quadratic | n=1,000: linear / quadratic |
|---:|---:|---:|---:|
| -2.0 | 0.0598 / 0.0671 | 0.0365 / 0.0442 | 0.0238 / 0.0182 |
| -1.2 | 0.0627 / 0.0743 | 0.0304 / 0.0404 | 0.0187 / 0.0188 |
| 0.0 | 0.0581 / 0.0766 | 0.0275 / 0.0418 | 0.0130 / 0.0158 |
| +1.2 | 0.0591 / 0.0712 | 0.0351 / 0.0449 | 0.0198 / 0.0189 |
| +2.0 | 0.0566 / 0.0705 | 0.0386 / 0.0421 | 0.0269 / 0.0191 |

Across these pilot cells, the quadratic fit's held-out RMSE was descriptively higher at n=100 and n=300, while the gap narrowed or sometimes reversed at n=1,000. At zero curvature, the simpler linear model also performed better in this pilot, consistent with the extra quadratic coefficient adding finite-sample variance when it is unnecessary. This is a useful candidate bias-variance pattern, not a reliable model-selection conclusion: each gamma × validation-size × retention × degree cell has only five replicates.

Because validation objects are excluded from the analysis catalogue, larger validation samples also leave fewer analysis objects. Absolute ESS therefore should not be compared without accounting for the analysis-catalogue size. The mean Kish ESS fraction (Kish ESS divided by analysis-catalogue size) was:

| Retention target | n=100 | n=300 | n=1,000 |
|---:|---:|---:|---:|
| 0.10 | 0.750 | 0.827 | 0.854 |
| 0.25 | 0.860 | 0.881 | 0.890 |
| 0.50 | 0.884 | 0.938 | 0.945 |

The corresponding mean analysis-catalogue sizes fell as validation size increased (for example, at target retention 0.10: 406, 384 and 310). Thus the rising ESS fraction indicates less concentrated estimated weights relative to the remaining catalogue in these cells, not a larger absolute effective sample. This is descriptive and still needs more replicates and held-out scenario validation.
