# Finite-validation probability pilot: independent analysis

**Date:** 2026-11-08  
**Branch:** `research/independent-posterior-audit-claim-fix`  
**Source run:** GitHub Actions workflow `42446224532`, commit `a013e00819de67e42d553c7ae95abc3a5bee0959`  
**Artifact:** `numerical-diagnostics-py3.12`, containing `finite_validation_probability_pilot/`

## Scope and audit status

The configured pilot contains 450 fit attempts: 5 curvature values, 3 target parent-retention fractions, 3 validation sizes, 2 polynomial degrees, and 5 replicates. Each attempt has 4,000 raw parent-frame records, for 1,800,000 raw records total. The independent auditor reports PASS after reconstructing all 450 experiments, the complete parent-frame partitions, the calibrated truth probabilities, fitted probabilities, score equations, observed-information covariance, held-out scores, and IPW diagnostics.

The run manifest records Python 3.12.16, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.1, a 37.03-second runtime, and SHA-256 hashes for the raw and summary files. The artifact is a five-replicate **implementation pilot**, not a confirmatory performance study.

## Fit outcomes

- 448/450 fits succeeded.
- 2/450 were marked FAILED because BFGS reported precision loss and the score infinity norm was slightly above the fitter's (10^{-6}) acceptance threshold.
- Both failures had well-conditioned observed information (condition numbers 2.36 and 8.37) and moderate coefficients; neither showed evidence of separation from the stored diagnostics.
- Re-running damped Newton iterations from the validation data for those two cases reached score infinity norms (6.81\times10^{-11}) and (9.70\times10^{-13}), respectively.

The last point motivates a numerical polishing fallback, now added to the fitter on this branch. A regression test simulating a BFGS precision-loss return was also added. **These two source changes have not yet been verified by a new GitHub Actions run.** They must not be described as CI-passed until a new run exists.

## Descriptive probability-stage results

Across successful fits, mean held-out probability RMSE was:

| Validation size | Mean held-out RMSE |
|---:|---:|
| 100 | 0.0656 |
| 300 | 0.0381 |
| 1,000 | 0.0192 |

This is descriptively consistent with better probability prediction at larger validation sizes. It is not an estimate of posterior calibration or H0 coverage.

Averaged over all curvature/retention/validation cells, the quadratic model had a higher held-out RMSE than the linear model (0.0444 versus 0.0378). This aggregate comparison is confounded by the deliberately balanced mix of correctly specified and misspecified scenarios and the small five-replicate pilot. Within the 1,000-unit validation cells, the quadratic model did better for the stronger curvature settings, while it generally paid a variance penalty at validation sizes 100 and 300. This is the expected bias-variance trade-off to investigate with more replicates, not evidence that linear models are generally preferable.

Mean estimated-weight Kish ESS was 956.6, 958.2, and 785.4 for validation sizes 100, 300, and 1,000, respectively. These pooled numbers are not monotonic because the fitted model and validation size affect weights, while the analysis catalogue and true retention conditions also vary. ESS is a concentration diagnostic, not a measure of whether the weighted density is correct. Estimated ESS exceeded oracle ESS in some individual paired cells; this is not evidence of estimator superiority.

Across retention targets, the average held-out RMSE was 0.0292, 0.0398, and 0.0542 for target retention 0.10, 0.25, and 0.50. These cross-target averages should not be interpreted as a causal effect of retention: the proper analysis is paired and stratified by curvature, validation size, model degree, and replicate.

## Interpretation limits

1. Five replicates per cell are insufficient for stable rankings, tail probabilities, or calibration claims.
2. Brier score and log loss score prediction of the binary retention label; probability RMSE scores recovery of the simulated inclusion probability. They answer related but different questions.
3. Kish/entropy ESS and top-weight mass describe the concentration of inverse weights. They do not measure bias in the reconstructed parent density.
4. The inverse observed-information matrix is a model-based asymptotic covariance estimate. It is **not** a posterior distribution for the logistic coefficients and is not yet justified as uncertainty propagation under this finite-population validation design.
5. The probability-stage benchmark does not yet build estimated-probability host KDEs, integrate coefficient uncertainty into the event likelihood, or assess H0 posterior coverage.
6. The independent audit establishes numerical reconstruction of saved records and summaries. It does not validate the scientific adequacy of the synthetic DGP or likelihood.

## Required next experiments

- Increase replication only after the current pilot is fully analyzed; preserve common random numbers and paired scenario identifiers.
- Report per-cell fit failures, held-out probability RMSE, Brier/log loss, probability calibration, coefficient bias/variance, Hessian eigenvalues, weight quantiles, ESS fractions, and top-1%/5% weight mass.
- Compare degree 1 and degree 2 within the same ((\gamma, r, n_{val}, replicate)) cells, not only pooled means.
- Add an explicit attempt ledger with SUCCESS, NON_EVALUABLE, FAILED, and NOT_RUN states before any long campaign.
- Derive a design-consistent uncertainty model for the validation coefficients. Do not treat the stored inverse-information covariance as a Bayesian posterior without an explicit approximation and validation.
- Only then implement a modular uncertainty-propagation experiment: for coefficient draws from a justified (p(\theta\mid\text{validation})), reconstruct the IPW host density, evaluate the selection-conditioned event likelihood for each draw, integrate the likelihood over coefficient uncertainty, and normalize the H0 posterior once. This is distinct from averaging separately normalized posteriors or treating a frequentist bootstrap distribution as a posterior.
- Run SBC and coverage on the estimated-inclusion arm and its paired oracle reference, with fixed and documented prior support and boundary-truth handling.

## Current decision

**GO:** probability-stage software and raw-record audit, subject to a fresh CI run for the newly added Newton polish.  
**NO-GO:** claims of estimated-inclusion benefit, H0 calibration degradation, calibrated coverage, or a validated uncertainty-propagation method.

The defensible project claim remains a planned synthetic benchmark of finite-validation uncertainty and misspecification. No posterior-calibration result has yet been established for estimated inclusion probabilities.


## Post-pilot Newton-polish verification (local, 2026-11-08)

After adding the damped-Newton fallback and its regression test:

- Ran the seven tests in the isolated local copy of `tests/test_estimated_inclusion.py`: **7 passed**.
- Refit all 450 validation groups from the downloaded raw pilot records using the updated fitter: **450 SUCCESS, 0 FAILED, 0 NON_EVALUABLE**.
- The two prior precision-loss failures were recovered. Their final score infinity norms were (2.61\times10^{-10}) and (3.29\times10^{-14}).
- For the other 448 groups, the maximum absolute difference between updated predictions and the predictions stored by the original pilot was (6.68\times10^{-9}).

This is a targeted local verification of the estimator module and the 450 pilot validation subsets. It is **not** a full repository test-suite run, a rerun of the benchmark/auditor pipeline, or GitHub Actions validation of the new commits. The original artifact's independent audit remains valid for the original BFGS outputs; a new full workflow artifact is still required to verify the modified source end to end.


## Exploratory reliability-bin calibration (post-hoc raw-record analysis)

I additionally calculated equal-frequency, 10-bin reliability diagnostics from the 448 successful groups in the downloaded raw artifact. The auditor's independent reconstruction algorithm was separately checked against the scoring implementation on all 448 groups.

| Validation size | Mean 10-bin ECE | Mean maximum bin gap |
|---:|---:|---:|
| 100 | 0.0562 | 0.1213 |
| 300 | 0.0370 | 0.0798 |
| 1,000 | 0.0238 | 0.0553 |

These descriptive values decrease with validation size in this pilot. ECE is binning-dependent and sample-dependent; five replicates per cell are not enough to establish a general calibration law. It is also not H0 posterior calibration.

The current branch now adds these reliability-bin metrics to the next probability-stage summary and independently recomputes them from raw held-out labels/predictions in the auditor. The corresponding diagnostics tests were added. **These latest calibration-metric source changes have not yet been run through the full repository suite or GitHub Actions.**
