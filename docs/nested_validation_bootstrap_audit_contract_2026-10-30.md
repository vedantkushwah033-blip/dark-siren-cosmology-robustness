# Nested validation-bootstrap comparison: audit contract

**Status:** exploratory implementation protocol; the integrated CI run is a smoke test, not confirmatory evidence.  
**Artifact schema:** `nested_validation_bootstrap_comparison_smoke_v1`.

## Scientific question and estimand

The benchmark asks whether resampling procedures provide useful uncertainty estimates for the *posterior mean functional* induced by finite validation sampling in a deliberately simplified dark-siren pipeline. For each synthetic parent scenario, the reference target is the posterior mean obtained by fitting the same inclusion/host-weighting pipeline to every record in that scenario's actual finite parent frame. It is not a superpopulation or pseudo-population target.

For scenario s, validation size n, and outer SRSWOR draw r, let θ_hat(s,n,r) be the posterior mean from the sampled validation records. The finite-frame target is θ_N(s). The empirical outer-draw spread estimates repeated-sampling variability conditional on that generated parent scenario. The inner resampling methods estimate uncertainty from one observed validation sample. Scenarios, not their nested outer draws, are the top-level independent simulation units.

## Compared procedures

1. **Paired record bootstrap:** resample validation redshift/outcome pairs with replacement and refit the inclusion model, weights, host density, and toy posterior on each replicate.
2. **Pseudo-population SRSWOR comparator:** form a multinomial pseudo-population of size N from the observed validation records, then sample without replacement. Its percentile summaries are exploratory only. They are not the studentized finite-population-mean procedure studied by Wang, Peng & Kim (2022), and that paper's theory does not validate this nonlinear posterior functional.
3. **Scalar FPC diagnostic:** multiply the paired-bootstrap standard error by sqrt(1-n/N), with a normal interval centered at the outer posterior mean. This is a diagnostic comparator, not a theoretically justified interval for the nonlinear estimator.

The bootstrap standard errors/intervals concern repeated-sampling variation of posterior means. They are not the conditional H0 posterior credible intervals produced within an individual fit.

## Independent artifact-audit requirements

The auditor must be a separate implementation. It may share data structures and stable public interfaces, but must not call the producer's scenario generator wrapper, posterior wrapper, or calibration summarizer as the authority for the checks it claims to independently reconstruct.

### A. Manifest and design integrity

- Validate schema/version, smoke-only status, model constants, grids, scenario count, outer-draw count, parent size, validation sizes, bootstrap replicate count, and declared method set.
- Reject missing, duplicate, reordered, or extra scenario/size/draw/method keys.
- Require all methods to share the same outer status, outer posterior mean, and finite-parent target for a given outer draw.
- Require the finite-parent target to be constant across all validation sizes and outer draws within a scenario.
- Recompute all seed derivations from the root seed and documented SeedSequence construction.

### B. Parent-frame and sampling reconstruction

- Regenerate the parent frame, observed GW-distance sample, and catalogue thinning from recorded seeds and model parameters.
- Independently recompute parent, catalogue, and distance fingerprints using a documented canonical byte representation; do not accept hashes merely because they are consistent across rows.
- Recompute the full-parent reference posterior mean from all actual parent records.
- Parse each stored outer validation index vector; verify integer type, length, bounds, uniqueness, and nested SRSWOR properties across sizes for the same draw.
- Independently refit the outer posterior and compare it with the stored value/status. Verify one-class and numerical/model-fit failures explicitly.

### C. Inner-resampling replay

For each method and attempted outer draw:
- Verify the deterministic bootstrap seed and exact number of attempts.
- Reconstruct every paired-bootstrap index vector or pseudo-population and selected-position vector from the saved seed and the declared algorithm.
- Verify dimensions, integer ranges, replacement/without-replacement rules, and the mapping from selected pseudo-population positions to original validation records.
- Independently recompute each successful replicate's posterior mean and reference mean; compare attempt-aligned arrays, success masks, failure classes/codes/messages, and successful counts.
- Require failed-attempt numerical outputs to be NaN and require all persisted indices and masks to align with the attempt count.
- Recompute standard errors, percentile endpoints, and the scalar FPC diagnostic from the attempt-level values. Do not infer missing attempt values from aggregate summaries.

### D. Summary reconstruction

Rebuild the calibration records from the attempt audit and independently calculate:
- empirical outer SRSWOR SD from unique successful outer draws, once per scenario/size (not separately by bootstrap method);
- mean inner-bootstrap SE, SE/outer-SD ratio, finite-frame target coverage, interval width, and separate outer/inner failure rates;
- scenario-balanced summaries using equal scenario weighting, with explicit denominators and missing-cell handling.

Compare every stored summary cell to the independent reconstruction. A summary table that merely reconciles row counts is not sufficient.

### E. Adversarial tests

The test suite must mutate one item at a time and confirm rejection for: root seed; scenario seed; model parameter; manifest grid; parent/catalogue/GW hash; outer index; non-nested sample; bootstrap seed; bootstrap index; pseudo-population mapping; success mask; attempt posterior mean; reference mean; failure class/code/message; finite-parent target; outer posterior mean; method coverage; duplicate/missing row; SE; interval endpoint; FPC factor; and each summary statistic.

## Interpretation and stopping rules

- The present smoke configuration (one scenario, two outer draws, three bootstrap replicates) checks execution and file contracts only. It cannot estimate coverage, calibration, or method ranking.
- Do not select a winning method from the smoke test.
- Before a confirmatory run, conduct a pilot to estimate failure rates and Monte Carlo variability. Choose scenario count and outer/inner replication counts using prespecified precision targets for the primary calibration contrasts, with scenario-level uncertainty and compute cost reported.
- Freeze the estimands, primary contrasts, handling of failures, and precision targets before the confirmatory run. Retain the full attempt-level audit trail.
- Report null, adverse, and inconclusive outcomes. Do not suppress failed fits or change the target after seeing results.

## Scope limitations

This is a synthetic low-redshift Euclidean toy model, not an LVK/GWTC analysis, not a new observational H0 measurement, and not a resolution of the Hubble tension. The observed catalogue and GW distances are held fixed within each generated parent scenario. Claims about real survey selection, galaxy catalogues, or astrophysical populations require a separate observationally grounded analysis and validation.

## Literature note

Wang, Zhonglei, Liuhua Peng, and Jae Kwang Kim (2022), “Bootstrap Inference for the Finite Population Mean under Complex Sampling Designs,” *Journal of the Royal Statistical Society: Series B* 84(4), 1150–1174, DOI: 10.1111/rssb.12506. The paper motivates careful design-aware resampling for finite-population means; it does not establish validity for the nonlinear logistic/KDE/selection-posterior functional used here.
