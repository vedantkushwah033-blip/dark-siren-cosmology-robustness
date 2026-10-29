# Preregistered follow-up: estimated catalogue inclusion probabilities from an external validation sample

**Date:** 29 October 2026  
**Status:** design proposal; NOT run; protocol requires implementation review before confirmatory execution  
**Branch context:** research/audit-posterior-grid-integrity

## 1. Why this is the next experiment

The existing host-density sensitivity screen already varies KDE bandwidth across 13 configurations, 40 paired replicates per configuration, and 1,560 fit rows, with independent reconstruction reported in the screening report. Repeating a four-arm bandwidth experiment at one truth value would largely duplicate existing work. More importantly, the existing known-completeness correction receives the true inclusion probabilities. That does not answer how weighting behaves when inclusion probabilities must be estimated.

This follow-up targets the more consequential gap:

> How much external catalogue-completeness validation data are needed for inverse-inclusion weighting to retain useful calibration, and when does estimation error or weak inclusion support make the correction unstable?

This is a synthetic statistical benchmark, not a realistic galaxy-survey model or observational H0 analysis.

## 2. Identifiability boundary

The inclusion function cannot generally be identified from the observed catalogue alone without additional assumptions or external information: the observed redshift distribution is proportional to the latent host distribution multiplied by the inclusion probability. Different latent-density/inclusion-function pairs can produce the same observed distribution.

Therefore this protocol does **not** claim to estimate completeness from the incomplete catalogue alone. It introduces an explicitly defined, independently sampled validation/audit subset with reference redshifts and known catalogue-retention indicators. The external validation sample supplies information unavailable from the incomplete catalogue by itself.

## 3. Data-generating process

Use the repository's controlled low-redshift parent-population and GW-distance toy as the starting point, after confirming the current implementation and preserving its limitations.

For each synthetic replicate:

1. Generate a parent population of galaxies with redshifts on the declared domain using the existing latent intensity and merger-weight convention.
2. Generate catalogue-retention indicators from a prespecified, smooth redshift-dependent inclusion law. For this experiment, use a logistic model whose form matches the calibration model: logit(pi(z)) = a + b * standardized_z. Standardize redshift using fixed, declared domain constants, not replicate-specific truth.
3. Generate the incomplete catalogue using those indicators.
4. Generate GW-like hosts, noisy distances and threshold-selected events using the fixed declared toy mechanism. Keep the event-generation and selection assumptions unchanged across estimator arms.
5. Draw an independent simple-random validation sample of parent-frame objects. For each validation object, the audit provides a reference redshift and the binary indicator of whether it appears in the catalogue. This represents an idealized external completeness audit; it is not assumed to exist in a real survey automatically.
6. Do not provide the full parent population, the true inclusion probabilities, or unobserved non-validation outcomes to non-oracle estimators. Save them only in a separate truth-only artifact for scoring.

The validation sample must be drawn independently of the catalogue-retention outcome, or its sampling design must be explicitly incorporated into the fit. Do not use a validation sample selected only from retained catalogue objects.

## 4. Estimators / analysis arms

- **A0 — Catalogue-only:** observed catalogue redshifts with equal weights.
- **A1 — Oracle diagnostic:** true inclusion probabilities supplied by the simulator. This is an oracle-information reference, not a deployable estimator and not an unconditional performance upper bound.
- **A2 — Estimated inclusion, small audit:** fit the prespecified logistic inclusion model using the smallest validation sample.
- **A3 — Estimated inclusion, medium audit:** same method with a larger validation sample.
- **A4 — Estimated inclusion, large audit:** same method with the largest validation sample.
- **A5 — Deliberately misspecified audit model:** use a prespecified wrong functional form or omit a known covariate in a separate, clearly labeled negative-control arm. Do not tune it after observing outcomes.

Use audit sizes of 100, 300 and 1,000 as initial candidate levels only. A non-confirmatory pilot must establish whether the smallest audit size provides enough retained and non-retained examples across the redshift range for stable fitting. If it does not, revise sizes before freezing the protocol. Do not silently discard failed logistic fits.

The estimator may use only the audit records and observed catalogue records. Define how fitted probabilities are applied to observed catalogue objects. Predefine a positivity guard, probability clipping or weight cap if used; report its threshold and sensitivity. Do not select the cap based on which one gives the best H0 answer.

## 5. Primary and secondary estimands

### Primary estimand

For each audit size m, estimate the paired replicate-level difference in absolute posterior-mean error between the estimated-probability arm and catalogue-only arm:

Delta_m = mean_r( |H0hat_catalogue_only,r - H0_true| - |H0hat_estimated_m,r - H0_true| ).

Positive values favor estimated-probability weighting under this synthetic model. Also report the paired difference between estimated-probability weighting and the oracle diagnostic, but do not treat the oracle contrast as a deployable comparison.

### Secondary metrics

- Bias and RMSE of the posterior mean.
- Empirical 68% and 95% credible-interval coverage, with covered counts and Wilson intervals.
- Median credible-interval width.
- Fitted-probability calibration against true probabilities in truth-only evaluation.
- Mean absolute probability-estimation error, reported only as a diagnostic.
- Effective sample size of inverse weights, maximum weight, weight quantiles and truncation rate.
- Logistic-fit failure rate, non-finite posterior rate and boundary posterior mass.
- Density L1 error and grid-convergence diagnostics.

The primary H0 estimand and probability-estimation diagnostics answer different questions; better probability prediction does not guarantee better H0 performance.

## 6. Pairing, randomization and seed management

Use the same generated parent population, catalogue, GW-like observations and validation-sample frame across all estimator arms within a replicate. This is an intentional paired design. Use a separate deterministic seed namespace for each replicate component (parent generation, catalogue thinning, GW generation, validation sampling and analysis bootstrap), or derive child seeds from a documented SeedSequence.

Record:
- replicate_id;
- parent_seed, catalogue_seed, gw_seed, audit_sample_seed;
- validation sample size;
- true and fitted inclusion-model parameters;
- estimator arm and weight-cap rule;
- software revision and configuration hash.

For paired uncertainty, resample replicate IDs as paired units. Do not interpret overlapping seed labels as a defect by themselves; verify whether records are paired or independent and analyze them accordingly. For comparisons across audit sizes, use nested validation samples within a replicate only if that nesting is explicitly intended and recorded; otherwise use independently sampled audit subsets and state the target contrast.

## 7. Precision and replicate count

Do not adopt 400 replicates merely because it was proposed in the previous bandwidth design. First run a small non-confirmatory pilot using seeds disjoint from confirmation. Estimate the standard deviation of the primary paired replicate-level contrast D_r. Choose the confirmatory replicate count from a prespecified target confidence-interval half-width w, using the planning approximation:

R >= (1.96 * SD(D_r) / w)^2.

The pilot is for variance/runtime planning, not hypothesis selection. Inflate the planned count for failed runs only by rerunning failures according to a documented rule; never remove difficult replicates silently. Empirical coverage uncertainty should be reported separately using covered counts and Wilson intervals. If the chosen budget cannot support meaningful precision, narrow the claim rather than imply adequate power.

## 8. Negative controls and deterministic checks

1. **No-completeness-gradient control:** set inclusion probability constant over redshift. The estimated model should recover a flat inclusion relationship within expected sampling uncertainty; weighting should not invent a strong redshift trend.
2. **Oracle recovery control:** with true inclusion probabilities supplied, reconstruct the expected weighted host-density target under the declared DGP. Verify normalization and the direct weighted KDE calculation.
3. **Audit-label permutation control:** permute retention labels among validation objects while preserving validation redshifts. The fitted redshift-inclusion association should be disrupted; specify expected metrics before execution.
4. **Misspecified audit model:** omit a prespecified nonlinear component or covariate. Measure, rather than assume, the direction of resulting H0 error.
5. **Positivity stress:** vary the minimum true inclusion probability in a small prespecified sensitivity analysis. Report maximum weights, ESS, fit failures and uncertainty; do not imply inverse weighting can recover regions with zero inclusion probability.
6. **Selection-likelihood control:** retain the repository's existing test that the selected log likelihood equals the unconditional log likelihood minus N times log detection efficiency for the implemented fixed-N conditional likelihood. A rate-based likelihood would require a different derivation.

## 9. Required implementation and release artifacts

Each run must emit:
- Versioned JSON configuration.
- Raw replicate-level results including failures.
- Validation-sample records with only the fields visible to the estimator.
- Separate truth-only records.
- Fit diagnostics and logs.
- SHA-256 manifest with file sizes.
- Source commit and dependency versions.
- Independent audit output that reconstructs the logistic fit inputs, inclusion probabilities, KDE densities, posterior grids, paired contrasts and coverage counts from saved records.
- A machine-readable summary linked to each manuscript table/figure.

The independent auditor should not import the production helper that fits inclusion probabilities or computes the primary aggregate. It may share a small declared data schema, but should independently reconstruct the central calculations.

## 10. GO / NO-GO criteria

**Pilot GO** only if:
- the inclusion model is identifiable under the declared audit sampling design;
- validation samples include sufficient retention variation to fit the model;
- estimated probabilities remain finite and diagnostics expose weak support;
- the implementation never leaks full-parent truth into non-oracle arms;
- all raw records and seeds are persisted.

**Confirmatory GO** only after the protocol, estimator definition, weight handling, primary estimand, replicate count, negative controls and independent audit are frozen.

**NO-GO** if:
- the estimator is described as learning inclusion probabilities from the incomplete catalogue alone without identifying assumptions or external data;
- zero-inclusion regions are presented as recoverable by inverse weighting;
- failed fits are silently dropped;
- oracle results are presented as realistic estimated-weight performance;
- or the implementation cannot be independently reconstructed.

## 11. Scope and interpretation

This experiment can test a narrow statistical question about external completeness calibration in a synthetic model. It cannot establish performance for actual galaxy surveys, realistic clustering, detector selection, redshift errors, or LVK cosmology. A positive result would support only the stated model and audit-sampling design. A negative result would still be useful by quantifying how validation sample size and positivity affect weighting stability.

**Decision:** design proposal only. Existing KDE bandwidth sensitivity work already exists; prioritize estimated-inclusion uncertainty as the next distinct question. No simulation has been launched by this document.
