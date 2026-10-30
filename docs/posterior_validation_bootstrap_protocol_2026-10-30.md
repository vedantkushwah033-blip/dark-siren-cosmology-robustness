# Posterior-level validation-bootstrap protocol (30 October 2026)

## Purpose and scope

The module `dark_siren_robustness.posterior_validation_bootstrap` propagates uncertainty in an external validation sample through the existing synthetic host-weighting and selection-aware H0 toy-inference pipeline.

For every resample it:
1. resamples paired validation redshift and catalogue-retention records;
2. refits the bounded logistic inclusion model;
3. recomputes inverse-inclusion weights and the reflected weighted host KDE;
4. recomputes the selection-aware toy likelihood for the fixed observed GW-distance sample;
5. normalizes the H0 posterior on the supplied H0 grid and stores its posterior mean, median, and 68% credible endpoints.

The observed incomplete catalogue and observed GW-distance values are held fixed. Thus this is a conditional uncertainty-propagation experiment, not total uncertainty for a real gravitational-wave analysis.

## Methods exposed

- `paired_bootstrap`: draws n paired validation records with replacement from the observed n-record validation sample.
- `design_aware_srswor`: generates a multinomial pseudo-population of size N from the n observed paired records, then draws n positions without replacement from that pseudo-population. It stores the pseudo-population mapping, sampled positions, and mapped paired-record indices.

The second method is an exploratory adaptation of the finite-population bootstrap construction for SRSWOR discussed by Wang, Peng & Kim (2022), DOI 10.1111/rssb.12506. Their theoretical results concern finite-population means; they do not establish validity for this logistic/KDE/selection-aware posterior functional.

## Estimands must not be conflated

For a given fixed catalogue and fixed GW-distance sample:

- The **true simulated H0** is the data-generating value. It is used for bias and coverage evaluation in repeated simulated datasets, not as a bootstrap center.
- The **full observed-validation plug-in posterior mean** is the reference used by the ordinary paired bootstrap. It is an estimator based on the observed validation sample, not the unknown finite-parent target.
- The **pseudo-population plug-in posterior mean** is recomputed in each design-aware replicate from that replicate's generated size-N pseudo-population. It is a replicate-specific bootstrap reference, not the true finite-parent posterior mean.
- The **posterior 68% credible interval** is conditional on a particular fitted host density and toy likelihood. It is not a confidence interval for the posterior mean across validation resamples.
- The **bootstrap percentile interval and bootstrap standard error** summarize the resampling distribution of posterior means. They must not be reported as H0 posterior credible intervals.

The finite-frame target (T_{s,N}) remains a separate quantity: it requires the same estimator fit using all N actual finite-parent validation records for scenario s. The pseudo-population reference is not a substitute for that target.

## Stored diagnostics

The function returns attempt-aligned posterior means, medians, credible endpoints, reference posterior means, inverse-weight ESS, maximum inverse weight, coefficient-boundary flags, success mask, exact paired-record indices, pseudo-population indices and selected pseudo-population positions. Failures are retained with attempt index, exception type and reason.

For successful draws, the full normalized posterior density on the supplied H0 grid is stored. The function checks that each stored density is finite, nonnegative, grid-aligned and integrates to one.

## FPC comparator

The returned (sqrt{1-n/N}) factor and FPC-adjusted standard error are diagnostic comparators only. Multiplying a nonlinear posterior-mean bootstrap standard error by the scalar SRSWOR FPC is not assumed to be valid. The experiment must compare these values against empirical repeated-SRSWOR variability before drawing conclusions.

## Confirmatory experiment still required

This module is implementation infrastructure, not a result. Before confirmatory use:

1. freeze the scenario-generating law, finite parent, catalogue thinning, GW-distance sample, validation sizes, H0 and redshift grids, seeds, KDE rule, and failure policy;
2. compute the actual finite-parent reference (T_{s,N}) for each independently generated scenario;
3. generate repeated independent SRSWOR validation samples from each fixed parent;
4. apply ordinary bootstrap, FPC comparator, and design-aware pseudo-population bootstrap to each outer sample;
5. compare estimated standard errors and interval coverage against the empirical outer-sample distribution and (T_{s,N});
6. independently regenerate all indices, refit every successful attempt, reconstruct posteriors and summaries, and replay failures;
7. report scenario-level uncertainty, failure rates, boundary-limited fits, ESS and extreme weights.

The pilot and confirmatory runs must be separated. No numerical conclusion, design-consistency claim, or general dark-siren inference claim follows from the implementation alone.

## Reproducibility test coverage

`posterior_validation_bootstrap_audit.py` independently replays the random index generation, refits the bounded logistic model, reconstructs the reflected weighted KDE and selection-aware posterior without calling the producer estimators, and checks posterior grids, summaries, diagnostics, FPC arithmetic, and failure labels. `tests/test_posterior_validation_bootstrap.py` checks deterministic replay for both methods, posterior normalization, paired-record mappings, independent reconstruction, and rejection of a tampered posterior mean. CI runs this test module on Python 3.10, 3.11 and 3.12.


## Independent failure-classification audit (30 October 2026)

The posterior-level auditor now reconstructs the exception class and exact
message for every failed resampling attempt, then compares both fields with
the persisted failure record. A failure label is therefore not accepted
merely because the independently replayed attempt also failed for some
possibly unrelated reason. The independent reconstruction also enforces the
same two-retention-class precondition as the producer's logistic inclusion
fit; this prevents a one-class resample from being misclassified as a
successful bounded logistic fit by the auditor.

Adversarial regression tests deliberately create one-class bootstrap
resamples, verify that the valid artifact audits, and then require the auditor
to reject a tampered exception class and a tampered failure message. Exact
exception messages can be sensitive to numerical-library versions, so these
checks are strongest when producer and auditor run in the same pinned
environment; a future artifact schema should consider stable structured
failure codes in addition to diagnostic messages.


## Failed-attempt missingness invariants (30 October 2026)

The independent posterior-bootstrap audit now also checks that attempt-aligned
posterior means, pseudo-population/reference means, credible endpoints, and
inverse-weight diagnostics remain NaN for failed attempts. This prevents a
failed replicate from carrying plausible-looking numerical values that could
be accidentally included in downstream summaries. A regression test tampers
with a failed attempt's posterior mean and requires the audit to reject the
artifact. The producer already initializes these fields as NaN; this change
makes that missingness contract independently enforceable.


## Failure record schema hardening (30 October 2026)

Every failed posterior-resampling attempt now stores both its exact exception class/message and a stable broad `failure_code`:
- `floating_point_failure`
- `arithmetic_failure`
- `invalid_validation_sample`
- `model_fit_failure`
- `unexpected_exception`

The independent auditor replays the attempt, derives the category from the independently reconstructed exception class, and requires an exact code match as well as the existing exact class/message match. The code is intended for stable aggregation across Python/library versions; the detailed reason remains the diagnostic record. Failed attempts remain in the attempted-run denominator.


## Input provenance fingerprints (30 October 2026)

Posterior-level resampling artifacts now include SHA-256 fingerprints for the canonical float64 representations (including array shape) of the redshift grid, H0 grid, observed catalogue redshifts, validation redshifts, binary validation outcomes, observed GW distances, and merger-weight array. When no merger-weight array is supplied, the fingerprint is computed from the explicit all-ones array used by the auditor.

The independent auditor reconstructs these fingerprints from its separately supplied inputs and rejects a mismatch before replaying RNG draws. This binds a result artifact to the exact numerical inputs used in the audit; it is a provenance/integrity check, not a cryptographic signature or a substitute for independently preserving the source data and code revision.


The addition of mandatory input fingerprints is recorded as artifact schema `validation_h0_posterior_resampling_v2`; schema-v1 artifacts do not satisfy the new provenance contract and must not be silently treated as v2.
