# Research-track and estimand separation decision

**Date:** 2026-11-05  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Status:** scope-control decision; not a preregistration and not a confirmatory result.

## 1. Why this record is necessary

The repository currently contains at least two related but statistically distinct research tracks. They share a synthetic dark-siren pipeline and some software components, but their targets, repeated-sampling units, uncertainty procedures, and scientific claims are not interchangeable.

The presence of both tracks in one repository must not be interpreted as evidence that one experiment validates the other.

## 2. Track A — inclusion-model misspecification and downstream H0 inference

**Working title:** *Finite-Validation and Model-Misspecification Effects in Inverse-Probability Host Reconstruction for Synthetic Dark-Siren Cosmology.*

**Primary question:** Under what conditions does estimating catalogue-inclusion probabilities from finite validation data improve host-population reconstruction and downstream H0 inference, and when do finite-sample uncertainty, model misspecification, and inverse-weight instability remove or reverse that benefit?

**Core estimands and outcomes (to be frozen before confirmatory simulation):**
- inclusion-probability prediction error against the known synthetic truth;
- host-density reconstruction error against the known parent-frame target;
- weight concentration and support diagnostics;
- posterior bias/error and interval coverage for the generating H0, with the exact repeated-sampling hierarchy stated;
- optimizer, boundary, and analysis failure rates, retaining all attempted replicates.

The estimator arms (for example oracle, estimated linear, and estimated nonlinear inclusion) must be defined precisely. Oracle probabilities are a diagnostic reference, not an implementable estimator and not a guaranteed downstream winner.

The current v1/v2 solver comparison is a software/estimator implementation comparison. A two-replicate smoke run can establish paired-input provenance and execution only; it cannot establish calibration, robustness, or superiority.

## 3. Track B — calibration of uncertainty for a finite-frame posterior-mean functional

The existing `manuscript/finite_validation_bootstrap_draft.md` describes a different primary target:
- construct a finite-frame reference posterior mean by applying the prespecified pipeline to all records in a generated finite parent frame;
- repeatedly draw validation samples by SRSWOR;
- evaluate inner resampling procedures for uncertainty in that finite-frame posterior-mean functional.

Its primary coverage target is the finite-frame reference posterior mean. It explicitly is **not** coverage of the generating H0 and is **not** the conditional credible interval from an individual posterior.

The pseudo-population/SRSWOR comparator and scalar FPC multiplier are exploratory/diagnostic comparators for this nonlinear functional. Finite-population-mean theory does not automatically validate them for the nonlinear posterior-mean functional.

## 4. Non-combinability rule

Do not pool Track A and Track B replicates, tables, coverage estimates, or conclusions.

In particular:
1. Coverage of a finite-frame posterior mean is not coverage of the generating H0.
2. A bootstrap standard error for the posterior mean is not a posterior credible interval.
3. A linear finite-population-mean benchmark does not validate nonlinear dark-siren posterior inference.
4. A passing numerical derivative or solver test does not establish repeated-sampling calibration.
5. Results conditional on a fixed catalogue and fixed detected-distance sample do not include all data-generation uncertainty.
6. A finite validation sample and a finite parent-frame sampling design must be described separately; do not imply SRSWOR unless that is the actual design.

## 5. Repository and manuscript decision

- Keep Track A and Track B as separate manuscript/result packages, even if they share tested library functions.
- The existing finite-validation bootstrap draft remains Track B unless a documented, scientifically justified scope change is approved.
- Track A needs its own manuscript file and its own preregistered estimand/design document before confirmatory execution.
- Shared code may be reused only with explicit tests that establish the relevant data-generating and estimand semantics.
- A result from one track may be cited as a software validation or methodological comparator in the other only when the precise scope of that evidence is stated.

## 6. Current release gate

As of this record, the latest known Track A branch head is `af95a843ef3f8d5f5dda28249ad88cc87c871de2`. Its push-triggered GitHub Actions matrix was observed in progress at the time of inspection; no green exact-head result is claimed here.

No confirmatory Track A experiment is authorized by this record. Before one is run, freeze:
- generating scenarios and validation-sample design;
- estimator arms and their exact parameterization/bounds;
- primary and secondary estimands;
- failure denominators and failure policy;
- positivity/support and weight-diagnostic reporting;
- posterior grid/prior policy and sensitivity plan;
- repetition counts justified by pilot-based Monte Carlo precision;
- independent artifact-audit requirements.

Historical optimizer incident D-01 remains unresolved unless original contemporaneous evidence is recovered. A new solver test cannot retrospectively identify the historical failure.

## 7. Permitted conclusion

The repository supports development of reproducible synthetic benchmarks and mathematical/software diagnostics. It does not yet support a claim that estimated inclusion weighting is calibrated or superior for dark-siren H0 inference, nor a claim of observational validation or a new cosmological measurement.
