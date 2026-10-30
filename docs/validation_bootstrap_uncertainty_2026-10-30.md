# Validation-sample bootstrap for estimated inclusion uncertainty

**Date:** 30 October 2026  
**Status:** implementation and unit tests added; full CI validation pending.  
**Scope:** conditional synthetic host-density uncertainty diagnostic.

## Motivation

The estimated-inclusion pilot fits a logistic retention model to an external validation sample, predicts inclusion probabilities for the retained catalogue, and uses inverse probabilities in a weighted KDE. A single fitted model suppresses the fact that a finite validation sample could produce materially different coefficients and host densities.

The new `bootstrap_validation_inclusion_density` function resamples the paired validation observations (redshift, retained/not-retained label) with replacement. For every resample it refits the logistic inclusion model and rebuilds the weighted host density, holding the observed catalogue fixed.

## What the bootstrap estimates

The output approximates repeated-validation-sample variability using an ordinary paired-record bootstrap with replacement. The simulated validation design samples without replacement from a finite parent frame; this implementation does **not** apply a finite-population correction. Its uncertainty estimate is therefore an iid-bootstrap approximation, not an exact design-based variance estimator. This distinction matters when the validation fraction is large. A future design-based replicate-weight method would be needed for calibrated finite-population uncertainty.

The output approximates the repeated-validation-sample variability of:
- the fitted logistic intercept and slope;
- the inverse-inclusion-weighted host-redshift density;
- boundary-limited fit frequency and bootstrap resample failure frequency.

The default bandwidth convention is recomputed in each resample using the existing weighted Scott rule. Supplying `bandwidth_factor` instead fixes the physical bandwidth multiplier while the weighted sample spread can still vary.

Every attempted resample is accounted for. A resample with only one retention class or an invalid fit is recorded in `failure_reasons`; it is not silently replaced. Coefficient-boundary-limited fits remain successful fits but are separately flagged.

## Interpretation limits

This is an **approximate nonparametric frequentist bootstrap conditional on the observed catalogue**,  the original validation sample, the assumed logistic selection family, the redshift domain, and the KDE construction. It does not apply a finite-population correction. It is not:
- a Bayesian posterior for the inclusion coefficients;
- a correction for logistic-model misspecification;
- uncertainty from a newly sampled galaxy catalogue or GW-event realization;
- uncertainty in GW distances, detection efficiency, cosmology priors, or the population model;
- a coverage guarantee for the final (H_0) interval.

The bootstrap density draws are ingredients for a downstream sensitivity analysis. If passed through the same fixed GW likelihood and posterior calculation, the resulting distribution of posterior summaries should be labelled a validation-bootstrap sensitivity distribution—not a posterior credible interval and not a replacement for a joint hierarchical model.

## Required reporting

For any study using this function, report:
1. validation design and sample size;
2. bootstrap seed and number attempted;
3. successful, failed, and boundary-limited resamples;
4. fitted-coefficient and host-density variability;
5. whether KDE bandwidth is fixed or recomputed;
6. the fixed catalogue/event assumptions;
7. the exact downstream transformation, if any, from density draws to posterior summaries.

Do not discard high failure rates. A large fraction of degenerate resamples is itself evidence that the available validation sample may be too small or poorly balanced for stable selection fitting.

## Current status and next gate

The reusable density-bootstrap utility and its deterministic/failure-accounting tests are implemented. A second helper, `run_validation_bootstrap_posterior_sensitivity`, passes each successful density draw through the fixed selection-conditioned toy GW likelihood and records the resulting posterior means, medians, and within-draw 68% intervals. Tests cover this end-to-end propagation and reject events below the detection threshold. These helpers are on branch `research/validation-bootstrap-uncertainty`. The next gate is a successful multi-version CI run, followed by persistence of bootstrap summaries in the pilot artifact schema and an independent auditor that reconstructs bootstrap refits and posterior summaries from saved primitives. No bootstrap result or cosmological uncertainty claim is asserted by this protocol alone.
