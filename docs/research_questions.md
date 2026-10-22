# Research questions and hypotheses

## Primary question

Under which controlled conditions do galaxy-catalogue incompleteness and host-population misspecification degrade the bias, uncertainty, or calibration of Bayesian dark-siren inference for the Hubble constant, and when can an explicit model of missing hosts improve inference?

## Questions

1. Does the inference implementation recover known answers in simplified cases?
2. How do random and distance-dependent catalogue inclusion mechanisms affect inference?
3. How do incorrect host-population weights affect bias, interval width, and coverage?
4. Do combinations of incompleteness and weight misspecification produce interactions beyond their separate effects?
5. Under what assumptions does a missing-population correction improve performance, and when does it fail?

## Working hypotheses (not results)

- H1: Some forms of catalogue incompleteness can change the effective host-redshift distribution and affect H0 inference if not modelled consistently.
- H2: The consequences of host-weight misspecification depend on how weights correlate with redshift, distance, and the simulated host population.
- H3: Explicitly modelling catalogue inclusion may improve calibration when the inclusion model is adequate, but can fail when that model is itself misspecified.

These hypotheses must be tested rather than assumed to be true.

## Scope and guardrails

- Begin with synthetic data and clearly stated assumptions.
- Distinguish catalogue selection from gravitational-wave detection selection.
- Separate posterior coverage, simulation-based calibration, bias, and precision.
- Do not describe synthetic results as observational measurements.
- Do not claim novelty until a literature review has established what is already known.
- Report negative and inconclusive results.
