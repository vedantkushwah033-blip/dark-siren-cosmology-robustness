# Paired degree-one / degree-two multi-replicate pilot

**Status:** preregistered engineering pilot; synthetic low-redshift Euclidean model only.

## Question

Across matched synthetic parent populations, detected GW events, validation frames and validation IDs, how do degree-one and degree-two catalogue-inclusion models differ in posterior error and weight stability?

This is a controlled sensitivity comparison. It is not a new estimator, a real-data H0 analysis, an LVK/GWTC analysis, or evidence that either degree is universally superior.

## Frozen pilot

- Truth: quadratic-logit catalogue inclusion, with curvature gamma in {-1, 0, +1}.
- Model arms: catalogue-only (A), full-parent oracle (B), known-probability IPW (C), MAP degree-specific IPW (D), and coefficient-integrated degree-specific IPW (E).
- Paired streams: 3 replicates; within each replicate and gamma, both model degrees receive identical parent, GW-event, validation-frame and validation-sample inputs.
- Parent size: 1,200; detected events: 30.
- Validation frame: 500; nested validation sizes: 100 and 500.
- H0 truth: 70 km/s/Mpc; prior support [30,110]; 281-point grid.
- Quadrature order: 5 for this multi-replicate pilot. The separate order 5/7/9 single-realization pilot has already tested numerical convergence, but this lower-order pilot is not a substitute for convergence checks in a future larger experiment.

## Pairing and audit gates

The orchestrator runs the independently audited multi-replicate workflow separately for degree 1 and degree 2 using the exact same master seed. It then:

1. Requires both child workflows to report lifecycle PASS and their independent audits to pass.
2. Compares parent, catalogue, GW events, event truth, validation-frame truth and every nested validation sample exactly across model degrees.
3. Requires the H0 grid and A/B/C likelihood curves to be bitwise identical across degrees. These arms should not depend on the fitted inclusion-model degree.
4. Reconstructs paired D2-minus-D1 posterior-error contrasts at replicate/gamma/validation-size/arm grain.
5. Reports paired Monte Carlo standard errors and t-based 95% intervals, alongside bias, RMSE and 68% coverage for each degree.
6. Hashes every outer artifact file and records all configurations, child audits, and logs.

An incomplete lifecycle cannot be promoted to a successful paired result. Failed/non-evaluable replicate records remain in the nested child ledgers; a future confirmatory study must predefine how these affect inference rather than silently dropping them.

## Critical coefficient-prior caveat

The generator uses a centered quadratic logit of the form

    eta(x) = alpha + beta*x + gamma*(x^2 - c_g).

A degree-two fit with raw basis [1,x,x^2] can represent this function, but independent Gaussian priors on raw coefficients (theta0, theta1, theta2) are not automatically the same prior as independent Gaussian priors on generator coefficients (alpha,beta,gamma). The degree-one and degree-two fits can therefore differ both in functional flexibility and in induced prior. The present result must be interpreted as a comparison of the complete fitted procedures under their implemented priors, not a clean isolated test of polynomial degree. A centered-basis, prior-matched sensitivity experiment is required before attributing differences solely to model flexibility.

## Interpretation limits

- Three replicates are an engineering smoke, not a calibration study. Coverage estimates are extremely noisy at this size.
- The H0 likelihood is a controlled synthetic construction. No observational measurement is claimed.
- The oracle arm uses simulation truth for diagnosis only.
- Galaxy-catalogue inclusion and GW detection selection remain distinct.
- IPW and an additional missing-host mixture correction must not be applied together.
- Posterior median shifts alone do not establish improved calibration.
- Quadrature order 5 is not evidence of order convergence; order sensitivity must be rerun for any apparent scientific effect.

## Promotion criteria

Do not describe either degree as superior unless a larger, independently audited prior-predictive experiment demonstrates stable calibration and paired improvement across curvature, validation size, event count, coefficient-prior parameterization and quadrature/grid sensitivity. Keep the current PR in draft state pending all workflow and artifact checks.
