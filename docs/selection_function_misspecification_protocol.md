# Protocol: selection-function misspecification in synthetic host-density inference

**Status:** design protocol; results must not be claimed until the implementation has been run and audited.

## Research question

How do catalogue-only estimation, correctly specified inverse-inclusion weighting, and misspecified inverse-inclusion weighting differ when a redshift-dependent galaxy selection process distorts the observed catalogue used in a simplified dark-siren H0 analysis?

## Data-generating process

For each replicate, generate one parent population with redshifts z_i from a documented parent density. Generate catalogue inclusion indicators C_i ~ Bernoulli(p_true(z_i)). The observed catalogue contains only objects with C_i = 1. Independently follow the existing event-generation mechanism for GW-like distances, preserving its limitations and its random-seed lineage.

The *true* selection function must be used only to generate inclusion indicators and to score the synthetic experiment. A non-oracle analyst receives either no selection model, the correct model, or an intentionally incorrect model. Do not allow the misspecified arm to access p_true.

## Arms

1. **Oracle reference:** estimate host density from the parent population. This is an unattainable reference in real incomplete-catalogue settings.
2. **Catalogue-only:** estimate from observed catalogue redshifts without inclusion correction.
3. **Correctly specified weighting:** use inverse p_true(z) weights, subject to positivity.
4. **Misspecified weighting:** use inverse p_assumed(z), where p_assumed differs from p_true by a declared parameterized perturbation.
5. **Optional stabilized weighting:** include only if the stabilization rule is specified before inspecting outcomes; report it as a bias-variance trade-off, not a free improvement.

The existing code may encode selection functions in a different form. Inspect it before implementing this protocol and preserve the exact existing scientific contract unless a change is justified and documented.

## Misspecification scenarios

Begin with a small set, not a large Cartesian grid:

- Correct model: p_assumed(z) = p_true(z).
- Shifted transition: the transition location is offset by a predeclared amount.
- Slope mismatch: the assumed redshift dependence is too shallow or too steep.
- Severe but positive mismatch: a stress test that still satisfies the positivity conditions required for inverse weighting.

Keep all probabilities within (0, 1]. Explicitly report regions where inclusion probability is too small for stable weighting. Do not silently clip probabilities; if clipping/stabilization is tested, treat it as a separate method with a declared rule.

## Primary outcomes

For each replicate and arm, preserve raw values for:

- host-density integrated squared error on a common grid, using a documented integration rule;
- H0 posterior mean and bias;
- H0 RMSE across replicates;
- empirical coverage of the declared 68% posterior interval;
- mean interval width;
- fit failure rate;
- weighted effective sample size, ESS = (sum_i w_i)^2 / sum_i w_i^2;
- maximum normalized weight and a high weight quantile;
- a clear record of any failed or non-finite fit.

Report Monte Carlo uncertainty. For coverage, report the covered count and denominator and a binomial interval. For between-arm contrasts, use paired replicate-level differences when the shared simulation design justifies pairing. Distinguish Monte Carlo uncertainty from posterior uncertainty within a replicate.

## Replicate and seed design

Use deterministic, recorded seeds. Within a replicate, share the parent population and underlying event realization across arms wherever this is scientifically valid, so the contrast isolates the catalogue-density method. Use independent replicate seeds. Preserve a replicate identifier and full configuration for every row. Do not discard failed fits without reporting them.

## Required validity checks

- Verify observed catalogue membership is generated from p_true, not p_assumed.
- Verify no non-oracle arm reads parent-population truth.
- Verify probabilities and weights are finite and in range.
- Test weight-scale invariance where applicable.
- Test expected behavior under the correctly specified model.
- Verify ESS bounds and behavior for equal weights.
- Verify density normalization and boundary behavior.
- Audit posterior normalization from saved grid values, not only a recorded scalar.
- Reconstruct posterior mean, quantiles, interval width, coverage, and aggregate metrics independently.
- Include adversarial tests for altered posterior arrays, altered aggregate summaries, missing arms, duplicate replicate-arm rows, and mismatched seeds.

## Pilot and scale-up

1. Run a deterministic unit-test-sized case.
2. Run 30–50 paired replicates for a smoke/pilot study if runtime permits.
3. Estimate runtime and Monte Carlo precision from the pilot.
4. Choose the main replicate count based on precision requirements and compute budget; document the rationale.
5. Run bandwidth factors 0.5, 0.8, and 1.2 only after verifying how the code interprets the bandwidth factor.
6. Assess H0-grid and redshift-grid sensitivity.
7. Preserve raw replicate outputs and a SHA-256 manifest.

A 30–50 replicate pilot is not by itself sufficient evidence for precise coverage claims.

## Interpretation guardrails

Do not presume that correct weighting always wins or that misspecification always harms H0. The experiment is intended to measure the behavior under its stated synthetic assumptions. A null or unfavorable result is valid.

Do not describe this simplified likelihood as a realistic LVK analysis. Do not claim a new H0 measurement or resolution of the Hubble tension. Generalize only to tested scenarios.

## Publication gate

Do not draft a numerical results section until runs are complete and independently audited. Before external review, ensure that the methods, random seeds, software versions, configurations, raw results, failure counts, audit report, figures, and limitations can all be reproduced from a clean checkout.
