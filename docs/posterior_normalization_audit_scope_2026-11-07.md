# Posterior normalization audit: claim boundary

## Scope

The host-density benchmark currently stores posterior summary statistics and a scalar `posterior_integral` in each replicate/arm record. The independent audit checks that this recorded scalar is finite and approximately one. It does not currently persist or reconstruct the complete posterior density vector.

## Correct interpretation

The audit therefore verifies **recorded posterior-integral claims**, not posterior normalization independently. The output key is named `recorded_posterior_integral_claims_checked` to avoid overstating what the artifact supports.

This is an audit-language correction. It does not establish that existing posterior grids are incorrect, and it does not independently validate the posterior mean, median, quantiles, or density shape.

## Work required for full independent verification

A later schema revision should persist, for each replicate and arm, either:

- the declared H0 grid and corresponding finite, non-negative posterior density vector; or
- a deterministic, versioned set of inputs sufficient for an independent auditor to reconstruct that vector.

The independent auditor should then recompute the integral, posterior mean, quantiles, interval width, and coverage indicator from the stored/reconstructed density. Tests should tamper with the density while leaving summary fields unchanged and verify that the audit fails. The artifact schema version and hashes should be recorded in the manifest.

## Scientific scope

The benchmark remains a synthetic exploratory toy study. It is not a real H0 measurement, does not use LVK event data, and does not validate a real detector or event-selection pipeline.
