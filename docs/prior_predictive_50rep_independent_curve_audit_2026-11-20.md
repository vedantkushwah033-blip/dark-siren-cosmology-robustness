# Independent posterior-rank and quantile reconstruction: 50-replicate SBC pilot

**Source workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44683807350  
**Source artifact:** 14613021797, five-arm-prior-predictive-50rep  
**Artifact SHA-256:** 1d8667e2f7391283237e27b0c6db26cac6e9dfd56d344649b6931ec047745c0c  
**Artifact size:** 168,798,378 bytes  
**Source commit:** 6cf444d25b0c422056c09082f8ea0e2120110a45

## Purpose

The existing independent orchestration audit verifies lifecycle completeness, child audit attestations, paired inputs, aggregate arithmetic, and artifact hashes. This additional audit independently reconstructs posterior ranks and quantiles directly from every saved H0 log-likelihood curve. It does not import the production rank or quantile helpers.

## Independent reconstruction

For each posterior curve, the audit reads the saved 561-point H0 grid and arm-specific log-likelihood, forms a stabilized density proportional to exp(log-likelihood), integrates over dH0 using the trapezoidal rule, evaluates the normalized CDF at the generating H0 truth stored in the child configuration, and inverts the CDF to reconstruct q2.5, q5, q10, q16, q25, q75, q84, q90, q95, and q97.5. It then compares all reconstructed values with the aggregate posterior-summary rows and rechecks all 2,858 outer-manifest entries and all 150 child audit PASS records.

## Result

**Independent curve audit: PASS**

| Check | Result |
|---|---:|
| Planned posterior curves | 2,250 |
| Posterior ranks independently reconstructed | 2,250 |
| Child audits rechecked | 150 / 150 PASS |
| Outer manifest entries verified | 2,858 / 2,858 |
| Maximum absolute posterior-rank difference | 2.0e-15 |
| Maximum absolute quantile difference | 1.6e-13 |

The maximum discrepancies are at floating-point roundoff scale. The audit also verifies each summary truth against its child configuration rather than trusting the truth column in the summary table.

## Descriptive rank pattern at validation size 500

Each gamma/arm stratum contains 50 independent replicate truths. Mean posterior ranks were:

| Curvature gamma | Full-parent oracle B | Known-pi IPW C | Degree-one MAP D | Degree-one integrated E |
|---:|---:|---:|---:|---:|
| -1 | 0.481 | 0.525 | 0.710 | 0.700 |
| 0 | 0.481 | 0.482 | 0.507 | 0.501 |
| +1 | 0.481 | 0.465 | 0.303 | 0.306 |

The catalogue-only arm's mean rank is 0.997, 0.993, and 0.983 at gamma=-1, 0, and +1, respectively. This is severe rank miscalibration in this synthetic design. The degree-one estimated-inclusion arms are close to 0.5 at zero curvature but show substantial curvature-dependent departures at gamma=+/-1. The full-parent oracle is the matched simulation reference, not a deployable estimator.

Ten-bin chi-square rank-histogram checks at 50 replicates per stratum are descriptive only. Some strata show strong departures, but no multiplicity-adjusted confirmatory claim is made. Gamma strata are paired through common parent/event realizations, and validation sizes are repeated within each replicate; these dependencies must be respected in pooled uncertainty analyses.

## What this audit does not establish

- It reconstructs ranks and quantiles on the saved 561-point grid. It does not recompute the underlying likelihood at a finer H0 resolution.
- It does not certify exact H0-grid convergence. The earlier 561-to-2,241 interpolation check remains interpolation-only.
- It does not validate the method on another DGP family or observational data.
- Fifty replicate truths per curvature/arm stratum remain a pilot-sized calibration sample.

**Interpretation:** the saved posterior ranks and central quantiles are reproducible from the saved likelihood curves to numerical precision. The calibration departures in the pilot are not caused by a mismatch between saved rank/quantile summaries and those curves. Exact H0-grid convergence remains an open numerical gate.

## Reproduction

    python scripts/audit_prior_predictive_posterior_curves.py path/to/five_arm_prior_predictive_50rep --output-prefix artifacts/prior_predictive_curve_audit
    pytest tests/test_audit_prior_predictive_posterior_curves.py -q
