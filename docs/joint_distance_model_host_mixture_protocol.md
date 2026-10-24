# Joint distance-model × host-mixture robustness protocol

**Status: protocol only; no implementation or results are claimed.** This is the next proposed experiment after the independently audited single-host distance-model diagnostic and the completed host-mixture/threshold factorial.

## Research question

Does the impact of using the low-redshift distance approximation instead of flat Lambda-CDM depend on host-redshift separation and intrinsic host weights in the simplified dark-siren likelihood? This interaction experiment distinguishes distance-relation mismatch from the additional error caused by collapsing a host mixture to one effective redshift.

This remains a synthetic Gaussian-distance diagnostic. It is not a real gravitational-wave analysis, a calibrated galaxy-catalogue model, or a claim about observational H0.

## Primary design

- Fixed true H0: 60, 70, and 80 km s^-1 Mpc^-1.
- Two-host redshift pairs: (0.02, 0.04), (0.02, 0.10), and (0.10, 0.40).
- Intrinsic weights: (0.5, 0.5), (0.7, 0.3), and (0.9, 0.1).
- 300 accepted simulated observations per cell.
- 3 × 3 × 3 = 27 primary cells; 8,100 accepted observations.
- Fixed Omega_m = 0.3 for all flat-Lambda-CDM calculations.
- Gaussian distance noise: sigma = max(20 Mpc, 0.05 times the true mixture-mean distance). Apply this rule identically and record it for every replicate.
- Primary inference prior H0 in [40, 100] km s^-1 Mpc^-1; use at least a 4,001-point grid and a prespecified 2,001-versus-4,001 grid check on selected cells.
- Seed: 20261025, with deterministic independent streams derived from a stable hash of the global seed and cell configuration.

Generate each host's distance using flat Lambda-CDM at the fixed true H0, draw a host according to intrinsic weights, then add the stated Gaussian measurement noise. The primary phase has no detection truncation or selection correction; omitting selection deliberately isolates the distance-model/host-mixture interaction. A later phase may add detection truncation after this baseline is audited.

## Paired inference arms

Analyze each identical observation with four model configurations:

1. Matched LCDM, two-host mixture: flat Lambda-CDM distance relation, both host redshifts, intrinsic host weights.
2. Low-z, two-host mixture: low-redshift relation, both host redshifts, same weights. This isolates distance-relation mismatch while retaining mixture structure.
3. Matched LCDM, collapsed host: flat Lambda-CDM at z_eff = sum_i w_i z_i. This isolates host-mixture compression under the matched distance relation.
4. Low-z, collapsed host: low-redshift relation at the same z_eff. This combines both misspecifications.

No row may be generated independently for the four arms; they must reuse the same observed distance, truth, cell, and replicate identifier. Do not treat intrinsic weights as detected-sample weights because this baseline has no detection conditioning.

## Primary estimands

Primary A: paired change in absolute H0 error between low-z and matched-LCDM inference, holding host representation fixed (two-host and collapsed analyzed separately).

Primary B: paired change in absolute H0 error between collapsed and two-host inference, holding distance relation fixed (matched-LCDM and low-z analyzed separately).

Primary C, interaction contrast: (MAE_low-z,collapsed - MAE_LCDM,collapsed) - (MAE_low-z,two-host - MAE_LCDM,two-host). Positive values mean distance-model mismatch increases MAE more under the collapsed approximation than under the two-host representation. Report paired bootstrap 95% intervals with at least 5,000 draws, a stable documented bootstrap seed, and a percentile convention.

For all four arms and every cell, report bias, MAE, RMSE, 68% and 95% interval coverage with Wilson intervals, interval widths, valid-run counts, and posterior medians near the prior limits. Treat the 27-cell equal-weight average as a design average, not a population-level estimate.

## Required validation

1. Exactly 8,100 raw observations and 27 unique cells, each with 300 replicates.
2. Exactly four inference arms for every observation key; same observed distance and true H0 across arms.
3. All raw numeric fields finite; posterior medians and interval endpoints lie inside the declared prior; intervals are ordered.
4. Recompute all summary metrics and paired contrasts independently from raw rows.
5. Independently recompute bootstrap intervals from paired replicate-level contrasts and the recorded seed.
6. Verify flat-LCDM distances against an independent high-accuracy numerical reference at selected redshifts and verify the low-z limit.
7. Compare posterior summaries across the prespecified grid resolutions; investigate any endpoint discrepancy above 0.05 km s^-1 Mpc^-1 before interpretation.
8. Run focused tests and the full repository test suite on supported Python versions. Archive raw rows, summaries, manifest, commit SHA, dependency versions, workflow URL, and artifact digest.
9. Do not publish a conclusion as a result until the workflow has succeeded and the raw artifact has been independently audited.

## Interpretation guardrails

- The experiment uses a single observed noisy distance, a two-host synthetic mixture, fixed Omega_m, and a finite H0 grid.
- It does not model galaxy-catalogue incompleteness, a detector selection function, redshift uncertainty, correlated host data, or a full hierarchical population likelihood.
- Fixed truths are intentionally chosen design points, not draws from a realistic cosmological population.
- Coverage is an empirical repeated-sampling property under the specified simulator; posterior calibration in this toy model does not establish real-world calibration.
- The interaction contrast is descriptive for this design; no astrophysical significance should be inferred from a bootstrap interval alone.
- The model with the smallest MAE in this finite design is not automatically the best model for real data.

## Planned sequence

1. Implement the four paired inference arms and focused tests.
2. Run a small smoke test with 2-5 replicates per cell to catch schema and numerical errors.
3. Run the full 300-replicate protocol.
4. Inspect raw artifacts, independently recalculate every summary and contrast, and apply grid checks.
5. Publish the results report and update the README only after the audit passes.
6. Only then design a second-stage experiment adding detection truncation and selection normalization, followed by a separate catalogue-completeness extension.