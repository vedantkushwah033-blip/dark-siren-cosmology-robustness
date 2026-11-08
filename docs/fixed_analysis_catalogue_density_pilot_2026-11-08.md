# Fixed-analysis-catalogue density-recovery pilot: corrected RNG run

**Run:** 42475919045  
**Artifact:** `fixed-analysis-catalogue-density-pilot-py3.12`  
**Artifact ID:** 13615247781  
**Commit:** `823773fb37555ccfb2c26c05c8ce62b9d4b586e5`  
**Disposition:** exploratory pilot; independent artifact audit PASS. This is not a confirmatory result.

## Design

- 5 independent simulation replicates.
- 4,000-object analysis parent frame per replicate/scenario.
- Separate independent validation frame of 1,000 objects.
- Nested validation sizes: 100, 300, 1,000.
- Curvature: -1.2, 0, +1.2.
- Target retention: 0.25.
- Logistic models: degree 1 and 2.
- Scott multiplier: 0.15.
- Five arms: catalogue-only, known-inclusion IPW diagnostic, parent-frame KDE diagnostic, estimated degree-1 IPW, estimated degree-2 IPW.
- The same analysis parent frame and retained analysis catalogue are used at all three validation sizes within each replicate and DGP scenario.

Independent random streams use `SeedSequence([seed, replicate, stream_id])`; analysis redshifts, validation redshifts, validation ordering, analysis retention uniforms and validation retention uniforms use distinct stream IDs.

## Execution and audit

| Check | Result |
|---|---:|
| Planned attempts | 225 |
| Successful attempts | 225 |
| Failed attempts | 0 |
| Non-evaluable attempts | 0 |
| Independently reconstructed density estimates | 225 |
| Frame/scenario partitions checked | 15 |
| Independent audit | PASS |

The raw validation-frame correlations were also inspected after extraction. Across the 15 replicate/curvature frame records, the maximum absolute correlation between validation redshift and retention uniform was 0.0165; the maximum absolute correlation between validation redshift and validation rank was 0.0406. These are diagnostic checks, not proofs of independence.

## Paired density-recovery contrasts

The table reports the replicate-paired mean difference
[
Delta=mathrm{IAE}_{mathrm{estimated IPW}}-mathrm{IAE}_{mathrm{catalogue-only}}.
]
Negative values favour estimated IPW on this density endpoint; positive values favour catalogue-only. Each cell has 5 paired replicates. The reported MC standard error is descriptive and should not be treated as a robust confidence interval with only five replicates.

| Curvature (gamma) | Validation n | Degree 1 mean (Delta) (MC SE) | Degree 2 mean (Delta) (MC SE) |
|---:|---:|---:|---:|
| -1.2 | 100 | -0.2322 (0.0261) | -0.1952 (0.0274) |
| -1.2 | 300 | -0.2776 (0.0069) | -0.2768 (0.0132) |
| -1.2 | 1,000 | -0.2946 (0.0077) | -0.2750 (0.0263) |
| 0.0 | 100 | -0.0292 (0.0422) | -0.0302 (0.0360) |
| 0.0 | 300 | -0.0725 (0.0221) | -0.0604 (0.0211) |
| 0.0 | 1,000 | -0.1143 (0.0078) | -0.0956 (0.0165) |
| +1.2 | 100 | +0.1240 (0.0396) | +0.0976 (0.0239) |
| +1.2 | 300 | +0.0424 (0.0159) | +0.0401 (0.0133) |
| +1.2 | 1,000 | -0.0027 (0.0033) | +0.0049 (0.0092) |

## Descriptive interpretation only

1. In the (gamma=-1.2) condition, estimated IPW has lower mean IAE than catalogue-only at all three validation sizes in this five-replicate pilot.
2. At (gamma=0), the estimated-IPW advantage is modest at n=100 and larger at n=1,000.
3. At (gamma=+1.2), estimated IPW performs worse than catalogue-only at n=100 and n=300; the paired mean difference is close to zero at n=1,000.
4. The oracle known-inclusion IPW arm has lower mean IAE than estimated IPW in all three curvature scenarios, consistent with estimation uncertainty being consequential.
5. The parent-frame KDE diagnostic has lower error than the retained-catalogue estimators, as expected for a benchmark with access to the complete analysis parent frame. It is a diagnostic ceiling, not a deployable estimator.

These patterns are consistent with a finite-validation information-size effect under this particular synthetic design. Five replicates are far too few for confirmatory inference, stable failure-rate estimation, multiplicity-adjusted claims, or choosing a preferred model. The bandwidth multiplier is fixed at 0.15 and is not established as universally optimal.

## What this does not establish

- No H0 likelihood or posterior is evaluated.
- No posterior calibration or coverage is evaluated.
- No GW detection-selection correction is evaluated.
- No real galaxy survey or LVK/GWTC data are used.
- No universal superiority of IPW, degree 1, or degree 2 is claimed.
- No evidence is provided about the operational disjoint finite-parent design; this FAC pilot isolates validation information by using an independent validation frame.

## Next gate

Scale the corrected FAC design to a predeclared, adequately powered replicate count. Freeze a primary contrast before the confirmatory run, report all secondary contrasts with multiplicity handling, and retain failure-aware summaries. Only after that density-recovery stage is stable should inclusion-coefficient uncertainty be propagated into the selected H0 likelihood.
