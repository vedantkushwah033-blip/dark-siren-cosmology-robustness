# Expanded absolute-bandwidth density pilot: independently audited readout

Date: 2026-11-15  
Workflow: `43617052462`  
Artifact commit label: `d94fdd62e514c81c20444196e166d2ebbb341a40`  
Status: exploratory synthetic density-stage pilot; not confirmatory and not bandwidth selection.

## Design and audit

The same 20 paired replicate seeds were evaluated at 11 common physical Gaussian-kernel standard deviations:
`0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080`.

Each floor cell used N=4,000 parent galaxies, a disjoint validation frame of 1,500, nested validation sizes 100/300/1,000, a 401-point redshift grid, and quadrature order 5. The three catalogue floors were 0, 0.10, and 0.25.

All three workflow jobs completed successfully. The downloaded artifacts' outer audit reports each say PASS:

| Catalogue floor | Successful / planned | Child raw audits | Manifest files verified | Metric rows reconstructed | Scott-width rows reconstructed |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 20/20 | 20 | 244 | 20,460 | 1,740 |
| 0.10 | 20/20 | 20 | 244 | 20,460 | 1,740 |
| 0.25 | 20/20 | 20 | 244 | 20,460 | 1,740 |

The earlier 19/20 floor-0.25 result is superseded for this expanded-grid run. The audit PASS verifies the artifact manifests, replicate lifecycle/seed ledger, raw-input reconstructions, and aggregate summaries for this run. It does not validate the scientific model against observations.

## Primary Arm E and five-arm IAE

Arm E's primary metric is now the loss of its posterior-mean density, `E_posterior_mean_density`. The posterior-weighted expected node-specific loss remains a distinct secondary metric; it is not substituted for the primary estimand.

Across the three floors, three validation sizes, five primary arms (A–E), and 20 paired replicates per cell, the descriptive equal-cell mean IAE over candidate widths was:

| Physical width | Equal-cell mean IAE |
|---:|---:|
| 0.006 | 0.273621 |
| 0.008 | 0.269198 |
| 0.010 | 0.266429 |
| 0.012 | 0.264526 |
| 0.016 | 0.262379 |
| 0.024 | 0.265097 |
| 0.032 | 0.278531 |
| 0.040 | 0.295289 |
| 0.050 | 0.314031 |
| 0.060 | 0.328680 |
| 0.080 | 0.341565 |

The lowest observed pooled value is at 0.016, but this is a same-pilot descriptive minimum, not a selected bandwidth. Arm-specific minima differ: catalogue-only A is still improving at the upper candidate edge for floors 0 and 0.10, whereas the oracle/IPW arms generally minimize around 0.016. This is evidence that the common-width trade-off is arm- and scenario-dependent, and that the current grid does not bracket every arm's minimum.

At width 0.016, pooled across the three validation sizes (60 paired replicate-size records per floor), selected paired contrasts in IAE were:

| Floor | D − A mean (MCSE) | D − C mean (MCSE) | E-primary − D mean (MCSE) |
|---:|---:|---:|---:|
| 0.00 | -0.304288 (0.011362) | +0.022297 (0.012944) | +0.013419 (0.002049) |
| 0.10 | -0.056417 (0.017678) | +0.135143 (0.017454) | +0.015930 (0.001707) |
| 0.25 | +0.063712 (0.012321) | +0.096548 (0.012626) | +0.011675 (0.001366) |

Negative D−A means the estimated-MAP IPW arm has lower IAE than catalogue-only; positive D−C means it has higher IAE than known-probability IPW. These are exploratory paired summaries, not adjusted confirmatory tests. The positive E−D values show that, in this pilot, averaging densities over coefficient uncertainty did not improve this primary IAE relative to the MAP plug-in density. The E posterior-mean loss is nevertheless correctly distinguished from posterior expected node loss.

## Interpretation

1. In the correctly specified floor-0 case, estimated inclusion correction improves density recovery substantially over catalogue-only, but does not reach the known-probability IPW/oracle benchmark at the selected descriptive width.
2. At floor 0.10, the estimated correction gives a smaller gain over catalogue-only and a large gap remains relative to known-probability IPW.
3. At floor 0.25, the MAP estimated-IPW arm is worse than catalogue-only at width 0.016, while known-probability IPW remains close to the parent oracle. This supports a controlled misspecification warning, not a general claim about all misspecified models.
4. The primary Arm E posterior-mean density is consistently slightly worse in IAE than D in this pilot. This is an empirical outcome, not a theorem that uncertainty integration is harmful.
5. A pooled width minimum must not be frozen from this same pilot. The A-arm boundary minima and different floor-specific optima require a predeclared objective, a grid that brackets relevant minima, and a disjoint training/holdout evaluation.

## Provenance caveat requiring cleanup

The artifact `config.json` contains a limitation string saying stored inclusion-posterior nodes are consumed from the production fitter and are not independently refitted. The audited child-auditor source on the active branch contains a separate reconstruction of the degree-one logistic posterior, Laplace proposal, quadrature nodes, and normalized importance weights from validation records. The runner source has now been corrected to state that the child auditor independently refits the degree-one logistic posterior and reconstructs its Laplace proposal/quadrature weights; a regression assertion protects this metadata. The already-produced artifact remains unchanged and retains the stale limitation string. The source correction is not yet validated by a new exact-head workflow, so a rerun is required before treating the new metadata as CI-verified. The PASS result above is reported as emitted by the original workflow; this note does not silently rewrite the artifact metadata.

## Explicit limitations and next gate

- Synthetic low-redshift Euclidean density recovery only.
- No LVK/GWTC/public gravitational-wave data, real galaxy-survey validation, observational H0 result, or H0 posterior is produced here.
- No bandwidth is frozen. No training/holdout split has yet been run.
- Density recovery does not establish posterior calibration, PIT, coverage, or estimator superiority.
- Before H0 propagation: reconcile the artifact metadata, freeze a training/holdout bandwidth protocol and objective, expand the grid where boundary minima remain, execute the disjoint holdout, and pass its independent audit.

