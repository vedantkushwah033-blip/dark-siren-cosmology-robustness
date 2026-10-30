# Estimated-inclusion pilot smoke results — 30 October 2026

**Status:** executed smoke test; non-confirmatory.  
**Source commit:** `829a27c6992624ade83aa85f37ff131a80f10d53`  
**GitHub Actions run:** `41001160389` (completed successfully on Python 3.10, 3.11, and 3.12).  
**Artifact:** `numerical-diagnostics-py3.12`, SHA-256 `ebadb5a2ed4f9168a9befd039fa2eab8aec2f73388d1cf1df0b34b89b7df09ad`.

## Configuration

- Seed: 20261030
- Parent population: 1,200 synthetic galaxies per replicate
- GW events: 8 per replicate
- Replicates: 2
- Validation sample sizes: 100 and 300, nested within each replicate
- Injected (H_0): 70 km s(^{-1}) Mpc(^{-1})
- Distance noise: 140 Mpc
- Detection threshold: 250 Mpc
- (H_0) grid: 41 points from 50 to 90
- Redshift grid: 201 points from 0.01 to 0.40
- True inclusion law: correctly specified linear logistic model with intercept 0.25 and slope -1.2

This was a pipeline smoke test, not a confirmatory simulation campaign.

## Results

| Arm | Successful / attempted | Mean bias | RMSE | Mean density ISE | Mean inclusion-probability MAE | Mean ESS fraction | Mean max inverse weight | Mean 68% interval width |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Catalogue-only | 2 / 2 | -3.666 | 3.669 | 0.16462 | — | 1.000 | 1.000 | 15.616 |
| Oracle | 2 / 2 | -0.727 | 0.728 | 0.07075 | 0.00000 | 0.920 | 3.582 | 14.917 |
| Estimated, (m=100) | 2 / 2 | -0.927 | 1.144 | 0.07779 | 0.02830 | 0.926 | 3.435 | 15.034 |
| Estimated, (m=300) | 2 / 2 | -0.862 | 0.900 | 0.06936 | 0.01976 | 0.927 | 3.433 | 14.985 |

The mean paired reduction in absolute posterior-mean error relative to catalogue-only was 2.739 for (m=100), 2.804 for (m=300), and 2.939 for oracle weighting. All contrasts were based on only two paired replicates.

Both estimated-inclusion validation sizes had lower probability MAE at (m=300) than at (m=100) in this tiny run. This is consistent with the expected benefit of a larger validation sample, but two replicates cannot establish that pattern reliably. Individual replicates were variable: for example, the (m=100) arm had absolute errors of 1.597 and 0.257 across the two seeds.

## Independent audit and failure accounting

The persisted artifact auditor reported:

- Overall status: **PASS**
- Replicate-arm rows expected / observed: 8 / 8
- Successful / failed rows: 8 / 0
- Boundary-limited successful fits: 0
- Paired successful replicates: 2 for each comparison arm
- Independent validation-logistic refit, weighted-KDE reconstruction, direct-quadrature posterior reconstruction, density-ISE reconstruction, and summary reconstruction all included in the audit scope.

The audit PASS establishes consistency of this smoke artifact with the implemented reconstruction checks. It does not establish that the assumed selection model is realistic or that the estimator is calibrated on real survey data.

## Interpretation

The smoke run demonstrates that the end-to-end synthetic pipeline executed reproducibly for this configuration and that its persisted results passed the independent artifact audit. In this particular correctly specified selection setting, estimated weighting moved the posterior-mean error and host-density ISE closer to the oracle benchmark than catalogue-only weighting.

**No inferential claim should be made from these two replicates.** The 2/2 coverage result has a very wide Wilson 95% interval (approximately 0.342 to 1.000). The normal intervals reported for paired contrasts are also extremely unstable at (n=2). The result does not establish nominal coverage, a general improvement from inverse-probability weighting, an observational (H_0) measurement, or performance on LVK/real galaxy catalogues.

## Next experiment

Use the prespecified 20-replicate smoke gate across the linear reference and nonlinear curvature conditions, inspect positivity and weight tails before scaling, and only then run the planned 200 replicates per cell. Add validation-sample bootstrap sensitivity as a separate uncertainty layer; do not conflate bootstrap variation of posterior summaries with posterior credible intervals.
