# Full grid-refinement audit: results and interpretation

## Run provenance

- Workflow: `Posterior grid-refinement audit`
- Run ID: `40088129105`
- Commit: `ebfb7ad22d1aabf90d5269c39dcbff17dc12efc4`
- Seed: `20261030`
- Replicates per design cell: 200
- Grid sizes: 1,001, 2,001, 4,001 points
- Artifact: `grid-refinement-audit-results` (contains `per_run.csv` and `summary.csv`)
- Artifact SHA-256: `df1e40bb421ab026fa49ef9f10aa8d3a0d6f46bb741f84097ab2995663092d31`
- Artifact retained until 2027-01-22 (per GitHub artifact metadata).

The run completed successfully and the workflow validated 21,600 posterior fits across 36 design cells, plus 36 cell-summary rows. This is a numerical audit of a synthetic two-host model, not observational inference.

## Reported grid discrepancies

All differences are absolute differences between summaries calculated on the same simulated observation at the indicated grid resolutions. Values are in the native units of the posterior summary (for (H_0)-like summaries, km s(^{-1}) Mpc(^{-1})).

| Diagnostic | Median of cell medians | Largest cell 95th percentile | Maximum observed |
|---|---:|---:|---:|
| Posterior median: 1,001 vs 4,001 | 0.00010945 | 0.00138315 | 0.00250182 |
| Lower 68% endpoint: 1,001 vs 4,001 | 0.00023138 | 0.00230339 | 0.00248676 |
| Upper 68% endpoint: 1,001 vs 4,001 | 0.00014267 | 0.00265184 | 0.00680109 |
| Lower 95% endpoint: 1,001 vs 4,001 | 0.00043612 | 0.00448120 | 0.00613489 |
| Upper 95% endpoint: 1,001 vs 4,001 | 0.00015498 | 0.00472141 | 0.13045331 |
| Posterior median: 2,001 vs 4,001 | 0.00002090 | 0.00027348 | 0.00057161 |

## Interpretation

Typical discrepancies and cell-level 95th-percentile discrepancies are small in this run. The posterior median differences decrease substantially from the 1,001-vs-4,001 comparison to the 2,001-vs-4,001 comparison, consistent with median estimates being numerically stable for this grid family.

However, the maximum discrepancy for the **upper 95% endpoint** is about 0.13045, much larger than the other reported maxima. This outlier must not be hidden by quoting only median differences. The existing log summary does not identify the responsible design cell or observation, so that row should be located in `per_run.csv` before claiming uniform convergence of all interval endpoints.

A likely mechanism to investigate is a discrete/rapidly changing posterior tail or quantile interpolation effect. This is a hypothesis, not a demonstrated explanation. Compare the affected observation's posterior density on all three grids, confirm normalization and quantile calculation, and test an even finer grid or a denser local grid around the tail.

## Important limitation in the current summary

The script reports each cell's median, 95th percentile, and maximum of the 1,001-vs-4,001 endpoint differences. The CI output therefore supports a descriptive convergence assessment, but it does not currently print the exact case responsible for the upper-95% maximum. The per-run artifact contains the individual endpoint values needed to recover it.

## Recommended prespecified acceptance criteria

Before using grid-based inference in the paper, choose a tolerance based on the precision needed for the scientific claim. A defensible operational rule would separately require:
1. The 95th percentile of the absolute discrepancy for every reported posterior summary to be below the prespecified tolerance in every design cell.
2. No unexplained maximum discrepancy above a separately stated guardrail.
3. The same conclusion under a targeted finer-grid check for any outlier.
4. Agreement of posterior normalization and quantile calculations with independent numerical checks.

Do not select the tolerance after inspecting which value makes this run pass. State the tolerance and rationale in the protocol, then evaluate against it.

## Scope

The simulation assumes a low-redshift two-host distance model, known redshifts and host weights, and a simplified noise model. Grid convergence establishes only that numerical discretization is controlled for this implementation and design. It does not establish astrophysical realism, correct selection modelling, calibrated coverage, or validity for real gravitational-wave observations.


## Outlier resolution and targeted follow-up

The archived per-run artifact was inspected and the exact case responsible for the largest upper-95% discrepancy was recovered:

- True synthetic H0: 52 km s^-1 Mpc^-1
- Second host redshift: 0.04 (first host redshift: 0.02)
- Noise fraction: 0.03
- Replicate index: 138 (zero-based in the simulation loop)
- Inference weighting: reversed (0.3, 0.7)
- Simulated observed distance: 119.47415214556634 Mpc
- Upper 95% endpoint at 1,001 points: 60.848373
- At 2,001 points: 60.742249
- At 4,001 points: 60.717920

The posterior median is already stable for this case (51.713915, 51.713591, and 51.713486 across the three grids), as are the 68% endpoints. The tail quantile is more sensitive: the 1,001-vs-4,001 difference is 0.130453, while the 2,001-vs-4,001 difference is 0.024330.

A targeted finer-grid check reconstructed this observation and recalculated the upper 95% endpoint:

| Grid points | Upper 95% endpoint |
|---:|---:|
| 1,001 | 60.848373 |
| 2,001 | 60.742249 |
| 4,001 | 60.717920 |
| 8,001 | 60.711905 |
| 16,001 | 60.710428 |
| 32,001 | 60.710055 |
| 64,001 | 60.709962 |

The sequence approaches a stable value near 60.710. The finer-grid changes shrink monotonically in this targeted check, supporting numerical discretization/quantile interpolation as the explanation for this particular discrepancy; this does not prove monotonic convergence in general. No defect in the likelihood or posterior normalization was demonstrated. The practical conclusion is that 1,001 points are not uniformly adequate for the upper 95% endpoint, and even 4,001 points leave a visible difference for this case.

The audit script now prints the identity and three-grid endpoint values of the largest upper-95% discrepancy. A targeted regression test preserves this observation and checks refinement through 64,001 points. The targeted diagnostic also records posterior area on the prior support; this area is recomputed with the same trapezoidal rule used by the inference function, so it is a consistency check rather than independent-quadrature validation.


## Follow-up: artifact-level checks of prior support and targeted convergence

The subsequent automated numerical-diagnostics workflow completed successfully on Python 3.10, 3.11, and 3.12 (run 40089986883). Its posterior_grid_stability.csv and targeted_grid_convergence.csv artifacts were downloaded and inspected directly.

### Normalization and numerical support

For the baseline case, the reported posterior area is 1.0. With the same prior interval [50, 90], extending the numerical grid from [50, 90] to [30, 110] changes the posterior median by about 1.4e-14 and leaves the reported credible-interval endpoints equal to numerical precision. This is the expected control: extending the computational grid without changing the prior should not change the posterior.

The diagnostic previously integrated over the zero-padded exterior grid as well as the prior support. It has now been corrected to integrate only over the prior support, matching the inference routine. Regression tests were added for this consistency check. The corrected diagnostic workflow for commit 123a99d6 passed on all three Python versions.

### Prior sensitivity is not grid sensitivity

A separate comparison widens the actual prior from [50, 90] to [30, 110], while retaining the wider numerical grid. The posterior changes substantially:

| Summary | Prior [50, 90] | Prior [30, 110] | Absolute change |
|---|---:|---:|---:|
| Median | 70.415477 | 68.200082 | 2.215394 |
| Lower 68% endpoint | 59.723136 | 40.935988 | 18.787148 |
| Upper 68% endpoint | 82.573676 | 89.985520 | 7.411844 |
| Lower 95% endpoint | 51.869120 | 31.720089 | 20.149031 |
| Upper 95% endpoint | 88.756256 | 105.078786 | 16.322531 |

These are outputs of the specified simplified synthetic model. They show that the inferred posterior is materially sensitive to the chosen prior interval in this diagnostic case. This must not be conflated with grid-discretization error: changing grid extent while holding the prior fixed is a numerical control; changing prior bounds changes the statistical model.

The baseline case also has nonzero posterior mass near its prior boundaries (approximately 1.08% in the lower 2%-wide prior band and 1.60% in the upper band). These values are worth reporting as a diagnostic, but by themselves do not establish that the prior is inappropriate.

### Targeted upper-tail convergence artifact

The saved seven-resolution output confirms the upper 95% endpoint for the previously identified case:

| Grid points | Upper 95% endpoint | Absolute change from previous grid |
|---:|---:|---:|
| 1,001 | 60.848373 | — |
| 2,001 | 60.742249 | 0.106124 |
| 4,001 | 60.717920 | 0.024330 |
| 8,001 | 60.711905 | 0.006015 |
| 16,001 | 60.710428 | 0.001477 |
| 32,001 | 60.710055 | 0.000373 |
| 64,001 | 60.709962 | 0.000093 |

This confirms that the targeted upper-tail sequence stabilizes as the grid is refined. It is a single-case result and cannot substitute for a tolerance-based assessment across all 36 design cells.

### Verification status

- The tests workflow for commit 123a99d6 completed successfully on Python 3.10, 3.11, and 3.12, including the corrected posterior-grid diagnostic and artifact upload.
- The full test suite for that same commit was still in progress at the time of the status check.
- The later full-suite run for commit 4c9ae127 had Python 3.10 completed successfully, while Python 3.11 and 3.12 were still in progress at the last status check.

Therefore the corrected diagnostic workflow is verified, but do not describe the latest full-suite run as completely finished until all matrix jobs have completed.
