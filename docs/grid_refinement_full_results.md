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
