# 100-replicate fixed-truth H₀ calibration readout

**Status:** completed exploratory fixed-truth repeated-sampling experiment; orchestration audit PASS. This is not joint simulation-based calibration (SBC), not an observational H₀ measurement, and not a confirmatory claim of estimator superiority.

## Run and artifact provenance

- Workflow: [44411667480](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44411667480)
- Runner commit: `652021943b31490302d586b5b9c77c7ccdfcc8f7`
- Artifact: `14485703115`, 475,742,896 bytes.
- ZIP SHA-256: `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`.
- Design: 100 paired replicate seeds × γ ∈ {−1, 0, +1} × validation n ∈ {100, 500}; 1,200 parent galaxies, four detected events, H₀ truth 70 km s⁻¹ Mpc⁻¹, H₀ grid 281 points, quadrature ladder 3/5/7/9.
- Status ledger: 600 planned cells; 599 SUCCESS, 1 NON_EVALUABLE, 0 FAILED, 0 NOT_RUN.
- Independent orchestration audit: PASS; 600 child audits verified; 300 nested validation blocks and 200 cross-curvature validation-pairing blocks checked; 14,106 aggregate-manifest entries verified; summary metrics recomputed.

The one NON_EVALUABLE cell is γ=−1, n=100. It remains in the planned denominator and is excluded from fitted-arm performance metrics for that cell. The known-π and full-parent oracle decomposition is being performed separately on the same raw realizations.

## Fixed-truth repeated-sampling results

Bias and RMSE are in km s⁻¹ Mpc⁻¹. Coverage denominators are successful/evaluable replicates.

| γ | Validation n | Success / planned | Mean median bias | RMSE | 68% coverage | 90% coverage | 95% coverage | Median minimum Kish ESS | Median maximum normalized weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −1 | 100 | 99 / 100 | +2.020 | 17.455 | 0.545 | 0.727 | 0.798 | 19.70 | 0.1552 |
| −1 | 500 | 100 / 100 | +2.255 | 17.215 | 0.510 | 0.740 | 0.780 | 60.89 | 0.0776 |
| 0 | 100 | 100 / 100 | +6.428 | 19.930 | 0.510 | 0.790 | 0.860 | 13.91 | 0.1976 |
| 0 | 500 | 100 / 100 | +6.768 | 19.843 | 0.510 | 0.770 | 0.850 | 42.41 | 0.0954 |
| +1 | 100 | 100 / 100 | +10.396 | 21.587 | 0.560 | 0.820 | 0.880 | 13.30 | 0.1623 |
| +1 | 500 | 100 / 100 | +11.114 | 22.071 | 0.540 | 0.720 | 0.880 | 34.53 | 0.0912 |

Wilson 95% intervals for the coverage estimates:

| γ | n | 68% coverage interval | 90% coverage interval | 95% coverage interval |
|---:|---:|---:|---:|---:|
| −1 | 100 | [0.448, 0.640] | [0.632, 0.805] | [0.708, 0.865] |
| −1 | 500 | [0.413, 0.606] | [0.646, 0.816] | [0.689, 0.850] |
| 0 | 100 | [0.413, 0.606] | [0.700, 0.858] | [0.779, 0.915] |
| 0 | 500 | [0.413, 0.606] | [0.678, 0.842] | [0.767, 0.907] |
| +1 | 100 | [0.462, 0.653] | [0.733, 0.883] | [0.802, 0.930] |
| +1 | 500 | [0.443, 0.634] | [0.625, 0.799] | [0.802, 0.930] |

## Paired contrasts

The γ=+1 minus γ=−1 mean median-bias contrast is:

- n=100: +8.273, MCSE 0.512, 99 paired evaluable replicates.
- n=500: +8.859, MCSE 0.466, 100 paired evaluable replicates.

The n=500 minus n=100 mean median-bias contrast is:

| γ | Paired replicates | Mean difference | MCSE |
|---:|---:|---:|---:|
| −1 | 99 | +0.152 | 0.229 |
| 0 | 100 | +0.340 | 0.266 |
| +1 | 100 | +0.718 | 0.258 |

The larger validation sample improved weight stability: median minimum Kish ESS increased and median maximum normalized weight decreased in every curvature cell. That improvement did not translate into a material reduction in H₀ median bias in this run.

## Interpretation

1. **Curvature dependence is pronounced in the fitted degree-one arm.** Mean median bias rises from about +2.0 to +2.3 at γ=−1, to +6.4 to +6.8 at γ=0, and +10.4 to +11.1 at γ=+1.
2. **Fixed-truth interval coverage is poor.** The nominal 68% intervals cover only about 51–56% of successful realizations. The nominal 95% intervals cover about 78–88%. These should not be described as calibrated intervals.
3. **Increasing validation n stabilizes IPW diagnostics but does not cure the H₀ bias.** The weight ESS improves substantially, while the paired bias changes are small and slightly positive.
4. **The γ=0 baseline is already biased.** The separate full-parent fixed-count baseline experiment found positive posterior-median bias even with the complete finite parent mixture and exact host weights. Therefore this 600-cell experiment alone cannot attribute the baseline bias to estimated catalogue inclusion.
5. **This is fixed-truth coverage, not SBC.** It does not test prior-averaged rank uniformity or joint calibration of all latent parameters. Four detected events per realization and the synthetic low-z setup sharply limit generalization.

## Next required analysis

The separate post-hoc oracle audit should reconstruct, on the same 100 replicate/event realizations:

- full-parent oracle posterior;
- known-π IPW posterior using the exact DGP inclusion probabilities;
- fitted degree-one IPW posterior from the pilot.

The contrasts separate event-sample/selected-likelihood baseline behavior, finite-catalogue realization effects, and the incremental effect of estimated/misspecified inclusion probabilities. These are diagnostic contrasts, not an additive causal decomposition. The oracle values are simulation-only and unavailable to a practical incomplete-catalogue analysis.

The oracle decomposition is implemented separately in draft PR #97; no generator or production likelihood changes are made by that analysis.
