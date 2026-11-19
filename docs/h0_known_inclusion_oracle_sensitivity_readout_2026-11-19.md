# Oracle-arm decomposition of the 90-cell fixed-truth H0 pilot

**Status:** completed exploratory post-hoc diagnostic on the corrected 90-cell artifact. This is not confirmatory calibration.

## Artifact and independent reconstruction

Input: GitHub Actions artifact 14476042847 from workflow 44408479569.

- ZIP SHA-256: `cadf1a8ad7c9d2052a73972338fef493cbbcb2b85b2da89527c7e2fe5e0805ba`.
- All 2,110 aggregate manifest entries independently re-hashed; zero mismatches.
- Ledger: 90/90 SUCCESS.
- The archived orchestration report records 90 child audits verified, 30 within-gamma nested validation blocks, 30 cross-gamma validation-pairing blocks, and all nine summary cells reconstructed.
- I independently recomputed the nine summary cells from the replicate ledger; maximum absolute discrepancy across the reported bias, RMSE, coverage, interval-width, and weight-diagnostic columns was below \(4\times10^{-15}\).

The artifact used a 141-point H0 grid (spacing 2 km s\(^{-1}\) Mpc\(^{-1}\)). The separate fixed-realization grid-refinement diagnostic is therefore relevant, and the 100-replicate follow-up uses 281 points.

## Additional diagnostic

I independently reconstructed the selected-event H0 likelihood from the archived observable event distances, catalogue, and truth-only simulation records, without calling the production shared-inclusion likelihood function. The selected likelihood is

\[
\log L(H_0)=\sum_{i=1}^{N_{\rm det}}
\log\!\left[\sum_j w_j\,\mathcal N\!\left(d_i;cz_j/H_0,\sigma^2\right)\right]
-N_{\rm det}\log\!\left[\sum_jw_j\,
\Phi\!\left(\frac{cz_j/H_0-d_{\rm th}}{\sigma}\right)\right].
\]

I computed two simulation-only oracle references:

1. **Known-inclusion-probability IPW:** observed catalogue hosts weighted by the exact DGP probability \(1/\pi(z_j)\), multiplied by the declared merger/base weight and normalized.
2. **Full-parent oracle:** every parent galaxy included with its declared merger/base weight. This is the exact finite parent mixture used to generate the GW hosts.

These are diagnostic benchmarks only. Neither oracle is available to a practical analysis with an incomplete catalogue.

## Results

Mean posterior-median bias (km s\(^{-1}\) Mpc\(^{-1}\)); 10 paired replicates per cell:

| Curvature \(\gamma\) | Validation n | Estimated-IPW bias | Known-\(\pi\) IPW bias | Full-parent oracle bias | Estimated minus known-\(\pi\) median |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +3.95 | +6.45 | +7.58 | -2.50 |
| -1 | 300 | +3.95 | +6.45 | +7.58 | -2.50 |
| -1 | 500 | +3.46 | +6.45 | +7.58 | -2.99 |
| 0 | 100 | +8.21 | +7.43 | +7.58 | +0.79 |
| 0 | 300 | +8.22 | +7.43 | +7.58 | +0.79 |
| 0 | 500 | +7.50 | +7.43 | +7.58 | +0.07 |
| +1 | 100 | +13.29 | +8.84 | +7.58 | +4.45 |
| +1 | 300 | +13.62 | +8.84 | +7.58 | +4.78 |
| +1 | 500 | +12.86 | +8.84 | +7.58 | +4.02 |

The full-parent oracle's mean median bias is +7.58 with Monte Carlo SE 3.15 across the 10 shared replicate/event realizations (paired t 95% interval approximately [+0.46,+14.70]). Its 68% interval coverage is 6/10. This indicates that the positive median bias seen in the pilot is not attributable solely to estimated catalogue inclusion: the same finite event realizations already yield a positive mean median bias under the complete-parent oracle. With only 10 replicates and four detected events, this is a noisy diagnostic, not proof of a systematic likelihood defect.

The estimated-minus-known-\(\pi\) contrast is directionally different under the two misspecified curvature signs:

- \(\gamma=-1\): mean contrast -2.50 at n=100 (paired t 95% interval [-4.27,-0.72]) and -2.99 at n=500 ([-4.73,-1.24]).
- \(\gamma=0\): +0.79 at n=100 ([-0.70,+2.27]) and +0.07 at n=500 ([-0.91,+1.06]); these intervals include zero.
- \(\gamma=+1\): +4.45 at n=100 ([+2.38,+6.53]) and +4.02 at n=500 ([+2.67,+5.38]).

The known-\(\pi\) IPW minus full-parent median contrast is small and uncertain in this 10-replicate sample: -1.13 (95% paired t interval [-5.30,+3.03]) at \(\gamma=-1\), -0.16 ([-2.97,+2.66]) at \(\gamma=0\), and +1.25 ([-0.16,+2.67]) at \(\gamma=+1\). This does not establish that finite-catalogue noise is negligible; the uncertainty is wide.

## Interpretation

This diagnostic supports a more precise decomposition of the next research question:

- **Event-sample / selected-likelihood baseline:** the full-parent oracle isolates the finite event-realization and selected-likelihood contribution for the exact parent mixture.
- **Finite catalogue realization:** known-\(\pi\) IPW versus full-parent oracle measures the catalogue finite-sample component.
- **Estimated inclusion model:** estimated-IPW versus known-\(\pi\) IPW measures the additional impact of learning the inclusion function from finite validation data, including degree-one misspecification and posterior uncertainty.

These posterior-median differences are diagnostic contrasts, not an additive causal decomposition. The sign reversal between \(\gamma=-1\) and \(\gamma=+1\) is consistent with a curvature-dependent misspecification effect, but the 10-replicate pilot is too small for a definitive conclusion.

The 100-replicate follow-up in PR #89 is intended to improve precision on the fixed-truth coverage and curvature contrasts. It does not replace a separate investigation of the full-parent baseline or a joint SBC study.

## Reproduction

The independent calculation is implemented in `scripts/audit_h0_known_inclusion_oracle_sensitivity.py`, with focused numerical tests in `tests/test_h0_known_inclusion_oracle_sensitivity.py`. The script verifies the aggregate artifact manifest, reconstructs the two oracle posteriors, checks that parent/catalogue/event inputs are invariant across validation sizes, and writes cell-level, summary, and paired-contrast CSVs.


## Refined-grid reanalysis (1,121 H0 points)

The archived 90-cell artifact was independently re-read and all 2,110 manifest entries were re-hashed with zero mismatches. The 90-cell data tables are byte-identical between this 141-point artifact and the later 281-point rerun artifact (only the inference grid metadata differs). I recomputed the known-\(\pi\) and full-parent oracle posteriors directly from the raw parent, catalogue, and event tables on a 1,121-point grid over [20,300]. Production posterior medians below come from the matching 281-point rerun artifact (run 44410603590, artifact 14477511669).

| Curvature γ | Validation n | Production median bias (281 grid) | Known-π IPW median bias (1121 grid) | Full-parent median bias (1121 grid) | Production − known-π (paired mean; 95% t CI) |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +3.95 | +6.45 | +7.59 | -2.50 [-4.27, -0.72] |
| -1 | 300 | +3.95 | +6.45 | +7.59 | -2.50 [-4.40, -0.61] |
| -1 | 500 | +3.46 | +6.45 | +7.59 | -2.99 [-4.74, -1.24] |
| 0 | 100 | +8.21 | +7.43 | +7.59 | +0.78 [-0.70, +2.27] |
| 0 | 300 | +8.22 | +7.43 | +7.59 | +0.79 [-0.41, +1.99] |
| 0 | 500 | +7.50 | +7.43 | +7.59 | +0.07 [-0.91, +1.06] |
| +1 | 100 | +13.30 | +8.84 | +7.59 | +4.46 [+2.38, +6.53] |
| +1 | 300 | +13.63 | +8.84 | +7.59 | +4.79 [+3.28, +6.29] |
| +1 | 500 | +12.87 | +8.84 | +7.59 | +4.03 [+2.66, +5.39] |

The full-parent oracle mean median bias is +7.585 km s^-1 Mpc^-1 (event-realization MCSE 3.147; paired t 95% interval approximately [+0.46,+14.71]) across 10 parent/event realizations. The known-\(\pi\) oracle mean biases are +6.451 (MCSE 3.688), +7.427 (MCSE 3.040), and +8.839 (MCSE 3.211) for γ=-1, 0, and +1, respectively. These remain noisy 10-replicate estimates.

The 1121-point oracle posteriors were also compared directly with coarser grids across all 30 replicate × curvature combinations:
- 561→1121: maximum change among the seven reported quantiles was 0.01041 for the full-parent oracle and 0.01126 for known-\(\pi\) IPW.
- 281→1121: corresponding maxima were 0.04847 and 0.05381.
- 141→1121: maxima were 0.21815 and 0.20268.

The oracle medians themselves are much less grid-sensitive than the tail quantiles: the maximum absolute median shift from 281 to 1121 was 0.00263 for the full-parent oracle and 0.00643 for known-\(\pi\) IPW. The earlier 141-grid oracle readout therefore gave broadly similar median-bias decomposition, but its reported interval endpoints were less numerically stable. The independent oracle auditor now accepts an explicit grid resolution and defaults to 1,121 points; CI runs it on the archived artifact.

The interpretation is unchanged but better resolved: positive median bias is present in the full-parent oracle; the estimated-minus-known-\(\pi\) contrast is negative under γ=-1, near zero under γ=0, and positive under γ=+1 in this small paired pilot. Known-\(\pi\) versus full-parent contrasts remain uncertain. These are paired diagnostic contrasts, not an additive causal decomposition or a claim of estimator superiority.
