# Independent oracle decomposition of the 100-replicate H₀ calibration

**Status:** completed independent post-hoc simulation diagnostic. Not joint SBC, not observational H₀ inference, and not a universal estimator comparison.

## Provenance and integrity

- Source calibration workflow: [44411667480](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44411667480).
- Source artifact: `14485703115`; SHA-256 `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`.
- Independent oracle workflow: [44434500153](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44434500153).
- Oracle artifact: `14488192032`; SHA-256 `8476e6e0f04892ba5b10967671b59a91141b27dd95d11357378f51ee1226db84`.
- Oracle H₀ grid: 1,121 points over [20, 300] km s⁻¹ Mpc⁻¹.
- Exact preregistered inventory: 600 planned cells; 599 SUCCESS, one NON_EVALUABLE, zero FAILED, zero NOT_RUN.
- Independent checks: 600 child audits; 300 nested-validation blocks; 200 cross-curvature validation-pairing blocks; 14,106 aggregate-manifest entries; 100 paired replicate blocks and 300 replicate/curvature blocks.
- The non-evaluable cell is `gm1_n100_r057` (replicate 57, γ=−1, n=100); the fitted inclusion quadrature did not converge through orders 3,5,7,9. It remains in the planned denominator.

## Arms and likelihood

The audit reconstructs three arms from the archived raw parent, catalogue, and event data:

1. **Estimated degree-one IPW:** production fitted-inclusion arm.
2. **Known-(pi) IPW oracle:** observed catalogue weighted by the exact DGP inclusion probability (1/pi(z)).
3. **Full-parent oracle:** complete finite parent population with the declared host-rate weight.

The oracle arms are available only because this is a simulation. They are not practical estimators for incomplete real catalogues. The selected-event likelihood is independently reconstructed from discrete host mixtures and the detection-truncation denominator.

## Mean posterior-median bias

Units: km s⁻¹ Mpc⁻¹. Parent/event oracle values use all 100 paired realizations in each cell; fitted-arm metrics use successful cells only.

| γ | Validation n | Estimated degree-one IPW (MCSE) | Known-(pi) IPW (MCSE) | Full-parent oracle (MCSE) |
|---:|---:|---:|---:|---:|
| −1 | 100 | +2.020 (1.751) | +5.771 (1.793) | +7.048 (1.809) |
| −1 | 500 | +2.255 (1.715) | +5.771 (1.793) | +7.048 (1.809) |
| 0 | 100 | +6.428 (1.896) | +6.474 (1.852) | +7.048 (1.809) |
| 0 | 500 | +6.768 (1.875) | +6.474 (1.852) | +7.048 (1.809) |
| +1 | 100 | +10.396 (1.901) | +6.852 (1.823) | +7.048 (1.809) |
| +1 | 500 | +11.114 (1.916) | +6.852 (1.823) | +7.048 (1.809) |

## Paired median-bias contrasts

Mean difference (MCSE; paired 95% t interval):

| γ | n | Estimated − known-(pi) | Known-(pi) − full parent |
|---:|---:|---:|---:|
| −1 | 100 | −3.595 (0.312; [−4.214, −2.975]) | −1.278 (0.502; [−2.274, −0.281]) |
| −1 | 500 | −3.516 (0.271; [−4.054, −2.978]) | −1.278 (0.502; [−2.274, −0.281]) |
| 0 | 100 | −0.046 (0.295; [−0.632, +0.539]) | −0.574 (0.355; [−1.277, +0.130]) |
| 0 | 500 | +0.294 (0.163; [−0.031, +0.618]) | −0.574 (0.355; [−1.277, +0.130]) |
| +1 | 100 | +3.545 (0.382; [+2.786, +4.303]) | −0.197 (0.252; [−0.697, +0.303]) |
| +1 | 500 | +4.262 (0.260; [+3.747, +4.778]) | −0.197 (0.252; [−0.697, +0.303]) |

The γ=−1,n=100 estimated-arm contrast uses 99 paired successful cells; all other fitted-arm contrasts use 100. Oracle contrasts use 100 paired parent/event realizations in every cell.

## Fixed-truth interval coverage

Coverage is success-only for the fitted arm; oracle coverage uses all 100 paired realizations.

| γ | n | Estimated IPW 68/90/95% | Known-(pi) oracle 68/90/95% | Full-parent oracle 68/90/95% |
|---:|---:|---:|---:|---:|
| −1 | 100 | 0.545 / 0.727 / 0.798 | 0.49 / 0.74 / 0.79 | 0.54 / 0.83 / 0.91 |
| −1 | 500 | 0.510 / 0.740 / 0.780 | 0.49 / 0.74 / 0.79 | 0.54 / 0.83 / 0.91 |
| 0 | 100 | 0.510 / 0.790 / 0.860 | 0.50 / 0.78 / 0.86 | 0.54 / 0.83 / 0.91 |
| 0 | 500 | 0.510 / 0.770 / 0.850 | 0.50 / 0.78 / 0.86 | 0.54 / 0.83 / 0.91 |
| +1 | 100 | 0.560 / 0.820 / 0.880 | 0.51 / 0.82 / 0.89 | 0.54 / 0.83 / 0.91 |
| +1 | 500 | 0.540 / 0.720 / 0.880 | 0.51 / 0.82 / 0.89 | 0.54 / 0.83 / 0.91 |

The one non-evaluable fitted cell yields all-attempt coverage bounds: at γ=−1,n=100, 68% coverage is [0.54,0.55], 90% is [0.72,0.73], and 95% is [0.79,0.80]. All other cells have no non-evaluable attempts, so their bounds equal the success-only estimates.

## Interpretation

1. **Positive baseline bias exists under the full-parent oracle.** The mean full-parent posterior-median bias is +7.048 km s⁻¹ Mpc⁻¹ (MCSE 1.809). This is already present with the complete finite parent mixture and the same four-event realizations; it cannot be attributed solely to catalogue incompleteness or fitted inclusion probabilities.
2. **The estimated-versus-known-(pi) contrast depends strongly on curvature.** It is about −3.59 under γ=−1, near zero under γ=0, and +3.54 to +4.26 under γ=+1. This is consistent with curvature-dependent misspecification; it is not evidence of universal estimator superiority.
3. **Known-(pi) versus full-parent differences are comparatively modest.** Their paired mean contrast ranges from −1.28 to −0.20, with uncertainty. This diagnoses finite-catalogue realization effects in this synthetic setup.
4. **Larger validation samples improve weight diagnostics but do not eliminate H₀ bias.** The fitted-arm bias changes only slightly from n=100 to n=500, while the previous readout showed improved Kish ESS and smaller maximum normalized weights.
5. **All three arms show undercoverage.** The full-parent oracle covers only 54% of truths with nominal 68% intervals, 83% with nominal 90%, and 91% with nominal 95% intervals. The baseline selected-event likelihood and four-event design therefore need separate investigation.

## Limits and next experiment

This is a synthetic low-redshift Euclidean fixed-truth study with four detected events per realization. It is not joint SBC, observational H₀ inference, public LVK/GWTC analysis, or validation against a real galaxy survey. The oracle arms are simulation-only. Posterior-median contrasts are diagnostics, not an additive causal decomposition.

A paired full-parent event-count experiment is underway to compare four detected events with a nested 30-event sample from the same realizations. That will test whether the baseline bias and undercoverage primarily reflect the very small event count before making claims about inclusion correction.
