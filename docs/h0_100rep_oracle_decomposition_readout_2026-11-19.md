# Independent oracle decomposition of the 100-replicate H₀ calibration

**Status:** completed post-hoc simulation diagnostic. This is not joint SBC, observational H₀ inference, or a universal estimator comparison.

## Provenance and design

- Calibration workflow: [44411667480](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44411667480).
- Calibration artifact: `14485703115`; SHA-256 `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`.
- Oracle audit workflow: [44433502099](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44433502099).
- Oracle readout artifact: `14487388966`; SHA-256 `126356fbe80d7ab4363485a8e865d0bda4ee79e606ca59b2acb98448fd199b29`.
- The independent oracle likelihood was evaluated on the same 281-point (H_0) grid as the fitted degree-one IPW arm, over [20, 300] km s⁻¹ Mpc⁻¹.
- The preregistered factorial contains 600 cells. The terminal ledger records 599 SUCCESS, 1 NON_EVALUABLE, 0 FAILED, and 0 NOT_RUN.
- The oracle audit independently verifies the complete planned-cell inventory and retains the non-evaluable cell in the planned denominator. Because the oracle posterior depends on the paired parent/catalogue/event realization rather than validation size, the known-π and full-parent oracle values are reconstructed for all 100 replicates in every gamma/validation-size cell. Fitted-arm metrics remain success-only (99 in γ=−1,n=100).

The three arms are:
1. **Estimated degree-one IPW:** the fitted procedure in the fixed-truth pilot.
2. **Known-(pi) IPW oracle:** the observed catalogue weighted by the exact DGP inclusion probabilities, (1/pi(z)).
3. **Full-parent oracle:** the complete finite parent population with its declared host-rate weight.

The two oracle arms are available only because this is a simulation. They are not operational estimators for an incomplete real catalogue.

## Mean posterior-median bias

Units: km s⁻¹ Mpc⁻¹. Estimated-arm values use successful fitted cells (99 or 100); both oracle-arm values use all 100 planned parent/event realizations per cell. MCSE is shown in parentheses.

| γ | Validation n | Estimated degree-one IPW | Known-(pi) IPW | Full-parent oracle |
|---:|---:|---:|---:|---:|
| −1 | 100 | +2.020 (1.751) | +5.769 (1.794) | +7.048 (1.809) |
| −1 | 500 | +2.255 (1.715) | +5.769 (1.794) | +7.048 (1.809) |
| 0 | 100 | +6.428 (1.896) | +6.474 (1.852) | +7.048 (1.809) |
| 0 | 500 | +6.768 (1.875) | +6.474 (1.852) | +7.048 (1.809) |
| +1 | 100 | +10.396 (1.901) | +6.851 (1.823) | +7.048 (1.809) |
| +1 | 500 | +11.114 (1.916) | +6.851 (1.823) | +7.048 (1.809) |

## Paired median-bias contrasts

The paired contrast is computed replicate-by-replicate. Values below are mean difference (MCSE; paired 95% t interval).

| γ | n | Estimated − known-(pi) | Known-(pi) − full parent |
|---:|---:|---:|---:|
| −1 | 100 | −3.593 (0.312; [−4.213, −2.974]) | −1.278 (0.502; [−2.274, −0.282]) |
| −1 | 500 | −3.514 (0.271; [−4.052, −2.977]) | −1.279 (0.502; [−2.274, −0.283]) |
| 0 | 100 | −0.046 (0.295; [−0.632, +0.540]) | −0.574 (0.354; [−1.277, +0.129]) |
| 0 | 500 | +0.294 (0.163; [−0.030, +0.618]) | −0.574 (0.354; [−1.277, +0.129]) |
| +1 | 100 | +3.545 (0.382; [+2.787, +4.304]) | −0.197 (0.252; [−0.696, +0.303]) |
| +1 | 500 | +4.263 (0.260; [+3.748, +4.778]) | −0.197 (0.252; [−0.696, +0.303]) |

The (gamma=-1,n=100) estimated-arm contrast uses 99 paired successful replicates. The oracle arms are reconstructed for all 100 planned parent/event realizations in every gamma/validation-size cell, including the one cell whose fitted arm is NON_EVALUABLE.

## 68% interval coverage

| γ | n | Estimated arm, success-only | Estimated arm, worst–best over all 100 attempts | Known-(pi) oracle | Full-parent oracle |
|---:|---:|---:|---:|---:|---:|
| −1 | 100 | 0.545 | [0.54, 0.55] | 0.490 | 0.540 |
| −1 | 500 | 0.510 | [0.51, 0.51] | 0.490 | 0.540 |
| 0 | 100 | 0.510 | [0.51, 0.51] | 0.500 | 0.540 |
| 0 | 500 | 0.510 | [0.51, 0.51] | 0.500 | 0.540 |
| +1 | 100 | 0.560 | [0.56, 0.56] | 0.510 | 0.540 |
| +1 | 500 | 0.540 | [0.54, 0.54] | 0.510 | 0.540 |

The worst–best interval treats the single non-evaluable attempt as uncovered or covered, respectively. This keeps the full planned denominator visible. The oracle coverage columns use all 100 paired parent/event realizations in every cell.

## Interpretation

1. **A positive baseline bias remains under the full-parent oracle.** The mean full-parent posterior-median bias is about +7.05 km s⁻¹ Mpc⁻¹ across the paired cells. This is already present with the complete finite parent mixture and the same four detected-event realizations. It cannot be attributed solely to estimated catalogue inclusion.
2. **Known-(pi) versus full-parent differences are smaller than the estimated-versus-known-(pi) contrast in the positive-curvature case.** The known-(pi) minus full-parent contrast ranges from about −1.31 to −0.20. It is a diagnostic of finite-catalogue realization effects, with uncertainty.
3. **The estimated-inclusion effect changes sign with curvature.** The degree-one fitted arm is about 3.5 km s⁻¹ Mpc⁻¹ less biased than known-(pi) IPW under (gamma=-1), essentially indistinguishable near (gamma=0), and about 3.5–4.3 more biased under (gamma=+1). This is consistent with curvature-dependent misspecification, not universal superiority of the fitted arm.
4. **Larger validation samples improve weight stability but do not eliminate the (H_0) bias.** The earlier calibration readout shows improved Kish ESS and smaller maximum normalized weights at n=500, while the paired bias differences from n=100 are small and slightly positive.
5. **Coverage is low in this 100-replicate inclusion sample, but this is not a general full-parent calibration conclusion.** The nominal 68% intervals cover about 0.49–0.56 of the 100 event/parent realizations in this particular pilot. A separate 500-replicate full-parent event-count experiment (documented below) found near-nominal 68/90/95% coverage for four events. Thus the low full-parent coverage in this pilot should be treated as sample-specific and not generalized; the fitted-versus-known-π curvature contrasts remain the central result.

## Limits and next steps

This is a synthetic low-redshift, fixed-truth benchmark with four detected events per realization. It is not a joint SBC study, an observational H₀ measurement, a public LVK/GWTC analysis, or validation against a real galaxy survey. The oracle arms are simulation-only. Posterior-median contrasts are diagnostic and should not be interpreted as an additive causal decomposition.

The next scientific step is a joint SBC study with prior-generated H₀ truths and a larger event-count/tail-coverage investigation. The present study remains a controlled synthetic sensitivity benchmark, not a universal calibration result.


## Refined-grid robustness check (561 points)

A separate independent reconstruction on a 561-point H₀ grid completed successfully (workflow 44433653290, artifact 14488165421, ZIP SHA-256 `84475df957aec989b73241a5555ade669fd2d60e693d2cc7d8091d01ffc540d5`). It reconstructed all 600 oracle rows, including the non-evaluable fitted cell by reusing that replicate/curvature's paired parent/catalogue/event realization from its successful validation-size cell.

The 561-point oracle mean median biases were:

| γ | n | Known-π IPW bias | Full-parent bias | Estimated − known-π |
|---:|---:|---:|---:|---:|
| −1 | 100 | +5.770 | +7.048 | −3.594 |
| −1 | 500 | +5.770 | +7.048 | −3.516 |
| 0 | 100 | +6.474 | +7.048 | −0.046 |
| 0 | 500 | +6.474 | +7.048 | +0.294 |
| +1 | 100 | +6.851 | +7.048 | +3.545 |
| +1 | 500 | +6.851 | +7.048 | +4.263 |

At both 281 and 561 points, the known-π biases are +5.769, +5.769, +6.474, +6.474, +6.851, and +6.851 in the same row order; the full-parent biases are +7.048 in all six cells. The updated 281-point workflow and the 561-point reconstruction both retain all 100 oracle realizations per cell, including the non-evaluable fitted-arm replicate.

The estimated-minus-known-π paired contrasts are nearly unchanged between grids: 281-point values −3.593, −3.514, −0.046, +0.294, +3.545, +4.263; 561-point values −3.594, −3.516, −0.046, +0.294, +3.545, +4.263. This shows the curvature-dependent contrast is not a coarse-grid artifact. The updated workflow now audits both 281 and 561 points while preserving all 600 planned oracle rows.



## Independent 500-replicate full-parent event-count cross-check

A separate paired experiment (PR #101; workflow 44434616291; artifact 14487983108; archive SHA-256 `306680d5547d3a2c793585334c44d609fdaa577836afa32c17a20887891d754d`) generated 500 independent parent populations and compared the first four detected events with a nested 30-event sample from each same realization. No catalogue thinning or inclusion fitting was used. The full-parent selected-event likelihood and uniform H₀ prior were held fixed on a 561-point grid.

| Detected events | Mean posterior-median bias | MCSE | RMSE | 68% coverage | 90% coverage | 95% coverage |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | +5.060 | 0.717 | 16.79 | 69.6% | 91.0% | 95.4% |
| 30 | +0.866 | 0.182 | 4.16 | 66.2% | 86.0% | 92.0% |

The paired 30-event minus 4-event median-bias difference is −4.194 km s⁻¹ Mpc⁻¹ (MCSE 0.699; 95% paired t interval [−5.567, −2.820]). Mean 68% interval width contracts from 23.44 to 7.35 km s⁻¹ Mpc⁻¹.

This independent cross-check shows that the four-event positive median bias is replicated (about +5 km s⁻¹ Mpc⁻¹), while four-event credible-interval coverage is close to nominal in the larger sample. Consequently, the 0.49–0.56 full-parent coverage observed in the 100-replicate inclusion pilot should be treated as a feature of that particular realization set, not as a general baseline calibration failure. At 30 events, 68% coverage is near nominal, but 90% and 95% coverage are below nominal in this fixed-truth experiment. This motivates a dedicated tail-calibration and joint-SBC study; it does not by itself establish a likelihood defect.

The event-count cross-check is a separate simulation with a different master seed and no catalogue thinning. It contextualizes the inclusion pilot but does not alter its paired curvature contrasts or turn the pilot into a five-arm A–E comparison.
