# Independent oracle decomposition of the 100-replicate H₀ calibration

**Status:** completed post-hoc simulation diagnostic. This is not joint SBC, observational H₀ inference, or a universal estimator comparison.

## Provenance and design

- Calibration workflow: [44411667480](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44411667480).
- Calibration artifact: `14485703115`; SHA-256 `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`.
- Oracle audit workflow: [44433502099](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44433502099).
- Oracle readout artifact: `14487388966`; SHA-256 `126356fbe80d7ab4363485a8e865d0bda4ee79e606ca59b2acb98448fd199b29`.
- The independent oracle likelihood was evaluated on the same 281-point (H_0) grid as the fitted degree-one IPW arm, over [20, 300] km s⁻¹ Mpc⁻¹.
- The preregistered factorial contains 600 cells. The terminal ledger records 599 SUCCESS, 1 NON_EVALUABLE, 0 FAILED, and 0 NOT_RUN.
- The oracle audit independently verifies the complete planned-cell inventory and retains the non-evaluable cell in the planned denominator. Oracle posterior reconstruction is based on the 599 successful cells; the one non-evaluable cell is not silently relabelled as an oracle success.

The three arms are:
1. **Estimated degree-one IPW:** the fitted procedure in the fixed-truth pilot.
2. **Known-(pi) IPW oracle:** the observed catalogue weighted by the exact DGP inclusion probabilities, (1/pi(z)).
3. **Full-parent oracle:** the complete finite parent population with its declared host-rate weight.

The two oracle arms are available only because this is a simulation. They are not operational estimators for an incomplete real catalogue.

## Mean posterior-median bias

Units: km s⁻¹ Mpc⁻¹. Each value is the mean across successful paired replicates; MCSE is shown in parentheses.

| γ | Validation n | Estimated degree-one IPW | Known-(pi) IPW | Full-parent oracle |
|---:|---:|---:|---:|---:|
| −1 | 100 | +2.020 (1.751) | +5.613 (1.805) | +6.924 (1.823) |
| −1 | 500 | +2.255 (1.715) | +5.769 (1.794) | +7.048 (1.809) |
| 0 | 100 | +6.428 (1.896) | +6.474 (1.852) | +7.048 (1.809) |
| 0 | 500 | +6.768 (1.875) | +6.474 (1.852) | +7.048 (1.809) |
| +1 | 100 | +10.396 (1.901) | +6.851 (1.823) | +7.048 (1.809) |
| +1 | 500 | +11.114 (1.916) | +6.851 (1.823) | +7.048 (1.809) |

## Paired median-bias contrasts

The paired contrast is computed replicate-by-replicate. Values below are mean difference (MCSE; paired 95% t interval).

| γ | n | Estimated − known-(pi) | Known-(pi) − full parent |
|---:|---:|---:|---:|
| −1 | 100 | −3.593 (0.312; [−4.213, −2.973]) | −1.311 (0.506; [−2.315, −0.308]) |
| −1 | 500 | −3.514 (0.271; [−4.052, −2.977]) | −1.279 (0.502; [−2.274, −0.283]) |
| 0 | 100 | −0.046 (0.295; [−0.632, +0.540]) | −0.574 (0.354; [−1.277, +0.129]) |
| 0 | 500 | +0.294 (0.163; [−0.030, +0.618]) | −0.574 (0.354; [−1.277, +0.129]) |
| +1 | 100 | +3.545 (0.382; [+2.787, +4.304]) | −0.197 (0.252; [−0.696, +0.303]) |
| +1 | 500 | +4.263 (0.260; [+3.748, +4.778]) | −0.197 (0.252; [−0.696, +0.303]) |

The (gamma=-1,n=100) estimated-arm contrast uses 99 paired successful replicates; all other fitted-arm contrasts use 100. Oracle contrast denominators are 99 for that same cell and 100 otherwise.

## 68% interval coverage

| γ | n | Estimated arm, success-only | Estimated arm, worst–best over all 100 attempts | Known-(pi) oracle | Full-parent oracle |
|---:|---:|---:|---:|---:|---:|
| −1 | 100 | 0.545 | [0.54, 0.55] | 0.495 | 0.545 |
| −1 | 500 | 0.510 | [0.51, 0.51] | 0.490 | 0.540 |
| 0 | 100 | 0.510 | [0.51, 0.51] | 0.500 | 0.540 |
| 0 | 500 | 0.510 | [0.51, 0.51] | 0.500 | 0.540 |
| +1 | 100 | 0.560 | [0.56, 0.56] | 0.510 | 0.540 |
| +1 | 500 | 0.540 | [0.54, 0.54] | 0.510 | 0.540 |

The worst–best interval treats the single non-evaluable attempt as uncovered or covered, respectively. This keeps the full planned denominator visible. The oracle coverage columns use the successful-cell denominators shown above.

## Interpretation

1. **A positive baseline bias remains under the full-parent oracle.** The mean full-parent posterior-median bias is about +7.05 km s⁻¹ Mpc⁻¹ across the paired cells. This is already present with the complete finite parent mixture and the same four detected-event realizations. It cannot be attributed solely to estimated catalogue inclusion.
2. **Known-(pi) versus full-parent differences are smaller than the estimated-versus-known-(pi) contrast in the positive-curvature case.** The known-(pi) minus full-parent contrast ranges from about −1.31 to −0.20. It is a diagnostic of finite-catalogue realization effects, with uncertainty.
3. **The estimated-inclusion effect changes sign with curvature.** The degree-one fitted arm is about 3.5 km s⁻¹ Mpc⁻¹ less biased than known-(pi) IPW under (gamma=-1), essentially indistinguishable near (gamma=0), and about 3.5–4.3 more biased under (gamma=+1). This is consistent with curvature-dependent misspecification, not universal superiority of the fitted arm.
4. **Larger validation samples improve weight stability but do not eliminate the (H_0) bias.** The earlier calibration readout shows improved Kish ESS and smaller maximum normalized weights at n=500, while the paired bias differences from n=100 are small and slightly positive.
5. **Coverage is poor in all arms.** The nominal 68% intervals cover only about 0.49–0.56 of successful realizations. Do not describe these posteriors as calibrated.

## Limits and next steps

This is a synthetic low-redshift, fixed-truth benchmark with four detected events per realization. It is not a joint SBC study, an observational H₀ measurement, a public LVK/GWTC analysis, or validation against a real galaxy survey. The oracle arms are simulation-only. Posterior-median contrasts are diagnostic and should not be interpreted as an additive causal decomposition.

The next scientific step is to investigate the positive full-parent baseline bias separately (selected-event likelihood and finite four-event realization effects), then design a joint SBC study with enough events and prior-generated truths. The present study is best reported as a controlled synthetic sensitivity benchmark.
