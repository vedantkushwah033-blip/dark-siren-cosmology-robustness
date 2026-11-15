# Hidden-host assumption sensitivity: audited synthetic results (2026-10-26)

## Executive summary

The hidden-host assumption factorial completed successfully and its independent reconstruction audit passed. It evaluates 105,600 posterior fits across 48 true-data cells, with each simulated observation reused under the omitted-host model, an oracle hidden-host model, and nine fixed assumed fraction/redshift combinations.

Two effects are clear in this simulator:

1. **Omitting a true hidden population progressively degrades performance.** At hidden fractions 0.15, 0.30, and 0.45, the equal-cell mean RMSE difference (omitted minus oracle) was +0.523, +1.402, and +2.286 km s^-1 Mpc^-1, respectively. The corresponding 68% coverage deficits were 3.3, 11.1, and 17.7 percentage points.
2. **An incorrect hidden-host model is not equivalent to the oracle.** When the true hidden fraction was zero, assuming a 45% hidden fraction at redshift 0.030 increased RMSE by 4.409 km s^-1 Mpc^-1 relative to the omitted/oracle model (95% stratified paired-bootstrap interval +4.160 to +4.662). When the true fraction was 0.30, using an assumed fraction of 0.30 and redshift 0.030 still increased RMSE by 0.888 (95% interval +0.660 to +1.128) and lowered 68% coverage by 6.2 percentage points. Matching the fraction alone does not remove redshift misspecification.

The oracle arm is supplied the true hidden fraction and true redshift. It is a simulation-only reference, not an estimator available to a real analysis. The findings are conditional on the declared toy model.

## Reproducibility record

- Successful workflow and independent audit: [GitHub Actions run 40284995493](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40284995493)
- Raw artifact: [hidden-host-assumption-sensitivity](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40284995493/artifacts/12632994472)
- Artifact SHA-256: `94b3cb3fa5167a680eeeeb52c18180ed244b04bf7178947be6911aecbc3e4356`
- Runner/auditor/summarizer revision: `49a4ab544798374516190d6cfc0b147c76635e56`
- Data seed: `20261032`; cell bootstrap seed: `20261033`
- 48 true-data cells × 200 replicates = 9,600 shared synthetic observations.
- 11 inference arms per observation = 105,600 posterior fits.
- 528 cell/arm summaries and 2,880 cell-level paired contrast rows.
- The independent auditor regenerated the observations from the seed and reconstructed all posterior outputs, cell summaries, and cell-level paired-bootstrap intervals. Audit: **PASS**.
- Equal-cell marginal summary script: [summarize_hidden_host_assumption_sensitivity.py](../scripts/summarize_hidden_host_assumption_sensitivity.py). It independently recomputes cell summaries from raw rows before calculating marginal contrasts, and deterministic/reconstruction tests pass.
- The full repository test matrix passed on Python 3.10, 3.11, and 3.12 in [run 40285005529](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40285005529).
- Prespecified design: [hidden-host assumption sensitivity protocol](hidden_host_assumption_sensitivity_protocol.md).

## Design

| Factor | Values |
|---|---|
| Injected H0 | 55, 70, 85 km s^-1 Mpc^-1 |
| Fractional distance-noise scale | 0.05 and 0.20 × c × 0.03/H0 |
| True hidden fraction | 0, 0.15, 0.30, 0.45 |
| True hidden redshift | 0.025, 0.035 |
| Assumed hidden fraction | 0.15, 0.30, 0.45 |
| Assumed hidden redshift | 0.025, 0.030, 0.035 |
| Inference arms | Hidden component omitted; oracle true fraction/redshift; nine fixed assumed combinations |

The prior is uniform on [50,90] and the posterior grid has 4,001 points from 25 to 125. The generator redraws nonpositive distances, and each host-specific Gaussian likelihood is conditioned on positive distance to match that generation rule.

Each of the nine fixed assumed models uses the listed-host weights 0.7/0.3 but inserts its own assumed hidden fraction and redshift. Contrasts are paired within cell and replicate. The marginal estimates below average equally across 12 truth × noise × true-redshift cells within each true hidden-fraction stratum. Their 95% intervals use 1,000 stratified bootstrap resamples with the same resampled replicate indices across arms within each cell.

## Baseline versus oracle across true hidden fractions

RMSE is the root mean squared posterior-median error. Coverage is the fraction of central posterior intervals containing the injected truth.

| True hidden fraction | Arm | Mean median error | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.00 | Omitted | +0.378 | 7.378 | 0.670 | 0.971 | 15.371 | 25.329 |
| 0.00 | Oracle | +0.378 | 7.378 | 0.670 | 0.971 | 15.371 | 25.329 |
| 0.15 | Omitted | +0.165 | 9.059 | 0.631 | 0.929 | 15.549 | 25.635 |
| 0.15 | Oracle | +0.406 | 8.536 | 0.664 | 0.982 | 17.879 | 28.378 |
| 0.30 | Omitted | -0.261 | 10.549 | 0.573 | 0.874 | 15.968 | 26.020 |
| 0.30 | Oracle | -0.100 | 9.147 | 0.684 | 0.978 | 18.929 | 29.192 |
| 0.45 | Omitted | -0.731 | 11.355 | 0.535 | 0.824 | 16.257 | 26.240 |
| 0.45 | Oracle | +0.038 | 9.070 | 0.712 | 0.984 | 18.932 | 29.631 |

### Paired omitted-minus-oracle contrasts

| True hidden fraction | RMSE difference (95% paired CI) | 68% coverage difference (95% CI) | 95% coverage difference (95% CI) |
|---:|---:|---:|---:|
| 0.00 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 0.15 | +0.523 [+0.336, +0.726] | -0.033 [-0.045, -0.020] | -0.053 [-0.062, -0.044] |
| 0.30 | +1.402 [+1.090, +1.716] | -0.111 [-0.128, -0.093] | -0.104 [-0.115, -0.092] |
| 0.45 | +2.286 [+1.938, +2.670] | -0.177 [-0.197, -0.157] | -0.160 [-0.174, -0.146] |

The zero-fraction negative control passed: the omitted and oracle arms produce identical posterior summaries when the true hidden fraction is zero.

## Selected fixed-assumption contrasts against the oracle

Each contrast is assumed model minus oracle. Positive RMSE differences mean worse point-estimation performance; negative coverage differences mean lower empirical coverage. Intervals are 95% stratified paired-bootstrap intervals across the 12 cells in each true-fraction stratum.

| True fraction | Assumed fraction | Assumed redshift | RMSE difference (95% CI) | 68% coverage difference (95% CI) | 95% coverage difference (95% CI) |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.15 | 0.030 | +0.159 [+0.067, +0.261] | -0.015 [-0.028, -0.003] | +0.005 [+0.003, +0.008] |
| 0.00 | 0.45 | 0.030 | +4.409 [+4.160, +4.662] | -0.036 [-0.052, -0.018] | +0.000 [-0.004, +0.005] |
| 0.15 | 0.15 | 0.030 | +0.108 [+0.010, +0.219] | -0.028 [-0.039, -0.017] | -0.018 [-0.025, -0.012] |
| 0.30 | 0.30 | 0.030 | +0.888 [+0.660, +1.128] | -0.062 [-0.078, -0.047] | -0.033 [-0.042, -0.025] |
| 0.45 | 0.45 | 0.030 | +2.154 [+1.907, +2.409] | -0.115 [-0.133, -0.096] | -0.055 [-0.065, -0.047] |

The assumed redshift 0.030 is intermediate between the two true redshifts 0.025 and 0.035. Even when the assumed fraction matches the true fraction, this redshift mismatch reduces coverage in the marginal averages. Larger assumed fractions can be especially damaging when the true fraction is small or zero. These are prespecified stress configurations; the report does not select a globally best assumption by scanning all outcomes.

## Interpretation and limitations

1. **Omission worsens with increasing hidden fraction in this setup.** RMSE rises and 68%/95% coverage falls for the omitted-host model across the tested strata.
2. **Wrong hidden-host assumptions can be worse than omission on RMSE.** Adding a spurious hidden population can shift the posterior away from the truth even when no hidden population exists.
3. **The assumed hidden fraction and redshift interact.** A fraction close to the truth does not guarantee calibration if the assumed redshift is wrong; a compromise redshift can also mask rather than eliminate misspecification.
4. **Coverage is not a single score of model quality.** Some incorrect assumptions increase one coverage measure while worsening RMSE or the other coverage level.
5. The bootstrap intervals quantify simulation variability conditional on this generator. They do not account for model uncertainty, prior uncertainty, the choice of factor levels, or the realism of the data-generating process.

The study remains a simplified Gaussian-distance, finite-host model. It omits realistic gravitational-wave likelihoods, detector selection, peculiar velocities, catalogue selection/completeness uncertainty, and full hierarchical population inference. The oracle arm is not available in real data. No observational H0 estimate or astrophysical validity claim follows from these simulations.
