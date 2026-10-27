> **Revalidation note — 2026-10-27:** This factorial was rerun on corrected revision `ecb260886716c54fdd2d5e47390701949f5a18d5`; the independent reconstruction passed. A row-by-row comparison of historical artifact `12632769624` and corrected artifact `12695694723` found identical categorical fields and summary/contrast values at displayed precision; the maximum absolute difference in any raw numeric field was `2.37e-12`, consistent with floating-point evaluation. Corrected artifact SHA-256: `0e82705ebc2cb6808bb320a0b3e7e9b94a4d8fa2bcb9dd51cf567366fdcc4fab`. See the [corrected rerun provenance report](corrected_hidden_host_rerun_provenance_2026-10-27.md). These remain conditional results from a synthetic model, not astrophysical validation.

# Prior-support × hidden-host interaction: audited synthetic results (2026-10-26)

## Executive summary

A new paired factorial experiment crossed prior support with hidden-host misspecification and passed independent reconstruction. The experiment reused each synthetic observation across all three priors and all four inference arms, then used stratified paired bootstrap resampling to estimate equal-cell marginal differences.

The central result is **not** that wider priors are always better. Wider supports markedly improved interval coverage when the injected truth was near the edge of the baseline [50,90] interval, but also increased RMSE in the equal-cell averages. The size of that trade-off depended on the hidden-host fraction and inference arm. In the [50,90] prior, omission of a hidden population at fractions 0.30 and 0.45 increased RMSE and sharply reduced 68% and 95% coverage. The correctly specified hidden-host arm improved RMSE and coverage in that narrow-prior setting, but this arm is a hybrid oracle-style diagnostic because it is supplied the true hidden fraction and redshift. Under wider priors, the coverage advantage was smaller or reversed for some metrics.

These findings apply only to the specified synthetic generator and finite design; they are not an observational H0 estimate or real-data validation.

## Reproducibility and audit record

- Full workflow and independent audit: [GitHub Actions run 40283661686](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40283661686)
- Raw and summarized artifact: [host-prior-support-factorial](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40283661686/artifacts/12632769624)
- Artifact SHA-256: `61f78104305fc7d3fb4992e0841b7aa8addc9b27864bd090a6bf719243279c16`
- **Corrected-revision rerun:** [workflow 40417072086](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072086), [artifact 12695694723](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072086/artifacts/12695694723); revision `ecb260886716c54fdd2d5e47390701949f5a18d5`; independent reconstruction passed.
- Runner/auditor/summarizer revision: `8f576dccb48ae7d5ae1b2106228f918fa0bc4cf1`
- Data seed: `20261030`; cell-bootstrap seed: `20261031`
- 36 shared-data cells × 200 replicates = 7,200 synthetic observations.
- 3 prior supports × 4 inference arms were applied to every observation: 86,400 posterior fits.
- 432 cell/prior/arm summaries, 1,296 within-prior arm contrasts, and 1,728 prior contrasts.
- The independent auditor reconstructed seeded observations, all 86,400 posterior summaries, all cell summaries, and both families of paired-bootstrap intervals. Audit status: **PASS**.
- The latest repository test matrix passed on Python 3.10, 3.11, and 3.12: [CI run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40283661740).
- Protocol: [prespecified factorial design](host_prior_support_factorial_protocol.md).
- Reproducible marginal summaries: [summary script](../scripts/summarize_host_prior_support_factorial.py).

## Frozen design

| Factor | Values |
|---|---|
| Injected H0 | 55, 70, 85 km s^-1 Mpc^-1 |
| Fractional distance-noise scale | 0.05 and 0.20 × c × 0.03/H0 |
| Hidden-host fraction | 0, 0.30, 0.45 |
| Hidden-host redshift | 0.025, 0.035 |
| Uniform prior support | [50,90], [40,100], [30,120] |
| Inference arms | Base weights 0.7/0.3; reversed 0.3/0.7; equal 0.5/0.5; hidden-host mixture |
| Posterior grid | 4,001 points from H0=25 to 125 |

The generator redraws nonpositive distances. Each host-specific likelihood is therefore conditioned on positive distance. The hidden-host arm uses the exact simulated hidden fraction and hidden redshift; it is a controlled misspecification benchmark, not a realistic treatment of unknown catalogue completeness.

The marginal estimates below average equally over the 12 truth × noise × hidden-redshift cells at each hidden fraction. Bootstrap intervals use 1,000 stratified resamples, with replicate indices paired across all priors and arms within each base cell. They quantify simulation uncertainty conditional on this design, not model uncertainty.

## Equal-cell results: base weights versus hidden host modelled

RMSE is the root mean squared posterior-median error in km s^-1 Mpc^-1. Coverage is the empirical fraction of central posterior intervals containing the injected truth.

| Hidden fraction | Prior | Arm | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---|---:|---:|---:|---:|---:|
| 0.00 | [50,90] | Base weights | 7.225 | 0.680 | 0.980 | 15.289 | 25.387 |
| 0.00 | [40,100] | Base weights | 10.686 | 0.777 | 0.979 | 23.213 | 37.598 |
| 0.00 | [30,120] | Base weights | 16.819 | 0.745 | 0.974 | 35.732 | 53.925 |
| 0.30 | [50,90] | Base weights | 10.702 | 0.573 | 0.872 | 15.793 | 25.902 |
| 0.30 | [40,100] | Base weights | 12.093 | 0.769 | 0.927 | 25.414 | 39.732 |
| 0.30 | [30,120] | Base weights | 16.755 | 0.790 | 0.960 | 38.222 | 55.898 |
| 0.45 | [50,90] | Base weights | 11.868 | 0.505 | 0.810 | 16.393 | 26.508 |
| 0.45 | [40,100] | Base weights | 12.957 | 0.753 | 0.891 | 26.865 | 41.080 |
| 0.45 | [30,120] | Base weights | 16.783 | 0.823 | 0.960 | 40.319 | 57.676 |
| 0.00 | [50,90] | Hidden host modelled | 7.225 | 0.680 | 0.980 | 15.289 | 25.387 |
| 0.00 | [40,100] | Hidden host modelled | 10.686 | 0.777 | 0.979 | 23.213 | 37.598 |
| 0.00 | [30,120] | Hidden host modelled | 16.819 | 0.745 | 0.974 | 35.732 | 53.925 |
| 0.30 | [50,90] | Hidden host modelled | 9.090 | 0.672 | 0.978 | 19.062 | 29.273 |
| 0.30 | [40,100] | Hidden host modelled | 10.572 | 0.785 | 0.977 | 27.204 | 41.412 |
| 0.30 | [30,120] | Hidden host modelled | 15.583 | 0.735 | 0.970 | 38.662 | 56.906 |
| 0.45 | [50,90] | Hidden host modelled | 9.235 | 0.712 | 0.983 | 19.229 | 29.953 |
| 0.45 | [40,100] | Hidden host modelled | 10.998 | 0.792 | 0.980 | 27.146 | 42.095 |
| 0.45 | [30,120] | Hidden host modelled | 15.577 | 0.760 | 0.974 | 37.900 | 57.355 |

## Paired contrasts: hidden host omitted minus modelled

A positive RMSE difference means omission performed worse. Negative coverage differences mean the omitted-host model covered the truth less often.

| Hidden fraction | Prior | RMSE difference (95% paired CI) | 68% coverage difference (95% CI) | 95% coverage difference (95% CI) |
|---:|---|---:|---:|---:|
| 0.00 | [50,90] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 0.00 | [40,100] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 0.00 | [30,120] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 0.30 | [50,90] | +1.612 [+1.320, +1.913] | -0.099 [-0.117, -0.083] | -0.106 [-0.118, -0.095] |
| 0.30 | [40,100] | +1.521 [+1.211, +1.836] | -0.016 [-0.031, -0.002] | -0.050 [-0.059, -0.041] |
| 0.30 | [30,120] | +1.172 [+0.836, +1.513] | +0.056 [+0.042, +0.068] | -0.010 [-0.016, -0.004] |
| 0.45 | [50,90] | +2.633 [+2.270, +2.997] | -0.207 [-0.227, -0.187] | -0.173 [-0.188, -0.160] |
| 0.45 | [40,100] | +1.959 [+1.587, +2.319] | -0.039 [-0.056, -0.022] | -0.089 [-0.100, -0.079] |
| 0.45 | [30,120] | +1.206 [+0.775, +1.634] | +0.063 [+0.050, +0.078] | -0.014 [-0.021, -0.007] |

The zero-fraction control passes exactly. For the narrow baseline prior, modelling the hidden component lowers RMSE and substantially improves both coverage levels at nonzero hidden fractions. Under the widest prior, however, the omitted-host arm has higher empirical 68% coverage than the modelled-host arm in these marginal averages, despite having higher RMSE. Coverage alone is therefore not an adequate ranking criterion, and the modelled arm's exact knowledge of hidden fraction/redshift makes it a diagnostic benchmark rather than a deployable estimator.

## Paired contrasts: prior support versus [50,90]

Each cell is an equal-cell marginal difference for the listed metric, with a 95% stratified paired-bootstrap interval. A positive RMSE difference means the wider prior increased RMSE; positive coverage differences mean higher empirical coverage.

| Hidden fraction | Arm | Prior comparison | RMSE difference (95% CI) | 68% coverage difference (95% CI) | 95% coverage difference (95% CI) |
|---:|---|---|---:|---:|---:|
| 0.00 | Base weights | [40,100] − [50,90] | +3.461 [+3.158, +3.764] | +0.098 [+0.083, +0.112] | -0.001 [-0.005, +0.002] |
| 0.00 | Base weights | [30,120] − [50,90] | +9.594 [+9.155, +10.029] | +0.066 [+0.047, +0.084] | -0.006 [-0.011, -0.002] |
| 0.30 | Base weights | [40,100] − [50,90] | +1.391 [+1.083, +1.697] | +0.196 [+0.180, +0.213] | +0.055 [+0.047, +0.064] |
| 0.30 | Base weights | [30,120] − [50,90] | +6.053 [+5.572, +6.496] | +0.217 [+0.197, +0.239] | +0.087 [+0.075, +0.100] |
| 0.45 | Base weights | [40,100] − [50,90] | +1.089 [+0.813, +1.376] | +0.247 [+0.230, +0.263] | +0.082 [+0.072, +0.092] |
| 0.45 | Base weights | [30,120] − [50,90] | +4.915 [+4.469, +5.333] | +0.318 [+0.298, +0.338] | +0.150 [+0.138, +0.164] |
| 0.30 | Hidden host modelled | [40,100] − [50,90] | +1.483 [+1.191, +1.747] | +0.113 [+0.098, +0.131] | -0.001 [-0.005, +0.003] |
| 0.30 | Hidden host modelled | [30,120] − [50,90] | +6.493 [+6.052, +6.900] | +0.062 [+0.043, +0.081] | -0.008 [-0.015, -0.002] |
| 0.45 | Hidden host modelled | [40,100] − [50,90] | +1.762 [+1.512, +2.028] | +0.080 [+0.066, +0.095] | -0.003 [-0.006, +0.001] |
| 0.45 | Hidden host modelled | [30,120] − [50,90] | +6.342 [+5.906, +6.751] | +0.048 [+0.030, +0.067] | -0.009 [-0.015, -0.004] |

The wider supports increase RMSE in the equal-cell averages for both base and hidden-modelled arms, while generally improving 68% coverage. At high hidden fractions, 95% coverage also improves for the omitted-host arms. This is a real trade-off in this generator: wider support can reduce boundary-driven undercoverage while allowing more posterior mass in regions that increase point-estimate error. It is not evidence that the narrow prior is scientifically preferable; prior choice must be justified independently of these results.

## Limits and next research priorities

1. The experiment uses one seeded run and 200 replicates per base cell. Bootstrap intervals quantify finite-simulation variability but not variation across generators or model structures.
2. The data model is a low-redshift Gaussian distance mixture. It omits realistic GW likelihoods, detector selection, peculiar velocities, catalogue selection/completeness uncertainty, and full hierarchical population inference.
3. The hidden-host arm knows the true hidden fraction and redshift; real analyses would need to infer or marginalize over those quantities.
4. The chosen prior intervals are sensitivity scenarios, not claims about the scientifically correct prior. The experiment does not identify a universally optimal prior.
5. The design selects truths 55, 70, and 85 to expose edge effects in the [50,90] baseline support. The findings should not be extrapolated to all H0 values or prior families.

Next, the most valuable work is to test robustness to uncertain hidden fraction/redshift (rather than oracle inputs), include multiple independent simulation seeds, and benchmark the likelihood/selection structure against a realistic simulated or public gravitational-wave/catalogue case before making astrophysical claims.
