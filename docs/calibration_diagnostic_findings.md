# Calibration diagnostic findings

## Scope

This note separates three diagnostics that answer different questions:

1. Fixed-`H0` repeated-sampling coverage at selected true values.
2. Prior-predictive coverage with `H0` drawn from the prior.
3. Simulation-based calibration (SBC), using posterior CDF ranks under a data-generating process matched to the likelihood.

All results are synthetic-model diagnostics. They are not results from real gravitational-wave observations.

## A. Fixed-`H0` coverage and prior sensitivity

Workflow: [Fixed-H0 coverage diagnostic](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39948705401)

Configuration: 500 replicates at each of seven true values (55, 60, 65, 70, 75, 80, 85 km/s/Mpc), with two uniform prior ranges, [50, 90] and [30, 120], for both known-host and unknown-host inference. The artifact contains 28 scenario/prior/H0 rows and 70 paired coverage contrasts.

Selected 68% coverage results:

| True H0 | Scenario | Prior [50, 90] | Prior [30, 120] |
|---:|---|---:|---:|
| 55 | Known host | 67.6% | 65.6% |
| 55 | Unknown host | 47.2% | 64.4% |
| 70 | Known host | 76.6% | 68.6% |
| 70 | Unknown host | 94.8% | 68.4% |
| 85 | Known host | 76.2% | 67.4% |
| 85 | Unknown host | 75.4% | 92.4% |

At true H0 = 70, changing the prior from [50, 90] to [30, 120] changes unknown-host 68% coverage from 94.8% to 68.4%. The paired difference (wide minus narrow) is -26.4 percentage points, with paired SE 1.97 percentage points. For the known-host control, the corresponding coverages are 76.6% and 68.6%, a paired difference of -8.0 percentage points (SE 1.21 pp).

The effect is not uniform across true H0 values: at H0 = 55 the unknown-host coverage rises from 47.2% to 64.4% with the wider prior, while at H0 = 85 it rises from 75.4% to 92.4%. This is strong evidence that conclusions from the narrow-prior run cannot be generalized without examining prior sensitivity and the mixture posterior's shape.

These are fixed-parameter frequentist coverage diagnostics. Bayesian credible intervals are not guaranteed to have nominal frequentist coverage at every fixed parameter value, so a departure alone is not proof of a coding bug. However, the magnitude and strong dependence on prior range make this a substantive robustness finding.

## B. Prior-predictive diagnostic: legacy noise design versus fixed-fiducial design

Workflow: [Prior-predictive comparison](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/39949631472)

Both runs draw 1,000 values of H0 from Uniform(50, 90). The legacy noise design sets measurement uncertainty to 15% of the true host distance. The fixed-fiducial variant sets it to 15% of the distance at fiducial H0 = 70 for the generated host redshift.

Posterior-rank diagnostics:

| Noise design | Scenario | Mean posterior CDF rank | KS statistic vs Uniform | Descriptive KS p-value |
|---|---|---:|---:|---:|
| Legacy fractional true-distance | Known host | 0.4846 | 0.0521 | 0.0085 |
| Legacy fractional true-distance | Unknown host | 0.4895 | 0.0510 | 0.0106 |
| Fixed-fiducial | Known host | 0.5026 | 0.0277 | 0.420 |
| Fixed-fiducial | Unknown host | 0.4949 | 0.0375 | 0.118 |

The legacy design has ranks less consistent with uniformity in these 1,000 draws. The fixed-fiducial design is more consistent with uniform ranks in this diagnostic, though these descriptive p-values are not proof of calibration. The fixed-fiducial variant still ties the reported uncertainty to host redshift; the independent-noise SBC in section C is the cleaner calibration check.

## C. Simulation-based calibration with a matched noise design

Workflow: [Bayesian simulation-based calibration](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39949582504)

This separate diagnostic draws H0 from Uniform(50, 90), draws measurement uncertainty independently of H0 and host identity, generates the observation from the same Gaussian likelihood used for inference, and evaluates the posterior CDF at the true H0. It uses 1,000 replicates per scenario.

| Scenario | Mean rank | KS statistic vs Uniform | Descriptive KS p-value | Ten-bin rank counts |
|---|---:|---:|---:|---|
| Known host | 0.5115 | 0.0315 | 0.268 | 93, 104, 90, 106, 96, 88, 104, 105, 103, 111 |
| Unknown host | 0.5032 | 0.0176 | 0.910 | 98, 96, 99, 114, 81, 105, 105, 108, 91, 103 |

These ranks are broadly consistent with Uniform(0, 1) for this matched synthetic model. That supports the numerical posterior/CDF implementation under this specific data-generating process; it does not validate the original heteroscedastic noise design or establish realism for actual GW data.

## Interpretation

- Do not collapse fixed-H0 frequentist coverage and prior-predictive SBC into one calibration number; they test different properties.
- The narrow prior [50, 90] strongly affects fixed-H0 unknown-host interval coverage in this experiment.
- The legacy noise model makes the reported measurement uncertainty depend on the true host distance, while the inference likelihood treats the supplied uncertainty as fixed across H0 and candidate hosts. That design can create a mismatch, especially when host identity is unknown.
- The matched-noise SBC gives a useful control: posterior ranks look uniform when the generative model and inference likelihood agree.
- The next experiment should vary catalogue incompleteness and host weighting using a clearly specified, matched observation model, and report both prior-predictive SBC and fixed-H0 coverage across a grid of true H0 values. Keep the original simulator pinned and unchanged for reproducibility; treat corrected noise designs as explicit sensitivity variants.

No claim is made that these results generalize beyond the specified synthetic model.
