# Prior-predictive calibration: legacy versus fixed-noise design

## Execution

- Workflow: [Prior-predictive calibration diagnostic](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39949631472)
- Replicates: 1,000 prior draws per noise design
- H0 prior and posterior support: Uniform(50, 90) km/s/Mpc
- Original simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Artifact: `prior-predictive-calibration-1000`

This experiment compares two data-generating designs while leaving the pinned inference code unchanged.

## Noise designs

1. **Legacy fractional-true-distance noise:** the simulator sets `sigma = 0.15 × true distance`. Since true distance depends on the drawn H0, the noise scale changes with the true parameter. The inference likelihood then conditions on the realized sigma as fixed for all candidate H0 values. This is a potential generative/inference mismatch.
2. **Fixed-fiducial noise:** sigma is set to `0.15 × d_L(z=0.04, H0=70)`, a constant independent of the drawn H0 and latent host identity. This matches the likelihood's fixed-sigma assumption more closely.

In both designs, H0 is drawn from the same prior used by the posterior. The host is generated uniformly from the catalogue. The unknown-host likelihood uses a uniform mixture over catalogue redshifts; the known-host control uses the true host alone.

## Posterior-CDF rank diagnostics

The posterior CDF evaluated at the true H0 should be approximately Uniform(0,1) under prior-predictive simulation when the generative model and inference model are consistent. KS tests below compare the 1,000 ranks with Uniform(0,1).

| Noise design | Scenario | Mean rank | Rank SD | KS statistic | KS p-value |
|---|---|---:|---:|---:|---:|
| Legacy fractional-true-distance | Known host | 0.4846 | 0.3008 | 0.0521 | 0.0085 |
| Legacy fractional-true-distance | Unknown host | 0.4895 | 0.3003 | 0.0510 | 0.0106 |
| Fixed-fiducial | Known host | 0.5026 | 0.3003 | 0.0277 | 0.4199 |
| Fixed-fiducial | Unknown host | 0.4949 | 0.3003 | 0.0375 | 0.1176 |

At this sample size, both legacy rank distributions show evidence against uniformity by the descriptive KS test; neither fixed-fiducial rank distribution shows a similarly clear departure. These p-values are diagnostics, not absolute pass/fail criteria, and should be interpreted alongside rank plots, Monte Carlo uncertainty and the exact data-generating assumptions.

## Prior-predictive interval coverage

Coverage averages across H0 values drawn from the prior. It is not the same as fixed-H0 coverage at a particular injected value.

| Scenario | Noise design | 50% coverage | 68% coverage | 80% coverage | 90% coverage | 95% coverage |
|---|---|---:|---:|---:|---:|---:|
| Unknown host | Legacy | 45.2% | 66.9% | 76.8% | 87.1% | 93.2% |
| Known host | Legacy | 45.7% | 65.5% | 77.3% | 88.5% | 94.0% |
| Unknown host | Fixed-fiducial | 46.2% | 66.0% | 76.8% | 86.7% | 93.1% |
| Known host | Fixed-fiducial | 46.0% | 65.0% | 77.4% | 88.0% | 94.0% |

Monte Carlo standard errors for coverage range from roughly 0.75 to 1.58 percentage points across these cells. The 95% coverage values are below 95% in this run, so the fixed-noise model should not be described as perfectly calibrated based on rank tests alone.

## Interpretation

1. The legacy noise design has non-uniform posterior-CDF ranks in both known-host and unknown-host cases, consistent with the suspected mismatch between how sigma is generated and how the likelihood treats it.
2. The fixed-fiducial design has ranks more consistent with uniformity in this 1,000-replicate run. This supports using a clearly specified fixed-sigma model for a controlled self-consistency experiment.
3. The interval coverage rates are broadly similar across the two prior-predictive designs, which illustrates why coverage summaries alone can miss rank-distribution deviations.
4. These results do not erase the fixed-H0 coverage findings. Prior-predictive diagnostics average over the H0 prior and can hide conditional behaviour at individual H0 values.
5. This is still a synthetic model. It does not validate real LVK event likelihoods, galaxy-catalogue selection functions, luminosity weighting, or observational inference.

## Recommended next step

Retain the legacy model as a documented baseline and use the fixed-fiducial noise design as a separate self-consistency control. Before proceeding to the main incompleteness/weighting experiment, compare fixed-H0 coverage and prior-predictive ranks under both models, and explicitly document the noise model in every results file.
