# Prior-predictive simulation-based calibration (SBC)

## Execution and provenance

- Workflow: [Bayesian simulation-based calibration](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39949582504)
- Original simulator pinned to `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Replicates: 1,000 prior draws per scenario
- Prior: H0 ~ Uniform(50, 90) km/s/Mpc
- Posterior grid: 1,601 equally spaced points from 50 to 90 km/s/Mpc
- Output artifact: `simulation-based-calibration`

## Generative design

For each replicate, H0 is drawn from the same uniform prior used by the inference model. A mock catalogue is generated, and the true host is selected uniformly from the catalogue. The observed distance is generated from a Gaussian likelihood with a measurement uncertainty drawn independently of H0 and host identity. This independence is intentional: it matches the likelihood's assumption that the supplied distance uncertainty is fixed when evaluating candidate H0 values and candidate hosts.

Two inference cases are evaluated:
- **Known host:** inference uses the true host redshift alone.
- **Unknown host:** inference marginalizes uniformly over all catalogue redshifts, matching the host-generation rule.

For each case, the diagnostic computes the posterior CDF evaluated at the true H0. Under a correctly specified prior-predictive model, these ranks should be approximately uniform. This is a check of internal simulation/inference consistency, not a test of whether the model describes real gravitational-wave observations.

## Results

| Scenario | Replicates | Mean rank | KS statistic vs Uniform(0,1) | Descriptive KS p-value |
|---|---:|---:|---:|---:|
| Known host | 1,000 | 0.5115 | 0.0315 | 0.268 |
| Unknown host | 1,000 | 0.5032 | 0.0176 | 0.910 |

Ten-bin rank counts:

- Known host: `93;104;90;106;96;88;104;105;103;111`
- Unknown host: `98;96;99;114;81;105;105;108;91;103`

The rank distributions show no conspicuous departure from uniformity in this 1,000-replicate run. The KS p-values are descriptive diagnostics, not pass/fail proof; they depend on finite Monte Carlo sampling and the chosen test.

## Interpretation

1. This prior-predictive SBC experiment is consistent with the inference code being self-consistent under this specific synthetic model and fixed-noise design.
2. It does **not** erase the fixed-H0 coverage findings. Fixed-H0 repeated-sampling coverage and prior-predictive SBC answer different questions: fixed-H0 coverage conditions on a particular true value, while SBC averages over the model's prior.
3. The original baseline's legacy simulator sets the distance-noise scale as a fraction of the true distance. In the inference likelihood, the realized scale is held fixed across candidate H0 values. That is a different repeated-sampling design from this SBC experiment, where sigma is drawn independently of H0 and host identity. The legacy design should be assessed separately rather than conflated with this result.
4. The unknown-host result checks only a uniform host prior over the simulated catalogue with no event-selection model. It does not validate luminosity weighting, catalogue incompleteness, real survey selection, or real LVK data.

## Next step

Compare the legacy fractional-distance-noise design with a deliberately specified fixed-noise design, documenting exactly what information the reported uncertainty carries. Then re-run fixed-H0 coverage and the planned incompleteness/weighting contrasts under each explicitly defined model. Keep the original pinned simulator immutable and treat alternatives as separate model variants.
