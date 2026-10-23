# Baseline calibration audit (1,000 simulated universes)

## Status and provenance

This report records an executed validation run, not an observational cosmology result.

- Original simulator repository: `vedantkushwah033-blip/bayesian-dark-siren-cosmology`
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Base seed: `2026`
- Replicates per scenario: `1000`
- GitHub Actions run: [Baseline calibration reproduction](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39939117327)
- Workflow artifact: `baseline-calibration-1000`
- Validation: the workflow's calibration-artifact validator passed.

The simulation uses a generated galaxy population, a Gaussian gravitational-wave distance likelihood, a finite H0 grid from 50 to 90 km/s/Mpc, and the simplified low-redshift relation d_L = c z / H0. These results are conditional on that model and configuration. They are not measurements from real LVK events.

## Repeated-sampling summaries

The injected value is H0 = 70 km/s/Mpc.

| Metric | Unknown-host mixture, complete catalogue | Known-host control |
|---|---:|---:|
| Replicates | 1,000 | 1,000 |
| Mean posterior-median error (reported as median_bias by legacy summarizer) | +1.986 km/s/Mpc | +1.475 km/s/Mpc |
| Mean posterior-mean error | +1.259 km/s/Mpc | +1.616 km/s/Mpc |
| RMSE of posterior median | 6.190 km/s/Mpc | 7.177 km/s/Mpc |
| Mean 68% interval width | 22.372 km/s/Mpc | 17.045 km/s/Mpc |
| Mean 95% interval width | 35.820 km/s/Mpc | 29.959 km/s/Mpc |
| Empirical 68% interval coverage | 94.4% (MC SE 0.73 percentage points) | 74.9% (MC SE 1.37 percentage points) |
| Empirical 95% interval coverage | 99.9% (MC SE 0.10 percentage points) | 98.6% (MC SE 0.37 percentage points) |
| Mean standardized posterior-median error (pull) | +0.143 | +0.200 |
| SD of standardized posterior-median error (pull) | 0.670 | 1.041 |

The legacy summary field named `median_bias` is actually the **mean** of replicate-level posterior-median errors. The factorial script explicitly recomputes the median error when it uses that name. Keep this naming caveat in mind when comparing reports.

## Coverage curve

Coverage is the fraction of simulated intervals containing the injected value. Monte Carlo standard errors (MC SE) quantify finite-replicate uncertainty only; they do not address model misspecification.

| Scenario | Nominal interval | Empirical coverage | MC SE | Mean interval width (km/s/Mpc) |
|---|---:|---:|---:|---:|
| Unknown-host | 50% | 88.3% | 1.02 pp | 15.697 |
| Unknown-host | 68% | 94.4% | 0.73 pp | 22.372 |
| Unknown-host | 80% | 97.4% | 0.50 pp | 27.561 |
| Unknown-host | 90% | 99.1% | 0.30 pp | 32.717 |
| Unknown-host | 95% | 99.9% | 0.10 pp | 35.820 |
| Known-host | 50% | 55.6% | 1.57 pp | 11.792 |
| Known-host | 68% | 74.9% | 1.37 pp | 17.045 |
| Known-host | 80% | 87.0% | 1.06 pp | 21.420 |
| Known-host | 90% | 95.4% | 0.66 pp | 26.345 |
| Known-host | 95% | 98.6% | 0.37 pp | 29.958 |

## Interpretation and next steps

1. The 68% and 95% equal-tail intervals over-cover in both controls in this repeated-sampling experiment; the effect is much larger for the unknown-host mixture.
2. The unknown-host mixture produces wider intervals and a substantially narrower pull distribution than the known-host control. This is a diagnostic of the current mixture model's uncertainty behaviour, not proof that either interval construction is generally calibrated or uncalibrated in all settings.
3. The posterior-median errors are positively biased in both controls under this configuration.
4. Do not treat the 20-replicate paired-seed smoke test as inferential evidence. It only established that the pipeline runs and emits structurally valid files.
5. Before interpreting catalogue-incompleteness effects, investigate why the unknown-host baseline's nominal 68% intervals cover 94.4% of repeated draws. Check the mixture likelihood, equal-tail interval calculation, fixed-H0 repeated-sampling design, and prior/grid effects. Keep this audit separate from changes to the original pinned simulator.
6. After that diagnostic, run the prespecified factorial experiment and compare its controls to this recorded baseline. Report negative or inconclusive results as well as any effects.

No claim is made here that the model is scientifically realistic or that the findings generalize beyond this synthetic setup.
