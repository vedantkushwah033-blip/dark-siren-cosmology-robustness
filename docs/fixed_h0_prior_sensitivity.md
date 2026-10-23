# Fixed-H0 coverage and prior-range sensitivity

## Purpose

This diagnostic checks whether fixed-H0 repeated-sampling coverage changes when the posterior's uniform prior support is widened. It uses the original synthetic simulator pinned to `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. It is a sensitivity check, not a real-event cosmology result.

- 500-replicate run: [GitHub Actions run](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/39947923471)
- Output artifact: `fixed-h0-coverage-diagnostic`
- A 100-replicate exploratory run is superseded by the 500-replicate results below.

## Design

- True injected values: 55, 60, 65, 70, 75, 80, 85 km/s/Mpc.
- Scenarios: known host and unknown host.
- Replicates: 500 per true-H0 × scenario × prior-range combination.
- Prior ranges: uniform [50, 90] and uniform [30, 120] km/s/Mpc.
- Random seeds are held constant across prior supports, so each prior comparison uses the same generated mock universes.
- Nominal central interval levels: 50%, 68%, 80%, 90%, 95%.
- The H0 grid spacing is 0.025 km/s/Mpc in both prior ranges.

## Results at H0 = 70

| Scenario | Prior support | 68% coverage | MC SE | 95% coverage | MC SE | Mean posterior-median error | Mean 68% width | Mean 95% width |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Known host | [50, 90] | 76.6% | 1.89 pp | 99.0% | 0.45 pp | +1.17 | 17.19 | 30.24 |
| Known host | [30, 120] | 68.6% | 2.08 pp | 96.4% | 0.83 pp | +4.38 | 23.46 | 45.63 |
| Unknown host | [50, 90] | 94.8% | 0.99 pp | 100.0% | 0.00 pp | +1.59 | 22.40 | 35.84 |
| Unknown host | [30, 120] | 68.4% | 2.08 pp | 99.8% | 0.20 pp | +8.53 | 35.24 | 66.72 |

Errors and widths are in km/s/Mpc; MC SE is the binomial Monte Carlo standard error in percentage points. A reported 0.00 pp MC SE at 100% empirical coverage means no misses occurred in 500 trials; it does not prove true coverage is exactly 100%.

## Findings

1. **Prior support materially changes the coverage estimate.** At H0=70, unknown-host 68% coverage changes from 94.8% with [50,90] to 68.4% with [30,120]. Known-host coverage changes from 76.6% to 68.6%.
2. **Better nominal coverage does not mean the inference is automatically better.** The wider prior also makes the mean 68% interval substantially wider and raises the mean posterior-median error in this simulation. For unknown-host, mean 68% width increases from 22.40 to 35.24 and mean posterior-median error from +1.59 to +8.53 km/s/Mpc.
3. Across the H0 grid, errors and coverage vary with the injected H0, and the patterns differ between known-host and unknown-host cases. The [50,90] prior can pull posterior summaries inward when the truth is close to its endpoints; the wider prior reduces that direct truncation but exposes long tails and substantial uncertainty in some runs.
4. A 95% interval need not literally touch a prior endpoint for prior support to affect posterior summaries. Posterior normalization and the shape of tails can change even when the interval endpoints do not coincide with the bounds.
5. The wider prior is **not** automatically the scientifically correct prior. Its support must be justified and sensitivity reported. These results show that a single fixed prior can conceal strong sensitivity; they do not select a preferred prior by themselves.
6. The unknown-host likelihood is a mixture over many candidate galaxies, so its behaviour must be diagnosed separately from the known-host control. Do not attribute all differences to catalogue incompleteness without a prespecified controlled contrast.

## Model specification question to investigate

The pinned simulator draws the observed distance with a standard deviation set to a fraction of the simulated true distance. During inference, that reported standard deviation is held fixed while candidate H0 values change the predicted distance. This is a possible generative/inference mismatch: in the simulator, the noise scale depends on the true distance, while the likelihood conditions on the realized scale as fixed. We should explicitly state whether the intended experiment treats the uncertainty as a fixed reported measurement uncertainty or as a fractional-error model whose scale changes with predicted distance. A corrected likelihood must be tested as a separate model variant, not silently substituted into the pinned baseline.

## Next steps

1. Add a paired per-universe prior-support comparison, reporting the coverage-indicator difference and its uncertainty rather than only comparing marginal coverage rates.
2. Audit the distance-noise generative assumptions and implement any alternative as a separate, versioned model.
3. Add Bayesian simulation-based calibration (SBC) as a distinct diagnostic: draw H0 from the stated prior, simulate data under the same generative model used by the likelihood, and examine posterior ranks.
4. Re-run the prespecified incompleteness/weighting contrasts under both the documented baseline and any explicitly justified corrected model.

All findings here are synthetic and conditional on the stated assumptions. No real LVK H0 measurement is claimed.
