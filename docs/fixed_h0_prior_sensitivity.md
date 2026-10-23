# Fixed-H0 coverage and prior-range sensitivity

## Purpose

This diagnostic asks whether the fixed-H0 coverage pattern changes when the posterior's uniform prior support is widened. It uses the original synthetic simulator pinned to `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. It is a sensitivity check, not real-event cosmology.

Workflow: [Fixed-H0 coverage diagnostic](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39947604821)

## Design

- True injected values: 55, 60, 65, 70, 75, 80, 85 km/s/Mpc.
- Scenarios: known host and unknown host.
- Replicates: 100 per true-H0 × scenario × prior-range combination.
- Prior ranges: uniform [50, 90] and uniform [30, 120] km/s/Mpc.
- The same generated mock universes are used for both prior-range comparisons at each H0/scenario combination.
- Nominal central interval levels: 50%, 68%, 80%, 90%, 95%.

## Initial results at H0 = 70

| Scenario | Prior support | 68% coverage | MC SE | 95% coverage | MC SE | Mean posterior-median error |
|---|---|---:|---:|---:|---:|---:|
| Known host | [50, 90] | 77% | 4.21 pp | 98% | 1.40 pp | +1.92 |
| Known host | [30, 120] | 64% | 4.80 pp | 96% | 1.96 pp | +5.59 |
| Unknown host | [50, 90] | 95% | 2.18 pp | 100% | 0.00 pp | +1.74 |
| Unknown host | [30, 120] | 70% | 4.58 pp | 100% | 0.00 pp | +8.36 |

Errors are in km/s/Mpc; MC SE is the binomial Monte Carlo standard error in percentage points.

These estimates are exploratory and noisy at n=100. A coverage estimate of 64% has MC SE about 4.8 percentage points; therefore this single value should not be treated as proof that 68% coverage is wrong. The unknown-host difference (95% versus 70%) is a strong warning that interval behaviour depends materially on prior support in this synthetic setup and deserves a larger paired comparison.

## Interpretation

1. The previous 1,000-replicate baseline used a uniform prior over [50, 90]. The fixed-H0 grid shows that coverage and bias vary with the true H0 and the prior support.
2. A 95% interval need not literally touch a prior endpoint for the prior to affect posterior summaries. The likelihood shape, truncation and normalization can still change posterior quantiles.
3. Widening the prior substantially increases interval widths, especially in the unknown-host mixture, and changes its 68% coverage at H0=70 from 95% to 70% in this exploratory run.
4. This is not evidence that the wider prior is automatically correct. Prior support must be scientifically justified and sensitivity reported. The next run should increase replicates and report uncertainty in coverage differences.
5. The unknown-host mixture uses a catalogue mixture likelihood with many candidate galaxies. Its behaviour should be investigated separately from the known-host control; do not attribute the full discrepancy to catalogue incompleteness without a controlled contrast.
6. The fixed-H0 repeated-sampling experiment is not the same as Bayesian simulation-based calibration (SBC). SBC should be implemented as a separate diagnostic that draws H0 from the prior, simulates data conditionally, and evaluates posterior ranks under a consistent generative model.

## Next step

Increase the prior-range comparison to at least 500 replicates per cell, compute paired per-universe coverage differences between prior supports, and report Monte Carlo uncertainty. Then add an independent SBC rank diagnostic. Keep all results labelled synthetic and conditional on the chosen model.
