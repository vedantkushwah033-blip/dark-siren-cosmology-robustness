# Single-host fixed-truth selection control: protocol

**Status:** protocol only as of 2026-10-24. No results are claimed by this document.

## Motivation

The existing single-host prior-predictive SBC experiment tests posterior ranks when the true parameter is drawn from the prior. The existing fixed-truth coverage maps use a two-host mixture. A fixed-truth single-host control helps determine whether edge undercoverage is already present with one host, before attributing it to mixture complexity.

## Prespecified design

- Host redshift: z = 0.03.
- Distance model: d_obs | H0 ~ Normal(c z / H0, sigma^2), sigma = 20 Mpc.
- Detection rule: 0 < d_obs < 180 Mpc.
- Analysis prior: uniform H0 in [50, 90] km s^-1 Mpc^-1.
- Fixed truths: 52, 55, 60, 65, 70, 75, 80, 85, 88.
- Replicates per truth: 300.
- Random seed: 20261030.
- Posterior grid: 4,001 points over [25, 125].
- Methods: uncorrected likelihood and likelihood divided by the analytically calculated detection efficiency beta(H0).
- Pairing: each detected observation is analysed by both methods.
- Primary metrics: bias and RMSE of posterior median; empirical equal-tailed 68% and 95% interval coverage with Wilson 95% intervals; mean interval widths; paired corrected-minus-uncorrected absolute-error differences with percentile-bootstrap 95% intervals.
- Expected per-run rows: 2,700 (9 truths × 300 observations); expected method summaries: 18 (9 truths × 2 methods).

The generator samples H0 as a fixed design point, draws the host (only one possible host), generates noisy distance, and rejection-samples until the stated detection condition is satisfied. It must not resample or discard the truth itself. The corrected posterior uses the same finite detection interval as the generator.

## Prespecified interpretation

1. This is a diagnostic control, not a test of real detector physics.
2. The correction is not presumed to win on fixed-truth RMSE or coverage. The aim is to quantify where it changes behavior.
3. Compare the single-host results with the existing two-host coverage-transition map, but do not pool the samples across experiments or treat different seeds/designs as one sample.
4. If edge undercoverage persists in the one-host control, it is not solely caused by the two-host mixture. It may still reflect the bounded prior, selection-conditioned sampling, likelihood geometry, finite-sample variation, or other assumptions.
5. A nominal 68% Bayesian credible interval does not guarantee 68% pointwise frequentist coverage. Interpret empirical coverage alongside Monte Carlo uncertainty and interval width.
6. No claims about catalogue incompleteness, realistic GW selection, LVK observations, or publication readiness follow from this control alone.

## Implementation acceptance criteria

- Validate the finite positive inputs and enforce that every generated observation lies strictly within (0, 180) Mpc.
- Verify paired methods use exactly the same observations.
- Verify the selected posterior is normalized on the prior support for every test threshold.
- Test deterministic reproducibility at a small replicate count and full-design output row counts.
- Report finite metrics and Wilson intervals for coverage.
- Store the seed, grid size, design, and package/runtime versions in a manifest.
- Run focused tests and the full project test suite before interpreting the artifact.
