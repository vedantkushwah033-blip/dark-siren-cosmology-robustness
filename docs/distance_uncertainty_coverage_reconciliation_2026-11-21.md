# Reconciliation of fixed-truth distance-uncertainty coverage (2026-11-21)

## Scope and immutable evidence

This note reconciles two coverage summaries from the auxiliary synthetic, one-known-host, low-redshift Gaussian distance-noise experiment. It is not the five-arm dark-siren benchmark and is not an observational or LVK analysis.

- Workflow: `44872071636`, Python 3.12 job.
- Code revision recorded in the artifact: `d252e0483dbf19379b924693ff280fcb51bf25e0`.
- Full run: 300 replicates per truth × fractional-error × likelihood-arm cell; 9 fixed truths (50, 55, ..., 90), 4 fractional errors, 3 arms.
- Raw rows: 32,400; summary cells: 108.
- Independent raw-record audit: PASS; 10,800 paired A/B replicate inputs verified.
- Raw CSV SHA-256: `62af121350d62b141f0f93917258e1967cd51b611234e6c4c9839f3088e2a4da`.
- Summary JSON SHA-256: `03000d3fbf00eb850016b4522b71c33f9a35a7ba3f4eba2b9c15888aea69a726`.

## The denominator correction

The earlier 57.0%–57.3% figure was calculated as 6,156–6,184 covered records divided by all 10,800 records per arm. That denominator includes the two exact prior-boundary truths, H0=50 and H0=90. For a continuous posterior on [50,90], the central equal-tail credible interval endpoints lie strictly inside the support, so these two exact boundary truths have zero coverage by construction.

That 57% figure is therefore the equally weighted average over the *nine-point fixed-truth grid*, which assigns 2/9 of its weight to two measure-zero endpoints. It must not be described as coverage among interior truths or as prior-predictive coverage under a continuous Uniform(50,90) truth distribution.

The independent fixed-truth interval audit explicitly separates the seven interior grid truths H0={55,60,65,70,75,80,85}. The corrected results are:

| Likelihood arm | 68% interval coverage | 95% interval coverage | Interior records |
|---|---:|---:|---:|
| A: fixed realized sigma | 73.29% (6,156/8,400) | 97.12% (8,158/8,400) | 8,400 |
| B: candidate-dependent fractional sigma(H0) | 73.62% (6,184/8,400) | 97.27% (8,171/8,400) | 8,400 |
| C: fixed reference sigma | 73.37% (6,163/8,400) | 97.00% (8,148/8,400) | 8,400 |

These are equal-weight aggregates over seven *fixed* truth values and four fractional-error settings. They are descriptive fixed-truth coverage, not a prior-predictive calibration estimate. Coverage is heterogeneous across truths: for Arm A, 68% coverage ranges from 66.67% at H0=55 to 79.08% at H0=80 (with intermediate truth-specific values retained in `fixed_truth_coverage_reference_interval_audit.json`). Report truth-stratified results rather than relying on one pooled number.

## Independent reference-interval comparison

For Arm A only, a separate audit inverted the central Gaussian distance interval using the recorded realized sigma. This is an oracle-conditional reference interval, not an analyst-available prescription when H0 is unknown.

Among the same 8,400 interior records, its coverage was:
- 68% reference interval: 67.62% (5,680/8,400), Wilson 95% interval [66.61%, 68.61%].
- 95% reference interval: 95.12% (7,990/8,400), Wilson 95% interval [94.64%, 95.56%].

This comparison shows that the Bayesian equal-tail intervals and the inverted Gaussian reference intervals have different fixed-truth coverage behavior under this parameterization and truth grid. It does not establish that either is universally preferable.

## Numerical stability

The independent 3,201-point reconstruction found:
- Maximum absolute posterior-quantile change across all 32,400 records: 0.0001764 H0.
- Maximum absolute changes by quantile: lower 68% 0.0001122; upper 68% 0.0001562; lower 95% 0.0001244; upper 95% 0.0001764; median 0.0001239.
- The 68% coverage indicator was unchanged for every arm under this grid refinement.

Thus the reported coverage pattern is not explained by the tested 1,601-to-3,201 grid refinement. This does not rule out other model or design effects.

## Interpretation and next experiment

1. Correct the prior statement: 57% is the nine-point all-truth-grid average including exact prior endpoints; the seven-point interior-grid averages are approximately 73.3%–73.6%.
2. Do not call the 73%–97% interior-grid aggregates nominal calibration. The seven truths are fixed and equally weighted, not draws from the continuous prior; the 55 and 85 cells also remain close enough to prior boundaries for truncation effects to matter.
3. The experiment's likelihoods are explicit and independently reconstructed. Arm A uses a truth-dependent realized sigma (oracle-conditional); Arm B uses candidate-dependent sigma with the Gaussian normalization; Arm C changes both the generator's noise scale and the likelihood's fixed reference scale.
4. A proper prior-predictive coverage study should draw continuous H0 truths from the stated prior using a dedicated stream and report the truth distribution, all attempted replicates, and boundary/support diagnostics. The main five-arm Stage 0 runner already adopts a per-replicate continuous Uniform prior draw, but its lifecycle ledger gate remains blocked.

No change to the dark-siren benchmark's scientific claims follows from this auxiliary known-host diagnostic.
