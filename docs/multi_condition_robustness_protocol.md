# Multi-condition robustness protocol

**Status:** prospective design; no results from this protocol are claimed yet.  
**Scope:** synthetic methodological study only. This is not a real LVK analysis or an empirical H0 measurement.

## Objective

Test whether the Stage B result is stable beyond the single fixed-truth scenario. The existing Stage B primary endpoint and output must remain immutable; this is a separate extension with its own configuration and seeds.

## Design principles

- Treat a replicate—not an event, galaxy, or grid point—as the independent sampling unit.
- Use paired seeds across the six arms within each scenario.
- Preserve truth-only inputs strictly for the analytic-truth/oracle arms and scoring; non-oracle fits must use only observable catalogue/event data and declared analyst assumptions.
- Report all prespecified scenarios, including failed runs and scenarios that reverse the ranking.
- Do not pool scenario-specific estimates into a single headline effect without a prespecified weighting scheme.
- Keep numerical verification separate from statistical conclusions.

## Minimum scenario grid

Use a staged design to keep runtime manageable.

### Stage R1: broad screening

For each scenario, use 100 paired replicates and a unique, non-overlapping seed range. Hold the baseline redshift domain and six arms fixed.

| Factor | Prespecified levels |
|---|---|
| True H0 (km s^-1 Mpc^-1) | 60, 70, 80 |
| Detected events | 15, 30, 60 |
| Distance-noise scale (Mpc) | 100, 140, 180 |
| Detection threshold (Mpc) | 200, 250, 300 |

The full factorial is 81 scenarios and is likely too expensive as a first pass. Use a transparent staged subset: (i) vary H0 alone around baseline; (ii) vary event count alone; (iii) vary noise and threshold in paired low/high stress cases. Record the exact scenario list before launching. Do not select scenarios based on preliminary outcomes.

### Stage R2: catalogue-selection stress test

The generator, runner, CLI, and independent auditor have now been extended to accept and record the true inclusion curve (floor, midpoint, width), ensure the correct-weighting arm uses the exact truth parameters, and reject inconsistent truth-parameter records during audit. The new `catalogue-selection-stress-screen.yml` workflow prespecifies three 50-replicate screening conditions with non-overlapping seed ranges. These runs are not considered complete until GitHub Actions artifacts and independent audits are available.

Use a predeclared grid such as: floor 0.05, 0.10, 0.25; midpoint 0.18, 0.22, 0.26; width 0.03, 0.05, 0.08. Start with 100 paired replicates per condition. After validating the generalization and auditor, expand to 300–500 replicates for a small number of predeclared representative stress conditions.

## Primary and secondary outcomes

For each scenario, retain the primary paired contrast:

D = abs(mean_H0(catalogue_only)-H0_true) - abs(mean_H0(correct_weighting)-H0_true).

Positive values favor correct weighting on absolute posterior-mean error.

Report mean D, paired SD/SE, percentile-bootstrap 95% interval; per-arm bias and RMSE; nominal 68% and 95% interval coverage with Wilson intervals; mean interval width; posterior-grid normalization and independently reconstructed posterior summaries; effective sample size and maximum normalized weight; runtime, package versions, seed range and manifest checksums.

The primary contrast is scenario-specific. Secondary contrasts (correct weighting versus analytic truth, oracle versus analytic truth, misspecified versus correct weighting) must be labelled secondary and interpreted as diagnostics, not as a leaderboard.

## Numerical acceptance gates

1. Repeat a fixed subset at H0-grid sizes 81, 161, and 321 and redshift-grid sizes 401, 801, and 1601.
2. Compare posterior means and 68% interval endpoints against tolerances locked before results are inspected.
3. Independently reconstruct posterior summaries from saved posterior grids.
4. Verify every scenario has all six arms exactly once per replicate, common seeds across arms, and no seed overlap with the original Stage A/B campaigns.
5. Validate true versus assumed inclusion parameters in the raw records and manifest.
6. Fail closed on missing rows, duplicate keys, non-finite values, invalid probabilities, or checksum mismatch.

## Stop/go interpretation

A scenario where correct weighting does not improve the chosen endpoint is a result to report, not a reason to discard the scenario. If rankings change across conditions, the conclusion must be conditional and describe the effect modifiers. If the analytic-truth arm remains biased or calibration fails materially, investigate the likelihood and selection normalization before making estimator-performance claims.

Passing this protocol would support a broader synthetic robustness statement only for the tested design. It would not establish realistic detector selection, event-rate inference, host-population inference, or a real-world H0 measurement.
