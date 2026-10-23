# Paired-seed experiment protocol

## Status

Design specification only. This document does not claim that the paired experiment has been implemented or run. The existing factorial experiment remains the baseline and should be preserved unchanged for reproducibility.

## Motivation

The current `scripts/run_factorial_sensitivity.py` assigns each scenario a separate seed block (`base_seed + scenario_index * 1_000_000`). Its prespecified contrasts therefore compare independent scenario samples. A paired design may improve precision by comparing conditions applied to the same underlying simulated universe, provided the simulator can expose and reuse the same generated population and random draws without changing the scientific meaning of the comparison.

## Primary question

For a fixed underlying mock universe, how much do catalogue selection and host-population/inference-weight mismatch change posterior error and uncertainty in the inferred Hubble constant?

## Design principles

1. Preserve the current baseline script and its pinned original simulator commit.
2. Do not call runs "paired" merely because they use the same integer seed. Pairing is valid only if the relevant underlying universe is shared intentionally and the scenario transformations are well-defined.
3. Separate random-number streams for universe generation, host selection, catalogue selection, and measurement noise where feasible. Record the stream strategy.
4. Apply each catalogue-selection mechanism to the same generated parent galaxy population for a given universe.
5. Keep the true generating value fixed at the current baseline, H0 = 70, unless a separately prespecified experiment changes it.
6. Record all scenario settings, seeds, model versions, and exclusions in machine-readable metadata.
7. Treat the experiment as synthetic; it is not an analysis of real LVK events or a measurement of the real Hubble constant.

## Proposed first paired contrasts

- Random 70% incompleteness versus complete catalogue, under uniform host generation and uniform inference weights.
- Redshift-dependent 70% incompleteness versus random 70% incompleteness, holding host-generation and inference-weighting settings fixed.
- Matched luminosity inference weights versus uniform inference weights when hosts are generated with luminosity weighting, at fixed catalogue condition.

The primary contrast list must be frozen before inspecting the results. Additional comparisons should be labeled exploratory.

## Required per-replicate outputs

For every universe and scenario, save at least:

- universe identifier and replicate index;
- base seed and any derived random-stream seeds;
- scenario name, completeness target and realized completeness;
- host-generation and inference-weighting settings;
- posterior median and mean;
- posterior standard deviation;
- 68% and 95% credible-interval bounds;
- whether each interval contains the true H0;
- convergence, numerical-validity, and exclusion flags with reasons.

A paired contrast file should contain one row per universe and contrast, with the replicate-level metric difference (scenario B minus scenario A). Never silently drop invalid runs; report counts and reasons, and prespecify the handling rule.

## Analysis

For continuous outcomes, report the mean paired difference and a confidence interval, and optionally the median paired difference as a robust companion. For coverage, report each scenario's empirical coverage and Monte Carlo uncertainty, plus the paired difference in per-universe coverage indicators with an interval suitable for paired binary data. For all metrics, show effect sizes and uncertainty; do not use a significance threshold as a substitute for scientific interpretation.

Report bias, RMSE, empirical 68% and 95% coverage, and credible-interval width together. Coverage is a repeated-sampling property, not a property established by one posterior.

## Implementation validation checklist

- Demonstrate in a small test that each paired scenario uses the same intended parent-universe identifier and parent-population data.
- Demonstrate that scenario-specific selection differs as expected while the parent population remains unchanged.
- Confirm deterministic reruns for the same configuration.
- Confirm that changing the base seed changes the simulated universes.
- Validate exact replicate counts, scenario labels, finite summaries, ordered intervals, and completeness values.
- Test paired-contrast code on hand-constructed examples with known differences.
- Check that paired and independent analyses agree in expectation in a simple null case, within Monte Carlo error.
- Archive raw per-replicate outputs, configuration, code commit, dependency versions, and analysis summaries.

## Go/no-go rule

Implement the paired design only after inspecting the pinned simulator's random-number generation and data flow. If the simulator cannot cleanly reuse a parent universe without invasive or scientifically ambiguous changes, retain the independent-scenario design and describe it honestly rather than claiming pairing.

## Reproducibility relationship to the baseline

Run and archive the current factorial workflow first. The paired experiment is a separate extension and must not overwrite the baseline outputs. Compare the paired design against the baseline only after each has passed its own validation checks.
