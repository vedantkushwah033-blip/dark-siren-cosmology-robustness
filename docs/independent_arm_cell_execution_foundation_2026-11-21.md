# Independent arm-cell execution foundation — 2026-11-21

## Purpose

The next execution unit for the synthetic low-redshift Euclidean dark-siren benchmark is one pair:

`(validation_size, arm)`

The current multi-replicate runner still executes one whole child per `(replicate, gamma)`, calculating all validation sizes and all five arms together. This module does **not** claim that the orchestrator has yet been refactored. It provides the numerical primitive required for that refactor.

## API

`dark_siren_robustness.arm_cell_execution.evaluate_arm_cell` accepts one arm and the shared immutable inputs needed by that arm. It returns an `ArmCellResult` containing:

- the arm's raw/component weights and explicit mixture weights;
- component KDE densities;
- component selected log likelihoods and detection efficiencies;
- the marginal selected log likelihood;
- posterior quantiles on the declared H0 support;
- weight-concentration diagnostics and realized KDE bandwidths;
- inclusion-posterior mode/nodes/weights only for D/E.

## Dependency isolation

| Arm | Host sample | Needs validation labels? | Needs true inclusion probabilities? | Fits logistic model? |
|---|---|---:|---:|---:|
| A catalogue-only | observed catalogue | No | No | No |
| B parent-oracle | full parent | No | No | No |
| C known-probability IPW | observed catalogue | No | Yes (synthetic truth) | No |
| D MAP plug-in IPW | observed catalogue | Yes | No | Yes |
| E posterior-integrated IPW | observed catalogue | Yes | No | Yes |

A/B/C must remain evaluable if the inclusion-model fit fails. A uses only catalogue inputs; B uses only parent inputs. D and E fit the validation posterior independently at the arm-cell boundary. E marginalizes the complete selected multi-event likelihood over coefficient-posterior nodes, rather than averaging weights or likelihoods event-by-event.

No pure-IPW arm is combined with a missing-host mixture correction.

## Current tests

The focused test module compares all five independently executed arm likelihoods against the established matched five-arm implementation to tight numerical tolerance. It also verifies:
- A/B/C do not call the validation inclusion fitter;
- A/B/C do not require validation data;
- A/B do not require synthetic true inclusion probabilities;
- A and B do not validate the other arm's host/base-weight inputs;
- C does not require validation labels;
- one arm's invalid input does not change another arm's deterministic result;
- D/E still require validation labels.

## Not yet delivered

This PR is not an end-to-end cell runner. Remaining gates:

1. Persist one immutable shared-input bundle per `(replicate, gamma)` before any arm is evaluated.
2. Create one checkpoint cell for every `(replicate, gamma, validation_size, arm)`.
3. Mark and execute each cell independently; atomically persist only that cell's terminal outcome.
4. Define arm-specific `NON_EVALUABLE` reason-code policy before using that status.
5. Write per-cell result artifacts and content hashes.
6. Build an independent auditor that reconstructs each arm's weights, KDE, selected likelihood, posterior quantiles, and lifecycle status without importing the execution function.
7. Prove restart equivalence, hard-kill recovery, and failure isolation.
8. Benchmark checkpoint persistence at the intended experiment scale.

No confirmatory prior-predictive calibration should start before these gates and the frozen numerical-convergence protocol pass. The research scope remains a synthetic sensitivity benchmark, not an observational H0 estimate, real survey validation, or LVK/GWTC data analysis.
