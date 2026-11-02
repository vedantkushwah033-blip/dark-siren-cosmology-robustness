# Multi-frame v3 manifest and frame generator

**Status:** implemented data-generation foundation; not a benchmark runner or confirmatory experiment.

Files:
- `src/dark_siren_robustness/multiframe_design.py`
- `tests/test_multiframe_design.py`

The manifest validator requires explicit scenario identifiers, distribution family and parameters, population size, frame count, sample sizes, outer repetitions, bootstrap repetitions, methods, and nominal interval level. It rejects duplicate identifiers, invalid integer/real values, unsupported distributions, and sample sizes that are not smaller than the finite population.

The frame generator uses the v3 domain-separated `frame_seed` helper. Each scenario × frame coordinate deterministically produces an immutable finite vector and metadata including its seed, exact realized arithmetic mean, population variance, and frame identifier. `generate_frame_grid` creates the balanced hierarchy in explicit manifest order.

The implementation currently supports Normal, Gamma, and Lognormal finite-frame generating families. These are controlled synthetic frames, not empirical populations.

## Deliberate scope boundary

This implementation does not yet generate outer SRSWOR samples, calculate intervals, run bootstrap attempts, persist artifacts, or independently audit results. In particular, no v3 method comparison or coverage estimate is available yet. Those require the next runner and auditor layers.
