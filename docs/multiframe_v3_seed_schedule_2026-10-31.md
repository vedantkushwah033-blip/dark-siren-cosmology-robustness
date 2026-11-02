# Multi-frame v3 seed schedule

**Status:** implemented deterministic utility; not a complete multi-frame experiment.

The v3 seed helper is `src/dark_siren_robustness/multiframe_seed_schedule.py`. It uses domain-separated `numpy.random.SeedSequence` inputs to derive seeds from the complete experimental coordinates.

| Stream | Coordinates |
|---|---|
| Frame generation | root seed, scenario ID, frame ID |
| Outer SRSWOR sample | root seed, scenario ID, frame ID, sample size, replicate ID |
| Method-specific random stream | all outer-sample coordinates plus stable method ID |
| Inner attempt | all method coordinates plus attempt ID |

The schedule is independent of loop order and method execution order. The outer sample seed is deliberately method-agnostic so every method can share the same sampled units for paired comparisons. Method-specific bootstrap randomness uses a separate domain.

The implementation rejects booleans, negative coordinates, and values above the supported 32-bit coordinate range. It records a schema identifier and domain constants through `seed_schedule_metadata()`.

Tests cover deterministic replay, coordinate sensitivity, domain separation, method-order invariance, inner-attempt uniqueness, invalid inputs, and SRSWOR reproducibility.

This module is only a building block. It does not yet implement a v3 manifest, frame generator, runner, artifact schema, auditor, or confirmatory design. The existing v2 benchmark and its historical outputs are unchanged.
