# Multi-frame v3 paired outer-sampling layer

**Status:** implemented deterministic SRSWOR sampling primitive; not an interval runner or confirmatory result.

Files:
- `src/dark_siren_robustness/multiframe_sampling.py`
- `tests/test_multiframe_sampling.py`

For every declared scenario, frame, sample size, and outer replicate, the sampler regenerates the finite frame from the manifest and selects units without replacement using the method-agnostic outer-sample seed. Consequently, every method applied to a given key receives the same immutable sampled values and indices. The returned record includes sample and frame seeds, hashes, sample mean and variance, and the exact finite-frame mean.

The grid generator emits the complete declared hierarchy in canonical manifest order. It rejects undeclared sample sizes, out-of-range frame or replicate identifiers, unknown scenarios, and boolean/non-integer coordinates.

This is a sampling primitive, not a complete simulation. It does not yet construct Wald, percentile, or studentized intervals, persist a dataset, or independently audit an artifact. The sample hash supports future artifact integrity checks but is not a substitute for independent reconstruction.
