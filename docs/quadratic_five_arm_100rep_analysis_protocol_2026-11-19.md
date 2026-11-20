# Analysis protocol for the 100-replicate five-arm 30-event runs

**Frozen before inspecting the 100-replicate output.** This is a precision/robustness follow-up to the 10-replicate smoke, not a confirmatory observational analysis.

## Run identities

Two Actions runs were launched with the same scientific configuration but disjoint master seeds:

- Run [44460222957](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44460222957): master seed 2026111929, 100 replicates. Its first 10 replicate seed pairs coincide with the earlier 10-replicate smoke and must not be double-counted if the smoke is pooled with this run.
- Run [44462487127](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44462487127): master seed 2026111930, 100 replicates. This is a disjoint independent seed stream.

Both use gamma = -1, 0, +1; validation sizes 100, 300, 500; 1,200 parent galaxies; 30 detected events; H0 truth 70; five arms A-E; 281 H0 points; Scott multiplier 0.7; and 25 coefficient-quadrature nodes. The second run was launched as a separate batch to add an independent seed stream, not to replace the first.

## Analysis units and dependence

The independent Monte Carlo unit is the parent/event replicate, not an individual arm-cell row. Within each replicate:
- all five arms share the same generated parent population and detected-event realization;
- gamma values are paired using common parent/event seeds;
- validation sizes 100/300/500 are nested samples from the same validation frame;
- the five inference arms are matched analyses of the same synthetic realization.

Therefore, all comparisons below are paired at replicate level. Do not treat the 45 rows per replicate as independent observations.

## Primary descriptive estimands

For every gamma × validation-size × arm cell, report:
1. mean posterior-median error (bias);
2. RMSE;
3. empirical 68% equal-tail interval coverage and Wilson 95% interval;
4. mean 68% interval width;
5. failure-inclusive RMSE and coverage bounds;
6. number of successful, failed, non-evaluable, and not-run planned replicates.

For the fitted methods, additionally report validation-fit and weight diagnostics at replicate level: optimizer status, inclusion-probability range, Kish ESS, maximum normalized weight, top-1% and top-5% weight mass, support gaps, and quadrature-node ESS / minimum node-wise Kish ESS.

## Paired contrasts

For each gamma and validation size, preserve replicate-level paired differences for:
- E minus D (effect of integrating inclusion-coefficient uncertainty versus MAP plug-in);
- E minus C (integrated estimated inclusion versus known-probability IPW);
- D minus C (MAP estimated inclusion versus known-probability IPW);
- B minus C (oracle reference versus known-probability IPW).

Report mean paired differences and 95% t intervals for:
- posterior-median error;
- absolute posterior-median error;
- squared error;
- 68% coverage indicator;
- 68% interval width.

The paired squared-error difference is the direct RMSE-comparison estimand; do not compare RMSEs as if independent. Treat the nine gamma-by-validation-size cells as a structured family, not nine unrelated experiments. Any formal p-values must be multiplicity-adjusted and labeled exploratory.

## Failure handling

The planned denominator is fixed at 100 replicates per gamma × validation-size × arm. Failed, non-evaluable, and not-run cells remain in the denominator. Report success-only metrics separately from failure-inclusive bounds. Under the frozen [50,90] H0 support and truth 70, the maximum absolute posterior-median error bound is 20. Do not silently discard a child whose independent audit fails.

## Precision and interpretation limits

With 100 replicates, a nominal 68% coverage estimate has binomial standard error about 4.7 percentage points (roughly ±9 percentage points at 95% before accounting for other uncertainty). Thus this run can reveal gross miscalibration and broad effects, but it is not a high-precision coverage certification. If both master-seed batches pass, their disjoint replicate streams may be combined after confirming identical configurations and disjoint seed pairs. Do not pool the earlier 10-replicate pilot separately with run 44460222957 because those seeds are already included in that 100-replicate batch.

The fixed Scott multiplier 0.7 and 281-point H0 grid remain engineering choices. No result from these runs alone establishes bandwidth optimality, grid convergence, universal estimator superiority, or observational validity.
