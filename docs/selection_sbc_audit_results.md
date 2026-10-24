# Selection-conditioned simulation-based calibration (SBC)

**Workflow:** [Toy selection-normalized SBC audit, run 40070612148](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40070612148)  
**Status:** workflow succeeded; output artifact retrieved and inspected.  
**Artifact:** `selection-sbc-audit-results`; seed 20261027.

## Design

- 2,000 prior-predictive replicates.
- Draw (H_0) uniformly from the analysis prior [50, 90] km s^-1 Mpc^-1.
- Conditional on each fixed draw of (H_0), repeatedly draw host and noisy distance until the event passes the toy detection rule (0<d_{m obs}<180) Mpc.
- Infer with the selection-conditioned likelihood, using the same host mixture, noise scale, and detection interval.
- Grid: 4,001 points spanning 25 to 125; prior support [50, 90].
- Compute the posterior CDF evaluated at the generating (H_0). Under a correctly implemented and calibrated continuous posterior, these ranks should be approximately uniform, subject to grid and Monte Carlo error.

A key design detail is that the parameter draw is **not** rejected when an event fails detection. We condition the data on detection at a fixed (H_0), which keeps the chosen analysis prior as the prior for the SBC experiment. This is distinct from resampling (H_0) itself until the event is detected, which would change the distribution of parameter values among accepted simulations.

## Results

| Diagnostic | Result |
|---|---:|
| Replicates | 2,000 |
| Mean posterior-CDF rank | 0.50092 |
| Rank variance | 0.08390 |
| KS statistic against Uniform(0,1) | 0.01115 |
| KS p-value | 0.96237 |
| Decile chi-square statistic | 3.39 |
| Decile chi-square p-value | 0.94681 |
| Minimum (eta(H_0)) on prior grid | 0.69948 |
| Maximum (eta(H_0)) on prior grid | 0.99679 |

Decile counts were 208, 191, 187, 213, 208, 196, 196, 193, 201, and 207.

The rank mean and variance are close to the uniform reference values 0.5 and 1/12 ≈ 0.08333. Neither the KS test nor the decile chi-square test detects a departure from uniform ranks in this run.

## Interpretation and limits

This is a useful positive calibration check for the *toy selection-conditioned model*. It does not prove the model is physically realistic, nor does a large p-value prove exact calibration. The test is conditional on the implementation, prior, host mixture, Gaussian noise model, selection interval, and finite posterior grid used here.

It is also not inconsistent with fixed-truth repeated-sampling coverage differing from nominal 68% at particular values of (H_0). Prior-predictive SBC averages over the prior distribution; fixed-truth coverage is a different diagnostic and can vary with parameter location, particularly near prior boundaries. Both should be reported.

The separate prior-sensitivity experiment shows that correction effects vary with true (H_0) and prior support. The right conclusion is therefore not that the method is universally superior, but that the selection-normalized toy implementation passes this particular prior-predictive rank check while fixed-truth performance remains heterogeneous.

## Next validation

1. Add a deliberate negative-control SBC branch that omits (eta(H_0)), using the exact same generated observations and seeds.
2. Compare calibrated and deliberately misspecified ranks and fixed-truth coverage side by side.
3. Repeat the calibration check with grid refinement and at least one nonzero lower detection threshold.
4. Keep catalogue incompleteness separate from event detection; this experiment only addresses the latter.

This remains a synthetic method-validation exercise, not an analysis of real LVK observations.
