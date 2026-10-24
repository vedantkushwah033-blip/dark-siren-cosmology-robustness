# Selection-conditioned simulation-based calibration (SBC)

**Baseline workflow:** [run 40070612148](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40070612148)  
**Paired negative-control workflow:** [run 40071738592](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40071738592)  
**Artifact:** `selection-sbc-audit-results`; seed 20261027. The negative-control run succeeded and its CSV artifacts were downloaded and inspected.

## Design

- 2,000 prior-predictive replicates.
- Draw (H_0) uniformly from the analysis prior [50, 90] km s⁻¹ Mpc⁻¹.
- At each fixed draw of (H_0), draw host and noisy distance jointly conditional on the toy detection rule (0<d_{obs}<180) Mpc.
- Analyse each identical observation twice: with the selection-conditioned likelihood and with the uncorrected likelihood.
- Grid: 4,001 points spanning 25 to 125; prior support [50, 90].
- Compute the posterior CDF evaluated at the generating (H_0). Under a correctly specified continuous posterior, ranks should be approximately uniform, subject to numerical and Monte Carlo error.

The parameter draw is not rejected when an event fails detection. Instead, the data are sampled conditional on detection at that fixed (H_0), preserving the specified prior for this SBC experiment. This differs from repeatedly drawing (H_0) until a detected event occurs, which would alter the distribution of accepted parameter values.

## Paired results

| Diagnostic | Selection-normalized | Uncorrected control |
|---|---:|---:|
| Replicates | 2,000 | 2,000 |
| Mean posterior-CDF rank | 0.50092 | 0.47034 |
| Rank variance | 0.08390 | 0.08480 |
| KS statistic vs Uniform(0,1) | 0.01115 | 0.05405 |
| KS p-value | 0.96237 | 0.0000161 |
| Decile chi-square statistic | 3.39 | 26.31 |
| Decile chi-square p-value | 0.94681 | 0.00182 |

Selection-normalized decile counts were 208, 191, 187, 213, 208, 196, 196, 193, 201, and 207. Uncorrected counts were 242, 218, 227, 213, 183, 205, 176, 175, 180, and 181.

For this particular toy setup, the selection-normalized ranks are consistent with a uniform distribution, while the uncorrected control shows a detectable departure under both the KS and decile tests. The uncorrected ranks are shifted downward on average (0.470 versus 0.501), consistent with a systematic rank distortion when the detection normalization is omitted.

These p-values are diagnostics, not posterior probabilities that either method is correct. The paired comparison is stronger than the earlier baseline-only check because the methods receive the same generated observations.

## Interpretation and limits

This is evidence that the selection normalization matters for prior-predictive calibration under this specific synthetic data-generating process. It does not prove exact calibration or physical realism. The experiment depends on the host mixture, Gaussian distance noise, selection interval, prior, and finite posterior grid.

It also does not resolve the fixed-truth coverage problems. Prior-predictive SBC averages over the prior; fixed-truth coverage evaluates performance at particular parameter values and may vary substantially, especially near prior boundaries. Both diagnostics should remain in the report.

The selection normalization here addresses **event detection**, not galaxy-catalogue incompleteness. Those mechanisms must be modelled separately. This is not an analysis of real LVK observations.

## Next validation

1. Confirm the updated test suite passes after fixing the negative-control test's undefined variable.
2. Repeat the paired SBC at a finer posterior grid and compare rank changes.
3. Check at least one nonzero lower detection threshold.
4. Continue investigating fixed-truth coverage and prior-boundary effects before expanding the simulation matrix.
