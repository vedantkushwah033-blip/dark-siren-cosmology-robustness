# Paired host-marginalization comparison (independent check of first audit artifact)

**Input:** artifact `matched-fixed-h0-coverage` from [run 39995536840](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39995536840).  
**Resampling:** 10,000 paired bootstrap replicates, seed 20261023.  
**Sign convention:** unknown-host minus known-host; 500 paired synthetic universes per fixed (H_0).  
**Status:** preliminary paired analysis of the first raw artifact. The workflow is being updated to generate this CSV automatically.

## Paired differences

| True H0 | Mean signed-error difference (95% bootstrap CI) | RMSE difference (95% bootstrap CI) | 68% coverage difference (95% bootstrap CI) |
|---:|---:|---:|---:|
| 55 | +7.478 [6.793, 8.151] | +8.032 [7.357, 8.718] | −28.6 pp [−33.8, −23.4] |
| 60 | +5.652 [5.016, 6.302] | +5.253 [4.729, 5.778] | −13.4 pp [−19.0, −7.8] |
| 65 | +3.744 [3.139, 4.365] | +2.685 [2.264, 3.106] | +23.0 pp [+18.4, +27.4] |
| 70 | +0.734 [0.180, 1.304] | +0.084 [−0.256, +0.432] | +22.2 pp [+17.8, +26.6] |
| 75 | −1.189 [−1.690, −0.687] | −0.979 [−1.456, −0.509] | +15.2 pp [+10.8, +19.6] |
| 80 | −2.696 [−3.192, −2.173] | −0.457 [−1.056, +0.135] | +12.4 pp [+8.4, +16.4] |
| 85 | −4.729 [−5.113, −4.329] | +2.946 [2.462, 3.432] | +1.0 pp [−3.4, +5.2] |

Here pp means percentage points. Bootstrap intervals quantify Monte Carlo uncertainty for this fixed simulator/design; they do not quantify uncertainty from misspecified physics or catalogue assumptions.

## Interpretation

- At true H0=55 and 60, unknown-host inference has markedly higher positive error and RMSE, with substantially worse 68% coverage.
- The pattern changes across H0: at 65–75 the unknown-host model has much higher 68% coverage, and its RMSE is lower at 75. At 80 the RMSE difference is uncertain under the paired bootstrap; at 85 it is higher.
- This H0-dependent pattern is not summarized well by a single pooled coverage statistic. It should be shown by truth value in any report.
- These comparisons do not yet establish that host marginalization alone causes the differences. The likelihood, catalogue distribution, prior boundaries, and conditioning on the catalogue need a line-by-line model audit before a causal claim is made.

## Validation next

The updated workflow now runs `scripts/analyze_matched_host_paired_differences.py` and writes `paired_comparison.csv`. After that run passes, test the implementation with deterministic cases: a one-galaxy catalogue, a two-galaxy catalogue with hand-computable likelihood terms, and a check that the host-generating process and the inference weights use the same distribution.
