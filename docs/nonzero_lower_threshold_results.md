# Nonzero lower-threshold selection-normalization audit

**Status:** Numerical audit completed and focused tests passed on 2026-10-24.  
**Workflow:** [Toy nonzero lower-threshold selection audit, run 40072191754](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40072191754)  
**Artifact:** `nonzero-lower-threshold-results`.  
**Scope:** deterministic synthetic two-host model; not a realistic detector selection function and not an LVK analysis.

## Question

Does the implementation apply the same finite detection interval in the selection efficiency and the observed-data validation when the lower threshold is nonzero?

## Design

- Detection interval: 20 < observed distance < 180 Mpc.
- Fixed observed distances: 30, 70, and 140 Mpc.
- Host redshifts: 0.02 and 0.04; mixture weights 0.70 and 0.30.
- Gaussian distance uncertainty: 20 Mpc.
- Uniform prior support: H0 in [50, 90] km s^-1 Mpc^-1.
- Posterior grid: 4,001 points spanning H0 = 25 to 125.
- Each observation is analyzed with the uncorrected likelihood and the likelihood divided by beta(H0).

The selection efficiency is

\[
\beta(H_0)=\sum_j w_j\left[
\Phi\left(\frac{180-cz_j/H_0}{20}\right)
-\Phi\left(\frac{20-cz_j/H_0}{20}\right)
\right].
\]

The selected likelihood is proportional to the host-mixture likelihood divided by \(\beta(H_0)\), then normalized over the stated prior support. The observed distance must be strictly inside the same interval.

## Results

Posterior medians and equal-tail intervals from the workflow output:

| Observed distance (Mpc) | Method | Median H0 | 68% interval | 95% interval |
|---:|---|---:|---:|---:|
| 30 | Uncorrected | 83.802 | [75.905, 88.314] | [66.785, 89.749] |
| 30 | Selection-normalized | 83.644 | [75.444, 88.277] | [65.957, 89.744] |
| 70 | Uncorrected | 76.957 | [65.597, 85.883] | [56.143, 89.353] |
| 70 | Selection-normalized | 75.965 | [64.196, 85.558] | [55.268, 89.301] |
| 140 | Uncorrected | 75.771 | [53.776, 85.731] | [50.480, 89.313] |
| 140 | Selection-normalized | 72.858 | [52.894, 85.124] | [50.386, 89.213] |

Across the full numerical grid H0 = 25 to 125, the workflow reported beta(H0) from 0.000971 to 0.990214. **This minimum is over the full computational grid, not just the prior support [50, 90].** Values outside the prior are not used in the posterior normalization. The low minimum still warns that the selection efficiency can become small in portions of the wider computational grid; follow-up sensitivity analysis should report the efficiency range restricted to the prior support.

## Interpretation

1. The nonzero lower-threshold path executes and returns finite posterior summaries for all three observations.
2. The correction shifts the posterior median modestly at 30 Mpc and more noticeably at 70 and 140 Mpc in this chosen toy setup. These are sensitivity demonstrations, not evidence that the corrected estimate is closer to truth: the observations here are fixed illustrative values, not repeated draws with a known true parameter.
3. Credible intervals remain broad and are influenced by the finite prior [50, 90]. In particular, the 140 Mpc example has a 95% interval close to the lower prior boundary.
4. The code and test enforce a shared lower/upper interval for beta(H0) and observed-data validity. This is a consistency check of the implementation, not validation of the assumed selection model.

## Validation status

- Workflow run 40072191754 completed successfully.
- The focused tests passed after fixing a test that had attempted to cast the string-valued `method` field to a number.
- The full test-suite workflows for the same commit also completed successfully (runs 40072191711 and 40072191718).

## Limitations and next checks

This model omits a real gravitational-wave detection statistic, network sensitivity, duty cycle, injection campaign, source-population distribution, cosmological volume factors, and galaxy-catalogue incompleteness. The finite interval in noisy distance is only a controlled testbed.

The next sensitivity experiment varies the lower threshold while holding the upper threshold, host mixture, noise, prior, and observed distances fixed. It reports the efficiency range on the prior support separately from the broader computational grid. A separate repeated-sampling study is still needed to diagnose fixed-truth interval coverage; this deterministic audit cannot answer that question.
