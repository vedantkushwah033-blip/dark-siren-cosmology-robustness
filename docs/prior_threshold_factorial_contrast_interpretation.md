# Paired-contrast interpretation audit for the compact factorial run

## Why this note exists

The compact factorial artifact contains paired comparisons because each simulated observation is reused across prior, threshold, and method configurations. This improves precision for comparisons, but it creates dependence between rows. Treating all 7,200 rows as independent would be pseudoreplication.

The artifact is a 100-replicate-per-truth CI smoke test, not the full protocol run. All values below are conditional on the toy model and this particular random sample.

## A useful implementation invariant

For the uncorrected method, the assumed detection threshold does not enter the likelihood. Holding truth, replicate, prior, and method fixed, posterior medians and interval endpoints are exactly identical across the four threshold labels. The audited artifact confirms this invariant (maximum absolute median difference across threshold labels: 0).

This is expected behavior, not evidence that the threshold factor has no scientific relevance. It means the uncorrected method is deliberately insensitive to the supplied threshold parameter; only the selection-corrected method uses it.

## Selected paired contrasts: corrected method, baseline prior [50,90]

The table gives comparison minus baseline (assumed threshold 120 Mpc). The bootstrap intervals are 95% percentile intervals over the 100 paired repetitions. A negative coverage difference means fewer intervals contained the truth in the comparison configuration.

| Truth | Assumed threshold | Median difference | 68% coverage difference | 95% coverage difference |
|---:|---:|---:|---:|---:|
| 52 | 150 | +1.079 | −0.05 | −0.04 |
| 52 | 180 | +0.529 | −0.05 | −0.04 |
| 52 | 210 | +0.926 | −0.05 | −0.04 |
| 70 | 150 | +0.307 | 0.00 | 0.00 |
| 70 | 180 | +0.109 | 0.00 | 0.00 |
| 70 | 210 | +0.786 | 0.00 | 0.00 |
| 88 | 150 | +0.067 | −0.01 | +0.02 |
| 88 | 180 | +0.124 | +0.01 | +0.03 |
| 88 | 210 | +0.660 | +0.01 | +0.06 |

For true (H_0=52), the compact run's corrected 68% coverage falls from 0.05 at the correctly specified threshold to 0.00 under all three larger assumed thresholds. This is a deterioration in an already badly calibrated baseline configuration. At true (H_0=70), empirical coverage happens to remain unchanged in this sample, even though posterior medians and interval widths shift. At true (H_0=88), the direction of coverage changes is small and not uniform across metrics.

Do not turn this table into a universal ranking of thresholds. It is a finite-sample, truth-specific sensitivity diagnostic.

## Why some bootstrap intervals are exactly [0,0]

For the central truth, the paired coverage indicators did not change between the baseline and the compared threshold for the selected configurations. Every paired difference was therefore zero, making the percentile-bootstrap interval exactly [0,0].

This does **not** mean the true coverage difference is known to be exactly zero. It means the empirical paired sample contains no variation in that contrast for those 100 repetitions. A bootstrap cannot invent variation absent from the observed paired differences. The full run should increase repetitions, and the report should distinguish Monte Carlo uncertainty in estimated coverage from uncertainty in a paired difference.

## Multiple comparisons and uncertainty

The artifact evaluates many configurations and six metrics. The bootstrap intervals are pointwise, not multiplicity-adjusted simultaneous intervals. Statements about individual intervals excluding zero should therefore remain exploratory, especially with only 100 repetitions per truth. The full experiment should report effect sizes, uncertainty, and consistency across truths rather than highlight isolated nominal intervals.

Coverage estimates are coarse at (n=100): each additional covered repetition changes empirical coverage by 0.01. The Wilson intervals in `summary.csv` quantify binomial uncertainty for each configuration's coverage estimate. They answer a different question from bootstrap intervals for paired differences; both should be retained.

## Prior-support contrasts: interpretation

At the correct 120 Mpc threshold, widening [50,90] to [40,100] or [30,110] changes both bias and interval width. Near truth values 52 and 88, the baseline prior places the truth close to a boundary, and the compact run shows severe 68% undercoverage. Widening the support improves empirical coverage in those cells but also increases the mean 68% interval width.

At truth 70, the narrow baseline prior yields higher empirical coverage and lower RMSE in this particular sample than the widest prior. This does not establish that the narrow prior is generally superior: the truth is centrally located in it, and the performance is conditional on the specified truth and simulation mechanism. Prior choice must be justified independently of which setting looks best after observing simulation results.

## Scientific limitations and next validation

- The generator is a low-redshift two-host Gaussian distance toy model with a hard observed-distance cutoff.
- The selection-correction denominator is analytically computed for that same toy host mixture and interval; it is not a realistic GW search selection function.
- The experiment has one synthetic event per repetition, not a hierarchical population of multiple GW events.
- Host redshifts and mixture weights are known exactly, and no realistic redshift, peculiar-velocity, waveform, or calibration uncertainties are included.
- The true values near the baseline prior edges are deliberate stress tests, but should be complemented with more truth values and an explicitly reported grid-convergence check.
- The compact CI workflow verifies code paths and artifact generation. It is not a substitute for the full 300-replicate/4,001-grid run, numerical convergence analysis, or independent scientific review.

Recommended next step: execute the full design, preserve raw outputs and software revision, then rerun the key summaries at a finer posterior grid on a smaller matched subset to quantify discretization error. Do not claim observational (H_0) constraints from this experiment.
