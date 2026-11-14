# Corrected three-floor density-pilot readout (2026-11-14)

## Scope and provenance

This report analyzes the exact-head artifacts from GitHub Actions workflow run `43617052487`, source commit `d94fdd62e514c81c20444196e166d2ebbb341a40`, for the three catalogue floors 0.00, 0.10, and 0.25.

All three jobs completed successfully. Each artifact has 20/20 successful replicates, zero failed replicates, 20 independently audited children, 244 verified manifest files, 9,300 reconstructed aggregate metric rows, and 1,740 reconstructed realized-Scott-width rows. Each outer audit reports PASS.

The run uses 4,000 parent objects, an independent 1,500-object validation frame, validation prefixes 100/300/1,000, a 401-point redshift grid, five candidate physical bandwidths (0.012, 0.016, 0.024, 0.032, 0.040), and 20 paired replicate seeds. The primary Arm E estimand is loss of the posterior-mean density. Posterior expected node-specific loss is retained as a separate secondary estimand.

This is a density-only exploratory pilot. It does not select a bandwidth, evaluate H0 posteriors, establish repeated-sampling calibration, or validate a real galaxy survey.

## Fixed-width results at h = 0.016

Mean IAE over 20 replicates at validation size 1,000:

| Catalogue floor | A: catalogue only | B: parent oracle | C: known-probability IPW | D: MAP estimated IPW | E: posterior-mean density |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.53049 | 0.18378 | 0.20391 | 0.20720 | 0.21032 |
| 0.10 | 0.37809 | 0.18378 | 0.18653 | 0.29618 | 0.30006 |
| 0.25 | 0.21349 | 0.18378 | 0.18065 | 0.27331 | 0.27605 |

The oracle arm is identical across floors because the parent population and target density are held fixed. The known-probability correction remains close to the parent oracle in all three cases.

## Paired contrasts at h = 0.016

The contrasts below are paired within replicate. Intervals are mean ± 1.96 Monte Carlo standard errors, and are descriptive normal approximations, not multiplicity-adjusted confirmatory intervals.

At validation size 1,000:

| Floor | Contrast | Mean | MCSE | Approx. 95% interval |
|---:|---|---:|---:|---:|
| 0.00 | D − A | -0.32330 | 0.01145 | [-0.34573, -0.30086] |
| 0.00 | D − C | 0.00329 | 0.01189 | [-0.02001, 0.02658] |
| 0.00 | C − B | 0.02013 | 0.00654 | [0.00730, 0.03295] |
| 0.00 | E − D | 0.00312 | 0.00035 | [0.00244, 0.00380] |
| 0.10 | D − A | -0.08192 | 0.01470 | [-0.11074, -0.05309] |
| 0.10 | D − C | 0.10964 | 0.01359 | [0.08301, 0.13628] |
| 0.10 | C − B | 0.00275 | 0.00401 | [-0.00512, 0.01061] |
| 0.10 | E − D | 0.00388 | 0.00012 | [0.00365, 0.00411] |
| 0.25 | D − A | 0.05983 | 0.00907 | [0.04205, 0.07761] |
| 0.25 | D − C | 0.09266 | 0.01012 | [0.07283, 0.11250] |
| 0.25 | C − B | -0.00313 | 0.00256 | [-0.00816, 0.00189] |
| 0.25 | E − D | 0.00274 | 0.00008 | [0.00258, 0.00290] |

The D − A contrast changes sign across the three floor conditions. In this synthetic setup, estimated IPW strongly improves over catalogue-only under the correctly specified floor-zero mechanism, retains a smaller benefit at floor 0.10, and is worse than catalogue-only at floor 0.25. At the latter floor, known-probability IPW remains close to the parent oracle, isolating estimation/model error from the correction's oracle potential.

This pattern is a hypothesis-generating result. It is not yet a general theorem or a claim about real catalogues.

## Validation-size trend

At h = 0.016, the paired D − C IAE contrast is:

| Floor | n=100 | n=300 | n=1,000 |
|---:|---:|---:|---:|
| 0.00 | 0.04621 | 0.01740 | 0.00329 |
| 0.10 | 0.16592 | 0.12987 | 0.10964 |
| 0.25 | 0.11890 | 0.07808 | 0.09266 |

Under correct specification, the MAP estimated correction approaches the known-probability benchmark as validation size grows. Under misspecification, increasing validation size does not eliminate the discrepancy. The non-monotonic floor-0.25 sequence also warns against reading a three-point sample-size trend as a convergence law.

## Post-hoc inclusion-probability and weight diagnostics

I additionally reconstructed posterior-mean inclusion probabilities from each replicate's stored coefficient nodes and quadrature weights, evaluated them on the parent redshifts, and compared them with the known synthetic inclusion curve. This is a post-hoc exploratory diagnostic, not a preregistered primary endpoint.

Mean parent-sample absolute probability error and mean posterior pointwise standard deviation:

| Floor | n=100 MAE / posterior SD | n=300 MAE / posterior SD | n=1,000 MAE / posterior SD |
|---:|---:|---:|---:|
| 0.00 | 0.03703 / 0.04775 | 0.02459 / 0.02814 | 0.01374 / 0.01542 |
| 0.10 | 0.04479 / 0.05045 | 0.03068 / 0.02977 | 0.02116 / 0.01637 |
| 0.25 | 0.05545 / 0.05292 | 0.04021 / 0.03130 | 0.03488 / 0.01705 |

The floor-zero MAE falls with validation size as expected under correct specification. At nonzero floors, posterior uncertainty narrows while absolute error remains larger, consistent with a residual functional-form discrepancy. At n=1,000, floor-0.25 MAE is about twice the floor-zero MAE, while posterior pointwise SD is similar. This is a diagnostic clue, not a formal model-checking theorem.

At n=1,000, the MAP-arm weight diagnostics also improve monotonically across floor conditions:

| Floor | Mean Kish ESS | ESS fraction | Maximum normalized weight | Top 1% mass | Top 5% mass |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 536.2 | 0.2445 | 0.01282 | 0.1466 | 0.3380 |
| 0.10 | 813.0 | 0.3429 | 0.00569 | 0.1086 | 0.2980 |
| 0.25 | 1,429.9 | 0.5421 | 0.00250 | 0.0608 | 0.2125 |

Thus the floor-0.25 estimated weights are less concentrated by these diagnostics, yet the density estimate is worse. The mechanism is not simply extreme inverse weights: the probability model is functionally wrong. The analysis should keep weight concentration and probability-model discrepancy as distinct diagnostic dimensions.

## Bandwidth-grid caution

The five-width grid does not support robust bandwidth selection. At n=1,000, the minimum among tested widths is at the upper endpoint 0.040 for catalogue-only at floor 0.00, and at the lower endpoint 0.012 for several arms at floors 0.10 and 0.25. Other arms have interior minima near 0.016. These are minima on a sparse exploratory grid, not validated optima. The expanded 11-width run must be audited before any selection protocol is frozen.

## Interpretation of Arm E

At h = 0.016, posterior-mean density loss is slightly larger than MAP plug-in loss in every floor and validation-size cell in this pilot. The E − D difference narrows as validation size increases. This is a descriptive finding, not evidence that posterior integration is generally inferior. Arm E's primary estimand is the loss of the posterior-mean density; its secondary posterior expected node-specific loss must not be substituted for it.

Aggregate Arm E summaries have no single weight vector and therefore no Kish ESS. Node-specific E quadrature rows retain their own ESS.

## Decision

1. Preserve the floor-family result as a completed exploratory density pilot.
2. Do not use it to choose a final bandwidth.
3. Complete and independently audit the expanded 11-width grid.
4. Run a separately preregistered quadratic-logit curvature family that preserves the parent-average inclusion rate while changing only curvature.
5. Keep H0 propagation blocked until density estimands, bandwidth selection, and independent audits are frozen.

