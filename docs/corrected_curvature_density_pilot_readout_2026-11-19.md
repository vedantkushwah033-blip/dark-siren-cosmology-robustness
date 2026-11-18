# Corrected seven-curvature density pilot: initial readout

Date: 2026-11-19

## Provenance and scope

This readout is based on the corrected 20-replicate grid workflow run
`44369630799` (commit `17ad47e3c4211abb7f0d3317c535d7fd4cbe6186`) and its
cross-gamma summary artifact
`corrected-quadratic-logit-curvature-cross-gamma-summary` (artifact
`14462504250`). The artifact was independently re-opened and its
`audit.json` and CSV summaries inspected.

The cross-gamma audit reports `PASS` for all seven curvature values
`{-2,-1,-0.5,0,0.5,1,2}`. Each cell contains 20 planned, 20 successful,
zero failed replicates, and 20 passing child raw-input audits. It verifies
paired seed schedules, common parent and validation samples, calibration
consistency, and per-cell raw inputs. The summary contains 91,140 metric rows,
4,557 summary rows, 2,205 paired contrast rows, and 756 gamma-versus-zero
paired contrast rows.

This is a density-only synthetic stress pilot. It does not include H0
inference, a bandwidth-selection exercise, real galaxy catalogues, or
observational gravitational-wave data. The 20 replicates per cell are
exploratory, not confirmatory.

## First empirical patterns

All comparisons below use paired replicate differences at each of the seven
predeclared absolute KDE bandwidths. The aggregate values below are simple
averages over those seven bandwidth-specific mean contrasts; they are not a
new pooled estimator or a bandwidth selected from the data.

### 1. Estimated IPW improves on the uncorrected catalogue baseline

For D (MAP plug-in estimated-probability IPW) minus A (catalogue-only), the
mean paired IAE difference is negative in every gamma/validation-size cell.
The paired 95% t interval excludes zero at all seven bandwidths in all 21
gamma-by-validation-size cells.

The average D-A difference across bandwidths ranges from about -0.32 to
-0.57 IAE units. This is strong evidence in this designed synthetic setting
that correcting the catalogue selection can materially improve density
recovery relative to leaving the catalogue uncorrected.

It is not evidence that estimated probabilities are always reliable: D is
still compared with a known-probability arm below, and the result is tied to
this DGP and its chosen bandwidth grid.

### 2. Estimating inclusion probabilities incurs a measurable penalty relative to known probabilities

For D minus C (MAP estimated IPW minus known-probability IPW), the mean
difference is positive in all 21 cells. The average penalty is particularly
large for positive curvature:

- gamma=+2: approximately +0.255 IAE across validation sizes;
- gamma=+1: approximately +0.103 to +0.132;
- gamma=+0.5: approximately +0.043 to +0.094.

At gamma=-2, the penalty is approximately +0.086 to +0.114. At gamma=0 it
falls from approximately +0.081 at n=100 to +0.012 at n=1000.

Bandwidth-specific paired intervals are uniformly above zero for all
bandwidths at gamma=-2, gamma=0, and positive curvature. At gamma=-1 and
gamma=-0.5, evidence weakens as validation size grows: several intervals
include zero. This is consistent with finite-validation uncertainty and
model misspecification interacting with the curvature design; it should not
be simplified to a monotone sample-size law.

### 3. The oracle parent reference remains substantially better

C (known-probability IPW) has higher IAE than B (the finite parent-sample
oracle reference) in all seven bandwidths across all gamma values. This is
not a failure of the inclusion correction: B uses a much larger complete
parent sample and is a Monte Carlo reference, not an estimator available
under catalogue incompleteness.

### 4. Posterior integration versus MAP remains to be audited

The original corrected cross-gamma summary predates the new explicit E-D
contrast. It contains E-vs-A/B/C but not E-vs-D. A separate posthoc workflow
has been added to calculate the paired posterior-mean-density minus MAP
plug-in contrast directly from the preserved corrected raw artifacts. Until
that workflow completes and passes, no conclusion should be drawn about
whether coefficient-posterior integration improves on the MAP plug-in
density estimate.

## Important interpretation limits

- The seven gamma values are not severity-matched. Negative and positive
  curvature can induce different ideal IPW weight instability even though
  the parent-average inclusion rate is calibrated to the same target.
- The D-A result addresses selection correction versus no correction; it
  does not establish that estimated inclusion probabilities are as good as
  known probabilities.
- The D-C gap is the more direct diagnostic of finite-validation/model-fit
  cost in this design.
- The B arm is a finite complete-parent-sample reference, not exact
  superpopulation truth.
- These are pilot estimates with 20 paired replicates per cell. Confidence
  intervals are exploratory paired t intervals, and the seven bandwidths
  are sensitivity points, not independent replications.
- Density IAE results do not establish H0 posterior calibration, coverage,
  or SBC performance.

## Next gates

1. Complete the independent posthoc E-D audit and regression test.
2. Inspect and independently validate the E-D contrast table, including
   whether the posterior-mean density endpoint is paired on identical
   replicate keys at every gamma, validation size, and bandwidth.
3. Inspect weight-tail/ESS and logistic-fit diagnostics alongside D-C,
   rather than relying on IAE alone.
4. Decide a predeclared bandwidth protocol. The current seven-bandwidth
   surface is sensitivity analysis; the parent-only Scott multiplier pilot
   must not be reused as a catalogue-IPW bandwidth selector.
5. Increase replication or use a predeclared MCSE stopping rule before
   confirmatory claims.
6. Keep H0 inference gated until density-stage diagnostics, bandwidth
   robustness, and raw-artifact audits are complete.

No claim of observational validation, universal estimator superiority, or
new estimator novelty is supported by this pilot.
