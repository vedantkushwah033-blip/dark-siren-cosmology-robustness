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


## Literature positioning check (2026-11-19)

A targeted literature check reinforces that the project must be framed narrowly. Galaxy-catalog incompleteness and its cosmological consequences are already active, rapidly developing topics; the present benchmark is not a new general incompleteness correction.

1. **Borghi et al. (2026), A&A 706, A199, “Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology.”** This work implements a self-consistent incompleteness and host-weighting framework in CHIMERA, studying 100 O5-like BBH events, multiple completeness cuts, and host-mass weighting. It reports H0 forecasts and robustness to host-weight misspecification. This is close prior art for end-to-end incompleteness/host-weight sensitivity, but its central design is not the same as our finite-validation estimation of an object-level inclusion function followed by IPW density recovery. DOI: 10.1051/0004-6361/202557354.

2. **Leyde, Baker & Enzi (2026), JCAP 01 (2026) 013, “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology.”** This work reconstructs the underlying three-dimensional galaxy field and magnitude distribution while accounting for catalogue selection, using spatial correlations and mock validation based on the Millennium simulation. It establishes that catalogue completion can require structure-aware reconstruction rather than a simple homogeneous missing-host component. DOI: 10.1088/1475-7516/2026/01/013.

3. **Tagliazucchi et al. (2026), PRD 114, 083523, “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs.”** This is especially relevant because it explicitly fits an underlying galaxy-redshift distribution to an incomplete catalogue through a hierarchical likelihood that accounts for the survey selection function, then propagates the reconstructed distribution into H0 inference. It reports improved H0 precision over homogeneous completion in its mock configurations. DOI: 10.1103/vmzf-1wxq; arXiv:2608.24992. The paper's approach is GP reconstruction of the latent redshift field, not our parametric logistic inclusion-probability estimation plus inverse-probability weighting.

4. **Katsumata (online 2024; journal 2025), “How should we estimate inverse probability weights with possibly misspecified propensity score models?”** This statistical-methodology paper directly studies the gap between fitting propensity probabilities by likelihood and estimating the inverse weights needed by an IPW target. It discusses nonlinear misspecification and finite-sample weight error, and proposes a different weight-targeted approach. This is close general statistical prior art for the probability-model/weight-instability component; it is not specific to galaxy catalogues or dark-siren density recovery. DOI: 10.1017/psrm.2024.23.

### Consequence for our contribution claim

The defensible contribution, if supported by full results, is a **controlled synthetic sensitivity benchmark of finite-validation estimation of a catalogue inclusion mechanism and its downstream IPW density-recovery consequences**, with paired oracle/estimated arms, curvature misspecification, validation-size variation, coefficient-uncertainty integration, and pre-posterior weight/support diagnostics.

The work should not claim:
- a new general galaxy-catalogue completeness correction;
- a new GP/latent-field reconstruction method;
- that IPW is preferable to hierarchical missing-host or field-reconstruction methods;
- that statistical propensity-score results alone establish novelty in cosmology;
- or that density-stage improvements imply H0 calibration.

A more complete literature review is still required before any novelty claim. In particular, compare the exact estimand, selection likelihood, use of validation data, treatment of uncertainty in the inclusion model, and downstream target between our setup and each cited framework. A possible next research question is whether diagnostics available before cosmological inference can identify when fitted inclusion weights will degrade the recovered host-redshift density relative to the known-probability benchmark. This remains a hypothesis to test, not an established result.


## Posthoc posterior-mean versus MAP audit (completed 2026-11-19)

The separate posthoc workflow run `44395771659` completed successfully. Its artifact
`corrected-curvature-posterior-vs-map-audit` (ID `14474192424`) reports
`status=PASS` for all seven gamma cells. Each cell has 20 planned, 20
successful, zero failed replicates, and 20 passing child raw-input audits.
The auditor reconstructed 91,140 metric rows and 2,940 paired contrast rows.
The paired E-D contrasts use identical replicate IDs at each gamma, validation
size, and absolute bandwidth; every contrast has 20 paired replicates.

### Primary E-D endpoint: IAE of the posterior-mean density

The contrast is posterior-mean density minus MAP plug-in density. Negative
values favour the posterior-mean density on this density-IAE endpoint; positive
values favour MAP. The following are equally weighted averages of the seven
predeclared bandwidth-specific paired mean differences, not a bandwidth
selected from the data:

| gamma | validation n=100 | n=300 | n=1000 | bandwidth-specific 95% paired intervals |
|---:|---:|---:|---:|---|
| -2.0 | -0.01988 | -0.00807 | -0.00257 | all 7 below zero at each n |
| -1.0 | -0.01570 | -0.00548 | -0.00189 | all 7 below zero at each n |
| -0.5 | -0.00576 | -0.00065 | -0.00070 | all 7 include zero at each n |
| 0.0 | +0.00203 | +0.00395 | +0.00101 | 7/7, 5/7, and 6/7 include zero |
| +0.5 | +0.01306 | +0.00902 | +0.00308 | 3/7 include zero at n=100; 0/7 at n=300 and 1000 |
| +1.0 | +0.01972 | +0.01270 | +0.00412 | all 7 above zero at each n |
| +2.0 | +0.03106 | +0.01420 | +0.00467 | all 7 above zero at each n |

This is a clear curvature-dependent pattern in this particular synthetic
design: posterior averaging improves the IAE relative to MAP under the two
negative-curvature settings, is indistinguishable at the tested bandwidth
resolution for gamma=-0.5 and mostly so at gamma=0, and has higher IAE than
MAP for positive curvature, with the penalty generally shrinking as validation
size increases. The result does **not** establish that posterior integration
is generally better or worse. In particular, integrating coefficient
uncertainty within the fitted degree-one model does not represent uncertainty
about the omitted quadratic functional form.

### Secondary endpoint: posterior-weighted node-wise IAE

The separately reported posterior-weighted node metric is not the IAE of the
posterior-mean density. Its average E-D contrast is positive in all 21
gamma-by-validation-size cells (about +0.0008 to +0.0609 across the seven
bandwidths). At gamma=-2, several bandwidth-specific intervals include zero;
for the other six gamma values, all seven intervals are above zero at each
validation size. This distinction is expected: averaging density curves and
averaging node-specific losses are different operations and answer different
questions. They must not be conflated.

### Interpretation

The two E endpoints move differently. The posterior-mean density can have
lower IAE than MAP even when the posterior-weighted node-wise IAE is higher,
because the former evaluates the averaged density curve while the latter
averages losses before density curves are averaged. This is a useful
demonstration of why the estimand must be specified before interpreting
uncertainty propagation.

The E-D results remain exploratory: 20 paired replicates per cell, seven
fixed bandwidth sensitivity points, and paired t intervals. They do not
validate H0 posterior calibration, coverage, SBC, observational inference, or
a universal uncertainty-integration benefit. The expanded 11-bandwidth
workflow remains separate and is not considered complete until its own
cross-gamma audit uploads a passing artifact.
