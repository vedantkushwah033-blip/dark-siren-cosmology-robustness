# Cross-experiment synthesis: selection normalization and prior boundaries

**Prepared:** 2026-10-24  
**Scope:** synthesis of existing synthetic experiments in this repository. This is not a new simulation run and does not establish an observational result.

## Executive conclusion

The current evidence supports a narrow, conditional statement: in the implemented toy model, dividing the likelihood by the modeled detection efficiency improves prior-predictive simulation-based calibration (SBC) relative to an uncorrected negative control. However, the correction does not uniformly improve fixed-truth point estimates or interval coverage. The strongest repeatable signal across the fixed-truth and prior-sensitivity experiments is that results near the edges of the bounded prior are highly sensitive to prior support and can have severe pointwise undercoverage.

These findings are not contradictory. Prior-predictive SBC averages over parameter values drawn from the prior and tests ranks under that joint prior-predictive procedure. Fixed-truth coverage conditions on a particular parameter value and measures repeated-sampling behavior there. A model can pass a prior-averaged rank diagnostic and still have poor coverage at selected parameter values.

## Evidence map

| Existing experiment | Design / evidence | Main defensible inference | Key limitation |
|---|---|---|---|
| Selection-conditioned SBC | 2,000 prior-predictive draws; two-host model; detection interval (0,180) Mpc; [50,90] prior | Corrected ranks were consistent with uniformity in the reported KS and decile tests; uncorrected ranks showed departures | A diagnostic for the stated toy model, not astrophysical validation |
| Single-host SBC control | 2,000 prior-predictive draws; one host at z=0.03 | Similar corrected-vs-uncorrected pattern persists in a simpler one-host control | Still uses Gaussian noise and a hard toy threshold |
| Fixed-truth coverage map | 500 replicates at each of seven truths; [50,90] prior | 68% coverage collapses near prior edges; correction helps some truths and hurts others | Fixed prior and two-host mixture limit generality |
| Prior-support sensitivity | 300 matched observations per truth at 52, 70, 88; priors [50,90], [40,100], [30,110] | Wider support greatly changes edge-truth coverage; correction remains truth-dependent | Wider support is not automatically better and does not select a scientifically correct prior |
| Threshold misspecification | 300 paired observations at truths 52, 60, 70, 80, 88; true threshold 120 Mpc, assumed 120–210 Mpc | Wrong thresholds can change bias, error, width and coverage, with effects depending on truth | Some metrics improve under misspecification by chance or compensating bias; this is not evidence that wrong thresholds are preferable |
| Full prior × threshold factorial | 21,600 rows; 300 replicates at truths 52, 70, 88; 3 priors × 4 thresholds × 2 methods | Interaction between prior support, assumed threshold, truth and method should be examined jointly | Synthetic model; the experiment does not simulate a calibrated detector or realistic galaxy catalogue |
| Selection-normalization quadrature check | Three observed distances; 4,001-point grid vs adaptive quadrature | Grid integration of the stated selected likelihood agreed with adaptive quadrature within the reported tolerance for those cases | Only three integrals; not a calibration or physical-model check |

## Quantitative anchors from prior reports

### Prior-predictive SBC

For the two-host selection-conditioned SBC audit, the corrected method had mean posterior-CDF rank 0.50092, KS statistic 0.01115 (p=0.96237), and decile chi-square p=0.94681. The uncorrected control had mean rank 0.47034, KS statistic 0.05405 (p=0.0000161), and decile chi-square p=0.00182.

The single-host control showed the same qualitative pattern: corrected mean rank 0.49834 and KS p=0.56946; uncorrected mean rank 0.47818 and KS p=0.0000680. These outcomes are consistent with the expected role of the selection normalization under the toy generative model, but non-significant tests do not prove exact uniformity.

### Fixed-truth coverage and prior boundaries

In the [50,90] prior-support experiment with 300 replicates, 68% coverage at truth 52 was 0.000 for both methods. At truth 88 it was 0.043 for both methods. In the same report, changing the prior to [40,100] raised coverage at truth 52 to 0.707 uncorrected / 0.727 corrected; at truth 88 it yielded 0.810 uncorrected / 0.587 corrected. The [30,110] prior yielded 0.667 / 0.707 at truth 52 and 0.923 / 0.873 at truth 88.

This is strong evidence that the narrow prior support interacts with edge truths in this setup. It is not a justification to select a prior after observing which one gives preferred coverage. A real analysis would need scientifically justified prior choices and sensitivity reporting.

### New full-factorial reference cell

The [50,90] prior and 120 Mpc assumed threshold reference cell in the 21,600-row experiment gave corrected 68% coverage of 0.067 at truth 52, 0.933 at truth 70, and 0.033 at truth 88. The exact values differ from earlier experiments because their detection thresholds, random seeds, and/or designs differ. Do not combine results across runs as though they were one sample or one common Monte Carlo estimate.

## What is and is not established

### Supported within the tested toy setup

1. Omitting an H0-dependent detection-normalization factor distorts prior-predictive ranks relative to the selection-normalized model in the reported SBC controls.
2. Fixed-truth performance is heterogeneous across H0 and prior support.
3. Near-edge truths can have severe 68% interval undercoverage under the [50,90] prior.
4. The correction is not a universal finite-sample improvement in RMSE, absolute error, or coverage at every truth.
5. Numerical normalization checks and artifact-integrity checks can pass while pointwise coverage remains poor.

### Not established

- That the model is a faithful approximation to LVK detection selection.
- That a real galaxy-catalogue incompleteness correction has been implemented.
- That the corrected method is uniformly superior for real dark-siren inference.
- That any one prior interval is scientifically correct.
- That the project is publication-ready or that the reported effects are novel.

## Recommended next experiment, in priority order

1. **A matched interior-and-edge coverage map under multiple priors.** Use a fixed grid of truths including interior values and values close to both prior boundaries. Reuse identical detected observations across methods and priors. Report Wilson intervals for coverage and paired bootstrap intervals for changes in absolute error, RMSE, and interval width.
2. **A single-host fixed-truth control.** Prior-predictive one-host SBC already exists; a fixed-truth one-host experiment would help isolate the interaction between the two-host mixture and the boundary effect.
3. **Predeclare numerical acceptance criteria and rerun at a finer grid.** Check all posterior interval endpoints across the complete design, not only a selected case.
4. **Only after these checks, extend model realism.** A realistic project needs a defensible source population, cosmological distance-redshift relation, measurement uncertainty, detector selection/injection efficiency, and a separate model of galaxy-catalogue completeness.

## Reporting recommendation

The current manuscript should frame this as a controlled methodological simulation study. The safest main claim is that selection conditioning matters for prior-predictive calibration in the stated toy model, while fixed-truth error and coverage depend strongly on truth and prior support. Avoid claims that the correction is universally beneficial, that the simulations measure H0 from gravitational-wave data, or that the toy detection rule represents an actual detector.
