# Expanded absolute-bandwidth density pilot: independent audit readout

Date: 2026-11-15 (UTC workflow execution: 2026-11-14)
Branch: `research/absolute-width-distribution-pilot-20261114`
Experiment: 20 paired replicates per catalogue floor; parent population size 4,000; independent validation frame size 1,500; validation prefixes 100, 300, and 1,000; 11 common absolute redshift bandwidths:
`0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080`.

## Audit result

The expanded-grid workflow `43617052462` completed successfully for all three floors. Each of the three jobs independently audited all 20 replicate children and reconstructed the aggregate metrics:

| Catalogue floor | Replicates | Child audits passed | Failed replicates | Manifest files verified | Reconstructed metric rows | Reconstructed Scott-width rows | Audit |
|---|---:|---:|---:|---:|---:|---:|---|
| 0.00 | 20/20 | 20 | 0 | 244 | 20,460 | 1,740 | PASS |
| 0.10 | 20/20 | 20 | 0 | 244 | 20,460 | 1,740 | PASS |
| 0.25 | 20/20 | 20 | 0 | 244 | 20,460 | 1,740 | PASS |

The independent auditor scope includes seed/ledger/manifest checks, child raw-input reconstruction, and aggregate reconstruction. It explicitly does **not** perform bandwidth selection or H0 inference.

Workflow artifacts:
- Run: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43617052462
- Artifact IDs: 14123894369 (floor 0.00), 14123867382 (floor 0.10), 14123634445 (floor 0.25).

A separate repeated five-width run, workflow `43619594762`, also completed with all 20/20 child audits and zero failed replicates for each floor. It reproduces the central pilot pattern; it is not pooled with the expanded-grid experiment.

## Main results at validation size 1,000

The table reports mean integrated absolute error (IAE) over 20 paired replicates at the common bandwidth (h=0.016). Lower is better.

| Arm | Floor 0.00 | Floor 0.10 | Floor 0.25 |
|---|---:|---:|---:|
| A: catalogue-only | 0.53049 | 0.37809 | 0.21349 |
| B: parent oracle | 0.18378 | 0.18378 | 0.18378 |
| C: known-probability IPW | 0.20391 | 0.18653 | 0.18065 |
| D: MAP estimated IPW | 0.20720 | 0.29618 | 0.27331 |
| E: posterior-mean density | 0.21032 | 0.30006 | 0.27605 |

Paired replicate-level contrasts at (h=0.016), reported as mean difference and Monte Carlo standard error (MCSE):

| Floor | D − A | D − C | E − D |
|---|---:|---:|---:|
| 0.00 | −0.32330 (0.01145) | +0.00329 (0.01189) | +0.00312 (0.00035) |
| 0.10 | −0.08192 (0.01470) | +0.10964 (0.01359) | +0.00388 (0.00012) |
| 0.25 | +0.05983 (0.00907) | +0.09266 (0.01012) | +0.00274 (0.00008) |

## Bandwidth-grid interpretation

The 11-width grid confirms that the key qualitative pattern is not an artifact of choosing only (h=0.016):

- At floor 0.00, D is substantially better than A across the common bandwidths. D and C remain close; the minimum mean IAE for D occurs near (h=0.016) (0.20720), with (h=0.024) almost identical (0.20800).
- At floor 0.10, D improves on A modestly at small/intermediate widths but remains far worse than C. D's best observed mean IAE is about 0.29618 at (h=0.016), while C is about 0.18653 there.
- At floor 0.25, D is worse than A throughout the evaluated grid. D's best observed mean IAE is about 0.27317 at (h=0.012), compared with A at about 0.21271 and C at about 0.18350 at that same width.
- The oracle arm B and known-probability arm C remain close over the useful central bandwidth region. This isolates the estimated inclusion model as a major source of error in the stressed floors.

These are exploratory minima over the same 20 replicates used to evaluate the grid. They are **not** a valid final bandwidth-selection procedure: selecting and evaluating on the same replicates is optimistic. The common grid is broad enough to expose the broad response shape, but a leakage-free training/holdout selection protocol remains required before fixing a production bandwidth.

## Important distinction: probability error versus weight concentration

The earlier three-floor analysis found that the estimated-probability model's parent-population MAE increased with floor/misspecification, while Kish ESS increased and normalized-weight tails became less concentrated. The expanded-grid results reinforce that these diagnostics must remain separate: less concentrated weights do not imply a better recovered density when the inclusion model is structurally wrong.

## Interpretation and limits

The result supports a narrow synthetic statement: in this controlled setup, estimated IPW can recover much of the catalogue-only density error under the correctly specified inclusion model, but model misspecification can erase that benefit or reverse it. Knowing the true inclusion probabilities performs much better in the stressed cases, so the degradation is not simply unavoidable variance from inverse weighting.

The 20-replicate design is a pilot, not a final uncertainty assessment. The floor experiments are separate data-generating regimes, not a realistic survey completeness calibration. No real galaxy catalogue, LVK/GWTC event, observational H0 measurement, or external cosmological data are used. This density-only experiment does not establish H0 posterior calibration or select an H0 likelihood.

## Next research gate

1. Build the seven-curvature DGP density runner using the calibrated parent-measure-centred quadratic-logit model already implemented.
2. Extend the five-arm weight builder so Arm C can consume externally supplied true inclusion probabilities, with strict shape/range validation and unchanged default behavior.
3. Reuse the same paired replicate seeds, 100/300/1000 validation prefixes, and independent raw-record audit design.
4. Freeze bandwidth selection on training replicates and score only held-out replicates.
5. Only after the density gate passes, propagate selected density estimates into the shared-theta selected H0 likelihood and perform separate grid/quadrature/calibration audits.

No H0 propagation is authorized by this pilot alone.
