# Parent-cluster uncertainty readout for the fixed-count H0 baseline

**Status:** independently replayed and audited; synthetic diagnostic only.

## Provenance and design

- Baseline workflow: `44417807075`
- Baseline artifact: `14479644393`
- Cluster-analysis workflow: `44433155472`
- Cluster-analysis artifact: `14486494486`
- The three baseline input files were SHA-256 verified by the analysis workflow.
- 10 independently generated parent populations, with 200 event realizations per parent (2,000 event realizations total).
- Four detected events per realization; fixed H0 truth 70 km s^-1 Mpc^-1.
- The independent sampling unit for unconditional uncertainty is the parent population, not the 2,000 event draws treated as independent.

## Results

| Quantity | Estimate | Parent-cluster MCSE | Parent-cluster t 95% interval | Cluster bootstrap 95% interval |
|---|---:|---:|---:|---:|
| Mean posterior-median bias (km s^-1 Mpc^-1) | +6.420 | 0.459 | [5.381, 7.460] | [5.693, 7.371] |
| 68% interval coverage | 0.672 | 0.0131 | [0.642, 0.702] | [0.648, 0.697] |
| 90% interval coverage | 0.884 | 0.0090 | [0.864, 0.904] | [0.866, 0.899] |
| 95% interval coverage | 0.943 | 0.0072 | [0.927, 0.959] | [0.929, 0.955] |
| Mean 68% interval width (km s^-1 Mpc^-1) | 24.198 | 0.294 | [23.534, 24.863] | [23.678, 24.765] |

For comparison only, treating all 2,000 event realizations as independent conditional on these ten archived parent populations gives bias MCSE 0.414 and t interval [5.609, 7.232]. This conditional event-level uncertainty is narrower than the parent-cluster uncertainty and is not an appropriate substitute for it when generalizing across parent populations.

## Interpretation and limitations

1. The fixed-count full-parent diagnostic exhibits positive posterior-median bias of about +6.42 km s^-1 Mpc^-1 even with the complete parent population supplied. Therefore the earlier positive bias in the catalogue/IPW pilot cannot be attributed solely to estimated catalogue inclusion probabilities or quadratic-logit misspecification.
2. Coverage estimates are below nominal at the 68% and 90% levels in this small parent-cluster sample; 95% coverage is closer to nominal. The intervals above quantify uncertainty across only ten independent parent populations, so they remain imprecise.
3. The 2,000 event realizations are nested within ten parent populations. The event-level MCSE is conditional on those populations and must not be reported as the unconditional precision.
4. This experiment conditions on exactly four detected events. It does not use event-count information, does not constitute joint prior-predictive SBC, and is not an observational H0 measurement.
5. This is a baseline diagnostic for the fixed-count toy likelihood, not a universal claim about dark-siren inference or a causal decomposition of all sources of bias.

The next scientifically useful step is to diagnose why the full-parent oracle itself has positive bias (distance likelihood, one-sided detection conditioning, fixed event count, prior measure, and finite H0 grid) before interpreting the catalogue-inclusion arms.
