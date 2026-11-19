# Cluster-aware uncertainty for the full-parent H0 baseline

## Why this audit is needed
The full-parent baseline repeats 200 event realizations within each of 10 archived parent-galaxy populations. Event realizations are independent conditional on a fixed parent, but they share that parent population. Treating all 2,000 event draws as 2,000 independent units understates uncertainty for conclusions that generalize over parent-population realizations.

The appropriate independent cluster for the unconditional finite-parent design is the parent population (10 clusters), not the event realization (2,000 nested draws).

## Analysis
Run:

```bash
python scripts/analyze_full_parent_h0_cluster_uncertainty.py <baseline-output> --output-dir <readout>
```

The script independently reconstructs parent-level means for posterior-median bias, 68/90/95% coverage, and 68% interval width from `replicates.csv`. It reports:
- parent-cluster means, standard deviations, Monte Carlo standard errors, and t intervals with (m-1) degrees of freedom;
- parent-cluster percentile-bootstrap intervals;
- the event-weighted median-bias summary and its event-level standard error, explicitly labeled conditional on the archived parent populations;
- the number of parent clusters and event draws per cluster.

For equal numbers of event draws per parent, the grand event-weighted mean equals the mean of parent means, but its event-level standard error is not the unconditional uncertainty across parent populations. The cluster-level readout is the relevant one for generalization over the parent-population distribution.

## Interpretation
The previously reported 2,000-event Wilson intervals and event-level MCSE are conditional on the 10 archived parent populations. The cluster-aware intervals quantify variation across those 10 parent realizations, but remain imprecise because there are only 10 independent clusters. Neither calculation is joint SBC. The analysis remains a fixed-count synthetic low-redshift Euclidean benchmark, not observational H0 inference.


## Audited 2,000-event result

The cluster analysis ran in workflow `44419277487` and produced artifact `14482031023`. The source baseline was artifact `14479644393`; all three entries in its aggregate manifest were verified before analysis. The result contains 10 independently generated parent populations and exactly 200 event realizations per parent (2,000 total).

| Endpoint | Estimate | Parent-cluster t 95% interval | Parent-cluster bootstrap 95% interval |
|---|---:|---:|---:|
| Mean posterior-median bias (km s^-1 Mpc^-1) | +6.420 | [+5.381, +7.460] | [+5.693, +7.371] |
| 68% interval coverage | 0.672 | [0.642, 0.702] | [0.648, 0.697] |
| 90% interval coverage | 0.884 | [0.864, 0.904] | [0.866, 0.899] |
| 95% interval coverage | 0.943 | [0.927, 0.959] | [0.929, 0.955] |
| Mean 68% interval width (km s^-1 Mpc^-1) | 24.198 | [23.534, 24.863] | [23.678, 24.765] |

The event-weighted mean median bias has a conditional-on-these-ten-parents MCSE of 0.414 and t 95% interval [+5.609,+7.232]. The parent-cluster MCSE is 0.459, with the intervals above. The distinction matters: 2,000 event draws do not constitute 2,000 independent parent-population draws.

The positive median bias is present within all ten parent populations (parent-specific means range approximately +4.90 to +9.89). This strengthens the conclusion that the positive bias in the 90-cell pilot is not solely an estimated-inclusion or curvature-misspecification effect. It does **not** establish a likelihood defect: the benchmark has only four detected events per realization, fixed-count selection, and a low-redshift Euclidean toy model. Coverage is close to nominal for 68% and 95% intervals at this cluster precision; 90% coverage is 0.884, with a bootstrap interval just below 0.90 but a t interval that includes it. With only ten parent clusters, treat this as a conditional baseline, not a definitive population-level calibration verdict.
