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
