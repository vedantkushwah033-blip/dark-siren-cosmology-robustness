# Preregistered next-stage protocol: quadratic-curvature H0 calibration

**Status:** protocol proposal; no multi-replicate H0 calibration results are claimed.

## Research question

In the synthetic low-redshift Euclidean dark-siren benchmark, does finite-validation estimation of the catalogue-inclusion model produce calibrated H0 posteriors when the inclusion logit is correctly specified, and how does degree-one misspecification change posterior calibration when the true inclusion logit has quadratic curvature?

This is a deliberately bounded simulation question. It is not an observational H0 measurement, LVK/GWTC analysis, real-survey validation, or a claim that degree two is universally superior.

## Design

### Data-generating conditions

- Parent size: (N=1,200).
- Redshift domain: (0.01 le z le 0.40), sampled from the declared two-overdensity parent intensity.
- Catalogue inclusion:
  [
  R_isimmathrm{Bernoulli}(pi_i),qquad
  operatorname{logit}(pi_i)=alpha_gamma+eta x_i+
  gamma(x_i^2-c_g),
  ]
  where (x=(z-0.205)/0.195), (eta=-3.9), and (c_g=E_g[x^2]) under the parent intensity.
- Curvature cells: (gammain{0,1}). The intercept is calibrated separately in each cell to preserve the same parent-averaged inclusion rate (approximately 0.510192).
- Independent validation frame: 1,000 objects per replicate; nested simple-random-sample prefixes (nin{100,300,1000}). The fitting data are validation redshifts and binary retention labels only.
- Event generation: four detected events in the first feasibility stage, drawn from the full parent population using the declared merger weight, with the common Gaussian distance noise and one-sided threshold already used by the smoke. The same parent and event seeds are paired across validation sizes and curvature cells where the generator permits common random numbers.
- Generating truth: (H_0=70) km s(^{-1}) Mpc(^{-1}).

### Inference arms

1. Catalogue-only: no estimated inclusion correction.
2. Degree-one estimated inclusion: validation posterior under the linear-logit model, integrated through the shared selected-event likelihood.
3. Degree-two estimated inclusion: validation posterior under the correctly specified quadratic-logit model, integrated through the same likelihood.
4. Oracle true-inclusion diagnostic: exact generating coefficients used only as a simulation reference, never as an operational estimator.

The same observed catalogue and event set must be used across validation sizes and inference arms within a replicate. The same inclusion-coefficient node is shared across all events when marginalizing the joint likelihood. Apply the merger/base-population weight once, and do not add a second missing-host mixture correction to IPW.

### Priors and numerical settings

- Fix coefficient coordinates to the parent domain: center 0.205 and scale 0.195.
- Record explicit coefficient-prior means and standard deviations in every artifact; do not rely on implicit library defaults.
- Primary quadrature order: 7. A predeclared numerical audit subset of at least 10 replicates per curvature cell compares orders 5 and 7 through the final H0 posterior.
- H0 support and grid must contain the generating truth in the interior. The initial feasibility stage may use 61 points over 55–85; the confirmatory run must repeat a grid-refinement audit on a fixed subset before locking the primary grid.

## Staged execution

### Stage 0 — feasibility and independent audit

- Run 2–5 paired replicates per curvature cell.
- Verify exact RNG replay, catalogue/validation truth joins, event thresholding, truth-host joins, posterior normalization, and independent likelihood reconstruction.
- Measure wall time and memory per replicate before scaling.
- Do not report these replicates as calibration evidence.

### Stage 1 — pilot calibration

- Run 20 paired replicates per curvature cell at all three nested validation sizes.
- Preserve every attempted cell and classify it as SUCCESS, FAILED, NON_EVALUABLE, or NOT_RUN.
- Publish raw per-replicate posterior summaries and diagnostics before aggregation.
- Report posterior-median bias, RMSE, PIT, central equal-tail interval coverage at 68%, 90%, and 95%, interval widths, and Monte Carlo uncertainty.
- Report Wilson intervals for coverage and evaluable fractions. Include all-attempt worst/best bounds when failed or non-evaluable cells make outcomes unknown.
- Treat pilot results as descriptive and use them to assess feasibility, not to make a final estimator ranking.

### Stage 2 — confirmatory calibration

Proceed only after Stage 1 passes its independent artifact audit and runtime budget is measured. Target at least 100 paired replicates per curvature cell, retaining the same fixed design and estimands. Do not change curvature cells, priors, H0 support, or success criteria in response to interim outcomes. If 100 replicates are computationally infeasible, report the achieved count and uncertainty honestly rather than silently reducing or selectively extending cells.

## Required independent audit

The auditor must reconstruct from raw artifacts, without importing the runner's aggregation outputs as truth:

1. Parent and independent validation frames from declared seeds and generator rules.
2. Bernoulli retention labels and nested validation prefixes.
3. Catalogue membership and observable/truth firewall.
4. Event host draws, distance noise, one-sided threshold, and truth-host joins.
5. Degree-one and degree-two logistic posterior targets, priors, quadrature nodes and normalized weights.
6. Node-specific inverse-inclusion host weights, Kish ESS, maximum normalized mass, and minimum log inclusion probability.
7. Each event mixture and the one-sided detection efficiency using the same host weights.
8. Shared-node marginalization across the joint event likelihood.
9. H0 posterior normalization, quantiles, PIT, coverage indicators, and grid/quadrature convergence.
10. Exact manifests, byte counts, checksums, complete failure ledger, and planned-cell accounting.

A smoke artifact may contain truth-only files for auditing, but non-oracle inference must consume only observable files.

## Primary estimands

For each ((gamma,n,	ext{arm})), report:

- Median bias: (E[operatorname{median}(H_0mid D)-H_{0,mathrm{true}}]).
- RMSE of posterior medians.
- Empirical PIT distribution and mean PIT, with a uniform-reference histogram and Monte Carlo uncertainty.
- Equal-tail 68%, 90%, and 95% interval coverage with Wilson intervals.
- Median interval width at each level.
- Evaluable fraction and all-attempt bounds.
- IPW stability diagnostics and their association with posterior failures/calibration.

The fixed-truth PIT is a repeated-sampling calibration diagnostic; it is not by itself simulation-based calibration under a prior-predictive distribution. Do not label this experiment SBC unless a separate prior-predictive truth-generation design is implemented.

## Stop/go criteria

**GO to Stage 2 only if:**
- Stage 0/1 independent raw-record audits pass.
- All planned cells are represented in the failure ledger.
- The H0 grid and coefficient quadrature checks meet predeclared tolerances on the numerical-audit subset.
- No unresolved double-counting of catalogue or GW selection remains.
- Runtime and memory permit the planned replicate count without changing the scientific design.

**NO-GO / revise before Stage 2 if:**
- Raw-record reconstruction fails.
- Event selection or host-population weighting is inconsistent between generator and likelihood.
- Posterior results are materially sensitive to H0 grid or quadrature order.
- Failure rates are large or differ systematically by model/validation size and are not transparently bounded.
- The model/prior parameterization cannot be independently reconstructed.

## Interpretation boundary

The degree-two model is correctly specified for the chosen synthetic inclusion mechanism; this does not guarantee it will outperform degree one in finite samples. Any observed advantage must be reported with paired uncertainty and the failure ledger. Oracle performance is diagnostic only. No result from this protocol supports claims about real galaxy catalogues or current gravitational-wave observations.
