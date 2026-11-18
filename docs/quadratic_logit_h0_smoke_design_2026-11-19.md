# Quadratic-logit catalogue inclusion propagated into selected H0 likelihood

**Status:** exact-head matched-DGP software smoke and independent artifact audit PASS; repeated-sampling calibration remains separate and not yet established.
**Scope:** synthetic low-redshift Euclidean benchmark only.

## Why this experiment exists

The corrected seven-curvature density grid uses a quadratic-logit catalogue-inclusion law, while the earlier H0 smoke used a different sigmoid-with-floor law. Those results could not be combined into a curvature-to-H0 claim. This experiment closes that generator mismatch for a small first smoke.

## Generating process

The catalogue and independent validation frame use the same calibrated inclusion law:

[
pi_gamma(z)=\operatorname{expit}\{\alpha_\gamma+\beta x+
\gamma(x^2-c_g)\},\qquad x=(z-0.205)/0.195,
]

with beta = -3.9, parent-intensity centering c_g, and alpha calibrated separately for each gamma to preserve the common parent-average inclusion rate. The generator primitive is the same one exercised by the seven-curvature density-grid regression tests.

The analysis parent redshifts are shared with the GW generator. GW events are generated from the complete parent population with the existing host/merger weight and detection-threshold mechanism; the catalogue retention draw is not used to generate GW events. The H0 smoke fixes the event realization under the analysis seed, so curvature changes catalogue inclusion without changing the paired GW event data.

The independent validation frame is separately seeded. Its observable table contains only validation redshift and retention labels. True inclusion probabilities remain in truth-only files and are used only by the independent auditor.

## Inference

The first smoke deliberately fits a degree-one logistic inclusion model to a quadratic-logit DGP at gamma = +1.0. It therefore exercises a deliberately misspecified estimated-inclusion path, rather than pretending that the fitted model knows the true quadratic curvature. The fitted coefficients are integrated using the existing posterior quadrature.

For each coefficient node, the shared-inclusion selected likelihood computes the event host-mixture likelihood and its detection-efficiency normalization; the selection correction is applied within each node before coefficient-node marginalization. The independent auditor regenerates the coefficient quadrature from the observable validation records and independently reconstructs the nodewise selected likelihood, marginal likelihood, H0 posterior normalization and archived outputs.

## What this smoke can and cannot establish

A passing smoke establishes that the matched quadratic-logit DGP can be carried through the existing H0 inference code and that the archived arithmetic and truth firewall are independently auditable for this one realization.

It does **not** establish:
- H0 posterior calibration, frequentist coverage, PIT uniformity or joint SBC;
- that estimated inclusion improves H0 inference;
- performance across gamma, validation sizes or repeated samples;
- realistic detector selection, real galaxy-survey completeness or observational H0;
- universal superiority of any estimator.

## Next gate

Only after this smoke and its independent audit pass should a small paired multi-replicate pilot be launched. The first repeated design should use gamma = 0 and at least one positive and one negative curvature, common analysis/event seeds across gamma, and nested validation sizes. It must preserve all attempted cells in a failure ledger and separately report posterior-median bias/RMSE, PIT, equal-tail coverage, interval widths, quadrature evaluability and IPW stability. Conditional H0 PIT and fixed-truth coverage are distinct estimands; neither should be labelled joint SBC.


## Exact-head smoke readout (run 44402206818; 2026-11-18)

The exact branch-head workflow completed successfully. Artifact `14472907117`
contains a passing independent audit:

- 1,200 parent galaxies; 576 observed catalogue galaxies.
- 500-object independent validation frame; 150 validation observations.
- Four detected events; 932 detected candidates reconstructed from 1,000 generation attempts.
- Gamma=+1 quadratic-logit truth; degree-one fitted inclusion model.
- Independent selected-likelihood reconstruction has zero reported posterior-density difference at quadrature orders 3 and 5. Manifest verification covers 16 files.
- Inclusion quadrature was declared converged at order 5, with one consecutive passing transition. This checks coefficient quadrature order stability, not H0-grid convergence or calibration.

The posterior is numerically well inside the expanded [20, 300] prior, but this single realization is far from the generating value H0=70:
- Posterior median: 130.03 km s^-1 Mpc^-1.
- Equal-tail 95% interval: [86.37, 181.06] km s^-1 Mpc^-1.
- The generating value is outside even this one-realization 95% interval.
- Prior-slope-scale sensitivity (multipliers 0.5, 1, 2) moves the posterior median from 127.24 to 130.90.

The catalogue weight diagnostics are nontrivial: minimum catalogue Kish ESS 43.20 and maximum normalized weight 0.0811 at the reference prior scale. The posterior inclusion predictions span approximately 0.0038 to 0.9951 at order 5.

This is not evidence of systematic bias or method failure: it is one realization with four detected events. It does, however, make the next step unambiguous: repeated fixed-truth coverage and posterior-median error summaries are required before any H0 performance statement. The 90-cell pilot in PR #74 is the next exploratory gate; it must not be promoted to confirmatory evidence. The current smoke does not establish PIT, SBC, coverage, or an observational H0 measurement.
