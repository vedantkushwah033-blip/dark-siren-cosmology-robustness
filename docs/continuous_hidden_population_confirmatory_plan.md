# Continuous hidden-population confirmatory simulation protocol

**Status:** preregistration proposal; no confirmatory simulation results have been generated.  
**Model implementation:** `src/dark_siren_robustness/hidden_population.py`.  
**Runner:** `scripts/run_continuous_hidden_population_recovery.py`.  
**Diagnostics:** `scripts/diagnose_hidden_population_identifiability.py` and `scripts/audit_hidden_population_grid_convergence.py`.

## 1. Model and fixed truth

Use the continuous log-distance model:
[
u_i=\log d_i,\qquad
u_i\mid z_i,H_0\sim\mathcal N(\log(cz_i/H_0),\sigma_u^2).
]

Listed hosts are (z=(0.02,0.04,0.06)) with equal relative weights. With probability (1-f), a listed host is drawn; with probability (f), the redshift is drawn from a Normal distribution with ((\mu_h,s_h)), truncated to fixed support ([0.005,0.10]).

Baseline truth:
- (H_0=70);
- (f=0.30);
- (\mu_h=0.075);
- (s_h=0.012);
- (\sigma_u=0.10).

The supplied analogue's other values should be introduced only as separately labelled sensitivity conditions.

## 2. Primary inference arms

All methods use the same simulated observations within a replicate.

1. **Oracle:** fix (f,\mu_h,s_h) at their generating values and infer H0. Non-deployable benchmark.
2. **Fixed wrong:** fix (f=0.10,\mu_h=0.060), with (s_h) fixed at truth.
3. **Plug-in:** estimate the hidden parameters jointly, then condition the H0 posterior on their joint-posterior marginal medians.
4. **Fully integrated:** marginalize the joint posterior over hidden parameters.

Run two parameter-recovery tiers:
- Tier A: condition (s_h) at truth; estimate (f,\mu_h,H_0).
- Tier B: estimate (s_h) jointly; this is the four-parameter recovery experiment.

Do not describe Tier A as full four-parameter recovery.

## 3. Confirmatory cell selection

Keep the initial confirmatory design small enough to be computationally reproducible:
- sample size (N\in\{50,200\});
- log-distance noise (\sigma_u\in\{0.10,0.20\});
- baseline truth above.

This creates four primary cells per recovery tier. Add (f=0), (f=0.45), and alternative hidden-redshift truths only as prespecified secondary stress tests, not after inspecting the primary results.

A minimum of 1,000 independent replicates per primary cell is recommended if the goal is a 95% Monte Carlo half-width of roughly 2.9 percentage points for coverage near 68%. If computation limits require fewer replicates, report the achieved Wilson intervals and label the results exploratory; do not silently treat 200 replicates as high precision.

## 4. Primary outcomes

For each method and cell:
- posterior median bias;
- RMSE;
- 68% and 95% fixed-truth coverage;
- mean 68% and 95% interval width;
- Monte Carlo standard errors for bias/RMSE/width;
- Wilson intervals for coverage;
- paired method-minus-oracle contrasts with bootstrap intervals.

Report posterior boundary accumulation for (f,\mu_h,s_h), and any multimodality diagnostic separately.

## 5. Local identifiability

Use the finite-difference negative Hessian of the mean log likelihood, averaged across independent datasets, as a Monte Carlo approximation to expected local information. Standardize parameter coordinates before reporting condition numbers.

Report:
- standardized information eigenvalues;
- condition number only when the matrix is positive definite;
- parameter correlations from the inverse information matrix;
- H0 variance inflation when nuisance parameters are free versus fixed.

Repeat with a second finite-difference step vector to check derivative stability. A finite positive-definite local matrix is not proof of global identifiability. If estimated information has negative eigenvalues, report that result rather than forcing positive definiteness.

## 6. Numerical convergence

For the same fixed dataset, compare:
- quadrature orders 24, 48, 96;
- coarse and refined H0/fraction/mu/sigma grids;
- posterior medians and central 68%/95% quantiles for all parameters.

Use paired data and report each parameter's maximum absolute quantile change. Choose tolerances in parameter units before confirmatory runs. A one-dataset convergence test is a numerical diagnostic, not a population-level guarantee.

## 7. Reproducibility

- Master seed and seed derivation must be written to the manifest.
- Preserve per-replicate raw posterior summaries and method-specific parameter estimates.
- Recompute all aggregate metrics and paired contrasts from raw rows in a separate auditor.
- Record source commit, package versions, grid dimensions, quadrature order, runtime and artifact SHA-256.
- Do not overwrite previous artifacts.
- No result is confirmatory until the exact source revision passes CI and an independent artifact audit.

## 8. Scope and interpretation

This is a low-dimensional synthetic statistical model. It is not a galaxy-field reconstruction, not a survey completeness estimator, not a realistic GW likelihood, and not a replacement for hierarchical dark-siren frameworks. The purpose is to study finite-sample posterior behaviour under transparent hidden-mixture assumptions.
