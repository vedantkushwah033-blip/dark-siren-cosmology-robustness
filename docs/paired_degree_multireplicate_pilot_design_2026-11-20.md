# Paired degree-one / degree-two multi-replicate experiment

**Status:** expanded paired experiment queued; synthetic low-redshift Euclidean model only.

## Question

Across matched synthetic parent populations, detected GW events, validation frames and validation IDs, how do degree-one and degree-two catalogue-inclusion models differ in posterior error and inverse-probability weight stability?

This is a controlled sensitivity comparison. It is not a new estimator, a real-data H0 analysis, an LVK/GWTC analysis, or evidence that either degree is universally superior.

## Current expanded design

- Truth: quadratic-logit catalogue inclusion, curvature gamma in {-1, 0, +1}.
- Arms: catalogue-only (A), full-parent oracle (B), known-probability IPW (C), MAP degree-specific IPW (D), and coefficient-integrated degree-specific IPW (E).
- Paired replicates: **20**, master seed `2026112020`.
- Within each replicate and gamma, D1 and D2 receive identical parent, GW-event, validation-frame and nested validation-sample inputs.
- Parent size: 1,200; detected events: 30.
- Validation frame: 500; nested validation sizes: 100 and 500.
- H0 truth: 70 km/s/Mpc; prior support [30,110]; 281-point grid.
- Quadrature order: 5 for this fixed-truth paired experiment. A separate single-realization order 5/7/9 experiment assesses numerical quadrature sensitivity, but does not establish calibration.

The previous three-replicate run was an engineering pilot on the pre-correction branch. It must not be pooled with the new prior-matched 20-replicate run.

## Prior-matched parameterization

The generator uses

    eta(x) = alpha + beta*x + gamma*(x^2 - c_g).

Degree one uses basis `[1,x]` and independent Gaussian priors on `(alpha,beta)`. Degree two now uses the centered basis `[1,x,x^2-c_g]`, the same marginal Gaussian priors on `(alpha,beta)`, and an independent Gaussian prior on `gamma`.

The centering constant is the known synthetic parent-measure moment `c_g=E_g[x^2]=0.2845005579191241`; it is fixed by the simulation design, not estimated from the validation or catalogue sample. It is a parameterization choice for this controlled benchmark, not a proposed operational source of oracle information. A real-data application would need a predeclared centering convention or would need to propagate uncertainty in any data-estimated center. With declared prior SDs `(5, 2.5, 2.5)`, the D1 covariance on `(alpha,beta)` and D2 marginal covariance on `(alpha,beta)` are both `[[25,0],[0,6.25]]`. The audit script `scripts/audit_degree_prior_parameterization.py` and its test check this equality and zero cross-covariance between `(alpha,beta)` and `gamma`. This removes the previously identified raw-versus-centered prior confound. It does not guarantee that D2 estimates inclusion probabilities accurately or produces stable IPW weights.

## Pairing and audit gates

The orchestrator runs the independently audited replicate workflow separately for degree 1 and degree 2 using the same master seed. It then:

1. Requires both child workflows to report lifecycle PASS and their independent audits to pass.
2. Compares parent, catalogue, GW events, event truth, validation-frame truth and each nested validation sample across model degrees.
3. Requires identical H0 grids and A/B/C likelihood curves within strict 1e-12 absolute/relative tolerance.
4. Reconstructs paired D2-minus-D1 posterior-error contrasts at replicate/gamma/validation-size/arm grain.
5. Reports paired Monte Carlo standard errors and t-based intervals alongside bias, RMSE, coverage and weight diagnostics.
6. Hashes outer artifacts and records configurations, child audits and logs.

Incomplete lifecycle cells cannot be promoted to successful paired results. Failed and non-evaluable records remain visible in lifecycle ledgers.

## Separate quadrature sensitivity pilot

A completed single-realization pilot used N_parent=400, four detected events, gamma=+1, validation frame=300, validation sizes 100/300, and quadrature orders 5/7/9. Its artifact is a numerical sensitivity check only, not the 20-replicate paired experiment.

All six child audits passed, all 114 entries in the outer manifest verified, and all six nested manifests verified. Maximum reconstructed likelihood discrepancy was approximately 1.1e-14. Integrated-IPW posterior-density L1 difference from order 7 to 9 was at most 1.26e-5.

However, quadrature convergence did not imply weight stability. At order 9, minimum node-level Kish ESS was 9.60 (D1, n=100) and 3.17 (D2, n=100); maximum node top-5%-weight mass was 79.6% and 92.3%, respectively. At n=300, corresponding minimum ESS values were 21.87 and 14.47. This is a material instability warning for both degrees, particularly D2.

## Interpretation limits

- The quadrature pilot is one realization with four events; posterior medians in it cannot rank model degrees.
- The 20-replicate fixed-truth experiment is exploratory and does not constitute prior-predictive SBC or high-precision calibration.
- The oracle arm uses simulation truth for diagnosis only.
- Galaxy-catalogue inclusion and GW detection selection remain distinct.
- IPW and an additional missing-host mixture correction must not be applied together.
- Posterior median shifts alone do not establish improved calibration.
- Quadrature order 5 in the paired experiment is not evidence of order convergence.

## Promotion criteria

Do not describe either degree as superior unless a larger, independently audited prior-predictive experiment demonstrates stable calibration and paired improvement across curvature, validation size, event count, prior parameterization, quadrature/grid sensitivity and node-level weight stability. Keep the PR in draft state pending workflow and artifact checks.
