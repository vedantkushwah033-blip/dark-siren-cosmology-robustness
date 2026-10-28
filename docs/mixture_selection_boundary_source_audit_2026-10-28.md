# Source-level audit: mixture normalization, selection efficiency, and prior boundaries

**Audit date:** 2026-10-28  
**Scope:** targeted source review of the current research/distance-uncertainty-likelihoods branch. This is not a full production run or an astrophysical validation.

## Executive finding

The inspected implementation has the expected mathematical structure for the repository's finite synthetic model:

- Catalogue and missing-host weights are normalized internally and multiplied by their respective population masses once.
- Likelihood mixtures are evaluated in log space with logsumexp and logaddexp.
- The selection-conditioned posterior uses p(d | H0, detected) proportional to p(d | H0) / beta(H0).
- beta(H0) is computed from the same finite host mixture and the same distance-detection interval as the event model.
- The selected missing-host API uses full-population weights in both the event likelihood and selection denominator.

The source review did not find evidence of the historical missing-mass double-counting error in the inspected current functions. This is a source-level conclusion, not a claim that every branch, workflow artifact, or scientific assumption has been validated.

## Files inspected

- src/dark_siren_robustness/hardened_inference.py, blob 0939f6a0244793c88d1eee9ca7310b5fff9e90cc
- src/dark_siren_robustness/selection.py, blob a69c6a16093c7c09af807a013ee7c6501fd70a9d
- src/dark_siren_robustness/cosmological_inference.py, blob 1b62073f9384a8c635541e8ddeb8b492848d2948
- src/dark_siren_robustness/uncertainty_likelihoods.py, blob 01b3282ccdefa083d5812f3a5cd19fc613b8eb89
- src/dark_siren_robustness/inference.py, blob c6160412be2e42a84b51b32c0e14a7deeec6b8a9
- tests/test_joint_selection_missing_hosts.py, blob 08dd64d876f5fe715ea598b1c318c505ae1185f1
- tests/test_cosmological_inference_weight_normalization.py, blob b71b9d381dd6b22066c0c6730ff1bd59c63d9d9e
- tests/test_hardened_inference.py, blob 968706ac04beb326486f026bd28aeedc1c82eb1a

## 1. Mixture weights and missing population mass

For normalized relative catalogue weights w_tilde_i = w_i / sum(w_j), the intended unselected model is

L(H0) = (1-f) sum_i w_tilde_i p(d | z_i,H0) + f sum_k v_tilde_k p(d | z_k_missing,H0).

In hardened_inference.posterior_h0, catalogue weights are normalized and multiplied by (1-f); missing-host weights are separately normalized, and the missing likelihood is multiplied by f in log space before logaddexp. This implements each population mass once. The selected API in selection.posterior_h0_selected_with_missing_hosts similarly scales the normalized catalogue and missing-host weights, concatenates the components, and sends the full population to the selected posterior.

The inspected tests include equivalence against an explicitly expanded full mixture, zero-missing-mass limits, rescaling invariance, and an independently constructed numerator/selection-denominator calculation. These are good targeted tests for the earlier double-counting failure mode.

**Important semantics:** this is correct if missing_weight=f means the prior population mass assigned to missing hosts. It is not interchangeable with a post-normalization likelihood fraction or with an empirical catalogue completeness fraction unless the model defines that equivalence.

## 2. Numerical stability

The inspected inference paths use scipy.special.logsumexp for weighted host mixtures and numpy.logaddexp for the two population components. They normalize extreme finite weights by dividing by their maximum before summing. Zero weights are represented as negative infinity in log-weight arrays rather than by taking log(0).

The selection code's _normal_interval_probability uses log-CDF and expm1 differences in the tails and an erf expression for intervals crossing zero. The implementation fails closed when a computed beta(H0) is not representable as a positive probability; it does not silently replace zero with machine-tiny. This is preferable to fabricating support, though sufficiently extreme parameter configurations require a revised numerical method or grid.

No inspected path uses naive sum(norm.pdf(...)) as the production host-mixture likelihood.

## 3. Selection normalization beta(H0)

For the declared toy detection rule d_min < d_obs < d_max, with d_obs | H0,z_j distributed as Normal(mu_j(H0), sigma^2), the code computes

beta(H0) = sum_j w_tilde_j [Phi((d_max - mu_j(H0))/sigma) - Phi((d_min - mu_j(H0))/sigma)].

The selected posterior is proportional to the host-mixture likelihood divided by beta(H0), followed by normalization over prior support. The missing-host variant includes the full catalogue-plus-missing mixture in both terms. The inspected tests check an independent direct calculation, zero-missing equivalence, a complete-population equivalence, extreme tail probabilities, and weight-rescaling invariance.

This establishes consistency with the toy finite-interval selection rule. It does not establish the correct selection normalization for a real GW search. A real analysis needs an appropriate detector/search selection function and population integration, including relevant distance, redshift, source parameters, and event-detection pipeline. The current beta must not be described as LVK-calibrated.

## 4. Gaussian likelihood and normalization

uncertainty_likelihoods.gaussian_log_likelihood includes -log(sigma) - 0.5 log(2 pi). This term is essential when sigma depends on candidate H0: omitting it would change the likelihood surface in Arm B. Arm A holds the realized sigma fixed over candidate H0; Arm B uses sigma(H0)=f c z/H0 including its normalization; Arm C uses a fixed reference sigma. The module documentation correctly describes this as a one-host low-redshift synthetic model.

hardened_inference.posterior_h0(..., positive_truncated=True) subtracts log Phi(mu/sigma), the normalizing probability for a Gaussian observation conditioned to be positive. This is appropriate only when the data generator's positive-distance truncation is part of the assumed sampling model; it should not be enabled automatically for a different observation model.

## 5. Boundary coverage

The posterior functions apply a uniform prior by zeroing likelihood density outside [H0_min,H0_max] and normalizing by trapezoidal integration over grid points inside the prior. posterior_summary trims zero-density exterior points before integrating, avoiding artificial trapezoid mass across the prior discontinuity. The tests explicitly exercise this boundary quadrature behavior.

For a continuous posterior supported on [50,90], a central equal-tail interval's endpoints are ordinarily strictly inside the support. Thus a fixed truth exactly at 50 or 90 can have zero coverage for such intervals by construction. This is a property of the interval-and-prior setup, not by itself proof of a broken likelihood. However, it is not a reason to suppress the result: report exact-boundary truths separately, and do not pool them with interior truths without clearly labeling the estimand.

The interior-truth Bayesian coverage above nominal in the archived simulation is a distinct question. It should not be called a likelihood bug solely because it differs from 68%/95%; Bayesian credible coverage at fixed truths need not equal frequentist confidence coverage. A matched reference interval and prior-sensitivity analysis are appropriate follow-ups. The prior-predictive PIT/SBC experiment is a different calibration estimand and cannot by itself settle fixed-truth coverage.

## 6. Remaining gaps and release gates

1. **Grid resolution:** an independent extreme single-host stress test has shown that a fixed 1,601-point grid can under-resolve deliberately narrow posteriors. Existing headline-design grid stability was good in its archived test, but the production grid should ideally adapt/refine based on posterior width or expose a convergence criterion, rather than rely only on a fixed grid.
2. **Zero-redshift edge case:** some APIs accept non-negative redshifts including zero, while uncertainty_likelihoods.candidate_distance_mpc requires strictly positive redshift. Clarify whether zero-redshift hosts are intentionally allowed in each public interface; otherwise align validation.
3. **Selection model scope:** finite distance-threshold selection is a controlled testbed, not a detector efficiency model. Realistic selection normalization remains an explicit external-validation gate.
4. **Independent run status:** this note records source inspection only. It does not claim tests were executed locally or that the newest GitHub Actions run passed. The new reference-interval audit workflow must be checked against its exact commit and artifact before release.
5. **Prior-sensitivity duplication:** reconcile the unintegrated nine-patch research bundle before starting another prior-support experiment.

## Bottom line

No source-level evidence of a remaining missing-mass double-counting error was found in the inspected functions. Weight normalization, log-space mixture arithmetic, and the toy selection denominator are internally consistent with the stated synthetic model and have targeted independent tests. The largest remaining threats are scientific scope, fixed-truth versus Bayesian calibration interpretation, prior-boundary reporting, posterior-grid convergence, and the gap between toy selection and a realistic GW detection pipeline.


## 7. API-contract details discovered during cross-module review

### Scalar versus host-specific measurement uncertainty

The conceptual likelihood sometimes written for this project uses host-specific uncertainties, sigma_{d,i}. The currently inspected selection APIs instead accept a single scalar gw_sigma shared by all hosts in a call. uncertainty_likelihoods.py also provides scalar/array-broadcasting functions for the three one-host sensitivity arms, but it is not a general heteroskedastic catalogue likelihood. This is not a bug for the explicitly declared shared-sigma toy model; it is a scope mismatch if the paper's equation claims arbitrary sigma_{d,i}. Either state the common-sigma assumption in the methods and notation, or implement and test host-specific uncertainties consistently in both event likelihood and beta(H0) before claiming that extension.

### Zero-redshift validation differs across APIs

hardened_inference.posterior_h0 and the selected mixture API accept finite non-negative host redshifts, while uncertainty_likelihoods.candidate_distance_mpc requires strictly positive redshifts. A zero-redshift entry maps to zero mean distance in the linear toy model, but is not a meaningful extragalactic host for this study. This interface difference should be resolved deliberately: either reject z=0 consistently in public cosmology-facing APIs or document why the lower-level generic likelihood permits it. Do not silently change the contract without a regression test.

### Selection assumes a shared noise model

The inspected beta(H0) calculation uses the same scalar Gaussian noise scale and fixed distance window for all mixture components. If catalogue and missing populations have different distance-uncertainty distributions, or detection depends on additional source properties, the selection denominator must integrate those variables and use the corresponding component-specific observation model. The present code does not support that generalization merely because it correctly sums over host weights.
