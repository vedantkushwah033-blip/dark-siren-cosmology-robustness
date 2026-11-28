# Mathematical audit remediation: likelihood, design variance, and selection contract

**Status:** remediation specification; not a claim that every PR #20 finding is closed.  
**Branch:** \`research/math-audit-remediation-20261128\`  
**Parent commit:** \`a81caef8c25978c2c49ac634d9802f5e2d217902\`  
**Scope:** synthetic finite-parent catalogue inclusion and fixed-count dark-siren inference only.

## Executive status

| PR #20 finding | Current disposition | Closure evidence still required |
|---|---|---|
| Parent-law quadratic centering | Generator computes the parent-intensity moment; a high-resolution regression test checks its value and rate calibration. | Cross-component test proving generator, fit, prediction, manifest, and independent auditor use the same declared measure and value. |
| Non-retention likelihood | Five-arm workflow uses a validation-only modular inclusion posterior. It must not be described as the full observed-data likelihood. | Freeze modular conditioning/estimand, or separately derive and implement the full latent non-retention likelihood. |
| SRSWOR finite-population uncertainty | A standalone first-order logistic sandwich covariance primitive and algebraic unit tests exist. It is not integrated or Monte-Carlo calibrated. | Integrate fixed-coordinate handling and validate against repeated SRSWOR draws on fixed finite frames. |
| GW detection direction | Current generator and KDE likelihood both use the upper-tail rule, measured distance >= threshold. | Exact analytic selection-efficiency regression test and generator trace checks on the settled exact head. |
| Complete attempt ledger | Prepared-bundle lifecycle work exists, but scenario-level reference failures still need an end-to-end independent ledger audit. | Reconstruct every planned key, including reference failures and not-run attempts. |

## 1. Parent-measure centering

The parent frame is sampled from the baseline galaxy intensity

    lambda_g(z) = 1 + 2.2 exp[-(z-.105)^2/(2*.025^2)]
                   + 1.6 exp[-(z-.275)^2/(2*.040^2)].

For x = (z - z_c)/s_z, the degree-two basis is

    [1, x, x^2 - c_g],
    c_g = integral[x(z)^2 lambda_g(z) dz] / integral[lambda_g(z) dz].

This is a property of the declared parent measure, not a validation-sample moment. The current generator exposes the calculation through \`quadratic_logit_parent_measure\` and uses it in quadratic inclusion-probability calibration. The five-arm inference representation carries an explicit \`quadratic_center\`.

The new regression test checks the high-resolution value c_g approximately 0.284500557016 and verifies the quadratic intercept calibration. It does not yet prove that every inference and independent-auditor component uses the same value. The independent auditor must recompute the integral rather than import the producer's result.

A finite-frame alternative c_N = N^(-1) sum_i x_i^2 is a different target and must be labelled separately.

## 2. Inclusion likelihood: choose and name the analysis

Let R_i be catalogue retention, Z_i parent redshift, and S_i membership in an SRSWOR validation sample of fixed size n, drawn independently of R conditional on the parent frame.

The validation-only likelihood is

    L_V(theta) = product over i:S_i=1 of
                 pi(Z_i;theta)^R_i [1-pi(Z_i;theta)]^(1-R_i).

Using p(theta | V) proportional to L_V(theta) p(theta), while keeping the GW/host likelihood outside this update, is a modular cut. It is a legitimate specified analysis, but it is not the full likelihood of all catalogue and validation observations.

If the intended analysis instead conditions on a complete finite parent frame, retention contributes

    L_R(theta | Z,R) = product over i=1,...,N of
                       pi(z_i;theta)^R_i [1-pi(z_i;theta)]^(1-R_i).

If redshifts of units outside the observed union A = C union V are latent under an iid superpopulation model, each such unit contributes

    integral over Z of [1-pi(z;theta)] g(z) dz.

Overlapping catalogue and validation units contribute each observed retention label once. The fixed-size SRSWOR design factor is parameter-free conditional on the frame, but belongs in the sampling-design statement.

**Decision for the current five-arm project:** retain the validation-only modular cut unless a separate, fully specified joint observed-data model is approved. All documentation and outputs must call it a modular validation posterior. Do not describe it as the complete observed-data likelihood. Do not insert latent non-retention terms into the modular cut; doing so changes the model and estimand.

## 3. SRSWOR design variance for a logistic estimating equation

The inverse observed-information matrix (X_V' W_V X_V)^(-1) is a model-based likelihood covariance. It is not automatically the design variance of a logistic estimator from an SRSWOR sample.

For a fixed finite frame with values (z_i,R_i), define the finite-frame pseudo-true parameter theta_N by

    mean_N u_i(theta_N) = 0,
    u_i(theta) = x_i [R_i - pi(z_i;theta)].

Let J_N be the finite-frame average derivative of u_i and let S_u,N be the finite-frame covariance of the individual score contributions. Under regularity conditions and a regular interior root, first-order SRSWOR linearization gives

    Var_d(theta_hat_V | Z,R)
      approximately J_N^(-1) [(1-n/N)/n * S_u,N] J_N^(-T).

The new \`design_based_logistic.py\` module implements a plug-in sample version using:
- sample-average logistic sensitivity A_hat = mean[x_i x_i' p_i(1-p_i)];
- sample covariance of individual score contributions;
- the explicit SRSWOR factor (1 - n/N)/n.

It is intentionally separate from the model-based inverse-Hessian covariance.

Cautions:
- The target is the finite-frame score root theta_N, not automatically the superpopulation coefficient.
- The design matrix must use fixed, predeclared coordinates. If centering/scaling is estimated from the validation sample, its uncertainty must be propagated or coordinates must be fixed.
- Separation, active coefficient bounds, rank deficiency, and small samples can invalidate the regular linearization.
- This conditional design variance excludes parent-frame generation and Bernoulli-retention variation. Those require outer repeated-scenario experiments.
- The modular Bayesian coefficient posterior and this design-based covariance are different uncertainty constructions. Do not multiply one by the other without a coherent derivation.

**Current disposition:** algebraic unit tests cover the explicit matrix formula, census zero-variance limit, and invalid inputs. The primitive is not integrated into the five-arm posterior and has not been calibrated against repeated SRSWOR draws on fixed frames. No design-based uncertainty-propagation claim is authorized yet.

## 4. GW detection-direction contract

The generator accepts a candidate when

    D_obs >= d_th,
    D_obs ~ Normal(mu(z,H0), sigma^2),
    mu(z,H0) = c z / H0.

Therefore

    P(detected | z,H0)
      = P(D_obs >= d_th)
      = Phi((mu(z,H0) - d_th)/sigma).

The implementation in \`five_arm_kde_likelihood.py\` evaluates this as \`exp(log_ndtr((means - threshold) / sigma))\`, the correct Gaussian upper-tail probability for the generator's rule. The opposite lower-tail indicator in the older PR #20 prose is stale for this current code path.

A new regression test compares the implemented efficiency against direct numerical integration of the Gaussian upper tail and checks that the opposite tail gives a materially different answer. This test still needs to pass on the settled exact head.

## 5. Closure gates

Task 1 may be marked complete only when:

1. The primary inclusion update is explicitly frozen as a validation-only modular cut, or a separately approved full observed-data likelihood is derived and implemented.
2. Parent-measure centering is independently recomputed and compared across generator, inference, manifest, and auditor.
3. The design-based variance primitive is integrated and validated against repeated SRSWOR Monte Carlo on fixed finite frames, or the scientific uncertainty claim is explicitly narrowed and any design-based FPC claim is removed.
4. The analytic detection-efficiency regression test passes on the exact head.
5. The full planned-attempt ledger is independently reconstructed, including reference failures and not-run attempts.
6. Exact-head Python 3.10-3.12 CI passes on the remediation branch.

No confirmatory freeze is authorized by this note.
