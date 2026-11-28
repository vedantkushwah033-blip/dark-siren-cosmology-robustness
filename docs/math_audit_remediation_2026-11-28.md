# Mathematical audit remediation: likelihood, design variance, and selection contract

**Status:** remediation specification; not a claim that every PR #20 finding is closed.  
**Branch:** `research/math-audit-remediation-20261128`  
**Parent commit:** `a81caef8c25978c2c49ac634d9802f5e2d217902`  
**Scope:** synthetic finite-parent catalogue inclusion and fixed-count dark-siren inference only.

## Executive status

This note separates findings that can be reconciled from findings that still require mathematical implementation and independent validation.

| PR #20 finding | Current disposition | Closure evidence still required |
|---|---|---|
| Parent-law quadratic centering | The current generator computes the parent-intensity moment using `quadratic_logit_parent_measure`; the later centered-basis implementation stores and propagates `quadratic_center`. | Exact-head test proving generator, fit, prediction, and independent auditor all use the same declared measure and value. |
| Non-retention likelihood | The five-arm workflow uses a validation-only modular inclusion posterior. It must be described as a cut/modular analysis, not as the full observed-data likelihood. | Freeze the modular estimand and document its conditioning; if a full joint observed-data analysis is desired, derive and implement the latent non-retention terms separately. |
| SRSWOR finite-population uncertainty | Not closed by the ordinary inverse-Hessian covariance or by a scalar FPC multiplier. | Implement and independently validate a design-linearized covariance estimator, or explicitly limit the analysis to model-based/modular posterior uncertainty and report design variation by repeated SRSWOR simulation. |
| GW detection direction | The current generator and KDE likelihood are consistent with acceptance when measured distance is greater than or equal to the threshold. The old PR #20 prose/formula is stale for this code path. | Exact analytic regression test for the survival probability and retain generator trace checks. |
| Complete attempt ledger | The prepared-bundle lifecycle has made progress, but the earlier scenario-level reference-failure gap must be checked in the active end-to-end coordinator. | Independent audit of every planned key, including scenario failures and not-run states. |

## 1. Parent-measure centering

The current parent frame is sampled from the baseline galaxy intensity
[
lambda_g(z)=1+2.2exp[-(z-.105)^2/(2cdot .025^2)]
+1.6exp[-(z-.275)^2/(2cdot .040^2)].
]
For (x=(z-z_c)/s_z), the degree-two basis is
[
[1, x, x^2-c_g],qquad
c_g=rac{int x(z)^2lambda_g(z),dz}{intlambda_g(z),dz}.
]
This is a property of the declared parent measure, not a sample moment to be re-estimated from each validation sample. The current generator exposes this calculation through `quadratic_logit_parent_measure` and uses the result when calibrating quadratic inclusion probabilities. The five-arm inference representation carries an explicit `quadratic_center`.

Required invariant: the scenario manifest's center, the generator's center, the fit's basis center, prediction basis, and independent auditor's reconstructed center must agree to a declared numerical tolerance. The auditor must calculate the integral independently, not import the producer's computed constant.

A finite-frame alternative (c_N=N^{-1}sum_i x_i^2) is a different estimand/design and must not be silently substituted for (c_g).

## 2. Inclusion likelihood: choose and name the analysis

Let (R_i) denote catalogue retention, (Z_i) parent redshift, and (S_i) membership in an SRSWOR validation sample of fixed size (n), drawn independently of (R) conditional on the parent frame.

The validation-only likelihood is
[
L_V(	heta)=prod_{i:S_i=1}
pi(Z_i;	heta)^{R_i}{1-pi(Z_i;	heta)}^{1-R_i}.
]
Using (p(	hetamid V)propto L_V(	heta)p(	heta)) while keeping the GW/host likelihood outside this update is a **modular cut**. It is a legitimate specified analysis, but it is not the full likelihood of all catalogue and validation observations.

If the intended analysis instead conditions on the complete finite parent frame, retention contributes
[
L_R(	hetamid Z,R)=prod_{i=1}^{N}
pi(z_i;	heta)^{R_i}{1-pi(z_i;	heta)}^{1-R_i}.
]
If redshifts of units outside the observed union (A=mathcal Ccupmathcal V) are latent under an iid superpopulation model, each such unit contributes the integrated factor
[
int_{mathcal Z}{1-pi(z;	heta)}g(z),dz.
]
Overlaps must contribute each observed retention label once. The fixed-size SRSWOR design factor is parameter-free conditional on the frame but belongs in the sampling-design statement.

**Decision for the current five-arm project:** retain the validation-only modular cut unless a separate, fully specified joint observed-data model is approved. All documentation and output labels must call it a modular validation posterior. Do not describe it as the complete observed-data likelihood. Do not add latent non-retention factors to the modular cut; doing so would change the estimand and model.

## 3. SRSWOR design variance for a logistic estimating equation

The inverse observed-information matrix
[
{X_V^	op W_V X_V}^{-1}
]
is a model-based likelihood covariance. It is not automatically the design variance of a logistic estimator from an SRSWOR sample.

For a fixed finite frame with values ((z_i,R_i)), define the finite-frame pseudo-true parameter (	heta_N) by
[
ar u_N(	heta_N)=N^{-1}sum_{i=1}^N u_i(	heta_N)=0,
]
where (u_i(	heta)=x_i{R_i-pi(z_i;	heta)}). Let
[
J_N=N^{-1}sum_{i=1}^Nrac{partial u_i(	heta_N)}{partial	heta^	op},
quad
S_{u,N}=rac{1}{N-1}sum_{i=1}^N
(u_i-ar u_N)(u_i-ar u_N)^	op.
]
Under regularity conditions and a regular interior root, the first-order SRSWOR linearization is
[
operatorname{Var}_d(hat	heta_Vmid Z,R)
approx
J_N^{-1}left{rac{1-n/N}{n}S_{u,N}ight}J_N^{-T}.
]
A practical estimator substitutes sample estimates of the sensitivity and score covariance. This is a sandwich design variance, not a scalar multiplier applied to the ordinary inverse Hessian.

Cautions:
- The target is the finite-frame score root (	heta_N), not automatically the superpopulation coefficient.
- If covariate centering/scaling is estimated from the validation sample, its randomness must be included or the coordinates must be fixed by the design.
- Separation, boundary solutions, rank deficiency, and small (n) can invalidate the regular linearization.
- This conditional design variance excludes parent-frame generation and Bernoulli-retention variation. Those require outer repeated-scenario experiments.
- The modular Bayesian coefficient posterior and this design-based covariance are different uncertainty constructions. Do not multiply one by the other without a coherent derivation.

**Current disposition:** no design-based FPC correction is claimed. Until an independently tested design-variance implementation exists, report the modular posterior as model-based coefficient uncertainty and evaluate SRSWOR variability through repeated finite-frame simulations. Keep design-based uncertainty propagation as an open gate.

## 4. GW detection-direction contract

The generator accepts a candidate when
[
D_{m obs}ge d_{m th},qquad
D_{m obs}sim N(mu(z,H_0),sigma^2),quad
mu(z,H_0)=cz/H_0.
]
Therefore
[
P(mathrm{detected}mid z,H_0)
=P(D_{m obs}ge d_{m th})
=Phi!left(rac{mu(z,H_0)-d_{m th}}{sigma}ight).
]
The implementation in `five_arm_kde_likelihood.py` evaluates this as
`exp(log_ndtr((means - threshold) / sigma))`, which is the correct Gaussian survival probability for the generator's upper-tail rule.

The opposite lower-tail indicator appearing in the older audit prose is not the current implementation contract. Preserve the generator's `>=` rule and test the analytic efficiency directly. If the DGP is changed later, change the generator, likelihood, auditor, protocol, and tests together.

## 5. Closure gates

Task 1 may be marked complete only when all of the following are true:

1. A frozen document states that the primary inclusion update is a validation-only modular cut, or a separately approved full observed-data likelihood has been derived.
2. Parent-measure centering is independently recomputed and compared across generator, inference, and auditor.
3. Design-based variance is implemented and validated against exact SRSWOR Monte Carlo on fixed finite frames, or the project explicitly narrows its uncertainty claim and removes any design-based FPC claim.
4. An analytic detection-efficiency regression test passes for multiple (z,H_0,sigma,d_{m th}) combinations.
5. The full planned-attempt ledger is independently reconstructed, including reference failures and not-run attempts.
6. Exact-head Python 3.10–3.12 CI passes on the remediation branch.

No confirmatory freeze is authorized by this note.
