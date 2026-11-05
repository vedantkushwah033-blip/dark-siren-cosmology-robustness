# Independent mathematical audit and corrected research design

**Status:** design audit only. No confirmatory campaign has been run.  
**Scope:** synthetic finite-parent host reconstruction; no real LVK data or observational H0 claim.

## Executive decision

The research question is defensible: quantify how finite validation samples, propensity misspecification, and inverse-weight instability affect reconstructed host densities and downstream fixed-count dark-siren inference. The supplied consortium specification is not approved unchanged. Four issues are critical: the parent-law centering constant, the observed-data likelihood, design-based uncertainty, and scenario-level attempt accounting. A further contract mismatch exists between the GW generator's detection direction and the supplied likelihood formula.

## 1. Target-measure hierarchy

The live generator in `src/dark_siren_robustness/generative_selection.py` samples parent redshifts from a baseline galaxy intensity
[
\lambda_g(z)=1+2.2e^{-\frac12((z-.105)/.025)^2}
+1.6e^{-\frac12((z-.275)/.040)^2},
\quad g(z)=\lambda_g(z)/\int_{\mathcal Z}\lambda_g(u)du.
]
It then applies a separate host/merger weight
[
r(z)=1+1.5(z-z_{\min})/(z_{\max}-z_{\min})
]
when drawing GW hosts. Therefore the superpopulation host density is
[
q(z)=\frac{g(z)r(z)}{\int_{\mathcal Z}g(u)r(u)du}.
]
The catalogue KDE/IPW stage estimates the baseline galaxy density (g); multiplying by (r) once constructs (q). This is not double weighting under the actual generator. The consortium's prose saying that the parent is already sampled proportional to comoving volume times merger rate describes a different DGP and must not be mixed into this code path.

Keep five functions distinct in code and manuscript: baseline galaxy density (g), catalogue inclusion (pi), host/merger weight (r), host density (q), and GW detection efficiency (p_{det}(z,H_0)).

## 2. Quadratic inclusion centering

The specification defines (c=E_g[x(Z)^2]), but asserts (c=1/3). That equality requires (x) uniform on ([-1,1]), which is false for the current non-uniform (g(z)).

For the current superpopulation DGP use the deterministic constant
[
c_g=\frac{\int_{z_{min}}^{z_{max}}x(z)^2\lambda_g(z)dz}
{\int_{z_{min}}^{z_{max}}\lambda_g(z)dz}.
]
Numerical integration of the current declared λ_g over [0.01, 0.40] gives approximately **0.28450** (the exact stored value should be generated deterministically at the declared quadrature resolution). Calculate it from the declared generator law, store it in each scenario manifest, and independently recompute it in the auditor. Do not estimate it from each validation sample. A finite-frame alternative (c_N=N^{-1}\sum_i x_i^2) is a different DGP and must be labelled separately.

The current producer and independent auditor both hard-code (1/3) in the quadratic term. Thus the document's claimed parent-law centering is not implemented.

## 3. Finite-frame and observed-data likelihood

Conditional on a complete finite frame (Z=(z_1,\ldots,z_N)), retention contributes
[
L_R(\theta)=\prod_{i=1}^N\pi(z_i;\theta)^{R_i}
[1-\pi(z_i;\theta)]^{1-R_i}.
]
An SRSWOR validation sample (S) of fixed size (n), independent of retention conditional on (Z), has design probability
[
P(S=s\mid Z,R)=\binom{N}{n}^{-1}1\{\sum_i s_i=n\}.
]
The design factor is parameter-free for the propensity likelihood but must be retained in the sampling-design description.

The validation-only likelihood
[
L_V(\theta)=\prod_{i:S_i=1}\pi(z_i;\theta)^{R_i}
[1-\pi(z_i;\theta)]^{1-R_i}
]
is the likelihood used by a **modular cut**. It is not the full likelihood of all observed catalogue and validation data.

Let (A=\mathcal C\cup\mathcal V) and (U=\mathcal P\setminus A). For (i\in\mathcal C\cap\mathcal V), retention (R_i=1) is observed and contributes once. For (i\in\mathcal V\setminus\mathcal C), (R_i=0) is observed. For (i\in\mathcal C\setminus\mathcal V), (R_i=1) is observed. For (i\in U), catalogue non-retention is known even though (z_i) is latent. Under independent superpopulation redshifts, each such latent unit contributes an integrated factor
[
\int_{\mathcal Z}[1-\pi(z;\theta)]g(z)dz.
]
The supplied observed-data expression omitted these non-retention factors for units outside both samples. Further conditioning on the observed union size/identities must be stated explicitly if used. The primary project can remain modular; it must not describe (L_V) as the complete observed-data likelihood.

## 4. Design-based uncertainty: no automatic scalar FPC

The SRSWOR sample-mean correction $1-n/N$ does not imply that logistic coefficient covariance equals the ordinary inverse-Hessian covariance multiplied by this scalar. Define the finite-frame pseudo-true coefficient $\theta_N$ by the full-frame score equation
\[
\bar u_N(\theta_N)=N^{-1}\sum_{i=1}^N u_i(\theta_N)=0,
\]
where $u_i$ is the Bernoulli-logistic score contribution for the realized pair $(z_i,R_i)$. Let $J_N=N^{-1}\sum_i \partial u_i(\theta_N)/\partial\theta^T$ and
\[
S_{u,N}=\frac{1}{N-1}\sum_{i=1}^N
[u_i(\theta_N)-\bar u_N][u_i(\theta_N)-\bar u_N]^T.
\]
For an SRSWOR validation sample of size $n$, the first-order design linearization of the sample-score root gives
\[
\operatorname{Var}_d(\hat\theta_V\mid Z,R)
\approx J_N^{-1}
\left(\frac{1-n/N}{n}S_{u,N}\right)
J_N^{-T}.
\]
This is conditional on the realized finite frame and retention indicators, and assumes a regular interior root and nonsingular sensitivity. With bounded coefficients, separation, or active constraints, the unconstrained linearization may fail and boundary-aware inference is needed. In practice $J_N$ and $S_{u,N}$ must be estimated or independently benchmarked; the formula is not a license to multiply an ordinary Hessian covariance by an FPC without checking score scaling and target definition.

This design variance addresses only validation-sampling variation conditional on $(Z,R)$. Parent-frame generation and Bernoulli-retention randomness are separate outer sources and must be assessed by repeated independent scenarios if the scientific estimand averages over them. Until this derivation is implemented and calibrated, Method C should be called a conditional sensitivity integration over a specified propensity-parameter distribution, not validated design-based uncertainty propagation.

## 5. IPW, KDE, and host-density perturbation

For (hat\pi(z)=\pi(z)[1+\delta(z)]), Hájek normalization gives the first-order baseline-density perturbation
[
\frac{\hat g(z)-g(z)}{g(z)}
\approx-[\delta(z)-E_g\{\delta(Z)\}].
]
This concerns (g) before applying (r). The finite-sample KDE adds catalogue variation, weight concentration, adaptive bandwidth, reflection, and normalization effects. The current estimator uses weighted Scott bandwidth (h=SD_w ESS^{-1/5}), reflects observations at both boundaries, applies (r) on the grid, then renormalizes. These exact choices define the estimator and must be audited.

Pure IPW must not also add a missing-host mixture fraction unless a separately defined residual missing-population component is introduced. For the current pure-IPW arm, enforce (f=0) and test it.

## 6. Fixed-count GW likelihood: detection-direction contract

The fixed-count likelihood is
[
L(H_0;\theta)=
\frac{\prod_{j=1}^{n_{det}}\int q_\theta(z)p(d_j\mid z,H_0)dz}
{\beta(H_0;\theta)^{n_{det}}},
\quad
\beta(H_0;\theta)=\int q_\theta(z)p_{det}(z,H_0)dz.
]
The modular likelihood average is (int L(H_0;\theta)p(\theta\mid\mathcal V)d\theta); do not average numerator and denominator separately.

**Critical mismatch:** the live generator accepts an event when measured distance is greater than or equal to the threshold. The supplied consortium formula states (1(d\le d_{thresh})). These are opposite selection directions. Production likelihood, independent auditor, and manuscript must match the generator or the DGP must be deliberately changed everywhere. This is a release-blocking contract test.

## 7. Optimizer oracle contract

Production propensity fitting currently uses bounded L-BFGS-B; the independent V2 audit uses a separately coded bounded SLSQP solve and projected-gradient check. This provides algorithmic independence but is not the proposed L-BFGS-B oracle. The reference must optimize the exact same objective, including ridge penalty (if any), penalty scaling, intercept treatment, bounds, and gradient. Record objective gap, feasibility, projected-gradient/KKT residual, active-bound agreement, rank/conditioning, optimizer status, and failure category.

Current production/auditor KKT tolerance is (10^{-5}). The proposed (10^{-6}) standard has not been established. Complete separation with finite coefficient bounds is a boundary-constrained optimum, not an unconstrained finite MLE.

## 8. Scenario lifecycle and complete attempt ledger

The live `finite_parent_sampling.py` calculates oracle, catalogue-only, and full-parent reference posteriors before entering the per-validation-draw exception handler. A scenario-level reference failure can abort before any planned draw rows are emitted.

Required scenario states:
- `SCENARIO_GENERATION_FAILED`
- `ORACLE_REFERENCE_FAILED`
- `CATALOGUE_REFERENCE_FAILED`
- `FINITE_PARENT_REFERENCE_FAILED`
- `SCENARIO_READY`
- `ARTIFACT_SERIALIZATION_FAILED`

Every planned key ((scenario,draw,validation\_size)) must appear exactly once. Attempt statuses must be `SUCCESS`, `FAILED`, `NON_EVALUABLE`, or `NOT_RUN). A required reference failure makes downstream comparisons non-evaluable; an attempt never started is not-run. Do not collapse both into failed.

Reconcile
[
N_{planned}=N_{success}+N_{failed}+N_{non\text{-}evaluable}+N_{not\text{-}run}.
]
The independent auditor must reconstruct the planned key set, reject missing/duplicate keys, verify scenario state and failure stage, and replay individual attempts whenever prerequisites exist.

## 9. Experimental hierarchy

Parent scenarios are the top-level independent units for superpopulation performance. Validation draws are nested within a scenario and share its catalogue and GW data. Report pooled descriptive summaries and scenario-balanced summaries, between-scenario variability, within-scenario validation variability, all four status fractions, and scenario-level Monte Carlo uncertainty.

The proposed (R=2000) is a future design target, not proof of adequacy. Sequence: deterministic tests; smoke runs; numerical stress tests; pilot with scenario-level precision estimates; then pre-registered confirmatory design. Factor inclusion curvature, validation size, fitted model family, and solver version. Oracle IPW is diagnostic only.

## 10. Release gates

| Gate | Required evidence | Status |
|---|---|---|
| Parent-law centering | deterministic (c_g), manifest, independent recomputation | open |
| Detection direction | generator and likelihood use same threshold inequality | open |
| Target measure | explicit (g,\pi,r,q,p_{det}) contracts and tests | partly clarified |
| Attempt ledger | complete keys under injected scenario-reference failures | open |
| Optimizer oracle | identical objective and independent KKT/active-set checks | open |
| Design uncertainty | derived SRSWOR covariance and calibrated coverage | open |
| Numerical convergence | posterior/density stability under grid and bandwidth refinement | open |
| Pilot precision | scenario-level MCSE and failure-rate precision | not started |
| Confirmatory preregistration | frozen hypotheses, estimands, exclusions, seeds, revision | not started |

## Permitted current claim

“We are developing a synthetic finite-parent benchmark to study the effects of validation-sample size, propensity-model misspecification, and inverse-weight instability on reconstructed host densities and downstream fixed-count dark-siren inference.”

Do not claim validated uncertainty propagation, optimizer passage at (10^{-6}), confirmatory power, universal estimator superiority, real-data validation, or an observational (H_0) measurement.

**Decision:** continue software and mathematical development; no 2,000-replicate confirmatory campaign yet.
