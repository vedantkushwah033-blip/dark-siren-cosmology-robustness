# Statistical estimands and paired analysis plan for the density stage

Status: working methods note, 2026-11-14. This document defines the density-stage estimands and analysis safeguards. It is not a report of confirmatory results or a frozen bandwidth-selection protocol.

## 1. Data-generating structure

For each synthetic replicate, generate a parent galaxy sample with redshifts (Z_i) from the declared parent redshift distribution. Catalogue membership is a separate Bernoulli thinning process:

[
R_imid Z_i sim operatorname{Bernoulli}{pi(Z_i)},qquad
pi(z)=c+(1-c)left[1+exp{(z-z_0)/s}ight]^{-1},
]

where (cin{0,0.10,0.25}), (z_0=0.22), and (s=0.05). The validation sample is an independently seeded sample from a finite validation frame; validation sizes 100, 300 and 1,000 are nested prefixes.

At (c=0), (operatorname{logit}pi(z)) is linear in redshift, so the fitted degree-one logistic model is correctly specified. At (c>0), the floor transformation generally makes (operatorname{logit}pi(z)) nonlinear, so the degree-one model is misspecified. This is a specific floor-based misspecification family. It is not the separate quadratic-logit curvature family (operatorname{expit}(alpha+eta x+gamma(x^2-c_x))).

Catalogue inclusion and gravitational-wave detection selection are distinct mechanisms. This density stage studies catalogue inclusion only; it does not estimate a GW detection efficiency or an (H_0) posterior.

## 2. Target density and base-measure accounting

The latent parent intensity used in this benchmark is

[
g(z)=1+2.2expleft[-rac{(z-0.105)^2}{2(0.025)^2}ight]
       +1.6expleft[-rac{(z-0.275)^2}{2(0.040)^2}ight],
]

and the base-population/merger measure is

[
b(z)=1+1.5rac{z-z_{min}}{z_{max}-z_{min}},
qquad zin[0.01,0.40].
]

The target is

[
f_0(z)=rac{g(z)b(z)}{int_{z_{min}}^{z_{max}}g(u)b(u),du}.
]

The factor (b(z)) is applied exactly once. In the density runner it is incorporated into the raw sample weights. The post-KDE multiplier is the identity. Applying (b(z)) again after KDE would instead target a density proportional to (g(z)b(z)^2), which is not the declared target.

## 3. Five matched density arms

All arms within a replicate use the same parent realization, observed catalogue, validation frame, and nested validation records.

Let (b_i=b(Z_i)) and let (widehatpi_i(	heta)) denote the degree-one logistic inclusion probability at a catalogue redshift.

- **A — catalogue-only baseline:** catalogue redshifts with raw weights (b_i). This intentionally omits the inverse-inclusion correction.
- **B — parent oracle:** complete parent redshift sample with raw weights (b_i).
- **C — known-probability IPW:** catalogue redshifts with raw weights (b_i/pi(Z_i)).
- **D — MAP plug-in IPW:** catalogue redshifts with raw weights (b_i/widehatpi_i(widehat	heta_{mathrm{MAP}})).
- **E — uncertainty-integrated IPW:** for posterior quadrature nodes (	heta_k) and normalized posterior weights (q_k), use node-specific raw weights (b_i/widehatpi_i(	heta_k)).

No clipping or winsorization is part of the primary arms. Numerical underflow/overflow or failed fitting is explicitly recorded, not silently repaired.

## 4. Density estimator and bandwidth estimands

For each arm and common physical redshift bandwidth (h), construct a reflected weighted Gaussian KDE on the declared redshift interval and normalize the resulting density after the base measure is applied. A single physical width (h) is used for A–D and every E quadrature node in a matched comparison.

The realized weighted-Scott widths are recorded as diagnostics to characterize how a common Scott multiplier would imply different physical smoothing scales across arms and validation sizes. They are not substituted for the common-physical-width sensitivity grid.

The candidate bandwidth set is a numerical sensitivity grid, not an estimator parameter selected separately for each arm. Any eventual common-width selection must use a frozen equal-cell objective and an independent holdout set. A boundary minimum triggers a documented grid expansion; it is not frozen as an optimum.

## 5. Arm E has two distinct loss estimands

For each posterior node (k), let (widehat f_{E,k,h}) be its normalized weighted KDE. The posterior-mean density is

[
overline f_{E,h}(z)=sum_k q_kwidehat f_{E,k,h}(z).
]

The **primary E estimand** is the loss of this posterior-mean density:

[
L(overline f_{E,h},f_0).
]

The **secondary E estimand** is the posterior expected node-specific loss:

[
sum_k q_k L(widehat f_{E,k,h},f_0).
]

These are not interchangeable. For convex losses such as integrated absolute error, integrated squared error, and boundary-restricted integrated absolute error, Jensen's inequality gives

[
L(overline f_{E,h},f_0)
leq
sum_k q_k L(widehat f_{E,k,h},f_0).
]

The posterior-mean density is a mixture of normalized densities and is itself normalized. It is not associated with a single sample-weight vector, so no Kish ESS is assigned to that mixture density or to the posterior expected node-loss summary. Node-specific E rows retain their own ESS.

## 6. Primary and secondary loss metrics

On a fixed redshift grid, the primary density metric is integrated absolute error:

[
operatorname{IAE}(widehat f,f_0)
=int|widehat f(z)-f_0(z)|,dz.
]

Secondary metrics are integrated squared error,

[
operatorname{ISE}(widehat f,f_0)
=int{widehat f(z)-f_0(z)}^2,dz,
]

and boundary IAE over the predeclared lower and upper edge regions. Grid-resolution sensitivity must be checked before treating these numerical integrals as stable.

Weight diagnostics include Kish ESS, entropy ESS, maximum normalized weight, and top 1%/5% weight mass. These diagnostics describe weight concentration; they do not by themselves establish bias, posterior calibration, or failure of an (H_0) posterior.

## 7. Paired repeated-sampling analysis

For each replicate, all arms and candidate widths are evaluated on matched inputs. Candidate-width contrasts are computed within replicate before averaging. Report:
- mean paired difference;
- replicate standard deviation and Monte Carlo standard error;
- uncertainty interval;
- successful/evaluable replicate count;
- explicit failure count and reason.

Do not treat candidate-width scores from the same training data as independent. Do not use an unpaired standard error when the design is paired.

For primary IAE, the distance between two normalized densities is bounded above by 2. If a preregistered training/holdout design uses failure-inclusive selection, a failed or non-evaluable primary cell can therefore receive IAE = 2 and remain in the fixed denominator. This penalty is a conservative design rule, not an imputed scientific estimate. It must not be transferred to ISE, which has no analogous universal bound here.

## 8. Inclusion-posterior numerical validation

The degree-one logistic model uses an explicit proper Gaussian prior with coefficient mean (0,0) and standard deviations (5,2.5) in the declared centered/scaled coordinate system. The Laplace Gaussian is an importance proposal; normalized tensor Gauss-Hermite importance weights target the proper posterior.

The production optimizer must finish with score infinity norm at most (10^{-6}) and positive-definite observed information. A deterministic exact-Hessian Newton refinement may follow a trust-exact stop, but it may only be accepted when the actual score improves and the final threshold is met. The independent raw-input auditor separately reconstructs the mode, covariance, nodes and normalized weights from observable validation records.

Quadrature-order stability must be checked on posterior summaries and downstream likelihoods separately. Posterior quadrature convergence does not automatically imply convergence of a downstream (H_0) posterior.

## 9. Gates before an H0 stage

The density stage must first pass:
1. exact-head unit and regression tests;
2. 20/20 replicate completion in each required floor cell, or a separately preregistered failure-inclusive design;
3. independent validation-posterior reconstruction for every validation sample;
4. independent raw KDE, target, IAE/ISE/boundary IAE, and weight-diagnostic reconstruction;
5. independent outer manifest, seed schedule, lifecycle ledger, and aggregate reconstruction;
6. a documented bandwidth-grid decision and a disjoint training/holdout selection design.

Only after those gates should a separate H0 stage be launched, with a frozen density-stage bandwidth, shared inclusion coefficient node across all GW events, the complete joint selected likelihood marginalized once, and separate GW detection-selection normalization. H0 grid convergence, quadrature convergence, repeated-sampling calibration, and SBC remain distinct validation requirements.

## 10. Scope and interpretation limits

This is a synthetic low-redshift Euclidean methodology benchmark. It is not LVK/GWTC/public-data analysis, an observational H0 measurement, validation on a real galaxy survey, or proof that one estimator universally dominates another. Density recovery alone cannot establish H0 posterior calibration.
