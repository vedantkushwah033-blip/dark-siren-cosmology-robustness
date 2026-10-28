# Mathematical specification for a minimal selection-aware toy model

**Status:** design specification; not implemented, benchmarked, or scientifically validated.  
**Date:** 2026-10-28  
**Purpose:** establish a small, falsifiable generative model and explicit correctness gates before implementing a catalogue-incompleteness benchmark.

## 1. Scope and non-claims

This is a controlled synthetic model for testing internal statistical consistency. It is not a gravitational-wave detector model, a reproduction of a published analysis, an analysis of LVK events, or a new measurement of $H_0$.

The first prototype uses the low-redshift approximation $D_L(z;H_0) = cz/H_0$. It must be labelled as such. A cosmological luminosity-distance relation, realistic detector selection, event-rate inference, and shared population-parameter uncertainty are separate future extensions.

## 2. Latent galaxy population and catalogue thinning

Let the redshift domain be $\mathcal Z=[z_{\min},z_{\max}]$, with $0<z_{\min}<z_{\max}$. Define a positive latent galaxy intensity

$$
\lambda(z;\psi)=A\left[1+\sum_{k=1}^{K} a_k
\exp\left(-\frac{(z-\mu_k)^2}{2\tau_k^2}\right)\right],
$$

where $A>0$, $a_k\geq 0$, and the feature locations and widths lie within declared bounds. This deliberately simple form produces a positive baseline plus overdensity features; an extension may add underdensities while constraining the total intensity to remain positive.

The catalogue selection probability $S_{\rm cat}(z)\in(0,1]$ is specified before simulation. For example, a logistic function can decrease smoothly with redshift. The observed catalogue is a thinned Poisson process with intensity

$$
\lambda_{\rm obs}(z)=\lambda(z;\psi)S_{\rm cat}(z).
$$

The simulation may retain the complete parent catalogue for scoring only. Non-oracle inference arms must not access the unseen redshifts or use the true latent intensity as an input. Catalogue selection and GW detection are different mechanisms and must be recorded separately.

For a fixed-size implementation, a conditional alternative is to draw a declared number of parent galaxies from normalized $p_{\rm gal}(z)\propto\lambda(z)$ and independently retain each with probability $S_{\rm cat}(z)$. The chosen Poisson or fixed-size convention must remain the same in the generating code, likelihood contract, and tests.

## 3. Merger-host redshift distribution

For this initial toy, let the merger-host density be

$$
p_{\rm host}(z\mid\psi,\rho)=
\frac{\lambda(z;\psi)\,\rho(z)}
{\int_{\mathcal Z}\lambda(u;\psi)\rho(u)\,du},
$$

where $\rho(z)>0$ is a declared merger-rate weight. Every simulation must state whether this weight is constant or redshift-dependent and record its normalization convention.

**Critical simplification:** in the first unit-tested toy, the same declared $p_{\rm host}$ governs both the event-host draw and the GW detection-efficiency integral. This is a deliberate shared-population assumption, not a claim that an event-specific candidate-host catalogue is generally the correct population distribution for selection normalization. If event-level host weights and population-level merger weights differ, the joint hierarchical model must be derived consistently before using a selection correction.

## 4. Event measurement and detection

For a generated host redshift $z_i$, draw a distance-like measurement

$$
d_i\mid z_i,H_0\sim\mathcal N\!\left(\frac{cz_i}{H_0},\sigma_i^2\right),
$$

where $c=299792.458\,\mathrm{km\,s^{-1}}$, $H_0$ is in $\mathrm{km\,s^{-1}\,Mpc^{-1}}$, $d_i$ is in Mpc, and $\sigma_i$ is in Mpc. This Gaussian low-redshift likelihood is a numerical toy, not a realistic GW distance posterior.

Use a declared threshold rule $d_i\geq d_{\min}$ for the first prototype. Conditional on redshift and $H_0$, the detection probability is

$$
P_{\rm det}(z,H_0)=
\Phi\!\left(\frac{cz/H_0-d_{\min}}{\sigma}\right),
$$

for the common $\sigma$ used by the initial prototype. The generating distribution and efficiency calculation must use the same noise convention. If the model truncates distances to $d>0$, both must include that truncation. They must never silently mix truncated and untruncated conventions.

Under the shared host-population assumption, the detection efficiency is

$$
\beta(H_0)=\int_{\mathcal Z}p_{\rm host}(z)P_{\rm det}(z,H_0)\,dz.
$$

The numerical integration domain, quadrature, normalization, and grid-refinement tolerances must be explicit. The efficiency must be in $(0,1]$ over the evaluated $H_0$ domain.

## 5. Conditional likelihood for detected events

For a detected event with observed $d_i\geq d_{\min}$, the conditional density in this toy is

$$
p(d_i\mid H_0,\mathrm{det})=
\frac{\int_{\mathcal Z}p_{\rm host}(z)
\,\mathcal N(d_i;cz/H_0,\sigma_i^2)\,dz}
{\beta(H_0)}.
$$

For $N$ conditionally independent detected events and fixed $N$,

$$
\log\mathcal L_{\rm det}(H_0)=
\sum_{i=1}^{N}\log\left[
\int_{\mathcal Z}p_{\rm host}(z)
\mathcal N(d_i;cz/H_0,\sigma^2)\,dz\right]
-N\log\beta(H_0)+C,
$$

where $C$ does not depend on $H_0$. This equation is the contract the toy implementation must satisfy. It assumes an unconditional event numerator and one consistent host-population distribution for the selection integral. Do not apply the $-N\log\beta$ term a second time if the event likelihood has already been conditioned on detection.

This is a fixed-detected-sample likelihood. It does not model event counts or infer a merger rate. If the number of detected events is part of the inference, a point-process/event-rate likelihood is required instead.

## 6. Catalogue-incompleteness inference arms

Implement these arms only after the generating model and likelihood contract above pass their correctness gates:

1. **Complete-catalogue oracle:** uses the full latent host distribution. It is a diagnostic reference only and not a deployable estimator.
2. **Homogeneous completion:** uses the observed catalogue plus a predeclared homogeneous missing-population model. It must not use the true missing redshifts or true latent intensity.
3. **Flexible reconstruction:** estimates a positive latent redshift density from the observed catalogue with the catalogue-selection model included. A plug-in estimate must be labelled as such; uncertainty in reconstructed density and selection must either be propagated or listed as limitations.
4. **Misspecification control:** deliberately perturb a stated selection or population assumption to measure sensitivity.

The same underlying parent population and paired GW observations should be used across methods within a replicate. Oracle results must be reported separately from non-oracle performance.

## 7. Required correctness gates before scaling

- **G1 — Positivity and normalization:** numerically verify that $p_{\rm host}$ integrates to one and remains nonnegative.
- **G2 — Catalogue thinning:** in a simple case, compare simulated retained counts and redshift distributions with the known thinning model, with Monte Carlo uncertainty stated.
- **G3 — Detection efficiency:** compare numerical $\beta(H_0)$ against an independent analytic result in a tractable special case and against empirical detection frequencies from repeated simulation.
- **G4 — Likelihood contract:** compare the implemented conditional likelihood against a separately coded reference calculation for a small deterministic example.
- **G5 — Selection sign and limiting behavior:** when $\beta(H_0)$ is constant over the tested domain, the selection correction must change the log likelihood only by an $H_0$-independent constant. A deliberately varying efficiency case must verify the sign and magnitude of the correction.
- **G6 — Oracle recovery:** test several interior generating values of $H_0$ and report bias, RMSE, interval coverage, interval width, and Monte Carlo uncertainty. Boundary truths must be reported separately.
- **G7 — Grid stability:** repeat selected cases with finer redshift and $H_0$ grids and independently reconstruct posterior summaries.
- **G8 — Reproducibility:** persist configuration, seeds, software versions, raw replicate-level data, checksums, failures, and a script that independently rebuilds summary tables.

If any gate fails, stop and diagnose it before a larger factorial experiment. Do not tune tolerances or remove difficult cells after inspecting favorable results.

## 8. Primary estimand for the later comparison

Before running the full benchmark, preregister a primary paired contrast. One option is

$$
\Delta_{\rm abs}=
\mathbb E\left[|\widehat H_{0,\rm flexible}-H_0|-
|\widehat H_{0,\rm homogeneous}-H_0|\right],
$$

averaged over a predeclared set of completeness and clustering conditions with an explicit cell-weighting rule. Negative values favor the flexible method on absolute error for that design; they do not prove general superiority. Report per-cell bias, RMSE, coverage, interval width, convergence failures, and compute cost alongside this contrast. A narrower interval alone is not evidence of a better estimator.

Select the replicate count from a pilot precision calculation for the primary contrast, not by repeatedly running until the result looks favorable.

## 9. Go/no-go decision

**Current decision: NO-GO for the full GP-versus-homogeneous performance benchmark.** This document defines the toy model and correctness gates; it does not mean those gates have passed. The next implementation step should be the smallest tractable numerical prototype and its independent tests, beginning with G1, G3, and G4. Only after those pass should the full catalogue-thinning experiment be built.
