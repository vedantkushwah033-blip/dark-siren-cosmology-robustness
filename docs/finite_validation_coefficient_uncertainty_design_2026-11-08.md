# Finite-validation coefficient-uncertainty propagation: design note

**Status:** design only; no implementation or numerical result is claimed.  
**Scope:** synthetic low-redshift Euclidean benchmark. This note does not authorize an observational H0 claim.

## 1. Why this stage is needed

The current finite-validation fitter returns an unpenalized logistic maximum-likelihood estimate and the inverse observed-information matrix. That matrix is a local asymptotic covariance estimate. It is **not** a Bayesian posterior, and drawing coefficients from a normal distribution with that covariance must not be described as exact Bayesian marginalization.

The density-recovery benchmark currently conditions on one fitted inclusion curve per validation sample. It can measure plug-in density error, but it does not yet propagate finite-validation coefficient uncertainty into the host-redshift density or an H0 likelihood.

The next question is therefore conditional and testable:

> How much additional density and posterior uncertainty arises when uncertainty in the fitted inclusion curve is propagated, and does that propagation improve calibration under correctly specified and misspecified inclusion models?

Uncertainty propagation cannot by itself repair a misspecified inclusion model. Under misspecification, a posterior may concentrate around a pseudo-true curve that remains systematically wrong.

## 2. Validation likelihood and proper coefficient posterior

For validation observations \\(V=\\{(z_i,R_i)\\}_{i=1}^{n_V}\\), use the declared logistic model

\\[
\\pi(z;\\theta)=\\operatorname{logit}^{-1}
\\left(\\sum_{j=0}^{d}\\theta_j x(z)^j\\right),
\\qquad
x(z)=\\frac{z-\\bar z_V}{s_V}.
\\]

Conditional on the validation redshifts and the predeclared SRS validation design, the Bernoulli-label likelihood is

\\[
L_V(\\theta)
=\\prod_{i\\in V}
\\pi(z_i;\\theta)^{R_i}
[1-\\pi(z_i;\\theta)]^{1-R_i}.
\\]

Use a proper, declared prior, for example a baseline independent
\\(\\theta_j\\sim N(0,2.5^2)\\), with prior-scale sensitivity runs
\\(1.5, 2.5, 5.0\\). These are candidate regularizing priors, not an assertion that one is uniquely correct. Report prior sensitivity, especially under near-separation and small validation sizes.

**Important parameterization control:** the current plug-in fitter standardizes redshift using each validation sample's own mean and standard deviation. Applying the same independent coefficient prior after that sample-dependent transformation changes the implied prior over the physical inclusion curve across validation sizes and replicates. For a fair prior-sensitivity experiment, either define the Bayesian model in a fixed, predeclared redshift coordinate based on the simulation support, or transform a single prior specified in that fixed coordinate into each sample-standardized parameterization (including its induced covariance). Do not reuse identical numerical prior scales across changing standardizations and call the prior invariant.

The coefficient posterior is

\\[
p(\\theta\\mid V)
\\propto L_V(\\theta)p(\\theta).
\\]

The validation sample must be drawn from the complete parent frame independently of retention outcomes, and its units must remain excluded from the analysis catalogue. The true simulation inclusion probabilities remain inaccessible to the fitting/posterior code; they are used only by the independent auditor and oracle diagnostic arm.

A finite, proper prior makes the posterior well-defined in cases where an unpenalized MLE is unstable or does not exist. It does **not** make a poorly identified fit informative: posterior concentration, prior sensitivity, and predictive diagnostics must still be reported.

## 3. Propagate each coefficient draw through the normalized density

For each posterior draw \\(\\theta^{(b)}\\), construct the IPW density from the observed analysis catalogue \\(C\\):

\\[
\\widetilde f_b(z)
=\\widehat g_{\\mathrm{IPW}}(z;\\theta^{(b)})
\\,m(z),
\\qquad
f_b(z)=
\\frac{\\widetilde f_b(z)}
{\\int_{z_{\\min}}^{z_{\\max}}\\widetilde f_b(u)\\,du}.
\\]

Here \\(\\widehat g_{\\mathrm{IPW}}\\) is the same predeclared reflected weighted-KDE operator used in the density benchmark, with object weights
\\(1/\\pi(z_i;\\theta^{(b)})\\). The merger-rate factor \\(m(z)\\) is applied once, exactly as in the current density estimator. Normalize **within each draw**; averaging unnormalized KDEs and normalizing only afterward is a different operation.

The implementation must preserve the current bandwidth protocol and use common random numbers / identical catalogue data across coefficient draws. Do not silently retune the Scott multiplier using the H0 outcomes. Keep bandwidth uncertainty separate from coefficient uncertainty.

Record, for every draw, at least:
- posterior draw identifier and chain/iteration provenance;
- minimum and maximum fitted inclusion probabilities over the analysis catalogue;
- maximum and upper-quantile inverse weights;
- Kish and entropy ESS of the normalized IPW weights;
- top-1% and top-5% normalized weight mass;
- density normalization and density quantiles on the fixed redshift grid.

A draw that produces numerically unusable probabilities or weights must be retained in a failure ledger, not silently discarded. Report posterior mass in the unusable region.

## 4. Correct likelihood marginalization

Let \\(\\mathcal L(H_0;f_b)\\) denote the existing event likelihood evaluated using density draw \\(f_b\\), with the gravitational-wave detection-selection factor kept as a separate, independently defined component. Then the modular/cut likelihood average is

\\[
\\mathcal L_{\\mathrm{marg}}(H_0)
=\\int \\mathcal L(H_0;f_\\theta)
\\,p(\\theta\\mid V)\\,d\\theta
\\approx
\\frac{1}{B}\\sum_{b=1}^{B}
\\mathcal L(H_0;f_b).
\\]

The resulting H0 posterior is

\\[
p(H_0\\mid d,C,V)
\\propto p(H_0)\\,
\\mathcal L_{\\mathrm{marg}}(H_0).
\\]

This is a **modular/cut construction**: the validation data update the inclusion coefficients, and the resulting uncertainty is propagated into the host-density likelihood. It is not automatically the same as a fully joint posterior for all latent galaxy, merger, and detector-selection quantities.

Numerical requirements:
1. Evaluate the likelihood average with a log-sum-exp calculation:
   \\(\\log\\mathcal L_{\\mathrm{marg}}(H_0)
   =\\operatorname{logsumexp}_b[\\log\\mathcal L(H_0;f_b)]-\\log B\\).
2. Reuse the same coefficient draws for every H0 grid point and every paired comparison. Do not draw a fresh coefficient sample at each H0.
3. Increase B until Monte Carlo error in the marginalized log likelihood and H0 summaries is below a predeclared tolerance.
4. Keep the GW detection-selection factor \\(\\beta_{\\mathrm{GW}}(H_0)\\) distinct from catalogue inclusion \\(\\pi(z;\\theta)\\). IPW corrects the catalogue sampling mechanism; it must not trigger a second missing-host mixture correction unless a separate likelihood derivation explicitly requires one.
5. Never average already-normalized H0 posteriors and call that the marginal likelihood. Average likelihoods over the coefficient posterior first, then normalize the H0 posterior.


## 4A. Repository-specific H0 likelihood adapter contract

The existing repository adapter for the intended synthetic fixed-count selection experiment is
`selection_aware_toy_likelihood` in `src/dark_siren_robustness/selection_toy.py`, not the minimal `dataset_likelihood.dataset_posterior_h0` prototype. The latter explicitly excludes selection conditioning and shared nuisance-parameter marginalization.

For each coefficient draw, pass:
- the fixed redshift grid;
- that draw's already normalized merger-weighted host density (f_b(z));
- the common H0 grid;
- the same detected-distance vector for every draw;
- the predeclared fixed Gaussian distance scale and detection threshold.

The adapter returns the unconditional joint event log likelihood, \\(\\beta_{\\rm GW}(H_0;f_b)\\), and the selected fixed-count log likelihood
\\[
\\ell_b(H_0)=
\\log p(d_{1:N}\\mid H_0,f_b)
-N\\log\\beta_{\\rm GW}(H_0;f_b).
\\]
It computes the event numerator and the detection efficiency using the **same supplied host density** and the same Gaussian noise/threshold assumptions. The generated event distances must all satisfy the declared detection rule. This is a toy fixed-count conditional likelihood; it does not model event rates or a realistic detector.

Therefore the coefficient-marginalized selected likelihood is
\\[
\\log \\mathcal L_{\\rm marg,sel}(H_0)
=\\operatorname{logsumexp}_{b=1}^{B}
[\\ell_b(H_0)]-\\log B.
\\]
Do not average the unconditional event likelihoods and detection efficiencies separately, and do not evaluate a single selection factor from the posterior-mean density. In general,
\\[
E_{\\theta\\mid V}
\\left[
\\frac{L_{\\rm uncond}(H_0;f_\\theta)}
{\\beta_{\\rm GW}(H_0;f_\\theta)^N}
\\right]
\\ne
\\frac{E[L_{\\rm uncond}(H_0;f_\\theta)]}
{E[\\beta_{\\rm GW}(H_0;f_\\theta)]^N}.
\\]
The nonlinear selection correction must remain inside each draw's likelihood before the likelihoods are averaged.

The H0 posterior is then formed by applying the declared H0 prior once to this marginalized selected likelihood and normalizing over the exact prior interval. Reuse the same coefficient draws at all H0 grid points. The event vector, noise scale, threshold, merger-weight convention, and H0/redshift grids remain fixed across paired arms so that this stage isolates coefficient uncertainty rather than changing the event model.

The existing host-density benchmark used a legacy absolute SciPy bandwidth factor of 0.8. That value is not to be reused as a Scott multiplier. The density-recovery benchmark's explicit `scott_multiplier` protocol is the current path; any H0 adapter must consume the resulting normalized density without silently applying another merger weight or smoothing operation.

## 5. Required comparison arms

For every paired simulated event realization, compare:

1. **Analytic-parent diagnostic:** known target density; isolates the event-likelihood and numerical baseline.
2. **Parent-sample KDE diagnostic:** finite parent sample without catalogue thinning; measures finite-parent/KDE error.
3. **Catalogue-only:** observed catalogue with no inclusion correction.
4. **Oracle IPW:** true inclusion probabilities used only as a simulation diagnostic.
5. **Estimated plug-in IPW:** one fitted MLE curve, matching the current benchmark.
6. **Estimated posterior-propagated IPW:** coefficient draws from \\(p(\\theta\\mid V)\\), with per-draw density normalization and likelihood marginalization.

The oracle arm is not a feasible estimator. It diagnoses the best-case effect of knowing the inclusion probabilities. No result should be phrased as proving that IPW universally dominates catalogue-only reconstruction.

## 6. Experimental sequence and stopping gates

### Gate A — posterior sampler validation without H0

- Validate posterior sampling on simulated logistic data with known coefficients.
- Check chain convergence, effective sample size, trace stability, and prior sensitivity.
- Compare posterior predictive retention probabilities on held-out parent units.
- Include zero curvature, positive/negative curvature, weak retention, and near-separation regimes.
- Test recovery against an independent numerical reference on small cases.
- Confirm that the posterior sampler never receives oracle probabilities.

### Gate B — density uncertainty

- Reuse the repository's frozen density-recovery design and paired seeds.
- Compare plug-in and posterior-propagated density IAE/ISE, pointwise bias, interval coverage of the analytic density, and weight-instability diagnostics.
- Separate conditional density uncertainty (given one parent/catalogue frame) from repeated-sampling variation across frames.
- Report all planned attempts, including failed/non-evaluable cases.
- Use a training/held-out split for any bandwidth choice; do not tune and evaluate on the same seeds.

### Gate C — H0 propagation

Only after Gates A and B pass:
- Connect density draws to the existing H0 event likelihood without changing the event model in the same experiment.
- First reproduce the analytic-parent and parent-sample baseline.
- Verify likelihood-grid convergence independently of coefficient-draw Monte Carlo convergence.
- Run simulation-based calibration under a correctly specified inclusion model.
- Run repeated-sampling interval coverage under both correct specification and declared misspecification.
- Report posterior bias, interval width, coverage, PIT/SBC diagnostics, numerical failures, and Monte Carlo error separately.

A fixed-truth coverage experiment and matched-model SBC answer different questions; neither substitutes for the other.

## 7. What would count as evidence

A useful result would not be “posterior propagation fixes IPW.” Evidence would be a reproducible map of conditions under which coefficient uncertainty:
- widens density and H0 uncertainty appropriately;
- exposes unstable weights that plug-in estimates conceal;
- improves or fails to improve calibration under correct specification;
- cannot compensate for systematic inclusion-model misspecification; and
- becomes dominated by prior sensitivity or weight concentration at small validation sizes.

The analysis must report paired contrasts and uncertainty across independent parent/event replicates, not treat coefficient draws as independent simulation replicates.

## 8. Explicit non-claims

This design does not establish:
- a new general solution to incomplete galaxy catalogues;
- superiority of estimated IPW over catalogue-only or existing completeness methods;
- validity for LVK/GWTC events or real galaxy surveys;
- a real observational H0 measurement;
- a universal optimal KDE bandwidth;
- posterior calibration from grid convergence alone; or
- a Bayesian interpretation for the existing inverse-information covariance.

**Current project gate:** repository-native density-recovery implementation and its independent artifact audit must pass on the exact branch head before this coefficient-posterior stage is implemented. The 100-replicate density run and its raw artifact have not been verified here.
