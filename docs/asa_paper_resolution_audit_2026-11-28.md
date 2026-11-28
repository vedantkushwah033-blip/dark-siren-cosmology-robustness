# ASA paper resolution: source audit and literature-gap assessment

**Date:** 28 November 2026  
**Repository:** `vedantkushwah033-blip/dark-siren-cosmology-robustness`  
**Branch:** `research/asa-paper-resolution-20261128`  
**Evidence base:** hidden-host source revision `ecb260886716c54fdd2d5e47390701949f5a18d5`, checkpoint `research/asa-paper-checkpoint-20261128`, and focused literature search completed 2026-11-28.

> This is a first-stage audit and decision record. It does not claim that all P0/P1/P2 tasks are complete. No confirmatory workflow was launched and no existing artifact was modified.

## 1. Executive decision

The existing hidden-host experiments are internally coherent for the stated *discrete, fixed-parameter toy model*, but they do not estimate the hidden-population parameters. The supplied continuous hidden-redshift analogue is a different model and is not present in the audited production path.

The paper should not claim novelty for showing that galaxy-catalogue incompleteness or missing hosts can affect dark-siren inference. Multiple recent papers already study catalogue incompleteness, host weighting, and reconstruction of missing galaxy distributions, including Borghi et al. (2026), Leyde et al. (2026), VanWyngarden et al. (2025 preprint), and the 2025 PRD study “Systematic bias in dark siren statistical methods…”.

A defensible possible contribution is narrower: a transparent, controlled statistical study of *fixed-truth posterior reliability under hidden-component omission and deliberate hidden-component misspecification*, with paired method comparisons, explicit Monte Carlo uncertainty, prior-support sensitivity, and a clear distinction between prior-predictive SBC and fixed-truth coverage. Whether this is sufficiently novel and useful for an ASA paper remains conditional on the full literature review and domain review.

## 2. Source-level audit of the existing hidden-host implementation

### 2.1 Data-generating model

Inspected:
- `scripts/run_host_misspecification_sensitivity_surface.py`
- `scripts/run_hidden_host_assumption_sensitivity.py`
- `scripts/run_host_prior_support_factorial.py`
- corresponding independent audit scripts and protocol documents
- `src/dark_siren_robustness/hardened_inference.py`

The generator uses listed redshifts 0.02 and 0.04, listed relative weights 0.7 and 0.3, and one hidden component at a fixed cell-specific redshift. For each observation it draws a host label with probabilities (((1-f)0.7,(1-f)0.3,f)), then draws a Gaussian distance with mean (cz/H_{0,true}) and cell-fixed (sigma), redrawing the same selected host until (d>0).

The cell noise is (sigma = r c(0.03)/H_{0,true}), with (r) the specified noise fraction. The same absolute (sigma) is supplied to inference for every candidate H0. This is a fixed-absolute-noise experiment conditional on the cell's injected truth; it is not a fractional-noise likelihood whose sigma changes with candidate H0.

### 2.2 Likelihood and truncation

For a host at redshift z, the implemented positive-truncated component density is
[
p(dmid H_0,z,d>0)=
\frac{\phi(d;cz/H_0,\sigma)}{\Phi(cz/(H_0\sigma))},\quad d>0.
]
The denominator is correct for a Gaussian conditioned on (d>0), because (P(D>0)=\Phi(\mu/\sigma)). Since the simulator selects the host first and then redraws the distance for that same host until positive, the host-mixture probabilities remain the stated probabilities; no additional host reweighting by the truncation acceptance rate is appropriate for this conditional generator.

The inference code evaluates the truncation correction in log space using `norm.logcdf(mean/sigma)`, and combines listed and hidden components with log-sum-exp/log-add-exp. The hidden fraction is fixed input to each inference arm, not a parameter inferred from the observations.

### 2.3 What is and is not estimated

**Implemented:** H0 posterior conditional on fixed listed-host weights, fixed hidden fraction/redshift when the hidden component is included, fixed event sigma, and a chosen uniform H0 prior.

**Not implemented in this experiment family:** posterior estimation of hidden fraction (f), hidden-redshift location (mu_h), hidden-redshift scale (s_h), or joint posterior (p(H_0,f,\mu_h,s_h\mid d)). The existing oracle arm is given the true hidden fraction and redshift. The fixed-assumption arms are stress tests, not fitted nuisance-parameter models.

The supplied log-distance continuous-mixture analogue is therefore not a simple parameter setting of the existing production model. It changes the observation model from Gaussian distance to Gaussian log-distance and replaces a point-mass hidden redshift with a truncated continuous distribution. It needs its own generative model, likelihood, tests, and audit; its results must not be mixed with the existing 278,400 posterior fits.

### 2.4 Existing implementation checks

The inspected source and protocol agree on:
- identical observations reused across inference arms within each replicate;
- the (f=0) negative control, where oracle/modelled and omitted components coincide;
- truths 55, 70, 85 (or 60, 70, 80 in the surface experiment) lying inside the relevant H0 prior supports;
- a 4,001-point H0 grid from 25 to 125;
- 200 replicates per base cell and 1,000 paired bootstrap resamples;
- separate independent posterior reconstruction in the auditor.

The current independent auditors reproduce the declared generator and posterior equations. This is strong internal reproducibility evidence, but the auditor shares the same mathematical model specification; it is not a substitute for a second domain-informed model review.

### 2.5 Items that remain open after source inspection

1. Full-grid numerical convergence on all predeclared primary cells is not established by the five-replicate-per-cell 4,001-vs-8,001 screen.
2. Coverage Monte Carlo uncertainty is material at 200 replicates. For a nominal 68% coverage rate, the binomial standard error is approximately 3.3 percentage points and a simple 95% half-width approximately 6.5 points. For nominal 95% coverage, only about ten misses are expected in 200 trials, so precision is also limited.
3. The continuous hidden-redshift estimability model has not been implemented or reproduced in the repository.
4. No domain expert review of the simplified host mixture, positive-distance truncation, or interpretation of the distance model is documented.
5. The literature gap is narrower than the earlier broad hidden-host framing suggested.

## 3. Focused literature review and novelty implications

This is a focused first-pass review, not yet an exhaustive paper-by-paper review of every relevant method.

| Work | What it already studies | Consequence for our claims |
|---|---|---|
| Borghi et al. (2026), “Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology,” A&A 706, A199, DOI: 10.1051/0004-6361/202557354 | Self-consistent catalogue-incompleteness and host-weighting framework in CHIMERA; 100 BBH events in an O5-like setup; varying completeness and weighting; H0 and population correlations. | We cannot claim that a simulation study of missing catalogue hosts and H0 robustness is itself new. Our simplified experiment must be positioned as a statistical diagnostic study, not a replacement astrophysical framework. |
| Leyde, Baker & Enzi (2026), “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology,” JCAP 01 (2026) 013, DOI: 10.1088/1475-7516/2026/01/013 | Reconstructs the underlying galaxy field and magnitude distribution, retaining large-scale spatial correlations; validates on Millennium-simulation mock data. | We must distinguish a single latent redshift-mixture stress test from spatially structured galaxy-field reconstruction. |
| VanWyngarden et al. (2025), “How Low Can You Go: Constraining the Effects of Catalog Incompleteness on Dark Siren Cosmology,” arXiv:2511.04786 | Studies catalogue incompleteness, host halo selection, clustering, and localization effects with mock galaxy catalogues. | Completeness thresholds and clustering-driven robustness are already active topics; avoid presenting simple incompleteness sensitivity as a novel physical result. |
| “Systematic bias in dark siren statistical methods and its impact on Hubble constant measurement,” Phys. Rev. D 112, 063561 (2025), DOI: 10.1103/vd36-3mys | Examines systematic biases from catalogue incompleteness and incorrect host weighting in dark-siren H0 inference. | The broad host-weighting/misspecification motivation is established; any novelty must be about a precise statistical estimand or diagnostic not already treated. |
| Talts et al. (2018), “Validating Bayesian Inference Algorithms with Simulation-Based Calibration,” arXiv:1804.06788 | Introduces SBC for checking Bayesian computation using prior-predictive simulations and posterior ranks. | The paper must cite SBC as a computational validation method and distinguish algorithm validation under the assumed model from robustness to a misspecified data-generating process. |
| Kleijn & van der Vaart (2012), “The Bernstein–von Mises theorem under misspecification,” Electronic Journal of Statistics 6, 354–381, DOI: 10.1214/12-EJS675 | Establishes posterior asymptotics under model misspecification and explains why Bayesian credible sets need not be valid frequentist confidence sets. | The manuscript should connect empirical coverage degradation to the established misspecification theory, without claiming a new theorem. |

### Novelty decision

Do **not** use these claims:
- “First study of missing hosts in dark-siren cosmology.”
- “First demonstration that catalogue incompleteness biases H0.”
- “A new Bayesian method for correcting galaxy-catalogue incompleteness.”
- “The proposed method is superior for real GW analyses.”

Potentially supportable claim, subject to a more exhaustive review:
> We provide a controlled, reproducible simulation study that separates hidden-component omission from fixed hidden-component misspecification and examines their effects on H0 point error and fixed-truth credible-interval coverage, alongside prior-support sensitivity and a separate prior-predictive calibration diagnostic.

Even this claim must be checked against all prior papers and should be framed as a scoped diagnostic contribution, not as a novel astrophysical completeness model.

## 4. Mathematical caution for the supplied continuous-mixture analogue

The analogue proposes (u=\ln d), (u\mid z,H_0\sim N(\ln(cz/H_0),\sigma^2)), and (z\sim TN(\mu_h,s_h^2; a,b)) with fixed truncation limits (a=0.005,b=0.10).

The claim that the all-hidden case (f=1) has an **exact** H0–hidden-redshift scale invariance needs to be corrected or proved carefully. If the truncated-normal support limits (a,b) remain fixed while (mu_h,s_h) vary, the transformation (H_0\mapsto kH_0), (z\mapsto kz) also transforms the truncation limits. It therefore does not generally leave the same fixed-support truncated-normal family invariant. An approximate weak direction may still occur away from the truncation boundaries, but exact invariance cannot be asserted from the current description.

The analogue's Fisher condition number must also be computed after a clearly specified parameter scaling; raw condition numbers depend on units. Local Fisher rank/conditioning is not a proof of global identifiability. The reported analogue results should remain exploratory until code, parameterisation, derivatives, priors, and finite-sample behavior are independently checked.

## 5. Recommended frozen research question and design boundary

**Candidate research question:**
> In a controlled low-redshift dark-siren toy model, how do omission and fixed misspecification of a latent host component affect H0 posterior point error and fixed-truth credible-interval coverage, and how do these effects interact with prior support?

**Primary outcomes:**
- posterior-median bias;
- RMSE;
- 68% and 95% fixed-truth coverage;
- mean 68% and 95% interval widths.

**Primary contrasts:**
- omitted component minus oracle component, labelled explicitly as a non-deployable benchmark;
- fixed assumed fraction/redshift minus oracle;
- selected prior support minus [50,90], treated as sensitivity rather than a search for the “best” prior.

**Separate diagnostic target:** prior-predictive SBC ranks under the older calibration model. Do not pool its results with the positive-truncated hidden-host experiments.

**Scope statement:** synthetic low-redshift Euclidean distance model; discrete listed hosts plus one point-mass hidden component; no real LVK/GWTC data, realistic detector-selection function, peculiar velocities, redshift errors, spatial galaxy clustering, or full hierarchical astrophysical population model.

## 6. Confirmatory simulation recommendation

Do not rerun every cell automatically. First select a small set of scientifically interpretable primary cells, then freeze:
- injected H0 values covering low, central, and high prior positions;
- a low and a high distance-noise condition;
- hidden fractions 0, 0.30, and 0.45;
- at least two hidden redshifts to expose redshift misspecification;
- omitted, oracle, and a small number of prespecified wrong-assumption arms;
- the three prior supports only for the designated prior-sensitivity subset.

For 68% coverage, 1,000 replicates give an approximate simple 95% binomial half-width of 2.9 percentage points near the nominal rate. This is a planning approximation; final uncertainty should use the paired structure and account for multiple cell/contrast reporting. Freeze a precision target before seeing the new results. Recompute all primary metrics and contrasts from raw run-level rows with an analysis script whose version and hash are archived.

## 7. Required next gates

| Gate | Required evidence | Status after this audit |
|---|---|---|
| Source/model contract | Generator, likelihood, truncation, fixed sigma, mixture weights independently reconciled | **Partially resolved** for current discrete toy model |
| Literature/novelty | Full overlap matrix including recent 2025–2026 dark-siren papers and statistical literature | **Partially resolved**; broad novelty claim rejected |
| Continuous hidden-parameter model | Separate production implementation, analytic/numerical likelihood checks, priors, identifiability study | **Not started in repository** |
| Monte Carlo precision | Frozen primary cells, prespecified precision target, increased replicates, raw-level reanalysis | **Not complete** |
| Numerical convergence | Full audit on primary cells with at least two finer grids and tail/quantile checks | **Not complete** |
| Domain review | Written review of toy DGP, host weights, selection/truncation interpretation, and limitations | **Not obtained** |
| Manuscript | Evidence freeze, complete methods/results, limitations, references, independent reproduction | **Not started** |
| Poster | Final paper figures and a concise conference narrative | **Not started** |

## 8. Reproducibility and safety of project operations

- Work is confined to this research branch.
- `main` is unchanged.
- No existing experiment artifacts were deleted or overwritten.
- No PR was opened or merged.
- No new confirmatory workflow was launched.
- This document records a first-stage audit, not completion of the full ASA paper.
