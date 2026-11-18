# Central scientific and novelty audit — 18 November 2026

## Scope and evidence boundary

This is a review of the GitHub-visible project state, not of any unpushed local files. The review branch starts from PR #55 head `95e8164bd1c17756e5416566ec3a28475b15f1b9`. No files on `main` were changed. No experiment was launched during this review.

The inspected material includes the repository README on `main`; PRs #44, #54, #55, #56, and #57; the quadratic-logit workflow; and the multi-replicate runner, outer auditor, cross-gamma auditor, five-arm weight builder, and synthetic selection generator at the PR #55 head.

Status labels:
- **Verified in GitHub metadata/source:** directly visible in the inspected repository/PR content.
- **Reported by PR description:** a claim recorded in the PR body, not independently rerun in this review.
- **Not verified here:** no current artifact/log/test evidence was independently inspected.
- **Not run:** no experiment was launched in this review.

## Executive decision

**Continue as a synthetic statistical-methodology benchmark, but do not describe the current work as a new dark-siren estimator or as publication-ready.**

The potentially defensible contribution is narrower:
1. a controlled comparison of catalogue-only, parent-oracle, known-probability IPW, MAP plug-in IPW, and validation-posterior-integrated IPW;
2. explicit separation of finite validation uncertainty from functional-form misspecification;
3. measurement of how density reconstruction errors and weight instability propagate into a correctly selection-conditioned dark-siren likelihood;
4. independent reconstruction of simulation artifacts and failure-inclusive reporting.

This is a research programme, not yet a demonstrated contribution. The current quadratic-logit workflow is density-only and has not been dispatched. The H0 stage is explicitly blocked.

## Repository and branch state

| Item | Current GitHub-visible status | Interpretation |
|---|---|---|
| PR #44 | Open draft; five-arm preregistration/implementation | Foundational design; not evidence of a completed confirmatory experiment |
| PR #54 | Open draft; batched likelihood/KDE and density-pilot infrastructure | Important computational work; PR body records historical pilot failure (19/20) and that a corrected full run/audit remains outstanding |
| PR #55 | Open draft; quadratic-logit stress DGP and paired density infrastructure | The seven-gamma workflow exists but is explicitly manual-only and not dispatched |
| PR #56 | Open draft; JSON newline/audit-output correction | Small audit serialization repair; exact-head CI is not established by the inspected PR body |
| PR #57 | Open draft; bandwidth holdout protocol | Design-only; no bandwidth selected, no holdout opened, no H0 propagation |
| `main` README | Describes earlier synthetic studies and limitations | README status is dated 2026-10-25 and does not fully describe the later five-arm/quadratic-logit programme |

PR #55 records earlier CI passes at commit `afb2b647...`, but says the paired gamma-vs-zero contrast code was added later and exact-head CI was pending. Its current head is `95e8164...`. The PR body explicitly says the seven-gamma experiment has not been dispatched. These historical CI claims are not treated as a current exact-head pass.

The workflow file `.github/workflows/quadratic-logit-curvature-density-pilot.yml` specifies seven gamma cells, 20 replicates per cell, N=4,000, validation frame 1,500, validation sizes 100/300/1,000, 11 absolute bandwidth candidates, 401 redshift-grid points, and a cross-gamma audit. It has only a `workflow_dispatch` trigger. No workflow dispatch action was available through the GitHub connector during this review, so it was not launched.

## Novelty review: material prior art found

The initial literature scan finds that catalogue incompleteness and reconstruction of missing galaxy structure are already active, mature topics. This materially narrows the novelty claim.

1. **Dalang & Baker (2024), “The clustering of dark sirens' invisible host galaxies,” JCAP 02 (2024) 024.** They study completion methods that use clustering information and compare them with homogeneous/multiplicative completion. The broad claim “catalogue incompleteness affects dark sirens and needs correction” is not novel.
   - DOI: https://doi.org/10.1088/1475-7516/2024/02/024

2. **Leyde, Baker & Enzi (2026), “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology,” JCAP 01 (2026) 013.** This develops a reconstruction method for the underlying galaxy field and magnitude distribution, respecting large-scale spatial correlations, and validates on Millennium-based simulations. A general claim of a new catalogue-completion method would overlap strongly.
   - DOI: https://doi.org/10.1088/1475-7516/2026/01/013
   - Preprint: https://arxiv.org/abs/2507.12171

3. **Borghi et al. (2026), “Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology,” A&A 706, A199.** This is a self-consistent completeness/host-weighting framework implemented in CHIMERA, with hierarchical cosmological and astrophysical population constraints and mock LVK-like studies. It further reduces the space for broad incompleteness/weighting novelty.
   - DOI: https://doi.org/10.1051/0004-6361/202557354
   - Preprint: https://arxiv.org/abs/2509.18243

4. **Tagliazucchi et al. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” Phys. Rev. D 114, 083523, published 13 October 2026.** This is particularly important and very recent: it uses a hierarchical Bayesian likelihood with an explicit survey-selection function to infer a non-parametric GP reconstruction of the LOS galaxy redshift distribution, then propagates it into H0 inference. Its mock catalogues extend to z<=0.4, directly overlapping the broad redshift regime of this project. The abstract reports improved precision relative to homogeneous completion in its simulated settings.
   - DOI: https://doi.org/10.1103/vmzf-1wxq
   - Preprint: https://arxiv.org/abs/2608.24992

5. **Datrier & Hendry (2026), “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” MNRAS 549, stag828.** A robust test for apparent-magnitude completeness limits has been applied to gwcosmo. This is another reason to avoid broad claims that the project introduces a general completeness diagnostic.
   - Article: https://academic.oup.com/mnras/article/549/1/stag828/8666390

6. **Gair et al. (2023), “The Hitchhiker's Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology,” AJ 166, 22.** This is essential background for likelihood structure, catalogue terms, and selection effects.
   - Preprint: https://arxiv.org/abs/2212.08694

### Novelty conclusion

The current evidence does **not** justify any of these claims:
- a new galaxy-catalogue completion estimator;
- a new general dark-siren likelihood;
- first use of Bayesian completeness modelling;
- first demonstration that catalogue incompleteness can bias H0;
- superiority over established reconstruction approaches.

A narrower possible contribution remains: a carefully controlled, independently audited study of **finite validation-sample uncertainty and propensity-model misspecification in IPW-style catalogue reconstruction, with downstream consequences for a toy dark-siren likelihood**. The literature scan has not yet established that this precise combination is novel. A full paper-by-paper comparison is still required, especially against the 2026 GP paper, CHIMERA completeness work, and statistical literature on estimated propensity weights.

## Technical and scientific audit observations

### 1. Five-arm estimands need one canonical mathematical statement

The inspected weight-builder docstring defines:
- A: catalogue-only;
- B: parent-population oracle;
- C: catalogue divided by true inclusion probability;
- D: catalogue divided by the validation MAP inclusion probability;
- E: node-specific inverse-probability weights integrated over a validation-posterior quadrature.

With a declared base weight b(z), the intended target appears to be q(z) proportional to g(z)b(z), where g is the parent redshift density and catalogue inclusion is pi(z). Under correct positivity and known pi, the population identity is:

  E_catalogue[f(Z)b(Z)/pi(Z)] / E_catalogue[b(Z)/pi(Z)]
  = E_parent[f(Z)b(Z)] / E_parent[b(Z)].

This identity should be stated in the methods and checked against every arm's implemented target. It is a density-recovery identity, not by itself a valid H0 likelihood or a guarantee of calibrated cosmological inference.

### 2. Validation-posterior uncertainty is conditional on the fitted model

Arm E integrates uncertainty over coefficients in the fitted degree-one logistic model. When gamma is nonzero, the data-generating inclusion mechanism is quadratic-logit while the fitted model remains degree one. Therefore E can integrate coefficient uncertainty while remaining structurally misspecified. It does not marginalize over the omitted quadratic term or model class.

The paper must distinguish:
- finite-sample coefficient uncertainty under a correctly specified degree-one model;
- pseudo-true parameter uncertainty under a misspecified degree-one model;
- uncertainty about the functional form itself;
- positivity/overlap failure and unstable inverse weights.

These are different sources of uncertainty and should not be compressed into “Bayesian uncertainty correction.”

### 3. The seven gamma values do not form a symmetric difficulty design

The independent calibration note at the PR #55 head reports equal parent-average inclusion rates across gamma, but markedly asymmetric propensity ranges and ideal-weight variability. At gamma=-2, the reported minimum pi is about 0.00379 and catalogue-law ideal-weight CV about 2.39; at gamma=+2, the corresponding values are about 0.0973 and 0.965. These are deterministic population calculations, not observed simulation results.

Consequently:
- do not describe +/-gamma as equally severe;
- report min/max pi, lower-tail pi quantiles, and population weight moments per cell;
- interpret gamma contrasts as a joint curvature-direction/positivity-stress experiment;
- consider a separate severity-matched exploratory design only after the current plan is preserved and reported.

### 4. Density and cosmological calibration are separate gates

The current workflow estimates density-reconstruction metrics (IAE, ISE, boundary IAE, weight diagnostics). It does not run the H0 likelihood or report bias, RMSE, coverage, PIT, or SBC. A density-stage improvement does not imply H0 calibration improvement.

The next H0 stage must explicitly:
- use one shared H0 parameter across events;
- use the full joint event likelihood and the correct GW-detection selection normalization;
- keep galaxy-catalogue inclusion and GW detection selection as distinct processes;
- define how catalogue inclusion uncertainty is conditioned/marginalized;
- compare posterior calibration and precision across repeated simulations;
- verify numerical convergence in H0/redshift grids and coefficient quadrature;
- avoid double-correcting missing hosts after IPW.

### 5. Validation sampling design must be made explicit

The methods section needs an exact statement of how the validation frame is sampled from the parent population, whether it is an independent probability sample, whether sampling fractions are known, and whether the same validation frame is nested across n=100/300/1000. The code's separate population and validation seeds are a good reproducibility feature, but the inferential assumptions must be stated, not inferred from the implementation.

If the validation data are an independent sample from the parent, say so and derive the likelihood for validation labels conditional on sampled redshifts. If not, specify the actual design and any design weights. Clarify whether the validation sample is external to the catalogue and whether overlap is possible.

### 6. Positivity and finite-sample failure accounting

The minimum inclusion probability in the negative-curvature stress cells is very small. A complete report should include:
- parent-measure and catalogue-measure pi quantiles;
- maximum theoretical inverse weight and finite-sample maximum realized weight;
- Kish ESS and entropy ESS;
- top 1% and 5% normalized weight mass;
- fraction of parent support with pi below prespecified thresholds;
- validation sample's coverage of low-pi regions;
- complete replicate lifecycle, including optimizer failures and non-evaluable cases.

Do not silently clip weights in the primary analysis. If truncation/stabilization is explored, make it a separately labelled sensitivity analysis with its changed estimand/bias implications.

### 7. Audit architecture: good direction, but exact-head execution still matters

The code has a useful separation between production and independent reconstruction. The outer auditor checks manifests, seeds, child audit results, and aggregate consistency; the cross-gamma auditor checks common configuration and matched keys before paired contrasts. These are valuable software controls.

However, the PR descriptions themselves state that the latest exact-head CI was pending/not verified and the seven-gamma experiment was not run. No test or experiment result is upgraded to PASS merely because the workflow YAML exists. The cross-gamma summary gate requires complete successful cells, so one failed replicate in any gamma cell blocks the paired summary by design; this is conservative but should be documented as a strict-complete-case gate and accompanied by failure-inclusive descriptive reporting.

## Required next steps (ordered)

1. **Freeze the scientific question and estimands.** Write one methods note that defines the target density, all five arms, validation sampling design, and the distinction between coefficient uncertainty and model misspecification.
2. **Complete exact-head software validation.** Run focused tests and full CI on the current PR #55 head; resolve any failures before expensive simulations. Review PR #56's serialization fix and ensure the outer auditor writes parseable JSON.
3. **Run the seven-gamma density pilot.** Only after exact-head CI is green, dispatch the manual workflow. Preserve all raw artifacts and all failed attempts; independently audit each cell and the cross-gamma result.
4. **Review density-grid adequacy.** The 11-width grid is broad, but confirm whether it brackets risk minima. Do not select a bandwidth from the same replicates used for final evaluation.
5. **Design a separate fresh-seed bandwidth training/holdout experiment.** Do not reuse the overlapping pilot seeds or the floor-DGP bandwidth choice as evidence for the quadratic-curvature family.
6. **Only then conduct a formal H0-stage design review.** Do not launch H0 propagation automatically after a density pass.
7. **Complete the literature review and novelty matrix.** Read full methods/supplementary material for the 2026 GP paper, Leyde et al., Borghi et al., Datrier & Hendry, and the core dark-siren likelihood review. Add statistical literature on estimated propensity scores, Bayesian propensity uncertainty, positivity, and finite-sample IPW.
8. **Choose journal pathway only after results.** The current work may fit an astrostatistics/computational methods venue better than a top cosmology journal unless the final H0-level results are both scientifically meaningful and robust. The statistical-methodology path requires a generalizable statistical insight beyond implementing standard IPW and logistic regression in a toy simulation.

## Current gate table

| Gate | Status |
|---|---|
| Core repo/PR state inspected | PASS for GitHub-visible subset |
| Local unpushed files checked | NOT CHECKED |
| Full literature/novelty review | INCOMPLETE |
| Five-arm estimands fully formalized in a single methods note | NOT VERIFIED |
| Latest PR #55 exact-head CI | NOT VERIFIED |
| PR #56 exact-head CI | NOT VERIFIED |
| Seven-gamma density pilot | NOT RUN |
| Independent cross-gamma artifact audit | NOT RUN |
| Bandwidth selected | NO |
| Fresh disjoint bandwidth holdout | NOT RUN |
| H0 propagation | BLOCKED |
| H0 calibration/SBC | NOT RUN for this new five-arm curvature programme |
| Publication readiness | NO |

## References identified in the initial scan

- Dalang, C. & Baker, T. (2024), JCAP 02, 024, DOI 10.1088/1475-7516/2024/02/024.
- Gair, J. R. et al. (2023), AJ 166, 22, arXiv:2212.08694.
- Leyde, K., Baker, T. & Enzi, W. (2026), JCAP 01, 013, arXiv:2507.12171.
- Borghi, N. et al. (2026), A&A 706, A199, arXiv:2509.18243.
- Datrier, L. E. H. & Hendry, M. A. (2026), MNRAS 549, stag828.
- Tagliazucchi, M. et al. (2026), Phys. Rev. D 114, 083523, arXiv:2608.24992.

This is an initial targeted scan, not an exhaustive systematic review.


## Statistical-methodology prior art added in the second scan

The statistical side also has substantial prior art on the exact failure modes being studied. This means “IPW is unstable under low inclusion probabilities” and “estimated propensity uncertainty matters” are established facts, not novel findings.

- **Lunceford-style finite overlap/misspecification comparison:** a simulation study of propensity-score weighting under limited overlap and model misspecification compares conventional IPW, trimmed IPW, overlap weights, matching weights, and entropy weights. It reports instability and bias/variance/coverage trade-offs under weak overlap. The project must cite this class of work and explain what is distinctive about the catalogue-reconstruction estimand and downstream cosmological likelihood.
  - DOI: https://doi.org/10.1177/0962280220940334
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/32693715/

- **Bayesian design-stage uncertainty:** a 2020 Statistics in Medicine paper explicitly considers uncertainty in estimated propensity scores and marginalizing over design-stage outputs in a two-stage Bayesian propensity-score analysis. Therefore Bayesian integration over propensity-model uncertainty is not itself a new statistical contribution.
  - DOI: https://doi.org/10.1002/sim.8486
  - Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC9170228/

- **Variance accounting for estimated weights:** a tutorial on variance estimation for inverse-probability weighting emphasizes that inference must account for the propensity-score estimation step. The project should compare its posterior integration target with established frequentist two-step variance accounting, rather than presenting estimated-weight uncertainty as an unexplored issue.
  - DOI: https://doi.org/10.1002/sim.10078

- **Positivity diagnostics and alternative estimands:** a 2025 Statistics in Medicine tutorial reviews estimand changes and diagnostics under positivity violations. Any proposed stabilized, truncated, overlap, or entropy-weight sensitivity arm must state whether it changes the target population; such methods cannot be compared as if they all estimate the same full-parent density.
  - DOI: https://doi.org/10.1002/sim.70329
  - Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12659692/

### Updated statistical novelty conclusion

The plausible contribution is not “Bayesian propensity uncertainty for IPW.” That is established. The defensible research question is whether the combination of (i) an external validation sample for catalogue retention, (ii) a finite validation design, (iii) explicit misspecification of the retention model, and (iv) downstream reconstruction of a host-redshift density used in a selection-conditioned cosmological likelihood creates distinctive calibration behavior that existing generic propensity-score studies do not address. This remains a hypothesis until the full comparative literature review and experiments are complete.

In particular, a strong statistical paper would need at least one result that generalizes beyond this one toy DGP—for example, a clearly defined decomposition or diagnostic with prospective predictive value for density/H0 calibration failure, supported by multiple parent laws, inclusion functions, validation designs, and independently held-out simulation conditions. Merely showing that a linear logistic model fails when the true model contains a quadratic term would be a standard misspecification demonstration.

