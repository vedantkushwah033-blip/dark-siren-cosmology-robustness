# Prospective diagnostic-comparison protocol (synthetic study)

**Status:** proposed protocol; not a report of completed results.  
**Date:** 2026-10-27  
**Scope:** controlled synthetic low-redshift host-mixture model only.

## 1. Motivation and boundary

Galaxy-catalogue incompleteness, host weighting, and selection effects are established topics in dark-siren cosmology. This project must not claim those topics as novel. The narrower question proposed here is:

> In a declared synthetic host-mixture model, how reliably do posterior predictive/inference diagnostics distinguish host-population misspecification from a misspecified distance-noise model, and which errors remain hard to detect at finite sample size?

This question is not yet established as novel. Before any manuscript claims originality, search and compare prior calibration, simulation-based calibration, posterior predictive checking, and dark-siren pipeline validation studies.

The current inference implementation uses the approximation $d_L \approx cz/H_0$. It is not a complete hierarchical GW likelihood and does not include detector-selection normalization. All findings from this protocol must therefore be labelled synthetic and model-conditional.

## 2. Primary estimand

For each prespecified data-generating condition $s$ and inference arm $a$, estimate the repeated-sampling coverage probability

\[
C_{a,s}(q)=P_{s}\{H_{0,\mathrm{true}}\in I_{a,s}(q)\},
\]

where $I(q)$ is the central posterior credible interval with nominal mass $q\in\{0.68,0.95\}$. Report the empirical coverage, its binomial uncertainty interval, posterior interval width, bias, and RMSE together. Coverage alone is not sufficient: a very broad interval can achieve coverage at the cost of weak information.

The primary contrast is the paired difference in a prespecified metric between two inference arms applied to the *same simulated observations*. Report the estimate and a paired uncertainty interval; do not infer a meaningful difference solely because one arm's interval excludes a target while another does not.

## 3. Data-generating factors

Use a small, declared factorial design with independently controlled factors:

1. True hidden-host fraction, including zero.
2. True hidden-host redshift/population distribution.
3. True distance-noise scale.
4. Assumed distance-noise scale.
5. Assumed hidden-host fraction and redshift distribution.
6. H0 truth and prior support, with explicit boundary cases.
7. Number of independent events/replicates.

Do not conflate the number of posterior fits with the number of independent simulated datasets. State the sampling unit for every uncertainty estimate. Keep seeds, cell definitions, and exclusions fixed before examining results.

## 4. Inference arms and controls

At minimum, compare:

- **Matched/oracle control:** likelihood matches the declared generator; this is a diagnostic ceiling, not an operational method because it uses truth-known settings.
- **Host misspecification only:** alter the assumed hidden fraction or redshift distribution while retaining the correct noise model.
- **Noise misspecification only:** alter the assumed distance-noise scale while retaining the correct host model.
- **Joint misspecification:** alter both.
- **Omission control:** fit without the hidden component when the generator includes it.
- **False-positive component control:** include a hidden component when the generator has zero hidden fraction.

The likelihood and generator must agree on whether observed distances are truncated/redrawn to remain positive. Test this choice explicitly; do not silently switch data-generation conventions between arms.

## 5. Diagnostics and reporting

Prespecify the following, where applicable:

- Empirical 68% and 95% interval coverage, with Wilson or another justified binomial interval.
- Bias and RMSE of posterior mean or another declared point estimator.
- Mean and distribution of interval widths.
- PIT or rank diagnostics only with a clear definition of the predictive quantity and reference distribution.
- Simulation-based calibration (SBC) ranks when parameters are drawn from the prior and data are generated from the same model; do not treat SBC as proof of robustness to misspecification.
- Paired contrasts on common observations, with the resampling unit stated.
- Sensitivity to Monte Carlo error, seed choice, grid resolution, prior boundaries, and finite replicate count.

A diagnostic that performs well for one misspecification class must not be described as a universal detector. If two misspecification classes produce indistinguishable diagnostic distributions under the declared design, report the non-identifiability rather than forcing a classification.

## 6. Statistical analysis plan

1. Freeze the factor grid, seeds, replicate counts, metrics, and pass/fail criteria in a version-controlled commit before the final run.
2. Include a pilot solely to check runtime and numerical stability; do not use the pilot to tune criteria and then report it as confirmatory.
3. Run independent simulations at each prespecified cell. Save one row per simulated dataset and inference arm, including seed/replicate identifier, generating parameters, assumed parameters, convergence/validity flags, posterior summaries, and diagnostic values.
4. Reconstruct summaries from raw rows in a separate audit script that does not import the runner's summary functions.
5. Report uncertainty intervals for every key estimate and contrast. Correct or otherwise account for multiplicity if many cells are used for confirmatory claims; alternatively designate the factorial as exploratory.
6. Repeat selected cells at larger replicate counts to assess Monte Carlo stability.
7. Keep failed fits in a failure ledger. Do not silently drop them or replace them without a documented rule.
8. Publish all scripts, dependency versions, source revision, seeds, raw synthetic records, checksums, and exact reproduction commands.

## 7. Go/no-go gates

Do **not** write a claim of a new diagnostic or a journal-ready conclusion unless all gates pass:

- [ ] Downstream call graph for the corrected missing-host likelihood is exhaustively inventoried.
- [ ] Unit tests cover limiting cases $f_{\mathrm{miss}}=0$, near-zero missing mass, and multiple missing-host weights.
- [ ] Data generator and likelihood conventions match for positive-distance truncation and all other sampling assumptions.
- [ ] Clean-checkout reproduction succeeds and CI passes on supported Python versions.
- [ ] Independent auditors reconstruct raw summaries and paired contrasts from raw records.
- [ ] The design includes enough independent replicates to make uncertainty intervals informative; the stopping rule is prespecified.
- [ ] Literature review compares the proposed diagnostic question with existing calibration and dark-siren validation methods.
- [ ] A reference-pipeline comparison is completed, or the manuscript explicitly limits itself to a toy-model methodological note.
- [ ] Every headline claim survives alternative prior support, grid resolution, and plausible generator variants.

A failed gate means narrow the claim, fix the problem, or stop the publication claim—not relax the criterion after seeing the result.

## 8. Related work to read before claiming novelty

- Gray, Messenger & Veitch (2022), *A pixelated approach to galaxy catalogue incompleteness*, MNRAS, DOI: 10.1093/mnras/stac366.
- Gair et al. (2023), *The Hitchhiker's Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology*, AJ, DOI: 10.3847/1538-3881/acca78.
- Finke et al. (2021), *Cosmology with LIGO/Virgo dark sirens*, JCAP, DOI: 10.1088/1475-7516/2021/08/026.
- Datrier & Hendry (2026), *Implementing a robust test of galaxy catalogue completeness for dark siren measurements of H0*, MNRAS, DOI: 10.1093/mnras/stag828.
- Tagliazucchi et al. (2026), *Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs*, Phys. Rev. D, DOI: 10.1103/vmzf-1wxq.
- Leyde, Baker & Enzi (2026), *Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology*, JCAP, DOI: 10.1088/1475-7516/2026/01/013.

This is a starting bibliography, not an exhaustive systematic review. Verify publication metadata and read the methods/results sections before using these papers to support manuscript claims.

## 9. Current status

This document defines a proposed next-stage protocol. It does not assert that the full design has been implemented, that any new diagnostic has been discovered, that external CHIMERA validation has passed, or that the project is ready for submission.
