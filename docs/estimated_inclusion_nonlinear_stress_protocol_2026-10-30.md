# Estimated-inclusion nonlinear-selection stress protocol

**Protocol date:** 31 October 2026  
**Status:** prospective design; the current code contains a non-confirmatory smoke path, not a completed factorial campaign.  
**Scope:** controlled synthetic, low-redshift dark-siren host-catalogue inference only.

## 1. Research question and estimands

The study asks how finite validation samples and inclusion-model misspecification interact when inverse-inclusion weights are used to reconstruct a host-redshift density.

The generating retention law is allowed to be quadratic in a standardized redshift coordinate. Two validation estimators are compared on the same replicate:

- **Linear logistic:** intercept plus linear redshift term.
- **Quadratic logistic:** intercept, linear term, and quadratic term.

The primary paired model contrast is
\[
D_r=|\widehat H_{0,r}^{\rm linear}-H_0^{\rm true}|
    -|\widehat H_{0,r}^{\rm quadratic}-H_0^{\rm true}|.
\]
A positive value means the quadratic estimator has smaller absolute posterior-mean error for that replicate. This estimand says nothing by itself about interval coverage, uncertainty calibration, or general estimator superiority.

A separate, secondary contrast compares each estimated-inclusion arm with the catalogue-only baseline. Oracle weighting is retained only as a simulation diagnostic using the known generating probabilities; it is not an implementable estimator.

All results are synthetic. They are not LVK/GWTC validation, a real-observation Hubble-constant measurement, or evidence that this toy selection law represents a particular survey.

## 2. Data-generating mechanism

For parent-frame redshift \(z\), define
\[
x=\frac{z-(z_{\min}+z_{\max})/2}{(z_{\max}-z_{\min})/2}.
\]
Catalogue retention is generated independently for each parent object:
\[
P(R=1\mid z)=\operatorname{logit}^{-1}
\{a+b x+g(x^2-1/3)\}.
\]

The quadratic coefficient \(g\) controls curvature. The subtraction of \(1/3\) centers \(x^2\) under a uniform reference measure on \([-1,1]\) only; it does not generally center the term under the simulated parent-redshift distribution. Consequently, \(a\) is not generally the average log-odds or midpoint log-odds.

The validation sample is a simple random sample without replacement from the complete parent frame and records both redshift and catalogue-retention status. This is a strong, explicit design assumption. It is not equivalent to validating only objects already present in the retained catalogue.

The current pilot keeps the intercept and slope fixed and varies \(g\). It uses the same parent population, retained catalogue, GW-event realization, and nested validation records for both fitted model families within each replicate. Only the fitted logistic degree changes.

## 3. Staged factorial design

Candidate confirmatory factors are:

| Factor | Candidate levels | Purpose |
|---|---|---|
| True curvature \(g\) | \(-1.2,-0.6,0,+0.6,+1.2\) | Correctly specified reference and symmetric nonlinear stress |
| Validation size \(m\) | 100, 300, 1,000 | Finite validation information |
| Parent population | 4,000 initially | Hold the population-generation scale fixed |
| Detected GW events | 30 initially | Hold event information fixed |
| Replicates | 20 smoke; then precision-justified confirmatory count | Separate pipeline checks from inferential precision |

These are candidate settings, not a frozen confirmatory design. First inspect the actual inclusion-probability range, fractions below 0.05 and 0.10, retention counts, logistic-fit failures, coefficient-boundary hits, and weight concentration for every curvature. If a cell approaches practical positivity failure, report it as a distinct severe-overlap regime rather than pooling it with moderate stress.

The runner accepts one curvature value per invocation. Each invocation must preserve its seed, exact CLI arguments, output files, and checksum. Combining conditions requires a later aggregation step that retains an explicit condition identifier; files must not simply be concatenated and treated as one homogeneous experiment.

## 4. Arms

1. **Catalogue-only:** unweighted retained-catalogue density.
2. **Oracle:** inverse weighting by the true retention probabilities; diagnostic benchmark only.
3. **Estimated linear:** validation-fitted degree-1 logistic inclusion model.
4. **Estimated quadratic:** validation-fitted degree-2 logistic inclusion model.

For each validation size, the linear and quadratic arms must use exactly the same validation records and the same replicate-level parent/catalogue/GW realization. Validation samples at different sizes must be nested within a replicate. Seed lineage and validation inclusion counts must be retained.

Do not silently clip probabilities or cap weights. Any stabilized estimator must be a separately named arm with a prespecified rule. A quadratic fit is more flexible, not automatically better: it can have higher variance, boundary instability, or worse out-of-sample probability estimates at small validation sizes.

## 5. Outcomes and reporting

### Primary outcome
- Paired linear-minus-quadratic difference in absolute posterior-mean \(H_0\) error, with paired replicate count, mean, standard error, and interval.
- Report the linear and quadratic arms' failure counts alongside the complete-case contrast. The paired contrast is conditional on both fits succeeding and must not conceal differential failures.

### Required secondary outcomes
- Signed posterior-mean bias and RMSE.
- Fixed-truth 68% interval coverage with Wilson uncertainty interval and mean interval width.
- Inclusion-probability MAE evaluated over the complete synthetic parent frame; truth is scoring-only and must never enter fitting or posterior construction.
- Host-density integrated squared error against the known generating target.
- Logistic-fit failure and coefficient-boundary-hit rates.
- Inverse-weight ESS and ESS fraction, maximum inverse weight, weight quantiles, and pre-posterior support/overlap diagnostics.

The current pilot writes per-arm summaries and paired absolute-error contrasts. It does **not yet** provide paired linear-versus-quadratic contrasts for coverage, interval width, probability MAE, or host-density error. Those comparisons must be added to the producer, independent auditor, and adversarial tests before they are treated as confirmatory analysis outputs. Do not infer these contrasts from the current absolute-error contrast.

### Monte Carlo precision
The 20-replicate stage is a pipeline/runtime smoke test only. Do not use it to rank estimators. Before confirmatory runs, use pilot variability and a stated compute budget to select replicate counts. For coverage near 0.68, the binomial approximation gives a 95% Monte Carlo half-width of about 2 percentage points at roughly 2,100 independent replicates; paired designs and failure rates require additional planning. This approximation is a planning reference, not a universal required sample size or evidence of achieved precision.

Treat the factorial cells as a family. Specify the primary contrast and multiplicity strategy before examining confirmatory outcomes. Report every cell, including failures and unfavorable results; do not select only the best curvature or validation size.

## 6. Analysis and audit discipline

- Freeze the generating coefficients, curvature levels, validation sizes, replicate seeds, primary estimand, and analysis code before the confirmatory stage.
- Use paired parent, catalogue, GW-event, and validation draws across estimator arms.
- Preserve every attempted replicate and explicit failure reason; never silently replace failed fits.
- Distinguish fixed-truth repeated-sampling coverage from simulation-based calibration. They answer different questions.
- Keep truth-only scoring payloads separate from inference payloads and independently verify that no truth enters the fit or posterior.
- Independently refit each logistic model from saved validation records, reconstruct inclusion predictions, weighted host density, selection-conditioned likelihood, posterior summaries, and all reported aggregate contrasts.
- Audit both successful and failed rows, seed lineage, nested validation samples, output schema, and checksums.
- A posterior reconstruction from persisted inputs does not by itself prove that the simulator regenerated the correct inputs. Confirmatory audit must independently regenerate the parent frame, retention process, GW observations, and validation samples from primitive seeds.

## 7. Stop/go gates

**Gate A — implementation integrity:** tests establish that \(g=0\) recovers the linear generating law, nonzero \(g\) is deterministic for fixed seeds, both estimator families are correctly labelled, and invalid coefficients are rejected. The independent auditor must reconstruct both degree-1 and degree-2 fits.

**Gate B — smoke validation:** run a small nonlinear condition; confirm all arms, seed/sample pairing, explicit failures, metadata, and independent artifact audit. A green CI job validates software checks only; it is not scientific evidence.

**Gate C — overlap and stability review:** inspect probability support, low-probability mass, retained counts, validation inclusion counts, fit-boundary rates, ESS, and extreme weights. Amend stress levels only through a dated protocol revision before the confirmatory campaign.

**Gate D — precision-based confirmatory design:** choose replicate counts from explicit Monte Carlo precision targets and available compute. Freeze contrasts and multiplicity handling before running.

**Gate E — confirmatory execution and release:** execute every declared cell, independently regenerate and audit raw artifacts, report all outcomes and deviations, and obtain domain-expert review before manuscript claims.

A failure at any gate blocks stronger conclusions. An oracle advantage does not establish that estimated inclusion weighting is useful when inclusion probabilities are unknown.

## 8. Interpretation limits

This design isolates one nonlinear inclusion-model misspecification mechanism under a known parent-frame validation design and a simplified low-redshift Euclidean GW likelihood. It does not cover covariate shift, redshift measurement error, validation nonresponse, spatially correlated survey masks, unobserved galaxy properties, realistic luminosity functions, or injection-calibrated detector selection. It also does not propagate uncertainty in fitted inclusion-model coefficients through a joint hierarchical posterior; the current plug-in procedure treats fitted probabilities as fixed.

**Current status:** paired linear/quadratic pilot implementation and independent artifact-audit code are present on the research branch. The CLI and nonlinear smoke-workflow declarations are source-consistent. The workflow now retains the nonlinear smoke output as an uploaded artifact. Exact-head CI, generated pilot artifacts, and an independent pilot artifact audit have not been verified in this session. No confirmatory factorial results are asserted.
