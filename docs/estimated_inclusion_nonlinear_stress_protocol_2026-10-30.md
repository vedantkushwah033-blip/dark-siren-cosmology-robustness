# Estimated-inclusion nonlinear-selection stress protocol

**Protocol date:** 30 October 2026  
**Status:** prespecified design proposal; not yet executed as a confirmatory campaign  
**Scope:** controlled synthetic dark-siren host-catalogue inference only.

## 1. Question and estimand

The experiment asks when inverse-inclusion weighting learned from an external validation sample stops being reliable when the true catalogue-retention law is nonlinear but the fitted estimator remains a linear logistic model.

The primary inferential estimand is the paired change in absolute error of the posterior mean of the injected synthetic $H_0$ between (i) validation-estimated inclusion weighting and (ii) catalogue-only inference, on the same simulated parent population and GW-event realization. Oracle weighting is a diagnostic upper-bound benchmark because it receives the true retention probabilities; it is not an implementable estimator.

No result from this protocol is a real-observation $H_0$ measurement, a claim about LVK catalogue completeness, or evidence that the synthetic selection law represents a particular survey.

## 2. Data-generating mechanism

For parent-frame redshift $z$, define the fixed standardized coordinate

$$x = \frac{z-(z_{\min}+z_{\max})/2}{(z_{\max}-z_{\min})/2}.$$

Generate retention indicators as independent Bernoulli draws with

$$P(R=1\mid z)=\operatorname{logit}^{-1}\{a+b x+g(x^2-1/3)\}.$$

The quadratic term is a deliberate stress perturbation. The subtraction of $1/3$ centers $x^2$ under a uniform distribution on $[-1,1]$ only; it does not generally center the term under the simulated parent-redshift distribution. The fitted model remains $\operatorname{logit}P(R=1\mid z)=\alpha+\beta x$. Thus $g=0$ is the correctly specified reference and $g\ne0$ is model misspecification.

The validation sample is sampled independently by simple random sampling from the complete parent frame and includes both reference redshift and retention indicator. This is an explicit design assumption. It is not equivalent to drawing validation objects only from the retained catalogue.

## 3. Factorial conditions

Use the following staged design rather than immediately launching a large grid:

| Factor | Prespecified levels | Purpose |
|---|---|---|
| Selection curvature $g$ | $0, -1.0, +1.0, -2.0, +2.0$ | Correct specification and symmetric curvature stress |
| Validation size $m$ | $100, 300, 1000$ | Audit-sample information |
| Parent population | $4,000$ initially | Match pilot scale |
| GW events | $30$ initially | Hold event information fixed |
| Replicates | 20 smoke; then 200 per cell only after runtime and failure review | Separate pipeline validation from inferential precision |
| Randomness | paired parent/GW draws within replicate; independent audit sampling stream | Make arm contrasts paired without leaking truth to fitted arms |

Keep the existing intercept and slope fixed for the first factorial screen. Before interpreting a curvature level, numerically tabulate the true probability range and the fraction of parent objects with probabilities below 0.05 and 0.10. If a cell creates near-positivity violations or almost no retained objects, report it as a separate severe-overlap stress condition rather than pooling it with moderate-curvature results.

The pilot CLI currently supports one curvature coefficient per invocation. Run separate invocations with distinct recorded seeds and coefficient values; do not concatenate outputs without adding condition identifiers and preserving each source manifest.

## 4. Arms

1. **Catalogue-only:** no inclusion correction.
2. **Oracle:** inverse weighting with the true retention probabilities; diagnostic benchmark only.
3. **Estimated-$m$:** fit the existing linear logistic model to the first $m$ independently sampled validation records, predict inclusion probabilities for retained catalogue objects, and apply the current untruncated inverse-probability weights.

Do not silently clip probabilities or cap weights. The estimator currently deliberately applies no probability floor or weight cap. Any stabilized estimator must be a separately named arm with a prespecified rule, not an after-the-fact repair.

## 5. Outcomes

### Primary
- Paired difference in absolute posterior-mean $H_0$ error: estimated arm versus catalogue-only.
- Report mean paired difference, standard error, and a paired bootstrap interval over simulation replicates. Positive reduction means the estimated arm has lower absolute error.

### Key secondary
- Signed posterior-mean bias and RMSE.
- Empirical 68% interval coverage with Wilson interval and interval width.
- Fit failure rate and coefficient-boundary-hit rate.
- Inclusion-probability MAE evaluated over the complete synthetic parent frame (truth used for scoring only).
- Inverse-weight ESS and ESS fraction; maximum, 95th and 99th percentile inverse weights.
- Host-density integrated squared error.

Report attempted, successful, and paired-success replicate counts for every contrast. A complete-case paired contrast is conditional on both arms fitting successfully; it must never conceal differential failure rates.

## 6. Analysis discipline

- Freeze seeds, factor levels, primary contrast, and analysis code before the 200-replicate stage.
- Use the same parent population and GW events across arms within each replicate.
- Keep validation samples nested by size within replicate, as in the current pilot.
- Preserve failed fits as explicit rows with reasons; do not replace failures or silently drop replicates.
- Report all factor cells, including unfavorable and non-evaluable cells.
- Treat the 20-replicate stage as a pipeline/runtime check only.
- For the 200-replicate stage, report uncertainty intervals and avoid declaring success from nominal significance alone.
- Treat the five curvature levels and three validation sizes as a factorial family. Label cell-wise contrasts exploratory unless a multiplicity strategy is specified before unblinding; report the full pattern, not only the best cell.
- Independently reconstruct persisted summaries from raw replicate rows and verify checksums before interpreting results.

## 7. Stop/go gates

**Gate A — implementation:** unit tests establish that $g=0$ reproduces the baseline design, nonzero $g$ is deterministic for fixed seeds, invalid coefficients are rejected, and metadata records the true data-generating family and coefficients.

**Gate B — smoke:** run 20 replicates per selected condition; verify all requested arms are present, failure rows are retained, validation sizes and seeds are correct, and the independent artifact auditor passes.

**Gate C — overlap review:** inspect the true probability range, low-probability mass, retained catalogue counts, fitted boundary hits, and weight-tail diagnostics. Revise the stress grid only in a documented protocol amendment, before the larger run.

**Gate D — confirmatory-scale simulation:** only after A–C, run 200 replicates per cell, independently audit raw records and summaries, and publish all cells including failures.

A failure at any gate blocks stronger claims. In particular, an oracle advantage does not establish that estimated inclusion weighting is useful when the inclusion law is unknown.

## 8. Interpretation limits

This experiment isolates one misspecification mechanism under a known parent-frame validation design. It does not cover covariate shift, redshift measurement error, validation nonresponse, spatial clustering, selection depending on unobserved galaxy properties, correlated survey masks, uncertain GW detection efficiency, or a realistic galaxy luminosity function. Those require separate studies and cannot be inferred from this factorial screen.

**Current status:** protocol only. No factorial results are asserted here.
