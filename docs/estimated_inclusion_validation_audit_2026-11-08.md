# Estimated catalogue-inclusion validation audit

**Status:** code-path audit and experiment protocol; no estimated-inclusion result is claimed.  
**Scope:** synthetic low-redshift Euclidean toy methodology only.

## Current benchmark limitation

The current host-density benchmark compares:
- `oracle`: parent-population redshifts;
- `catalogue_only`: observed catalogue redshifts without correction;
- `known_completeness`: observed catalogue redshifts weighted by the exact declared inclusion law.

The last arm is a known-probability oracle reference, **not an estimated inclusion model**. It does not fit coefficients from a finite validation sample, propagate coefficient uncertainty, or test misspecification. The current benchmark therefore cannot answer whether estimated inclusion probabilities improve inference.

The KDE routine accepts supplied probabilities and forms inverse weights. It does not estimate probabilities or quantify their uncertainty.

## Required experiment design

For every paired replicate:

1. Generate parent redshifts and catalogue-retention indicators using the declared inclusion law.
2. Draw a validation subset independently of retention outcomes. Validation records contain redshift and the binary retained indicator.
3. Explicitly declare whether validation objects are removed from the inference catalogue. Prefer disjoint validation and analysis samples for the primary design.
4. Fit the prespecified logistic model using validation data only.
5. Predict inclusion probabilities for analysis-catalogue objects and form IPW host-density estimates.
6. Keep parent truth and true probabilities out of non-oracle inference; use them only for scoring and oracle diagnostics.
7. Reuse identical parent/catalogue/GW realizations across all arms.

Track parent size, validation size, analysis-catalogue size, and retained/non-retained validation counts separately.

## Required arms

- Catalogue-only, uncorrected.
- Known-probability oracle, diagnostic reference only.
- Estimated-probability plug-in from finite validation data.
- Estimated-probability uncertainty-aware method, with coefficient uncertainty propagated by a derived marginalization—not by casually averaging bootstrap likelihoods.
- Deliberately misspecified fitted model when the generating law has curvature.
- Optional stabilized-weight diagnostic, reported as a distinct procedure.

Do not apply a second missing-host mixture correction on top of an already IPW-corrected host density unless a separately derived likelihood requires it.

## Fit diagnostics and failure ledger

For every replicate, record validation sample size and class counts; design-matrix rank and condition number; optimizer status, iterations, gradient norm and termination reason; coefficient/Hessian diagnostics; predicted-probability quantiles and boundary fraction; weight concentration; Kish and entropy ESS; separation warnings; and support gaps.

Use mutually exclusive statuses: `SUCCESS`, `NON_EVALUABLE`, `FAILED`, `NOT_RUN`. Retain every attempted replicate in the ledger and report unconditional failure rates. A positive probability floor does not repair structural lack of positivity. Declare any clipping or regularization in advance.

## Outcomes

Assess probability estimation separately from downstream inference:
- probability calibration, Brier score and log loss on held-out labels;
- inclusion-probability error against synthetic truth;
- parent redshift-density reconstruction error;
- normalized weight concentration and ESS;
- H0 posterior bias, interval width and repeated-sampling coverage;
- coverage conditional on successful fitting, alongside unconditional failure rates.

Report Monte Carlo uncertainty for bias/RMSE and binomial uncertainty for coverage. Numerical posterior stability is not calibration.

## Historical D-01 optimizer issue: unresolved

The original D-01 failure records and exact optimizer conditions have not been recovered in the inspected branch artifacts. No cause is assigned. Separation, rank deficiency, poor conditioning, boundary constraints and convergence failure remain hypotheses.

Before changing the historical implementation, recover its exact code revision, input data/seed, design matrix, optimizer settings, exception/status and logs. Reproduce the failure unchanged, then compare a revised fit on identical inputs. Do not describe a new implementation as a historical fix without that reproduction.

## Go/no-go

Do not proceed to posterior claims until the validation sampling design is reproducible; ordinary, imbalanced, rank-deficient and separated cases are tested; every replicate has a terminal status and diagnostics; probability quality is assessed out-of-sample; the inclusion correction is applied exactly once; and all arms are paired on identical realizations.

**Conclusion:** the existing known-completeness arm is a useful oracle reference, but it is not evidence that finite-sample estimated inclusion probabilities work. D-01 remains unresolved, and no estimator-superiority claim is supported.


## Implementation added after the code-path audit

The research branch now contains an initial probability-stage implementation:

- `estimated_inclusion.py`: unpenalized logistic MLE with linear or quadratic standardized-redshift terms, analytic score, observed-information covariance, and explicit diagnostics. It returns `SUCCESS`, `NON_EVALUABLE`, or `FAILED`; it does not silently regularize separated fits.
- `validation_design.py`: simple random sampling without replacement from the entire parent frame, independent of retention labels. Validation objects are removed from the inference catalogue. The returned analysis catalogue contains only retained object IDs and redshifts; hidden inclusion probabilities are not passed to the fitter.
- `inclusion_diagnostics.py`: held-out Brier score and log loss, plus unmodified inverse-probability weight concentration summaries (Kish ESS, entropy ESS, maximum normalized weight, and top-1%/5% weight mass).
- `estimated_inclusion_benchmark.py` and `scripts/run_finite_validation_probability_benchmark.py`: paired first-stage probability benchmark over validation sizes 100, 300, and 1,000; zero/nonzero quadratic curvature; and linear/quadratic fitted models. Every planned cell produces an attempt row, including fit failures. Scoring uses held-out parent labels and true probabilities only after fitting, for synthetic evaluation.

### Interpretation and limitations of this first implementation

1. This is a probability-estimation and weight-diagnostic benchmark only. It does not yet propagate estimated inclusion probabilities into the host KDE, event likelihood, or H0 posterior.
2. The primary validation design is disjoint: validation objects are excluded from the analysis catalogue. Held-out labels from the non-validation parent complement are used only for simulation scoring, never for fitting.
3. The fitted model is an unpenalized logistic MLE. The separation coefficient threshold and Hessian condition threshold are diagnostic guardrails, not universal statistical tests for separation.
4. The reported inverse-Hessian covariance is a local plug-in covariance under the fitted model. It is not finite-population corrected and is not yet used for coefficient uncertainty propagation.
5. The quadratic model is correctly specified for the chosen synthetic logistic DGP; a linear model under nonzero curvature is deliberately misspecified. This is a controlled design, not evidence about a real catalogue.
6. The benchmark uses paired population seeds and validation splits across model degrees within each scenario. Validation sizes are nested in the same generated parent scenario but are independently sampled by size-specific seeds; they are not claimed to be nested validation samples.
7. The smoke benchmark is not confirmatory. No performance ranking, H0 calibration, posterior coverage, or estimator-superiority claim is supported until the probability-stage CI and independent artifact audit pass and downstream likelihood propagation is separately derived and validated.

The new code and tests are included in the research branch's CI workflow. Exact-head CI must complete successfully before these implementation tests are described as passing.
