# Finite-parent validation-sampling benchmark: preregistered protocol

**Status: protocol only. No experiment or scientific result is claimed by this document.**

## Motivation and estimand

The existing validation-bootstrap workflow resamples paired validation records with replacement. The synthetic data-generating design, however, draws the validation sample without replacement from a finite parent frame. The ordinary bootstrap is therefore an approximation whose calibration must be measured rather than assumed. Do not attach an arbitrary finite-population correction (FPC) to the nonlinear estimator: the logistic fit, inverse-probability weights, KDE, and posterior map do not generally admit a single scalar FPC.

The target is the **conditional repeated-sampling distribution** of the complete estimator under simple random sampling without replacement (SRSWOR) from one fixed synthetic parent frame, conditional on a fixed observed catalogue, fixed GW-distance data, fixed merger-rate specification, fixed likelihood, and fixed analysis grid. This is a design-based benchmark, not a claim about real survey sampling or real LVK data.

## Frozen design

For each independently generated parent-frame scenario:

1. Generate and freeze a finite parent frame with redshift and binary catalogue-retention indicator. Preserve the full parent frame in the *simulation truth* artifact only.
2. Generate and freeze the observed incomplete catalogue and GW-distance data once. These must not change across validation draws within that scenario.
3. Draw validation records by exact SRSWOR from the parent frame at prespecified sizes (initial development: 100, 300, 1,000; confirmatory sizes must be justified by parent size and class counts). Do not resample with replacement in the design reference.
4. For each validation draw, fit the same prespecified bounded logistic inclusion model, compute inverse-inclusion weights on the fixed observed catalogue, construct the same reflected weighted KDE, and propagate the density through the same fixed toy likelihood.
5. Record every attempted draw, including fit failures, one-class samples, boundary-limited fits, extreme weights, and numerical failures. Never discard failures from denominators.
6. Separately, for each validation draw, run the ordinary paired-record bootstrap with a prespecified number of replicates and seed. The bootstrap is nested within a design draw and is an uncertainty estimator, not another independent parent-frame sample.

## Comparators

- **SRSWOR reference:** empirical distribution across repeated independent validation samples drawn without replacement.
- **Ordinary bootstrap:** percentile uncertainty summaries from paired-record bootstrap draws, conditional on each observed validation sample.
- **Oracle diagnostic (simulation-only):** optional density/posterior obtained using the true parent inclusion law. It is a benchmark, never an inference input.
- **No-validation baseline:** catalogue-only estimate, where implemented, to contextualize the value and limitations of validation data.

## Primary outcomes

At the validation-sample level, retain raw records and estimate:

- bias and RMSE of posterior mean and median for fixed known synthetic truth;
- empirical coverage and mean width of prespecified nominal 68% and 95% intervals, with Wilson confidence intervals;
- density integrated squared error against the parent-frame target density;
- inclusion-probability MAE and calibration diagnostics over a frozen redshift grid;
- inverse-weight ESS, maximum, p95 and p99 weights;
- fit failure, boundary-limited, and posterior-computation failure rates;
- paired differences between bootstrap interval width/uncertainty summaries and the SRSWOR repeated-sampling reference, including interval inclusion/failure rates where a comparison is mathematically defined.

The bootstrap's interval for the distribution of posterior means must not be called an H0 posterior credible interval. Distinguish (a) within-draw H0 posterior intervals from (b) across-validation-sample sampling variability and (c) bootstrap estimates of that variability.

## Monte Carlo precision and staged execution

- **Stage 0 (unit tests):** exact SRSWOR uniqueness and range; deterministic seed replay; no parent-truth leakage into inference inputs; fixed catalogue/GW data hashes across draws; explicit failure accounting.
- **Stage 1 (smoke):** at least 20 parent scenarios × 20 validation draws per scenario for pipeline and artifact validation only. No scientific claims.
- **Stage 2 (pilot):** at least 100 parent scenarios × 100 validation draws per scenario, profiling runtime and failure rates. Revise the confirmatory design only before examining confirmatory outcome contrasts.
- **Stage 3 (confirmatory):** choose the number of parent scenarios and draws by a Monte Carlo precision calculation for coverage and paired contrasts; target Wilson half-width no larger than 0.02 for key coverage estimates where computationally feasible. Report actual intervals and failures, not only nominal target. Freeze seeds, code commit, scenario grid, outcomes, and analysis before execution.
- Stop or redesign if class counts make validation fits non-identifiable, positivity is materially violated, failure rates differ across arms, or Monte Carlo precision is inadequate. Do not tune caps, model complexity, or thresholds against H0 error and then report the same runs as confirmatory.

## Dependence and inference boundaries

Repeated validation samples within a parent scenario are dependent through their shared finite frame and fixed catalogue/GW data. Uncertainty summaries must respect the hierarchy: parent scenario is the top-level independent unit; validation draws are nested within scenario. Do not treat all nested draws as independent parent scenarios. Use scenario-clustered resampling or a hierarchical analysis for uncertainty in aggregate performance. Report both within-scenario and between-scenario variation.

The fixed observed catalogue means this experiment isolates validation-sampling uncertainty. It does **not** quantify total uncertainty from catalogue realization, GW-event realization, detector selection, cosmology, distance likelihood misspecification, or uncertainty in the merger-rate model. A separate outer simulation layer would be needed for those targets.

## Artifact and independent-audit requirements

Persist:
- protocol version, code commit, environment, all seeds, scenario identifiers, and parameter grid;
- parent-frame truth separately from inference payload;
- hashes of fixed observed catalogue and GW-distance arrays, verified identical across nested validation draws;
- exact validation sample indices, with uniqueness and bounds checks;
- all fitted coefficients, convergence/boundary status, probabilities, weights, ESS, density draws, posterior arrays, and failure reasons;
- raw scenario-level records and machine-readable summaries.

The auditor must independently regenerate SRSWOR index sets from declared seeds, verify uniqueness and exact replay, verify no parent-truth fields appear in inference payloads, reconstruct selected model/density/posterior quantities independently, recompute all aggregate summaries from raw records, and test tamper detection. A PASS means artifact integrity and numerical reconstruction under the stated implementation; it is not external scientific validation.

## Decision rules and interpretation

The ordinary bootstrap is considered adequately calibrated only if its uncertainty summaries reproduce the SRSWOR repeated-sampling variability within prespecified Monte Carlo tolerances across the prespecified scenario grid, without unacceptable failure or tail-weight behavior. If it fails, report the failure and use the SRSWOR reference for this synthetic design; do not invent a correction factor post hoc.

All results must be labeled synthetic and conditional on the declared design. No claim of real-catalogue validity, real LVK performance, a new H0 measurement, or superiority to established dark-siren pipelines follows from this benchmark.
