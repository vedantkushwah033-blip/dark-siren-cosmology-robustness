# Confirmatory multi-frame design for the linear finite-population bootstrap benchmark

**Status:** prespecified design proposal; not executed; not a result.  
**Applies to:** the linear finite-population arithmetic-mean benchmark on the active research branch.  
**Does not establish:** calibration of nonlinear dark-siren posteriors, a new bootstrap theorem, or a preferred interval method.

## 1. Why a multi-frame design is necessary

The v1 pilot generated one finite population for each distribution regime and repeatedly sampled that same frame. Its repeated-sampling coverage therefore describes performance conditional on those three realized frames. It does not estimate how performance varies across independently generated finite populations from the same superpopulation regime.

The next experiment must distinguish two sources of variation:

1. **Within-frame sampling variation:** repeated SRSWOR samples from one fixed finite population.
2. **Between-frame variation:** changes in the finite-population shape, tail realization, skewness, and finite-population moments across independent frames.

These are different inferential targets. Pooling all replicate rows and treating them as independent would understate uncertainty when replicates share a frame.

## 2. Target estimands

Let frame j be generated independently under scenario s, with exact finite-population mean theta_js. For method k, nominal level q, and outer sample replicate r, define I_jrks as the indicator that the interval is successfully constructed and covers theta_js.

### Conditional-on-frame coverage

C_jks(q) = sum_r I_jrks / R_jks_success,

with the denominator and interval-construction failures reported alongside it. This is conditional coverage among successfully constructed intervals, not unconditional operating success.

### Unconditional successful-and-covering probability

U_jks(q) = sum_r [interval_success_jrks AND coverage_jrks] / R_attempted.

This estimand penalizes methods that fail to construct intervals. Report interval-construction success probability separately so the decomposition is transparent.

### Across-frame average

For a clearly defined frame-generating distribution, estimate expected frame-level performance by averaging frame-level estimates with equal frame weights, unless a different weighting target is explicitly justified before analysis. Also report the distribution of frame-level performance (median, quantiles, and range). Do not let frames with more successful intervals silently receive greater weight.

### Paired method contrasts

For two methods k and l evaluated on the same outer sample replicate, define paired differences in:
- coverage indicators;
- successful-and-covering indicators;
- interval widths when both intervals exist;
- absolute estimation error (where relevant).

The same generated frame and outer SRSWOR sample must be reused across methods. Pairing reduces Monte Carlo noise and makes the comparison interpretable.

## 3. Design hierarchy

The confirmatory data structure should be:

scenario -> independent finite frame -> sample size -> outer replicate -> method -> inner bootstrap attempt.

Every level must have an explicit stable identifier and reproducible seed lineage.

### Scenario dimensions

Retain the v1 baseline regimes:
- Normal;
- Gamma with shape 0.7;
- Lognormal.

Add prespecified stress regimes only as separate labelled scenarios, not as post-hoc replacements:
- stronger right-tail/heavy-tail finite frames;
- contamination or rare-large-value mixtures;
- finite-frame skewness strata based on realized frame diagnostics.

Do not call a finite frame "heavy-tailed" solely because it has a high sample skewness. Define the generating distribution and report realized finite-frame summaries.

### Frame generation

Generate at least 50 independently seeded frames per distribution scenario for the confirmatory stage, subject to a computational pilot and precision-based refinement. The number 50 is a design starting point, not a universal guarantee of adequate power. Save:
- exact finite-frame values or a deterministic generator specification;
- frame seed and frame hash;
- N, exact finite-frame mean, variance and prespecified shape diagnostics;
- generating scenario and code version.

Do not discard inconvenient frames after seeing method performance. Any frame exclusions must be based on prespecified validity checks independent of interval results and must be fully logged.

### Within-frame sampling

For each frame and sample size, use repeated SRSWOR draws. Preserve the same outer sample indices and seed across all interval methods for a replicate. The two v1 sampling fractions (n/N = 0.10 and 0.40) can remain as baseline conditions, with any added fractions clearly marked as additional factors.

A staged planning target of 1,000 outer draws per frame and cell may be used for workload estimation, but final R must be selected from a documented Monte Carlo precision analysis. The total run size is M x R x number_of_cells x number_of_methods and must be estimated before launching.

## 4. Monte Carlo precision and decision rules

For a simple independent Bernoulli coverage estimate near p, the approximate Monte Carlo standard error is sqrt[p(1-p)/R]. For example, at p=0.95:
- R=100 gives MCSE about 0.0218;
- R=1,000 gives MCSE about 0.0069;
- R=2,000 gives MCSE about 0.0049.

These calculations describe within-cell binomial noise only. They do not account for between-frame heterogeneity. The multi-frame analysis must estimate both components.

Before confirmatory execution:
1. Set a target precision for the across-frame average coverage and paired coverage differences.
2. Run a separate design simulation to estimate how many independent frames and within-frame replicates are needed.
3. Use frame-clustered uncertainty intervals or a hierarchical model that preserves the two-level design.
4. Report uncertainty for all headline estimates; do not rank methods by point estimates alone.
5. If the precision target is not met, label the result inconclusive rather than selecting the most favourable method.

For method comparisons, predefine the primary contrasts and apply a multiplicity strategy (for example, a prespecified family-wise adjustment or simultaneous intervals). Secondary distributional and sample-size analyses should be labelled accordingly.

## 5. Required metrics

For each method, scenario, sample size, and nominal level, report:
- attempted outer intervals;
- successful and failed outer intervals;
- interval-construction success probability;
- conditional coverage among successful intervals;
- unconditional successful-and-covering probability;
- coverage uncertainty, with frame-level clustering respected;
- mean, median and quantiles of interval width;
- median absolute error and RMSE of the point estimator;
- paired coverage and width contrasts against prespecified references;
- inner bootstrap attempts, successes, failures and failure-code counts;
- runtime and memory diagnostics where measured.

Do not interpret narrower intervals as better unless coverage and failure behaviour are considered jointly. A narrow interval with poor coverage is not an improvement.

## 6. V2 artifact and audit requirements

The v2 schema must preserve the full hierarchy and all attempt-level information. For the studentized pseudo-population method, each outer replicate must include:
- inner attempts, successes and failures;
- canonical failure-code counts;
- deterministic seed lineage;
- enough replay information for the auditor to regenerate every inner draw;
- a stable hash for any serialized attempt-level record or canonical failure summary.

For methods without an inner bootstrap, store explicit NOT_APPLICABLE values rather than zeroes that could be mistaken for measured attempts.

The independent auditor must:
1. regenerate every finite frame from its seed and verify its hash;
2. regenerate every outer SRSWOR sample;
3. regenerate every inner multinomial pseudo-population and SRSWOR draw;
4. independently classify each inner replicate as success or failure;
5. reconstruct studentized statistics, quantiles and final endpoints;
6. reconcile attempted = successful + failed and failure-code totals;
7. reconstruct every summary from raw replicate records;
8. reject missing, duplicate, altered or impossible records even if checksums are refreshed.

The auditor must not import the benchmark runner's interval-construction or summarization functions as its source of truth. Shared low-level deterministic primitives should be minimized and documented.

## 7. Interpretation boundaries

- The v1 pilot remains a historical exploratory artifact and must not be silently overwritten.
- The v2 multi-frame study remains a simulation benchmark for a linear finite-population mean.
- The result cannot establish validity for nonlinear estimators, Bayesian posterior means, credible intervals, hierarchical GW likelihoods, or real gravitational-wave data.
- The theoretical results of Wang, Peng and Kim (2022) are conditional on their stated regularity conditions. An empirical implementation check does not prove those conditions or establish a new theorem.
- No method winner should be declared unless the confirmatory precision and audit gates pass.

## 8. Release gates

A confirmatory release requires all of the following:
1. Frozen protocol and scenario definitions before confirmatory results are inspected.
2. V2 artifact schema and complete internal-attempt logging.
3. Independent adversarial tests for tampered counts, failure histograms, seeds, endpoints, frame identifiers and summaries.
4. Multi-frame generation and independent audit.
5. Exact-head CI passing on all supported Python versions.
6. Reproducibility instructions and immutable artifact hashes.
7. Manuscript tables generated from audited artifacts, not manually transcribed.
8. Clear separation of conditional coverage, unconditional success-and-covering probability, and between-frame variability.

## 9. Relation to the published theoretical method

Wang, Peng and Kim (2022) develop a studentized pseudo-population bootstrap for several complex sampling designs, including SRSWOR, and establish second-order accuracy under explicit regularity conditions using Edgeworth expansions. This repository's benchmark is an implementation and empirical diagnostic exercise, not a proof of those results. Its purpose is to expose finite-sample operating behaviour, implementation failures, and reproducibility limits in a controlled design.

Reference: Wang, Z., Peng, L. & Kim, J. K. (2022). “Bootstrap inference for the finite population mean under complex sampling designs.” *Journal of the Royal Statistical Society: Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506
