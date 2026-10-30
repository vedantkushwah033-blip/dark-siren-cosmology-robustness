# Pilot and confirmatory design plan for finite-validation bootstrap calibration

**Status:** design plan, not a declaration that a confirmatory experiment has been run.  
**Applies to:** nested validation SRSWOR / posterior-level resampling benchmark.

## 1. Primary scientific estimands

The analysis is conditional on each generated finite parent frame. For scenario s and validation size n:

- Finite-frame reference: \(T_s\), the posterior mean obtained by applying the prespecified host-inclusion and selection-aware toy pipeline to all N actual parent records.
- Outer estimator: \(\hat T_{s,n,r}\), the same posterior-mean functional fitted to outer SRSWOR validation draw r.
- Repeated-sampling SD: \(\sigma_{s,n}=SD_r(\hat T_{s,n,r})\), estimated using successful outer draws with failed outer draws reported separately.
- Inner uncertainty estimate: \(\hat\sigma_{s,n,r,m}\), the method-m-specific bootstrap standard error conditional on outer draw r.
- Target coverage: the fraction of successful method-specific intervals containing \(T_s\), with the number of attempted, successful and failed intervals always reported.

The outer SD is computed once per scenario/size from unique outer draws and is shared across inner methods. It must never be recalculated after filtering to a method's successful bootstrap fits.

## 2. Prespecified comparisons

The three methods answer different diagnostic questions and are not interchangeable:

1. Paired nonparametric bootstrap: baseline resampling of observed validation redshift/outcome pairs.
2. Pseudo-population SRSWOR bootstrap: exploratory design-aware comparator. The implementation's percentile interval is not the studentized finite-population-mean method in Wang, Peng & Kim (2022), and their theoretical results are not transferred to this nonlinear posterior functional.
3. Scalar finite-population correction: a diagnostic transformation of the paired-bootstrap SE and normal interval. It is not assumed valid for the nonlinear posterior-mean estimator.

Primary descriptive contrasts:
- method-specific mean inner SE divided by empirical outer SRSWOR SD;
- method-specific finite-parent target interval coverage;
- mean interval width;
- outer fit failure fraction and inner bootstrap failure fraction.

No single metric determines a winner. A narrow interval with poor coverage is not an improvement. A method with low failure rate but systematic undercoverage is not calibrated.

## 3. Pilot phase (non-confirmatory)

Run a pilot only after the independent artifact auditor passes. A practical initial pilot may use approximately 20 independently generated parent scenarios, 30 outer draws per validation size and 100 inner replicates per method/draw. These are workload-planning values, not final sample-size claims; reduce the pilot if runtime or failure rates make it impractical.

The pilot must estimate:
- outer and inner failure rates, including reasons and whether failures cluster by scenario or validation size;
- distribution of scenario-level SE/outer-SD ratios;
- scenario-to-scenario variability in interval coverage and width;
- runtime and memory per scenario/outer draw;
- frequency of boundary-limited inclusion fits, low effective sample size and extreme inverse weights;
- stability of conclusions under a larger bootstrap-replicate count on a small, prespecified subset.

Do not use the pilot to declare a method superior or to silently tune the estimator. Any algorithmic change prompted by the pilot requires a versioned protocol amendment and a fresh pilot for affected quantities.

## 4. Precision-based confirmatory sizing

Choose the number of independent parent scenarios from precision requirements, not from a desire for a large run count. The scenario is the top-level independent unit; outer draws nested within a scenario do not create independent parent scenarios.

For coverage near 0.95, a simple binomial calculation gives MCSE approximately sqrt(0.95*0.05/S) only under independent Bernoulli trials. Here interval outcomes are clustered within parent scenarios, so this formula is a lower-bound planning approximation, not the final uncertainty calculation. Use pilot-estimated between-scenario variance and scenario-level resampling or a justified hierarchical model to select S such that the confidence interval for the scenario-balanced coverage contrast has the prespecified half-width.

For the SE/outer-SD ratio and interval-width contrasts, use the pilot's scenario-level variance and a scenario-level confidence interval to determine S for a prespecified half-width. Report the resulting planned precision and sensitivity to pilot variance uncertainty.

Choose outer draws per scenario to estimate each scenario's repeated-sampling SD with adequate stability. Choose inner replicates to make Monte Carlo noise in each bootstrap SE and percentile endpoint small relative to the between-outer-draw variation. Validate this with repeated inner seeds on a prespecified subset.

The confirmatory design must state:
- scenario count and outer draws per scenario/size;
- bootstrap replicates per method/outer draw;
- exact random seed schedule;
- primary estimands and contrasts;
- target confidence level and precision half-widths;
- failure handling and denominators;
- stopping rule and compute budget;
- all deviations from the frozen protocol.

## 5. Confirmatory analysis and reporting

- Freeze code commit, environment, protocol, and seeds before execution.
- Preserve all raw attempt records and hashes.
- Audit the complete artifact independently before interpreting results.
- Report scenario-level values and scenario-balanced summaries; do not treat outer draws as independent scenarios.
- Provide uncertainty intervals for scenario-balanced method contrasts.
- Include failure rates alongside every calibration result.
- Report null, adverse and inconclusive results without selective filtering.
- Do not select a method using the same confirmatory data and then report its apparent performance as if it were prespecified.

## 6. Go/no-go criteria

Proceed to confirmatory execution only if:
1. the smoke experiment passes on all supported Python versions;
2. the structural and summary auditor passes and has adversarial tests;
3. a separate replay auditor can reconstruct parent frames, outer samples, inner resampling and posterior calculations;
4. the pilot indicates a feasible design with acceptable, explicitly reported failure rates;
5. precision targets and analysis decisions are frozen before confirmatory data are generated.

If any criterion fails, fix the pipeline and rerun the relevant smoke/pilot stage. Do not bypass an audit failure to reach a publication deadline.

## 7. Claim boundary

The experiment is a controlled synthetic, low-redshift Euclidean toy study. It is not an analysis of LVK/GWTC observations, not a new observational H0 result, and not evidence that the method corrects real galaxy-catalogue incompleteness. Any generalization to real dark-siren analyses requires a separate observational likelihood, realistic selection modeling, and validation against appropriate external benchmarks.

## 8. Scenario-level uncertainty intervals in the summary artifact

The calibration summarizer now reports nominal 95% Student-t intervals across independent parent scenarios for scenario-balanced mean SE/outer-SD ratio, finite-frame target coverage, interval width, and outer/inner failure fractions. The number of scenarios contributing to each interval is reported separately. Outer draws and inner replicates are not used as independent degrees of freedom.

These intervals are descriptive small-sample intervals under an independent, approximately well-behaved scenario-level sampling distribution. They are not exact binomial intervals for coverage and do not correct a poorly designed simulation. Confirmatory precision planning must still use pilot-estimated scenario-level variability, assess skewness/outliers, and may require a scenario bootstrap or other justified hierarchical method. If fewer than two independent scenarios contribute, interval endpoints are unavailable rather than fabricated.
