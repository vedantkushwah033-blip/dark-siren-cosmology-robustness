# Repeated-sampling benchmark for finite-population mean intervals

**Status:** pilot infrastructure added; results are not confirmatory until exact-head CI, a pilot run, and independent artifact audit complete.

## Purpose and estimand

This experiment isolates a linear finite-population estimand before the project's nonlinear dark-siren posterior-mean functional. Each scenario creates one fixed finite frame of size N. Repeated samples are drawn by simple random sampling without replacement (SRSWOR), and every interval is evaluated against that frame's exact arithmetic mean. The frame mean is not the generating-distribution mean.

Three interval procedures are compared:

1. **Studentized pseudo-population SRSWOR bootstrap.** Generate a multinomial pseudo-frame of size N from the observed n sample records, draw n units without replacement from it, studentize the bootstrap sample mean around the pseudo-frame mean, and invert the bootstrap-t quantiles.
2. **SRSWOR Wald interval.** Use the conventional variance estimate (1 - n/N) s²/n and a standard-normal 97.5th percentile.
3. **Ordinary percentile bootstrap.** Resample the observed sample with replacement and take percentile endpoints. This is a diagnostic comparator; it does not reproduce the finite-population sampling design.

## Population regimes

The runner cycles through fixed frames generated from normal, skewed Gamma(shape=0.7), and lognormal distributions. It varies sample fractions through prespecified sample sizes. Each scenario is a distinct fixed frame; repeated samples within a frame are not independent new populations. The current defaults are a pilot configuration, not a power/precision calculation.

## Reproducibility and audit

The runner writes:
- finite_population_frames.csv: every finite-frame value, unit index, and frame SHA-256;
- linear_finite_population_bootstrap_replicates.csv: one row per method and repeated sample, exact SRSWOR indices, sample and bootstrap seeds, target, estimate, interval endpoints, coverage flag, and failure classification;
- linear_finite_population_bootstrap.csv: scenario/sample-size/method summaries;
- manifest.json: design, status, and checksums.

The independent auditor reconstructs each finite-frame hash and exact mean, replays SRSWOR sample indices, independently regenerates ordinary-bootstrap and studentized pseudo-population intervals, checks interval endpoints and coverage flags, and re-derives summary coverage, widths, failure accounting, and conditional Monte Carlo standard errors.

Coverage is reported conditional on successful interval construction, alongside the failure count. This distinction matters: silently dropping failures without reporting them can make a procedure look better than it is. The conditional coverage MCSE is a binomial approximation and does not include uncertainty from selecting a finite frame from a superpopulation.

## Interpretation limits

- A pilot's apparent coverage near 95% is not proof of nominal coverage. Monte Carlo precision must be adequate and reported.
- The three fixed frames do not represent all possible finite populations. Heavy-tail, lattice, bounded, highly skewed, and near-degenerate cases may need separate prespecified regimes.
- The Wald interval is a useful baseline, not a universal gold standard.
- The ordinary percentile bootstrap is deliberately a with-replacement diagnostic, not an equivalent SRSWOR design bootstrap.
- The studentized pseudo-population method is based on established finite-population bootstrap literature. The implementation is a computational benchmark, not a new theorem.
- Even successful calibration for a linear finite-population mean does not establish validity for the dark-siren posterior mean, which includes estimated inclusion probabilities, inverse weighting, weighted KDE smoothing, selection normalization, and posterior integration.

## Literature

Wang, Z., Peng, L. & Kim, J. K. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 84(4), 1150–1174. DOI: 10.1111/rssb.12506. The paper's SRSWOR method generates a bootstrap finite population, reapplies the sampling design, and studentizes; its theoretical conclusions are for finite-population means under stated regularity conditions.

## Execution

Pilot example:

    python scripts/run_linear_finite_population_bootstrap_benchmark.py --scenarios 3 --repetitions 100 --population-size 1000 --sample-sizes 100 400 --bootstrap-replicates 199 --seed 20261031 --output-dir artifacts/linear_finite_population_bootstrap_pilot
    python scripts/audit_linear_finite_population_bootstrap_benchmark.py --output-dir artifacts/linear_finite_population_bootstrap_pilot --report artifacts/linear_finite_population_bootstrap_pilot/audit_report.json

Do not use these defaults as a confirmatory design without first estimating Monte Carlo error and choosing replication counts to meet a declared precision target.
