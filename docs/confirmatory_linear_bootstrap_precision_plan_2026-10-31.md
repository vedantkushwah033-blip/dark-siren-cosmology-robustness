# Confirmatory design and Monte Carlo precision plan

**Status:** prospective protocol; no confirmatory results have been generated.

## 1. Inferential target

The benchmark targets repeated-sampling coverage of nominal 95% intervals for the exact arithmetic mean of each realized finite population under simple random sampling without replacement (SRSWOR). The primary finite-population estimand is

\[
\bar Y_N = N^{-1}\sum_{i=1}^{N}Y_i.
\]

This is not the mean of the distribution used to generate a frame. Repeated samples drawn from one fixed frame quantify design-based operating characteristics conditional on that frame. Independently generated frames are a second, distinct source of variation and are needed before making claims that generalize beyond the selected frames.

## 2. Two-level replication

Use a nested design:

1. Generate independent finite frames from each prespecified population regime.
2. Within every frame, draw independent SRSWOR samples at each prespecified sample fraction.
3. Apply every interval method to the same realized sample so method comparisons remain paired.
4. Retain the frame identifier, sample identifier, exact sample indices, method-specific random seeds, endpoints, failures, and coverage indicators.

Do not treat all sample replicates from one frame as independent evidence about superpopulation variability. Report (a) within-frame design-based results and (b) between-frame variability separately. Aggregate across frames with frame-level summaries or a hierarchical analysis whose estimand and assumptions are stated in advance.

## 3. Monte Carlo precision for coverage

For a fixed frame and method, let \(\hat p\) be the empirical coverage over \(R\) successful, independent repeated samples. The simple binomial approximation is

\[
\operatorname{MCSE}(\hat p) \approx
\sqrt{\hat p(1-\hat p)/R}.
\]

At nominal coverage 0.95, an MCSE target of 0.005 requires approximately 1,900 successful repetitions per frame-cell. A conservative worst-case bound using \(p(1-p)\leq 1/4\) requires 10,000 repetitions for MCSE at most 0.005. These are precision calculations, not evidence that any interval method achieves nominal coverage.

A practical alternative is to prespecify a 95% Monte Carlo interval for coverage and require its half-width to be no greater than 0.01 in each primary cell. Use a Wilson or exact binomial interval rather than relying only on a plug-in normal approximation, especially when coverage is near zero or one. Determine the required repetitions from the chosen interval and target precision before confirmatory execution.

Coverage conditional on successful interval construction must always be accompanied by the failure rate. A method with a high conditional coverage among a selected subset of successful fits is not acceptable if its failure rate is material. Report unconditional failure frequencies and do not silently discard failed attempts.

## 4. Primary and secondary comparisons

Prespecify a small primary family:

- nominal 95% coverage against the exact finite-frame mean;
- interval width;
- interval-construction failure rate.

Treat paired differences in coverage indicators and paired interval widths as secondary comparisons. Preserve pairing by comparing methods on the same sample replicate. If inferential intervals are reported for method contrasts, resample or otherwise account for the frame and sample nesting; do not use a naive independent-samples standard error.

The Wald SRSWOR interval is a conventional design-based comparator. The ordinary percentile bootstrap is a with-replacement diagnostic comparator, not a method that reproduces the SRSWOR design. The studentized pseudo-population method is an established method from the finite-population bootstrap literature, not a newly derived procedure.

## 5. Design coverage

The current pilot includes normal, Gamma(shape 0.7), and lognormal finite frames and sample fractions 0.1 and 0.4. Before confirmatory execution, decide whether to add prespecified bounded, bimodal, heavy-tailed, or near-degenerate frames and a larger sampling fraction. Add regimes only with a clear scientific rationale and document them before examining confirmatory outcomes.

Use multiple independent frames per regime. The number of frames should be selected using pilot estimates of between-frame variability and a declared precision target for regime-level summaries. Three frames are insufficient to characterize a broad superpopulation of possible finite populations.

## 6. Analysis and reporting rules

- Freeze the protocol, code commit, seeds, population regimes, sample fractions, methods, and stopping rule before the confirmatory run.
- Use pilot results only to estimate computational cost and Monte Carlo variability; do not recycle pilot replicates into confirmatory estimates unless a combined analysis was explicitly prespecified.
- Report exact numerators and denominators for coverage, not percentages alone.
- Report Monte Carlo intervals, frame-level dispersion, and failure counts.
- Preserve raw replicate records and finite-frame values, with checksums and independent reconstruction.
- Have the independent auditor verify the design, all expected frame/sample/method cells, attempt uniqueness and completeness, sample-index replay, interval endpoints, coverage, failures, and summaries.
- If the audit fails, do not repair the data by hand. Fix the code, regenerate from the frozen protocol, and rerun the audit.

## 7. Boundary to the dark-siren analysis

Success in this linear benchmark validates neither estimated inverse-probability weights nor the nonlinear dark-siren posterior functional. The latter also contains an estimated inclusion model, weight normalization, weighted KDE smoothing, selection normalization, and posterior integration. It needs a separate simulation-based calibration and repeated-sampling study with its own estimand, truth construction, design, precision targets, and independent audit.

The current finite-validation bootstrap work must therefore remain a separate, explicitly exploratory line until its own operating characteristics are established.

## 8. Literature basis

Wang, Z., Peng, L. & Kim, J. K. (2022). “Bootstrap inference for the finite population mean under complex sampling designs.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 84(4), 1150–1174. DOI: 10.1111/rssb.12506. The paper develops a studentized pseudo-population bootstrap for several sampling designs, including SRSWOR, and establishes its properties under stated regularity conditions. Those results concern finite-population mean inference and should not be extrapolated to the nonlinear dark-siren posterior.
