# A reproducible pilot comparison of intervals for a finite-population mean under simple random sampling without replacement

**Manuscript status:** exploratory simulation note; not confirmatory; not submitted.  
**Analysis source commit:** `74e6dacd6690875fdeb6c6d6fc0a1fcaa40fde65`.  
**Artifact:** Python 3.12 `numerical-diagnostics-py3.12` from workflow run `41126411371`.  
**Audit:** independent outer-artifact audit PASS (3 frames, 1,800 method attempts, 18 summary cells).  
**Important audit qualification:** the v1 artifact does not serialize row-level internal bootstrap-attempt failure histograms. See the v2 remediation protocol. This note does not claim full internal-attempt accounting.

## Abstract

Finite-population inference under sampling without replacement requires attention to the sampling design and to the finite-population correction. We report a reproducible pilot comparison of three nominal 95% intervals for the arithmetic mean of a realized finite frame: a studentized multinomial pseudo-population bootstrap followed by simple random sampling without replacement (SRSWOR), a normal/Wald interval using the SRSWOR variance estimator, and an ordinary percentile bootstrap. Three synthetic frames of size 1,000 were generated from normal, Gamma(shape 0.7), and lognormal distributions. For each frame, 100 independent SRSWOR samples were drawn at sample sizes 100 and 400, and the methods were compared on the same outer samples. The pilot produced 1,800 method-level interval attempts across 18 cells. Empirical coverage ranged from 0.89 to 0.98. The studentized pseudo-population and Wald intervals both had 0.89–0.90 coverage in the normal-frame, 40%-sampling cell; the ordinary percentile interval had 0.98 coverage there and was wider. These estimates are imprecise: each cell contains only 100 outer repetitions, and only one finite frame represents each distributional regime. An independent artifact audit reconstructed the three finite frames, the complete outer attempt grid, and all reported interval results. The artifact does not yet preserve independently auditable inner-bootstrap failure histograms, so an attempt-level v2 release is required. The results are descriptive and do not identify a generally superior method. They provide a pilot foundation for a multi-frame, precision-designed study.

**Keywords:** finite population; SRSWOR; bootstrap; studentization; coverage; Monte Carlo precision; reproducibility.

## 1. Introduction

Sampling without replacement from a finite frame differs from independent sampling from an unlimited population. In simple random sampling without replacement, the variance of the sample mean contains a finite-population correction (FPC). Bootstrap procedures designed for independent and identically distributed observations do not automatically reproduce this design feature.

A pseudo-population approach offers one way to represent finite-population sampling variability. In a common construction, a bootstrap finite population is generated from the observed sample and the original sampling design is reapplied. Wang, Peng and Kim (2022) develop bootstrap inference for finite-population means under complex sampling designs, including design-aware constructions and studentization. Their theoretical results are tied to the estimand and regularity conditions studied in that work; they do not establish validity for arbitrary nonlinear estimators.

This note reports a small computational pilot designed to check an implementation of a studentized pseudo-population SRSWOR interval and compare it with two transparent benchmarks. The study has three limited aims:

1. establish a deterministic, reproducible repeated-sampling implementation for a finite-frame arithmetic mean;
2. inspect pilot coverage and interval-width behaviour under three contrasting synthetic frame shapes and two sampling fractions; and
3. identify design, audit, and Monte Carlo limitations that must be addressed before confirmatory work.

The study is not a new bootstrap theorem, does not claim a new general method, and is not a dark-siren posterior experiment.

## 2. Estimand and design

For each realized finite frame (U_s=(y_{s1},\ldots,y_{sN})), the estimand is its exact arithmetic mean,
\[
\bar Y_s=N^{-1}\sum_{i=1}^{N}y_{si}.
\]
The target is conditional on the realized frame. It is not the superpopulation mean used to generate the frame.

A simple random sample without replacement of size (n) is drawn from each frame. The two sample sizes are (n=100) and (n=400), with (N=1000), giving sampling fractions 0.10 and 0.40. For each frame/sample-size pair, 100 repeated samples are generated. All three methods use the same outer sample in each repetition, creating a paired comparison.

Three frames are generated using deterministic child seeds from root seed 20261031:
- normal: independent standard normal values;
- Gamma: shape 0.7 and scale (1/0.7);
- lognormal: log-scale mean 0 and log-scale standard deviation 1.

There is one realized frame per regime. Therefore, the pilot cannot estimate variability across independently generated frames within a regime.

## 3. Interval procedures

### 3.1 Wald SRSWOR interval

Let (s^2) be the sample variance. The estimated design standard error is
\[
\widehat{SE}_{W}=\sqrt{(1-n/N)s^2/n}.
\]
The interval is the sample mean plus or minus 1.9599639845 estimated standard errors. This is a normal approximation, not an exact finite-sample interval for arbitrary finite frames.

### 3.2 Ordinary percentile bootstrap

The ordinary comparator resamples (n) observations with replacement from the outer sample, computes the bootstrap sample mean, and uses the 2.5th and 97.5th percentiles. It does not explicitly reapply SRSWOR to a size-(N) pseudo-population and is included as a conventional contrast.

### 3.3 Studentized pseudo-population SRSWOR bootstrap

For each bootstrap replicate, multinomial counts of total (N) are generated over the (n) observed sample records. Expanding the sample according to those counts creates a size-(N) pseudo-population. A size-(n) SRSWOR sample is then drawn from that pseudo-population. The bootstrap statistic is centered at the pseudo-population mean and studentized by the replicate SRSWOR standard-error estimate. The interval inverts the empirical bootstrap-t quantiles around the observed sample mean.

The pilot uses 199 inner bootstrap replicates for each outer sample. The implementation records deterministic seeds and replay arrays in the core function. The benchmark's v1 CSV artifact, however, does not serialize an independently checkable inner failure-code histogram; this limitation is addressed in Section 7.

## 4. Reproducibility and independent audit

The prespecified pilot used:
- root seed: 20261031;
- finite-frame size: 1,000;
- frames: 3;
- repeated outer samples: 100 per frame and sample size;
- sample sizes: 100 and 400;
- inner bootstrap replicates: 199;
- methods: 3.

This yields 3 × 2 × 100 × 3 = 1,800 outer method-level attempts and 3 × 2 × 3 = 18 summary cells.

The independent auditor checked the manifest schema and pilot status, regenerated the finite frames from the root seed, verified frame hashes and unit indices, reconstructed the SRSWOR sample seeds and indices, independently recomputed interval endpoints and coverage, and checked the complete method-by-frame-by-sample-size summary grid. The audit report returned PASS for three finite frames, 1,800 reconstructed outer attempts, and 18 summary cells.

The exact workflow run and artifact are linked from the repository's pilot-results report. The source commit and checksums are recorded there. Reproducibility of the stored outer results is distinct from scientific validity of the method.

## 5. Results

All 18 cells had 100 successfully constructed outer intervals and zero outer interval-construction failures. The following table reports empirical coverage and mean interval width. Widths are in the units of each frame and should not be compared across distributional regimes as though the scales were identical.

| Frame | Sampling fraction | Studentized pseudo-population coverage (mean width) | Wald coverage (mean width) | Ordinary percentile coverage (mean width) |
|---|---:|---:|---:|---:|
| Normal | 0.10 | 0.96 (0.366) | 0.97 (0.372) | 0.94 (0.380) |
| Normal | 0.40 | 0.89 (0.147) | 0.90 (0.151) | 0.98 (0.190) |
| Gamma(0.7) | 0.10 | 0.93 (0.449) | 0.96 (0.434) | 0.96 (0.444) |
| Gamma(0.7) | 0.40 | 0.94 (0.172) | 0.94 (0.175) | 0.98 (0.222) |
| Lognormal | 0.10 | 0.95 (0.788) | 0.93 (0.720) | 0.95 (0.735) |
| Lognormal | 0.40 | 0.95 (0.284) | 0.95 (0.284) | 0.97 (0.354) |

The clearest pilot deviation from nominal 95% coverage occurs in the normal frame at sampling fraction 0.40: the studentized pseudo-population and Wald procedures have empirical coverage 0.89 and 0.90, respectively. The ordinary percentile comparator has coverage 0.98 and a larger mean width in that cell. The ordinary percentile procedure also produces wider intervals and higher empirical coverage in both Gamma cells at sampling fraction 0.40.

These observations are descriptive. For a coverage probability near 0.95 and 100 independent repetitions, the binomial Monte Carlo standard error is approximately 0.022. Thus, the apparent differences of a few percentage points are not precise. The 0.89–0.90 observations warrant targeted follow-up, but do not establish systematic undercoverage.

## 6. Interpretation

The pilot does not support declaring one method superior. It demonstrates that the methods can behave differently under a fixed finite-frame design, and that sampling fraction and frame shape should be represented in a larger evaluation. It also shows why interval width must be reported alongside coverage: higher empirical coverage in a small pilot may be associated with wider intervals rather than better-calibrated efficiency.

The Wald interval provides a useful design-aware baseline but relies on a normal approximation. The ordinary percentile bootstrap is not a design-respecting finite-population bootstrap and should not be treated as such. The studentized pseudo-population construction explicitly generates a finite pseudo-frame and reapplies SRSWOR, but this pilot is not a theoretical validation of its coverage.

The normal-frame result at 40% sampling should be investigated with more independent frames and repetitions. A confirmatory study should distinguish variation across frames from repeated-sampling variation within a frame. Treating all repeated samples from one frame as independent superpopulation replications would overstate the breadth of evidence.

## 7. Audit limitations and required v2 work

The v1 independent audit verifies the outer artifact: frames, sample seeds and indices, outer method attempts, interval endpoints, coverage indicators, and summary cells. However, the outer CSV does not store the studentized method's internal bootstrap attempt/success/failure totals and failure-code histogram. The auditor replays the interval but does not compare a serialized row-level accounting of every inner bootstrap replicate.

Accordingly, the statement that all 18 cells had 100 successful intervals and zero outer interval-construction failures must not be interpreted as evidence that every inner bootstrap replicate succeeded. The inner bootstrap primitive does retain success masks and failure codes in memory, but v1 does not provide the full independently auditable row-level record.

A v2 release must serialize inner attempts, successes, failures, and stable failure-code counts for each studentized outer attempt; retain explicit not-applicable values for methods without inner bootstrap; independently reconstruct and compare every count; and include checksum-refreshed adversarial tampering tests. The existing v1 pilot should remain immutable as a historical pilot rather than be retroactively relabelled as v2.

## 8. Confirmatory design requirements

A larger study should be frozen before execution. At minimum it should:
1. generate multiple independent finite frames per population regime;
2. select outer repetitions using explicit Monte Carlo precision targets;
3. use paired outer samples across methods;
4. report exact coverage numerators and denominators;
5. separate between-frame and within-frame variation;
6. record and audit every inner bootstrap attempt and failure;
7. include additional finite-frame shapes and controlled tail/heavy-skew regimes;
8. evaluate sensitivity to confidence level and sampling fraction;
9. report interval width and coverage jointly; and
10. publish seeds, code commit, software environment, raw artifacts, manifest, and independent audit report.

A coverage target near 95% with an absolute binomial Monte Carlo standard error around 0.005 requires approximately 1,900 successful independent repetitions per cell under a simple binomial approximation. This is a planning calculation, not a guarantee of precision for a hierarchical multi-frame design. The final design must account for between-frame variation and expected failures.

## 9. Relation to the dark-siren project

This benchmark addresses only a linear finite-population arithmetic mean. It does not validate inverse-probability host weighting, estimated catalogue completeness, a nonlinear dark-siren posterior mean, posterior credible intervals, or an estimate of the Hubble constant. A separate nonlinear simulation must define its own finite-frame target, repeated-sampling design, uncertainty estimators, failure accounting, precision targets, and independent reconstruction.

The distinction is central: design-aware resampling for a linear mean cannot simply be transferred to a nonlinear fitted posterior functional.

## 10. Limitations

The study has only three realized finite frames, one per distributional regime; 100 outer repetitions per cell; 199 inner replicates; and no confirmatory protocol. The normal approximation may be imperfect for skewed or heavy-tailed finite frames. The ordinary percentile bootstrap does not reproduce the without-replacement design. The pseudo-population procedure is evaluated empirically rather than justified by a new theorem. Finally, the v1 artifact's inner-attempt serialization limitation prevents a full attempt-level audit claim.

## References

Wang, Z., Peng, L., & Kim, J. K. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *Journal of the Royal Statistical Society: Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506

## Reproducibility record

- Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`
- Source commit for pilot: `74e6dacd6690875fdeb6c6d6fc0a1fcaa40fde65`
- Workflow run: `41126411371`
- Audit status: PASS for the outer artifact checks described above.
- Artifact fingerprints and exact audit scope: `docs/audited_linear_finite_population_bootstrap_pilot_results_2026-10-31.md`.
- Required follow-up: `docs/inner_bootstrap_attempt_audit_gap_and_v2_protocol_2026-10-31.md`.


## 11. Subsequent v2 audit-engineering update (not a new result)

After the v1 pilot described above, the benchmark runner and independent auditor were extended to a v2 artifact contract. For each studentized outer attempt, the runner now serializes the number of inner pseudo-population draws attempted, successful, and failed, together with a canonical sorted failure-code histogram. Methods without an inner bootstrap carry the explicit value `NA`; an outer studentized failure retains the inner accounting when the primitive returned a result.

The auditor independently replays the multinomial pseudo-population counts and SRSWOR selections from the recorded deterministic seed, recomputes the replicate variance/standard-error success classification, and compares the reconstructed counts and failure histogram with the CSV fields. The benchmark test suite also includes checksum-refreshed tampering tests for altered inner counts and non-`NA` fields on non-studentized methods.

This is an implementation and provenance improvement, not a scientific result. The v1 pilot artifact and numerical table above remain unchanged and must not be described as a v2 run. The v2 workflow associated with the latest code revision is still pending at the time of this manuscript update. No v2 artifact audit pass, confirmatory coverage estimate, or multi-frame result is claimed. The next release gate is a successful exact-head test workflow followed by generation and independent audit of a separately labelled v2 pilot artifact.


The v2 auditor also retains a versioned compatibility path for historical v1 artifacts. A regression test converts a generated test artifact to the v1 schema (removing the v2-only columns and refreshing its checksum) and verifies that the legacy audit still succeeds. This preserves the distinction between historical v1 evidence and future v2 attempt-audited evidence; it does not upgrade any v1 artifact to v2.


## 12. Separating conditional coverage from unconditional operational performance

The v2 summary schema adds two quantities alongside the historical conditional coverage estimate. Let (S=1) denote successful interval construction and (C=1) denote that the constructed interval covers the finite-frame target. The existing coverage field estimates
[
P(C=1mid S=1),
]
which describes coverage among successfully constructed intervals. The new interval-success rate estimates (P(S=1)), while `coverage_unconditional_success_and_cover` estimates
[
P(S=1,C=1),
]
counting an interval-construction failure as not delivering a covering interval. The corresponding Monte Carlo standard error is calculated over all prespecified outer repetitions, not only successful intervals.

These are different estimands and should not be conflated. The unconditional success-and-cover probability is an operational performance measure, not a replacement for conditional coverage or a claim that failures are ordinary noncoverage events in every inferential framework. Both quantities and their denominators must be reported together, with failure counts and reasons.

The independent auditor reconstructs these summary quantities from the replicate-level records and checks the stored values. This is a v2 schema addition; historical v1 summaries remain interpreted using their original conditional-coverage fields. No numerical v2 result is claimed until the exact-head workflow passes and a v2 artifact is generated and independently audited.


## 13. Interpretation and decision limits

The pilot is useful for implementation diagnosis, not for selecting a universally preferred interval. Its 18 method-by-frame-by-sample-size cells are based on only three realized finite frames and 100 outer samples per cell. The nominal 95% coverage estimates consequently have substantial Monte Carlo uncertainty, while the one-frame-per-regime design leaves between-frame heterogeneity unidentified. The observed 0.89–0.90 coverage in the normal frame at a 40% sampling fraction is a diagnostic signal for replication, not evidence of systematic undercoverage.

Interval width is not a stand-alone measure of inferential quality. A narrower interval is not preferable if its coverage or construction reliability is inadequate. Conversely, the unconditional success-and-cover probability introduced in v2 is an operational summary and must be interpreted jointly with conditional coverage and interval success. No method ranking is justified by the pilot.

The ordinary percentile comparator does not explicitly reproduce the without-replacement design. Its inclusion is intended to illustrate the consequences of using a conventional with-replacement resampling benchmark, not to define a fair competitor under every inferential objective.

## 14. Reproducibility, software, and artifact status

The benchmark is implemented in Python using NumPy. Deterministic seed derivation is used for frame generation, outer SRSWOR samples, and method-specific bootstrap resampling. The artifact design stores finite-frame values, replicate-level method outcomes, summary statistics, and a manifest containing schema and file checksums. The independent auditor reconstructs the finite frames and outer sample sequence, recomputes interval endpoints, and verifies summary quantities rather than relying only on stored checksums.

The v1 pilot artifact and its audit remain the only reported numerical benchmark in this note. The v2 code path adds serialized inner-bootstrap attempt accounting and independent replay checks, plus the distinction between conditional coverage and unconditional successful-and-covering probability. A v2 result must not be treated as released until an exact-commit continuous-integration run succeeds, a v2 artifact is generated, and the independent auditor passes on that artifact. The existence of code and tests alone is not evidence that a v2 experiment has run successfully.

The confirmatory multi-frame design is documented separately and remains prospective. Its frame count and repetition count must be selected using an explicit precision analysis before execution; proposed planning values are not experimental results.

## 15. Relation to dark-siren cosmology

The finite-population arithmetic mean is a deliberately simple estimand. Dark-siren inference instead combines gravitational-wave distance likelihoods, galaxy redshift information, catalogue selection, host-population assumptions, and cosmological parameters. The resulting posterior summaries are nonlinear functionals of data and model components. Good performance for a finite-population mean does not establish calibrated uncertainty for a dark-siren posterior.

The linear benchmark can contribute a narrow software-engineering and statistical-validation lesson: resampling procedures should reproduce the sampling design, expose failed internal attempts, preserve seed lineage, and be independently auditable. Any transfer to dark-siren inference requires separate simulation experiments that target posterior bias, credible-interval coverage, simulation-based calibration, selection-function misspecification, and uncertainty propagation. Those experiments must use the nonlinear inference pipeline itself and must be independently reconstructed.

## 16. Conclusion

This exploratory benchmark provides a reproducible pilot comparison of three interval constructions for the arithmetic mean of fixed finite frames under SRSWOR. The initial audit reconstructed all 1,800 outer method attempts and 18 summary cells, while also identifying a gap in serialized inner-bootstrap failure accounting. The v2 implementation addresses that provenance gap in code and adds explicit reporting of both conditional coverage and unconditional successful-and-covering probability. These are improvements to auditability and interpretation, not new scientific findings.

The pilot is too small and too frame-limited to support a general performance ranking. The appropriate next step is a precision-designed multi-frame experiment with paired method comparisons, frame-clustered uncertainty, complete failure accounting, and an independently audited release. No conclusion about nonlinear dark-siren posterior calibration follows from this linear study.

## References

Wang, S., Peng, L., & Kim, J. K. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *Journal of the Royal Statistical Society: Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506
