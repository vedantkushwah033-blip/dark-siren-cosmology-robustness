# Calibrating finite-validation uncertainty in inverse-probability host weighting: a reproducible dark-siren simulation benchmark

**Manuscript status:** research draft / pre-results. This file is not a submission-ready paper. All result-bearing statements and numerical fields are deliberately marked as pending until the pilot and confirmatory experiments have been executed and independently audited.

## Abstract (to be completed after confirmatory analysis)

Dark-siren inference combines gravitational-wave distance information with galaxy-catalogue information to constrain cosmological parameters without an electromagnetic counterpart. When catalogue inclusion probabilities are estimated from a finite validation sample, uncertainty in the inclusion model propagates through inverse-probability weights, host-density reconstruction, and the resulting cosmological posterior. This study develops a controlled synthetic benchmark for evaluating that propagation in a simplified low-redshift Euclidean model. For each generated finite parent population, we define a finite-frame reference posterior-mean functional by fitting the same prespecified pipeline to all parent records. We then repeatedly draw nested simple-random samples without replacement (SRSWOR) and compare three uncertainty procedures: an ordinary paired-record bootstrap, an exploratory multinomial pseudo-population/SRSWOR percentile procedure, and a scalar finite-population-correction diagnostic. The primary evaluation concerns calibration of intervals for the finite-frame posterior mean, not coverage of the generating Hubble constant and not the conditional credible interval of an individual posterior. Results will be reported only after the pilot establishes feasible replication and a precision-based confirmatory design is frozen. **Confirmatory results: pending.** The benchmark is a methodological simulation study; it is not an analysis of LVK/GWTC observations or a new observational measurement of the Hubble constant.

**Keywords:** dark sirens; finite-population sampling; catalogue incompleteness; inverse-probability weighting; bootstrap calibration; SRSWOR; uncertainty propagation; simulation study.

---

## 1. Introduction

Gravitational-wave standard sirens provide luminosity-distance information that can be combined with redshift information to study cosmology. In dark-siren analyses, where no unique electromagnetic counterpart identifies the host, galaxy catalogues supply a probabilistic representation of possible host galaxies. Catalogue incompleteness and uncertain host-population descriptions can therefore affect the inference.

A practical pipeline may estimate galaxy-catalogue inclusion probabilities from validation data and use inverse-probability weights to correct the observed host population. The resulting uncertainty is not limited to the final cosmological likelihood: it is propagated through the estimated inclusion model, the weights, the reconstructed host-redshift density, and the posterior. When validation data are sampled without replacement from a finite parent population, the sampling design also matters.

This paper isolates one narrow question: **how well do common resampling procedures estimate repeated-sampling uncertainty in a nonlinear posterior-mean functional when validation records are sampled without replacement from a finite synthetic parent frame?**

The aim is not to propose a new general dark-siren likelihood. It is to build an independently auditable benchmark in which the target, sampling design, resampling algorithms, failure handling, and calibration summaries are explicit.

### 1.1 Scope and contribution

The planned contribution is an empirical comparison of uncertainty procedures passed through a specific nonlinear synthetic pipeline. The benchmark distinguishes:
1. the finite-frame reference posterior mean;
2. the outer repeated-SRSWOR distribution of posterior means;
3. inner-bootstrap uncertainty estimates conditional on an outer validation sample;
4. conditional credible intervals for the Hubble constant within an individual posterior.

These are different statistical objects and are not interchangeable.

The contribution is deliberately narrower than the established literature on catalogue incompleteness, host weighting, selection effects, or finite-population bootstrap theory. The pseudo-population comparator is exploratory and is not claimed to implement the studentized finite-population-mean procedure of Wang, Peng and Kim (2022).

## 2. Statistical target and estimands

Let scenario (s) contain a generated finite parent frame of (N) galaxies. Each parent record has redshift (z_i) and a binary catalogue-retention indicator (R_i). A fixed observed catalogue and a fixed set of detected gravitational-wave distance measurements are generated once per scenario.

Let (mathcal{P}(cdot)) denote the prespecified pipeline that:
1. fits a logistic model for catalogue inclusion using a validation sample;
2. estimates inverse-probability weights;
3. reconstructs a weighted host-redshift density using the specified reflected weighted KDE;
4. evaluates the simplified selection-aware distance likelihood;
5. normalizes the posterior over the declared (H_0) grid; and
6. returns its posterior mean.

### 2.1 Finite-frame reference

The finite-frame reference is
[
T_s = mathcal{P}({(z_i,R_i):i=1,ldots,N}).
]
It is calculated from all actual records in the generated parent frame. It is not the generating value (H_{0,mathrm{true}}), a superpopulation target, or a multinomial pseudo-population reference.

### 2.2 Outer repeated-sampling distribution

For validation size (n), let (V_{s,n,r}) be outer draw (r), sampled by SRSWOR from the (N)-record parent frame. The corresponding posterior mean is
[
widehat T_{s,n,r}=mathcal{P}(V_{s,n,r}).
]
The scenario-specific repeated-sampling standard deviation is
[
sigma_{s,n}=operatorname{SD}_r(widehat T_{s,n,r}),
]
estimated from unique successful outer draws. Outer failures are retained and reported; they are not silently removed from the attempted denominator.

### 2.3 Inner uncertainty procedures

For each successful outer draw, method (m) produces an estimated standard error (widehatsigma_{s,n,r,m}) and a percentile or diagnostic interval for the posterior mean. The principal standard-error diagnostic is the scenario-level ratio
[
Q_{s,n,m}=
rac{operatorname{mean}_r(widehatsigma_{s,n,r,m})}
{sigma_{s,n}}.
]
A ratio near one is a calibration diagnostic, not by itself proof of valid interval coverage.

The primary coverage outcome is whether the method-specific interval contains (T_s). It is not whether the interval contains (H_{0,mathrm{true}}).

### 2.4 Hierarchical unit of replication

The generated parent scenario is the top-level independent simulation unit. Outer validation draws are nested within scenarios; inner bootstrap replicates are nested within outer draws. They must not be treated as independent parent scenarios. Scenario-balanced summaries give each scenario equal weight. The current analysis code also reports nominal 95% Student-t intervals across independent scenario-level summaries for the calibration ratio, target coverage, interval width, and failure fractions. These intervals use parent scenarios—not outer draws or inner replicates—as the independent units. Their small-sample and distributional assumptions will be assessed in the pilot; they are descriptive uncertainty summaries, not exact coverage intervals.

## 3. Synthetic data-generating process

The model is a low-redshift Euclidean toy model. The parent redshift distribution is generated from a smooth positive intensity with two broad overdensity features over the declared redshift interval. Catalogue retention is generated using a known smooth redshift-dependent probability. Gravitational-wave hosts are drawn from the complete parent frame with a specified redshift-dependent merger weight. Measured distances are generated from a Gaussian perturbation to the Euclidean approximation (d simeq cz/H_{0,mathrm{true}}), and events below the declared detection threshold are rejected.

The benchmark fixes the observed catalogue and detected distance sample within each generated parent scenario while repeatedly resampling validation records. This isolates uncertainty propagation from finite validation sampling conditional on those scenario-level observations. It does not capture all sources of uncertainty in a real gravitational-wave analysis.

All constants, grids, random seeds, software versions, and input hashes must be recorded in the machine-readable experiment manifest.

## 4. Compared procedures

### 4.1 Ordinary paired-record bootstrap

Within each outer validation sample, resample the observed ((z_i,R_i)) pairs with replacement. On each replicate, refit the inclusion model, recompute inverse weights and the host density, and rerun the toy posterior. The standard deviation and percentile interval of successful replicate posterior means estimate resampling variability conditional on that outer sample.

### 4.2 Multinomial pseudo-population followed by SRSWOR

Generate a pseudo-population of size (N) by multinomial resampling from the observed validation records, then draw a size-(n) sample without replacement and rerun the complete nonlinear pipeline. This is an exploratory comparator. Its percentile interval is not the studentized finite-population-mean procedure of Wang, Peng and Kim (2022), and theory for finite-population means is not assumed to validate this nonlinear posterior functional.

### 4.3 Scalar finite-population correction diagnostic

Multiply the ordinary bootstrap standard error by
[
sqrt{1-n/N}
]
and form a normal interval centered on the outer posterior mean. This is included as a transparent diagnostic comparator, not as a theoretically justified interval for the nonlinear estimator.

### 4.4 Interpretation distinction

The intervals compared in this paper concern repeated-sampling uncertainty in posterior means. They are not the conditional (H_0) credible intervals generated inside each individual fit. Those posterior credible intervals are stored as diagnostics but are not substituted for the bootstrap intervals.

## 5. Failure handling and auditability

Every attempted outer draw and every attempted inner replicate must remain represented in the audit trail. The analysis reports outer failures separately from inner-bootstrap failures. Failed attempts retain their indices, exception class, stable failure category, and diagnostic message where available. Attempt-aligned numerical outputs for failed fits must be NaN.

The independent audit is designed to reconstruct:
- parent redshifts and detected distance samples from declared seeds;
- catalogue thinning and canonical input fingerprints;
- the finite-frame reference posterior;
- outer nested SRSWOR index vectors and outer posterior means;
- inner bootstrap indices, pseudo-population mappings, posterior calculations, and failures;
- FPC arithmetic and all scenario-level and scenario-balanced summaries.

The audit must be run on the complete pilot and confirmatory artifacts before scientific interpretation. A passing smoke test verifies execution and file contracts only; it is not evidence of calibration.

## 6. Pilot and confirmatory design

The pilot is used to estimate failure rates, runtime, scenario-to-scenario variability, outer-SD stability, bootstrap Monte Carlo variability, boundary-limited fits, effective sample sizes, and extreme inverse weights. Pilot outcomes may inform a versioned protocol amendment but must not be presented as confirmatory evidence.

The confirmatory number of independent parent scenarios is selected from prespecified precision targets using pilot-estimated between-scenario variability. Outer draws per scenario are selected to estimate scenario-specific repeated-sampling SDs with adequate stability. Inner replicates are selected to make Monte Carlo error in bootstrap SEs and interval endpoints acceptably small.

Before confirmatory execution, freeze:
- code commit and software environment;
- scenario count, outer draws, validation sizes, and inner replicates;
- random-seed schedule;
- primary estimands and contrasts;
- interval confidence level and target precision;
- failure denominators and handling;
- stopping rule and compute budget.

No method is selected as the winner using the same confirmatory data on which its reported performance is evaluated.

## 7. Results (pending)

No confirmatory numerical results are available in this draft. Populate the following only from an independently audited confirmatory artifact.

### 7.1 Execution and failure accounting

| Quantity | Confirmatory value |
|---|---:|
| Independent parent scenarios attempted | Pending |
| Parent scenarios successfully reconstructed | Pending |
| Outer SRSWOR draws attempted / successful | Pending |
| Inner bootstrap attempts by method | Pending |
| Inner bootstrap failures by method | Pending |
| Boundary-limited inclusion fits | Pending |
| Low-ESS / extreme-weight frequency | Pending |
| Runtime and compute environment | Pending |

### 7.2 Calibration summaries

| Validation size | Method | Scenario-balanced SE/outer-SD ratio | Finite-frame target coverage | Mean interval width | Outer failure fraction | Inner failure fraction |
|---:|---|---:|---:|---:|---:|---:|
| Pending | Paired bootstrap | Pending | Pending | Pending | Pending | Pending |
| Pending | Pseudo-population SRSWOR | Pending | Pending | Pending | Pending | Pending |
| Pending | Scalar FPC diagnostic | Pending | Pending | Pending | Pending | Pending |

All estimates must include scenario-level uncertainty intervals and explicit denominators. Outer draws and bootstrap replicates must not be treated as independent top-level units.

## 8. Discussion (pending)

The discussion will distinguish:
- calibration of uncertainty for the finite-frame posterior-mean target;
- stability and failure behaviour of the inclusion model and inverse weights;
- the extent to which the scalar FPC diagnostic tracks the empirical outer-SRSWOR spread;
- differences between the paired bootstrap and exploratory pseudo-population comparator;
- sensitivity to validation size and scenario heterogeneity.

No conclusion that a method is superior will be made from smoke or pilot data. If the confirmatory study is null, adverse, or inconclusive, that result will be reported rather than reframed as success.

## 9. Limitations

1. The data-generating model is synthetic, low-redshift, and Euclidean.
2. The observed catalogue and detected GW distances are fixed within each parent scenario during outer validation resampling.
3. The benchmark does not use LVK/GWTC observations and does not provide a new observational (H_0) measurement.
4. The benchmark does not resolve the Hubble tension and does not establish a general dark-siren likelihood.
5. The pseudo-population method is an exploratory nonlinear adaptation, not the published studentized finite-population-mean procedure.
6. The scalar FPC comparator is diagnostic only.
7. Results from this benchmark cannot be generalized to real survey selection or astrophysical populations without additional observationally grounded validation.

## 10. Reproducibility and data availability

The code, protocol, seeds, manifests, raw attempt-level artifacts, independent audit outputs, and environment metadata should accompany any public release. The final paper should identify the exact commit used for the confirmatory run and provide checksums for all released artifacts. No result table should be manually transcribed without an automated check against the audited output.

## References

Wang, Z., Peng, L., & Kim, J. K. (2022). Bootstrap Inference for the Finite Population Mean under Complex Sampling Designs. *Journal of the Royal Statistical Society: Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506

**Bibliography status:** incomplete. Add only verified references relevant to dark-siren inference, galaxy-catalogue incompleteness, selection effects, and finite-population resampling before submission.
