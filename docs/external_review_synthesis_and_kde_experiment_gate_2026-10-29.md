# External review synthesis and corrected next-experiment gate — 29 October 2026

## Status and provenance

This document synthesizes the external review text supplied to the project owner on 29 October 2026. The supplied Perplexity report explicitly says it could not inspect the live repository and that no repository changes, tests, data downloads, or workflow dispatches were performed. The supplied statistical review provides recommendations but no source files, executable code, raw replicate records, or independently reproduced numerical calculations.

These reports are review inputs, not independently verified repository findings. This document does not claim that the cited papers, headline estimates, or alleged defects have been re-verified here.

## Executive decision

- Retain the narrow synthetic calibration/robustness framing. Do not claim a new astrophysical method, observational H0 measurement, general LVK validation, or CHIMERA v2.2 validation.
- The literature review is useful but not exhaustive. Verify bibliographic records and every claim against primary sources before manuscript use.
- Some statistical-review statements are stronger than the evidence permits. Overlapping seed ranges alone do not establish a confirmed defect.
- Do not launch the proposed KDE-bandwidth campaign unchanged. Its estimand and the mechanism by which KDE estimates inclusion probabilities are insufficiently specified.
- CHIMERA v2.2 remains NO-GO until pinned inputs, environment, end-to-end execution, and output validation are established.

## 1. Literature report: accepted conclusions and remaining checks

The following conclusions are appropriately conservative:

1. Existing work already covers galaxy-catalogue incompleteness, hidden-host mixtures, host weighting, selection normalization, clustering-aware completion, GP reconstruction, and robustness to host-population assumptions.
2. The defensible project scope is a controlled synthetic calibration benchmark, not a new general incompleteness-correction method.
3. Oracle inclusion probabilities are a special, favorable assumption and must not be presented as equivalent to estimated probabilities in a real survey.
4. Documentation or API discovery alone is not a CHIMERA benchmark.

Before citing the submitted matrix, verify each entry against its primary publication and record DOI/arXiv identifier, version, publication status, and the exact source supporting the described result. Treat unreviewed leads (including Finke et al., Dalang & Baker, Cross-Parkin et al., and Perna et al.) as search leads, not confirmed evidence. Do not carry opaque bracket tokens such as [3266] into the manuscript unless mapped to retrievable source identifiers.

## 2. Statistical red-team: corrections to overstrong statements

### 2.1 Seed overlap is not automatically a defect

Overlapping seed ranges across conditions do not, by themselves, prove dependence or invalid uncertainty estimates. The relevant issue is whether simulated records are dependent and whether the estimand and resampling scheme account for that dependence.

- If the same simulated dataset is deliberately used across bandwidth arms, this is a paired design and may improve precision for within-dataset contrasts.
- If arms are intended to be independent but reuse random streams in a way that induces dependence, uncertainty calculations must account for it.
- Numeric seed overlap alone does not reveal whether random streams or generated records overlap.
- For independent conditions, use separate seed namespaces and persist replicate IDs. For paired comparisons, deliberately share the exogenous dataset realization across arms, label the pairing, and bootstrap replicate IDs as pairs.

Therefore, the submitted label “confirmed High defect” is not established without examining seed-generation code, records, and bootstrap implementation. Reclassify it as an audit question pending evidence.

### 2.2 Boundary coverage claims require raw counts and truth grid

The supplied review mentions pooled coverage near 57% and interior coverage near 73%, but provides no source artifacts or raw counts. Those figures are not verified by the supplied text alone.

Required check: report covered/total counts and Wilson intervals for every prespecified truth stratum, especially prior-boundary, near-boundary, and interior truths. Define strata before analysis. Do not use pooled coverage to conceal heterogeneity, but do not call boundary failure “catastrophic” without evidence.

### 2.3 Oracle weighting is a limitation, not automatically an upper bound

Oracle inclusion probabilities remove one source of estimation uncertainty, but do not guarantee best performance for every metric or under every misspecified model. Describe them as an oracle-information benchmark, not an unconditional upper bound. Estimated weights can introduce variance and bias; the direction and magnitude must be measured.

### 2.4 Selection normalization must match the likelihood

For a detected-event likelihood conditioned on observing N detections, a per-event selection normalization commonly contributes beta(H0)^(-N), if beta is the relevant detection probability under the modeled population. A rate-based Poisson likelihood has a different structure. Check the implementation against its exact data-generating process; do not apply a universal formula without defining the likelihood.

### 2.5 Negative controls need deterministic expected behavior

- **Uniform selection:** when epsilon(d)=1 everywhere, the selection probability is one under a properly normalized population density. The selection correction should reduce to the unconditioned likelihood, subject to the model's definitions. This is a valuable deterministic unit test.
- **Shuffled inclusion probabilities:** permuting probabilities across hosts is useful only if host covariates and true inclusion probabilities have a meaningful association under the DGP. Define the expected direction and magnitude in advance. “Catastrophic bias” is not guaranteed.
- **Random coordinates/noise:** replacing galaxy coordinates with random noise does not necessarily force infinite variance or total inferential collapse. Replace this vague control with a well-defined permutation or DGP-mismatch control and specify its measurable expected signature.

## 3. Proposed KDE bandwidth experiment: design blockers

The submitted proposal uses four arms (oracle, KDE at h_opt, 0.5 h_opt, 2 h_opt), 400 replicates per arm, and 300 events per dataset, with H0=70 and hidden-host fraction f=0.2. It is a useful candidate but is not ready to run.

### Blocker A — What exactly is being estimated by KDE?

The protocol says KDE estimates inclusion probabilities, but does not specify the mapping from the observed catalogue to each host's inclusion probability. A KDE of observed galaxy positions estimates an observed spatial density; it does not automatically estimate the probability that a particular galaxy entered the catalogue.

The protocol must define:
- the complete underlying host population available to the simulator;
- the catalogue inclusion mechanism and true pi(x) or pi(z, properties);
- exactly which data the estimator may observe;
- how the KDE estimate maps to estimated inclusion probability;
- whether a separate validation/training sample or external completeness information is available;
- how weights are stabilized or truncated, and whether that rule is part of the estimator.

Without this, “KDE-estimated inclusion probability” is not a reproducible procedure.

### Blocker B — Separate bandwidth effects from oracle-vs-estimated information

The four proposed arms compare oracle inclusion probabilities with three estimated-weight arms. The oracle-versus-KDE contrast changes more than bandwidth. Define the primary comparison among estimated-weight procedures; use the oracle arm as a diagnostic reference. Make data reuse in the standard-bandwidth KDE explicit.

### Blocker C — The null hypothesis and precision claim are not aligned

The proposed null (“±50% bandwidth perturbations do not significantly alter coverage”) is not a complete estimand or an equivalence test. A nonsignificant test is not evidence of equivalence.

Choose a primary endpoint, such as the paired difference in empirical 68% coverage between each perturbed bandwidth and the prespecified standard bandwidth. Define a scientifically meaningful equivalence or noninferiority margin before simulation if the goal is to establish similarity.

For nominal coverage p=0.68 and R=400 independent replicates, the binomial Monte Carlo standard error is sqrt(0.68*0.32/400) ≈ 0.0233 (2.33 percentage points) for one arm. A rough normal 95% margin is about ±4.6 percentage points for a single coverage estimate. This does not by itself establish power to detect a 7-percentage-point difference, especially for paired contrasts, multiple arms, and stratified coverage. Calculate power or precision for the actual planned comparison before locking R.

### Blocker D — One truth value cannot establish general calibration

A single truth H0=70 and one hidden-host fraction can estimate fixed-truth performance only at that condition. It cannot support broad calibration claims across prior boundaries, hidden-host fractions, or selection regimes. Either narrow the claims explicitly to that operating point or use a prespecified truth grid including interior and boundary/near-boundary values. Broad coverage claims require sufficient replicates per truth stratum.

### Blocker E — Computational estimate is unverified

The proposed “approximately four minutes” is a forecast, not an observed runtime. Measure a small non-confirmatory pilot on the actual environment, then estimate runtime and memory. Keep pilot seeds separate from confirmatory seeds.

## 4. Revised preregistration gate

Do not run the confirmatory simulation until the protocol answers all of the following:

1. **Question:** How does KDE bandwidth misspecification affect bias, RMSE, and credible-interval coverage of a fully specified inverse-inclusion-weighted estimator under a defined synthetic catalogue-inclusion mechanism?
2. **DGP:** Define complete hosts, inclusion mechanism, missing-host population, GW event/distance generation, measurement noise, selection, and detection conditioning.
3. **Observed information:** State exactly what the estimator can observe and what is oracle-only.
4. **Estimator:** Define KDE inputs, bandwidth rule, transformation/boundary treatment, pi-hat mapping, weight cap/truncation, and posterior construction.
5. **Arms:** Oracle probabilities as a diagnostic; standard KDE; prespecified under- and over-smoothing. Ensure the principal contrast isolates bandwidth where possible.
6. **Primary estimand:** Define a paired replicate-level contrast in 68% coverage or a prespecified calibration-loss metric. Name one primary endpoint and specify direction.
7. **Secondary estimands:** Bias, RMSE, 95% coverage, interval width, effective sample size of weights, maximum weight, non-finite fit rate, and posterior boundary mass.
8. **Pairing and seeds:** State whether the same DGP realization is used across arms. Persist seed, replicate ID, and shared-realization ID. Use paired resampling for paired arms and independent resampling only when justified.
9. **Truth grid:** Either restrict conclusions to H0=70, f=0.2 or justify a small prespecified grid. Do not claim general calibration from one point.
10. **Precision/power:** Justify R from the primary estimand's uncertainty, not just per-arm binomial SE. Define equivalence margins if equivalence is the objective.
11. **Negative controls:** Uniform-selection normalization test; inclusion-probability permutation/mismatch control with defined expected signature; boundary/grid-convergence controls.
12. **Reproducibility:** Raw replicate-level output, config, source revision, dependency versions, seeds, hashes, failure records, and independent audit script.
13. **Stopping rule:** Define when numerical failures pause a run, but preserve failed replicates and report their rate. Never silently drop difficult replicates.
14. **Pilot:** Use a small non-confirmatory pilot to test implementation and estimate runtime; freeze the protocol and use new seeds for confirmatory runs.

## 5. Recommended sequence

1. Resolve the manuscript-to-artifact map and validate existing headline numbers from immutable files.
2. Inspect actual seed code and bootstrap code before classifying seed overlap as a defect.
3. Add the uniform-selection identity as a deterministic unit test if it is not already covered.
4. Write and review the KDE-to-inclusion-probability estimator definition.
5. Freeze estimand, truth grid, pairing plan, precision target, and failure policy.
6. Run a small pilot, then conduct confirmatory simulations only after independent protocol review.
7. Keep CHIMERA v2.2 NO-GO until external inputs and execution are verified.

## Release decision

**Current decision: CONDITIONAL for internal synthetic-method development; NO-GO for manuscript submission or claims of external validation.**

The external reports strengthen the case for conservative positioning and a more carefully defined calibration question. They do not, by themselves, establish new repository defects, verify headline estimates, or authorize a large new simulation campaign.
