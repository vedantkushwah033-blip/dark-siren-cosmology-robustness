# Confirmatory protocol: finite-validation IPW host-density recovery

**Protocol date:** 2026-11-08  
**Branch:** `research/independent-posterior-audit-claim-fix`  
**Status:** preregistration/design specification; no confirmatory repository-native density experiment has yet been executed  
**Scope:** synthetic low-redshift Euclidean host-density recovery only; no H0 posterior or observational claim

## 1. Purpose and decision boundary

The probability-stage benchmark now has an independently audited 450-attempt pilot. A separate 20-replicate local density pilot suggests that estimated inclusion weights can improve or worsen parent-host density recovery depending on curvature and validation size. That local pilot was an isolated implementation mirror, not a repository-native benchmark, and its results are exploratory.

The next experiment must answer a narrower question before any H0 likelihood is changed:

> Under a declared synthetic parent population and catalogue-inclusion mechanism, when does fitting an inclusion model on a finite, independently sampled validation set improve or degrade reconstruction of the merger-weighted parent-host redshift density relative to using the observed catalogue without correction?

This protocol does **not** test posterior calibration, H0 coverage, or the validity of a real survey selection function. It is a prerequisite for those later tasks.

## 2. Estimand and data-generating process

Let the parent redshift density be proportional to the declared intensity

\[
\lambda(z)=1+
2.2\exp\!\left[-\frac12\left(\frac{z-0.105}{0.025}\right)^2\right]+
1.6\exp\!\left[-\frac12\left(\frac{z-0.275}{0.040}\right)^2\right],
\quad z\in[0.01,0.40].
\]

The target for host-density recovery includes the declared merger weight

\[
m(z)=1+1.5x,\qquad x=(z-0.01)/0.39,
\]

and is the normalized density

\[
g_0(z)=\frac{\lambda(z)m(z)}
{\int_{0.01}^{0.40}\lambda(u)m(u)\,du}.
\]

Catalogue inclusion is generated separately:

\[
R_i\mid z_i\sim\operatorname{Bernoulli}(\pi(z_i)),
\qquad
\pi(z)=\operatorname{expit}\{\alpha-x+\gamma(x^2-c)\},
\]

where

\[
c=E_{\lambda}[x^2]
=\frac{\int\lambda(z)x(z)^2\,dz}{\int\lambda(z)\,dz}.
\]

For each curvature and retention target, solve \(\alpha\) by numerical quadrature so that

\[
\frac{\int\lambda(z)\pi(z)\,dz}{\int\lambda(z)\,dz}=r.
\]

The centering constant must be computed under the declared nonuniform parent law. Do not replace it with 1/3.

The parent frame has \(N=4,000\) units. The validation set is selected by simple random sampling without replacement from the complete parent frame, independently of the retention labels. Validation units are excluded from the analysis catalogue. The retained analysis catalogue is

\[
C=\{i:i\notin V, R_i=1\}.
\]

The true inclusion probabilities are available to the simulator and auditor only. They must never be supplied to the fitted logistic estimator.

## 3. Confirmatory design: staged rather than an indiscriminate full factorial

### Primary confirmatory block

Use:

- Parent size: 4,000.
- Curvature: \(\gamma\in\{-1.2,0,+1.2\}\).
- Target parent retention: \(r=0.25\).
- Validation sizes: \(n_V\in\{100,300,1000\}\).
- Fitted inclusion models: polynomial degree 1 and 2 in standardized redshift.
- Replicates: 100 per curvature condition.
- Random seeds: deterministic and recorded in the run manifest.

This yields 1,800 model-fit attempts:

\[
3\times 3\times 2\times100=1,800.
\]

The three validation sizes should be **nested within each replicate** using one random permutation of parent IDs: the first 100 validation IDs are contained in the first 300, which are contained in the first 1,000. This preserves pairing across validation sizes. Use one parent redshift frame and one vector of retention uniforms per replicate, shared across curvature scenarios where common random numbers are intended. Record the exact seed and ID mapping.

### Secondary stress block

After the primary block passes its audit, add:

- Curvature: \(\gamma\in\{-2,+2\}\).
- Retention targets: \(r\in\{0.10,0.50\}\).
- Validation sizes: 100, 300, 1,000.
- Degrees: 1 and 2.
- At least 50 replicates per cell.

Treat this as a stress/heterogeneity analysis, not a replacement for the primary block. If computational or artifact limits require reducing it, report the exact omitted cells and rationale.

### Why not immediately run every possible combination?

A large factorial campaign is only useful if raw artifacts can be independently reconstructed. The primary block is designed to test the candidate asymmetry at positive curvature, the correctly specified linear case at zero curvature, and the opposite curvature direction. The secondary block tests whether the result survives stronger curvature and different retention levels.

## 4. Frozen density-estimation protocol

Use the existing `estimate_host_density` implementation and its explicit `scott_multiplier` API.

Primary bandwidth:

\[
h=0.15\times h_{\mathrm{Scott,weighted}},
\]

where the weighted Scott factor is computed from the original, unreflected analysis sample and its effective sample size. This value is a **candidate frozen protocol**, motivated by a separate training/held-out pilot; it is not claimed to be universally optimal. The confirmatory seed set must not be used to tune it.

Run a predeclared sensitivity analysis at Scott multipliers 0.40 and 1.00. Report all arms at all three multipliers in the sensitivity output, but designate 0.15 as the primary analysis. Do not select the multiplier that gives the most favorable IPW contrast or H0 coverage.

The implementation must preserve the legacy meaning of `bandwidth_factor`: an absolute SciPy covariance factor. The confirmatory runner must use `scott_multiplier`, not the legacy argument.

All non-analytic KDE arms must use the same grid, reflection rule, merger weight and bandwidth multiplier. Normalize the final merger-weighted density by numerical integration on the declared bounded support.

## 5. Comparison arms

For each shared parent frame, retention realization and validation split, construct:

1. **Catalogue-only:** KDE from retained non-validation catalogue redshifts, with unit weights.
2. **Known-inclusion IPW:** KDE from the same analysis catalogue with weights \(1/\pi(z_i)\). This is an oracle-probability diagnostic, not an operational estimator.
3. **Estimated linear IPW:** KDE from the same analysis catalogue with weights \(1/\widehat\pi_1(z_i)\), where the degree-1 logistic fit uses validation labels only.
4. **Estimated quadratic IPW:** corresponding degree-2 fit.
5. **Parent-sample KDE diagnostic:** KDE using all parent-frame redshifts. This is an oracle sample diagnostic, not the analytic target.
6. **Analytic target:** \(g_0(z)\), used as the exact scoring reference, not as an estimator arm.

Do not conflate the parent-sample KDE with the analytic target. The former has finite-sample KDE error; the latter is the declared population target.

The validation sample must not be added back to the analysis catalogue. The oracle arm may use the true probabilities only inside its own diagnostic calculation.

## 6. Primary and secondary outcomes

### Primary outcome

Integrated absolute error against the analytic merger-weighted host density:

\[
\operatorname{IAE}=\int_{0.01}^{0.40}|\widehat g(z)-g_0(z)|\,dz.
\]

Primary contrasts are paired replicate-level differences:

\[
\Delta_{d,\gamma,n}
=\operatorname{IAE}_{\mathrm{estimated},d,\gamma,n}
-\operatorname{IAE}_{\mathrm{catalogue-only},\gamma,n}.
\]

Negative values favor estimated IPW; positive values mean estimated IPW is worse.

### Secondary outcomes

- Integrated squared error (ISE).
- Boundary-region IAE, with boundary regions declared before execution.
- Density normalization error before final normalization.
- Maximum absolute density error.
- Density error by redshift subregion.
- Probability RMSE, Brier score and log loss from the validation model.
- Equal-frequency reliability-bin ECE and maximum bin gap.
- Kish ESS and entropy ESS as fractions of analysis-catalogue size.
- Maximum normalized weight.
- Top 1% and top 5% normalized weight mass.
- Weight quantiles and fraction of fitted probabilities below declared numerical thresholds.
- Fit failure and non-evaluable rates.

ESS is a weight-concentration diagnostic, not a density-accuracy measure. Probability calibration metrics are not posterior calibration metrics.

## 7. Statistical analysis plan

For every scenario and fitted degree, report the replicate-level paired differences, their mean, Monte Carlo standard error, and a 95% uncertainty interval. Use a paired bootstrap over replicate IDs as the primary interval procedure, with a paired t-based interval as a sensitivity check. Resample complete paired replicate records, never individual galaxies independently.

Predeclare the main confirmatory contrast family as estimated degree-1 and degree-2 IPW versus catalogue-only at \(\gamma=+1.2\), for all three validation sizes. Control family-wise error for this six-contrast family if formal hypothesis tests are reported. The zero- and negative-curvature blocks are supporting contrasts and must be labelled accordingly.

Do not interpret an interval excluding zero as evidence of universal estimator superiority. Report effect sizes, Monte Carlo precision, and the exact scenario scope.

Report both:
- performance conditional on successful fitting; and
- attempt-level fit status/failure rates.

Never silently discard failed or non-evaluable fits. Density outcomes are undefined for failed fits and must remain explicitly missing with their status retained.

## 8. Raw artifact and manifest requirements

The experiment output must preserve enough information for independent reconstruction. At minimum, save:

- Immutable configuration and configuration SHA-256.
- Git commit SHA and branch/ref.
- Python, NumPy, SciPy and pandas versions.
- Parent-frame seed and replicate ID.
- Parent object IDs and redshifts.
- Retention uniforms or sufficient exact data to reconstruct retained labels.
- True inclusion probabilities and DGP parameters, clearly marked simulation-truth/auditor-only.
- Validation IDs and split roles.
- Retained labels.
- Fitted coefficients, standardization center/scale, fit status and optimizer diagnostics.
- Fitted probabilities for analysis-catalogue units.
- Raw inverse weights or sufficient probabilities to reconstruct them.
- KDE grid, bandwidth multiplier, density vectors and their checksums.
- Analytic target density and target checksum.
- IAE, ISE and all secondary metrics.
- Per-attempt status: SUCCESS, NON_EVALUABLE, FAILED or NOT_RUN.
- Runtime, artifact sizes and SHA-256 hashes.

Avoid needless duplication of the same parent frame across polynomial degrees. If raw records are shared between model degrees, make the experiment keys and separate fitted-probability columns explicit so the auditor can reconstruct each arm without ambiguity.

The run manifest must state that the output is a density-stage synthetic benchmark and contains no H0 inference.

## 9. Independent auditor specification

The auditor must not import the benchmark runner, estimator or scoring functions. It should independently:

1. Validate required columns, types, unique experiment keys and complete parent ID ranges.
2. Reconstruct disjoint validation, analysis-catalogue and non-retained partitions.
3. Reconstruct the parent-law curvature center and calibrated intercept from the declared DGP.
4. Recompute the truth probabilities and verify stored retention labels against the stored randomization record.
5. Reconstruct fitted probabilities from stored coefficients and standardization metadata.
6. Independently evaluate logistic score and observed information; verify covariance only as a model-based covariance estimate.
7. Reconstruct catalogue-only and IPW weights.
8. Independently reconstruct weighted KDEs using the declared reflection and bandwidth semantics.
9. Verify density normalization and recompute the analytic target.
10. Recompute IAE, ISE, boundary metrics, weight diagnostics and paired replicate contrasts.
11. Verify every summary against raw records within predeclared numerical tolerances.
12. Verify configuration and artifact hashes.
13. Verify failed/non-evaluable attempts have not disappeared from the ledger.

Include adversarial tests that alter a probability, ID, split role, coefficient, bandwidth, density vector, summary metric or checksum and confirm that the auditor rejects the artifact.

The audit report must distinguish numerical reconstruction from scientific model validation.

## 10. Predeclared interpretation rules

- If known-inclusion IPW does not improve density recovery over catalogue-only in a cell, first investigate the DGP, target density, finite-parent noise and bandwidth; do not blame estimated probabilities.
- If known-inclusion IPW improves but estimated IPW does not, quantify the gap between known and estimated probabilities and investigate finite-validation uncertainty/misspecification.
- If estimated IPW performs worse without extreme weight concentration, do not force an ESS-based explanation. Probability-model bias or variance may dominate.
- If the result changes materially across the frozen bandwidth sensitivity grid, report the interaction and do not make a bandwidth-invariant claim.
- If the parent-sample KDE itself has substantial IAE, retain the analytic density as the scoring target and report the finite-sample oracle-KDE gap.
- Do not use H0 coverage to select a density bandwidth.
- Do not claim that pre-posterior diagnostics predict posterior coverage failure from this density-only experiment.

## 11. Exit criteria for moving to H0 inference

The density stage is ready to inform the next stage only when:

1. The repository-native runner completes the primary block.
2. The independent auditor passes on the full raw artifact.
3. Adversarial auditor tests pass.
4. The manifest and hashes verify.
5. The primary results are stable to reasonable numerical grid refinement.
6. Bandwidth sensitivity is reported and not hidden.
7. Fit failures are fully accounted for.
8. Known-inclusion IPW, catalogue-only and estimated-IPW arms are all compared against the same analytic target.

Passing these criteria validates the density-stage computational experiment only. It does not establish a calibrated H0 posterior.

## 12. Explicit non-goals

This protocol does not authorize:

- An observational H0 measurement.
- LVK/GWTC public event analysis.
- Claims about real survey completeness.
- Universal IPW superiority.
- A new estimator claim.
- A claim that the inverse observed-information covariance is a Bayesian posterior.
- H0 coverage conclusions.
- A second missing-host mixture correction on top of IPW.

The subsequent H0 stage must separately derive the selection-conditioned likelihood and a justified uncertainty model for inclusion coefficients. Only after that should posterior calibration and repeated-sampling coverage be studied.

## 13. Audit follow-up: CSV round-trip tolerance and failure-aware contrasts (2026-11-08)

A historical CI failure in the independent retention-stream audit was traced to an overly strict comparison, not to a demonstrated RNG-stream reuse:

- At failed commit `758ba999443133a9703b2bc1d8206592612202a9`, the auditor compared the stored CSV `retention_uniform` values with regenerated NumPy doubles using `np.array_equal`.
- The same auditor already used an absolute tolerance for the stored parent-redshift values, but required bitwise equality for retention uniforms.
- The runner serializes these floating-point values to CSV and the auditor reads them back. Decimal serialization/parsing need not preserve the exact binary floating-point bit pattern, even when the stored value is numerically equivalent to the generated value.
- The current auditor uses `np.allclose(stored, expected, rtol=0.0, atol=5e-15)` for retention uniforms, matching the explicit round-trip tolerance policy used for parent redshifts. The adversarial test still changes a uniform by a material amount and must be rejected.

This is a numerical serialization-tolerance correction. It does not establish that the old smoke artifact was scientifically valid, nor does it excuse a mismatch larger than the declared tolerance. The exact current-head CI run must still pass before this fix is marked validated.

The paired density-contrast summary has also been extended to report attempt-level lifecycle outcomes alongside conditional paired IAE:

- successful, failed and non-evaluable attempts for each arm;
- replicates with both attempts present;
- both-success pairs;
- cases where estimated IPW is not successful but catalogue-only succeeds, and the reverse;
- cases where neither arm succeeds; and
- the evaluable-pair fraction relative to the planned replicate count.

The independent auditor reconstructs these counts from the raw attempt ledger and compares them with the stored summary. The conditional mean contrast and its intervals remain explicitly conditional on both arms being successful; failure counts are not folded into an arbitrary numerical IAE penalty.

This extension improves transparency but does not yet define a single composite utility for comparing accuracy and failure probability. Both dimensions must be reported separately.
