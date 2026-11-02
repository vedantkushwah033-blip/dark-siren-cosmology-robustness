# Audited linear finite-population bootstrap pilot: results and interpretation

**Artifact source:** GitHub Actions run 41126411371, Python 3.12 job  
**Source commit:** `74e6dacd6690875fdeb6c6d6fc0a1fcaa40fde65`  
**Artifact status:** independent audit PASS  
**Scientific status:** pilot only; not confirmatory  
**Scope:** finite-population arithmetic mean under SRSWOR. This is not a dark-siren posterior experiment.

## 1. Prespecified design

The workflow executed the pilot with seed 20261031, three fixed finite frames of size 1,000, 100 repeated SRSWOR samples at each sample size (100 and 400), and 199 internal bootstrap replicates. The three frames were generated from normal, Gamma(shape 0.7), and lognormal regimes. Each outer sample was analyzed using (i) a studentized pseudo-population bootstrap with SRSWOR, (ii) a conventional Wald interval with finite-population correction, and (iii) an ordinary percentile bootstrap.

The estimand in each frame is its exact arithmetic mean. The 18 method-by-frame-by-sample-size cells contain 1,800 replicate/method attempts.

## 2. Independent audit

The machine-readable audit report records:

- overall status: **PASS**;
- finite frames independently verified: 3;
- replicate attempts independently reconstructed: 1,800;
- summary cells verified: 18;
- frame-seed reconstruction: PASS;
- design completeness: PASS.

The audit verifies the recorded frames, seed schedule, sample indices, interval endpoints, coverage indicators, failures, and summary grid. The artifact manifest labels the experiment `PILOT_NOT_CONFIRMATORY`.

## 3. Pilot results

Coverage below is conditional on successful interval construction. All 18 cells had 100 successful intervals and zero failed interval constructions in this pilot. With only 100 repetitions per cell, Monte Carlo uncertainty is substantial.

| Frame | n/N | Studentized pseudo-population | Wald SRSWOR | Ordinary percentile |
|---|---:|---:|---:|---:|
| Normal | 0.10 | 0.96 (mean width 0.366) | 0.97 (0.372) | 0.94 (0.380) |
| Normal | 0.40 | 0.89 (0.147) | 0.90 (0.151) | 0.98 (0.190) |
| Gamma(0.7) | 0.10 | 0.93 (0.449) | 0.96 (0.434) | 0.96 (0.444) |
| Gamma(0.7) | 0.40 | 0.94 (0.172) | 0.94 (0.175) | 0.98 (0.222) |
| Lognormal | 0.10 | 0.95 (0.788) | 0.93 (0.720) | 0.95 (0.735) |
| Lognormal | 0.40 | 0.95 (0.284) | 0.95 (0.284) | 0.97 (0.354) |

Entries are empirical coverage; parenthesized values are mean interval widths in the units of the frame variable. Widths are not directly comparable across different frames because their scales differ.

## 4. Interpretation

The pilot is descriptive and does not establish that any method attains nominal 95% coverage. The most visible pattern is that the ordinary percentile comparator produces wider intervals in several 40%-sampling cells and correspondingly higher empirical coverage in this small pilot. The studentized pseudo-population and Wald methods are close in some cells but show deviations from nominal coverage, including 0.89–0.90 in the normal frame at n/N=0.40.

These are hypotheses for a confirmatory design, not stable conclusions. A 100-repetition cell has a coverage Monte Carlo standard error of roughly 0.02 near 95% coverage, so differences of a few percentage points are not persuasive. Only three frames are available, so the pilot also cannot characterize between-frame variability within each population regime.

No method is declared a winner. The next study should increase independent frames and repeated samples according to a frozen precision plan, report exact coverage numerators and denominators, and preserve the pilot separately from confirmatory data.

## 5. Strict boundary to the dark-siren work

This experiment validates only the execution and empirical behavior of interval procedures for a linear finite-population mean in these three synthetic frames. It does **not** validate inverse-probability host weights, estimated catalogue completeness, a nonlinear dark-siren posterior mean, posterior credible intervals, or an H0 measurement. Those require separate nonlinear simulation experiments with their own estimands, data-generating processes, precision targets, failure accounting, and independent audits.

## 6. Reproducibility

The raw CSVs, manifest, and independent audit report are retained in the workflow artifact `numerical-diagnostics-py3.12` from run 41126411371. The pilot can be reproduced using the workflow command and parameters recorded in the manifest and repository workflow.


## 7. Artifact fingerprints and audit boundary

The Python 3.12 workflow artifact manifest records these SHA-256 fingerprints:

- finite-frame CSV: `4dd68eb0c45565130d431910d75443a96503b7a32d4020025e08950de8c95e15`;
- replicate CSV: `792119eaa9cfcd7012066d19b13c893caea5c1873527e0380b9d866ff095097e`;
- summary CSV: `4249e3d37d2d66271741e6a6d36b5c42eaecdf4aaa60cdccbcb5ff44103ebc8e`.

The audit report independently records 3 verified frames, 1,800 reconstructed outer method attempts, 18 verified summary cells, and PASS for design completeness and frame-seed reconstruction.

**Important audit limitation:** the v1 replicate artifact does not serialize the studentized method's internal bootstrap-attempt counts and failure-code histogram. The primitive retains these in memory, and the auditor replays the final interval, but the artifact-level audit does not yet independently compare a stored row-level accounting of every inner attempt. Thus “zero failed intervals” must not be read as “zero inner bootstrap failures.” This gap is explicitly documented in [the v2 remediation protocol](inner_bootstrap_attempt_audit_gap_and_v2_protocol_2026-10-31.md). The v1 pilot remains a valid historical pilot for the outer interval results, but a v2 attempt-level audit is required before claiming complete internal failure accounting.
