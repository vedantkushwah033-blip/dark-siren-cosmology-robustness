# Corrected density smoke: artifact and independent audit

Date: 18 November 2026

Status: **ONE ENGINEERING SMOKE EXECUTED; BUILT-IN AND SECONDARY RAW-ARTIFACT CHECKS PASS. NOT A SCIENTIFIC RESULT.**

## Code provenance

- Source correction PR: [#59](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/59), still draft and unmerged at time of audit.
- Exact PR head used: `390ba806c42376b1a2d4b8d083b92b8d206c8a88`.
- Validation branch: `research/corrected-density-smoke-validation-20261118`.
- Validation workflow commit: `3b527b888435b0a446bb60a1d28c4534b5c6f5b5`.
- Dedicated Actions run: `44363575233`, workflow `corrected-density-single-replicate-audit`.

The validation branch was created from the exact PR #59 head. It contains no changes to `main`, and no PR was merged.

## Executed design

This was deliberately a small engineering smoke, not a pilot or confirmatory design:

- Parent population: 700.
- Independent validation frame: 300.
- Nested validation sizes: 80 and 160.
- Curvature: gamma = 0.5; beta = -3.9.
- Quadrature order: 3.
- Redshift grid: 101 points.
- Absolute KDE bandwidth candidates: 0.012 and 0.024.
- Seeds: 2026111801 (population) and 2026111802 (validation).

The runner used the corrected convention:
- A and B have unit sample weights.
- C, D, and E carry inclusion-only inverse-probability weights.
- The common host-rate factor is applied once after KDE.

## Workflow execution

The dedicated workflow completed the following steps successfully:
1. Install package and test dependencies.
2. Run corrected one-replicate density smoke.
3. Independently audit raw smoke artifacts.
4. Upload raw inputs, outputs, manifests, and the independent audit.

The built-in independent audit reported:
- `status: PASS`
- 11 manifest-listed files verified.
- 60 density metric rows reconstructed.
- 26 realized-Scott-width rows reconstructed.
- Validation sizes 80 and 160.
- No H0 inference in scope.

## Secondary audit performed on downloaded artifact

Downloaded artifact:
- Artifact ID: `14455017523`.
- Name: `corrected-density-single-replicate-3b527b888435b0a446bb60a1d28c4534b5c6f5b5`.
- Size: 166,324 bytes.
- SHA-256 archive digest reported by GitHub: `6492cd9583e318687c15f1d1e470cd2d558f3a37d7f2d9a5bed48fe567a0325b`.
- Expiry shown by GitHub: 2 December 2026.

I downloaded and inspected the ZIP independently:
- All 11 manifest entries matched their declared byte counts and SHA-256 hashes; zero hash failures.
- The archive contained the expected config, manifest, observable records, truth records, density inputs, metrics, weight diagnostics, Scott widths, and independent-audit JSON.
- A and B sample-weight arrays were exactly all ones for both validation sizes.
- C weights were independently reconstructed from the configured quadratic-logit inclusion function. Maximum absolute difference from stored C weights was exactly zero.
- In this realized catalogue, minimum true inclusion probability was 0.04330314, maximum was 0.98672014, and maximum oracle inverse weight was 23.093012.
- The calibrated parent-average inclusion rate was 0.5101919962291072 against target 0.5101919962291073.

The built-in auditor independently refitted the validation logistic posterior/quadrature and reconstructed the density curves and metrics. The secondary audit above separately checked the archive manifest, A/B unit weights, and C oracle weights.

## Descriptive smoke values — not inferential evidence

For orientation only, the single replicate produced the following IAE values:

| Validation n | Bandwidth | A catalogue | B parent oracle | C known pi | D MAP | E posterior-mean density |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 0.012 | 0.6430 | 0.0627 | 0.1525 | 0.2179 | 0.2557 |
| 80 | 0.024 | 0.6395 | 0.0828 | 0.0979 | 0.1838 | 0.2329 |
| 160 | 0.012 | 0.6430 | 0.0627 | 0.1525 | 0.3105 | 0.3305 |
| 160 | 0.024 | 0.6395 | 0.0828 | 0.0979 | 0.2784 | 0.3046 |

These values are a one-seed, two-bandwidth engineering smoke. They must not be interpreted as estimator rankings, evidence of benefit/harm, a selected bandwidth, or a curvature effect. In particular, the fact that D/E change with validation size while A/B/C do not is a feature of this single nested validation realization and fitted model; repeated simulations are required.

The config labels the run `ENGINEERING_SMOKE_NOT_SELECTION_NOT_CONFIRMATORY`. No H0 likelihood, posterior, coverage, PIT, or SBC result was generated.

## Exact-head PR #59 CI status at last poll

PR #59 exact head: `390ba806c42376b1a2d4b8d083b92b8d206c8a88`.

PR-triggered full test workflow: `44362767807`.
- Python 3.10 job: completed successfully (113 steps; no failed steps).
- Python 3.11 job: still running at the last poll (93 completed steps; no failures observed).
- Python 3.12 job: still running at the last poll, in the 100-replicate bandwidth-grid sensitivity step (37 completed steps; no failures observed).
- The targeted five-arm selected-H0, reflected-KDE selected-H0, and absolute-bandwidth density smoke/auditor test steps had completed successfully in all three Python jobs at the last poll.
- Therefore the full exact-head CI was **not yet green as a whole** at that poll.

## Decision and next gate

The corrected one-replicate runner and independent raw-input auditor have now executed successfully and the downloaded artifact passed a second integrity/weight check.

This clears the single-replicate smoke gate only. It does not clear:
1. Full exact-head PR #59 CI.
2. Corrected multi-replicate orchestration and outer audit.
3. Seven-gamma density pilot.
4. Bandwidth training/holdout separation.
5. Any H0 propagation or calibration.

Do not select a bandwidth or make scientific claims from this smoke. The next valid computational step, after full PR #59 CI completes, is a small corrected multi-replicate orchestration smoke with independent outer audit. Only after that should the seven-gamma density-only experiment be considered.
