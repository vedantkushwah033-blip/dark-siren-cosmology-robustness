# Corrected quadratic-logit curvature: two-replicate smoke audit

Date: 18 November 2026

Status: **TWO-REPLICATE CURVATURE-SPECIFIC ORCHESTRATION SMOKE PASSED; RAW ARTIFACTS INDEPENDENTLY INSPECTED. DEVELOPMENT EVIDENCE ONLY.**

## Provenance

- Source code is the exact PR #59 correction head: `390ba806c42376b1a2d4b8d083b92b8d206c8a88`.
- Validation branch: `research/curvature-multireplicate-smoke-20261118`.
- Workflow commit: `41330c66ecc53712aeaf2522ba08139f3adbf138`.
- Dedicated workflow run: [44366177286](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44366177286).
- Workflow result: **completed / success**.
- Artifact ID: `14456104491`.
- Artifact name: `corrected-curvature-multireplicate-smoke-41330c66ecc53712aeaf2522ba08139f3adbf138`.
- Artifact size: 340,964 bytes.
- GitHub archive SHA-256: `531dfaeba4dfcb22de8215200fce8e6c080fe567dbffd9a94745a14dc02e42a4`.
- Artifact expiry shown by GitHub: 2 December 2026.

This branch was created from the exact PR #59 head. No code on `main` was changed, no PR was merged, and no other research branch was modified.

## Executed configuration

| Parameter | Value |
|---|---:|
| Replicates | 2 |
| Parent population per replicate | 700 |
| Independent validation frame | 300 |
| Nested validation sizes | 80, 160 |
| Catalogue inclusion model | Quadratic logistic |
| Curvature | `gamma = 0.5` |
| Linear coefficient | `beta = -3.9` |
| Quadrature order | 3 |
| Redshift grid | 101 points |
| Absolute KDE widths | 0.012, 0.024 |
| Scott multiplier | 0.7 |
| Master seed | 2026111811 |
| Seed stride | 10 |

The parent-rate calibration was `alpha = 0.06312269530994806`, target inclusion rate `0.5101919962291073`, achieved rate `0.5101919962291072`, and centering constant `E_g[x^2] = 0.2845005579191241`.

Replicate seeds:
- Replicate 0: population 2026111811; validation 2026111812.
- Replicate 1: population 2026111821; validation 2026111822.

## Workflow checks

All dedicated workflow steps succeeded:
1. Checkout.
2. Python 3.12 setup.
3. Package and development-dependency installation.
4. Two-replicate quadratic-logit curvature-density run.
5. Independent outer orchestration audit, including independent raw-input audit of each successful child.
6. Artifact upload.

The outer auditor reported:
- `status: PASS`.
- 2 planned, 2 successful, 0 failed replicates.
- 2 child raw audits passed.
- 26 outer-manifest files verified.
- 120 aggregate density-metric rows reconstructed.
- 52 Scott-width rows reconstructed.
- No H0 inference or bandwidth selection.

## Independent inspection of the downloaded ZIP

I downloaded artifact ID `14456104491` and inspected the archive independently.

- Archive contained 30 files.
- All 26 outer-manifest entries matched their declared byte counts and SHA-256 digests.
- Both child manifests contained 11 entries each; all child file hashes and byte counts matched.
- Replicate ledger enumerated the two planned seeds with no failures.
- Both child configurations independently recorded `catalogue_selection_model = quadratic_logit`, `gamma = 0.5`, `beta = -3.9`, and the calibrated inclusion rate.

The corrected weight convention was verified directly from the saved raw arrays:
- Arm A catalogue-only sample weights were exactly one.
- Arm B parent-oracle sample weights were exactly one.
- Arm C weights matched `1 / pi_true` reconstructed from the retained parent-truth rows; maximum absolute discrepancy was (2.85\times10^{-14}) in replicate 0 and (1.43\times10^{-14}) in replicate 1.
- The common host-rate factor is applied once downstream of the KDE; the duplicate sample-level base weighting is absent from the corrected runner.

Realized catalogue and retained inclusion-probability summaries:
| Replicate | Catalogue size | Minimum retained pi | Maximum retained pi | Maximum oracle 1/pi |
|---:|---:|---:|---:|---:|
| 0 | 362 | 0.05565 | 0.98683 | 17.97 |
| 1 | 375 | 0.04810 | 0.98645 | 20.79 |

## Descriptive density-recovery values

The table gives mean IAE across only these two replicates for selected arms. These values are included to preserve the smoke output, not to rank estimators.

| Validation n | Absolute width | A catalogue | B parent oracle | C known pi | D degree-one MAP | E posterior-mean density |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 0.012 | 0.7916 | 0.1054 | 0.2208 | 0.3991 | 0.3784 |
| 80 | 0.024 | 0.7860 | 0.1254 | 0.1865 | 0.3907 | 0.3646 |
| 160 | 0.012 | 0.7916 | 0.1054 | 0.2208 | 0.3529 | 0.3400 |
| 160 | 0.024 | 0.7860 | 0.1254 | 0.1865 | 0.3406 | 0.3245 |

These are tiny-sample development diagnostics with only two independent population/validation streams, two bandwidth values, a 101-point grid, and quadrature order 3. They do not establish estimator superiority, a stable ranking, bandwidth optimality, calibration, or generalization. The oracle arm is a simulation diagnostic, not an implementable estimator. Degree-one D/E remain misspecified at nonzero curvature; E integrates coefficient uncertainty within the degree-one model and does not average over omitted quadratic model uncertainty.

## What this clears and what remains blocked

**Cleared:** the corrected quadratic-logit child runner, two-replicate orchestration, seed schedule, outer manifest, aggregate reconstruction, and per-child independent raw-input audits have executed successfully for this small `gamma=0.5` smoke.

**Not cleared:**
1. Full exact-head PR #59 CI (at the time of this note, Python 3.12 was still completing the last numerical-diagnostic steps).
2. Multi-replicate curvature pilot at production-like `N=4000`, validation frame 1500, sizes 100/300/1000, and adequate replicate count.
3. Seven-gamma curvature experiment.
4. Bandwidth training/holdout protocol.
5. Degree-two D2/E2 comparison.
6. H0 propagation, posterior calibration, PIT, coverage, or SBC.

The correct next step is to finish exact-head PR #59 CI, then run a larger but still clearly labelled development pilot at one or more curvature values. Do not promote this two-replicate smoke into a scientific conclusion.
