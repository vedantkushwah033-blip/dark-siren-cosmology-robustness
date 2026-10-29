# Independent Review Response: Posterior-Grid Convergence Audit

**Reviewed snapshot:** `78f3b5419569a21ed4ef0395d69e2c1a56157856`  
**Response branch:** `research/audit-posterior-grid-integrity`  
**Response status:** changes committed; final CI verification pending  
**Scope:** synthetic known-host low-redshift distance-uncertainty sensitivity diagnostic only

## Reviewer's release recommendation

The independent reviewer assigned a **conditional** release decision for the synthetic method study and explicitly stated that real-data validation is not established. The review reported no Critical or High issues, but identified Medium concerns in convergence tolerances, verdict logic, and stale CI evidence, plus Low concerns in boundary-mass masking and test strength.

This response does not reinterpret the review as approval for real-data claims.

## Findings and response

### Medium: convergence tolerance was not posterior-width-aware

**Finding:** endpoint changes were normalized by the endpoint value, which can make changes look negligible even when they are material relative to a narrow posterior interval.

**Response:** the convergence audit now includes both an absolute endpoint tolerance (default 0.01 in H0 units) and a maximum endpoint change relative to the fine-grid 95% posterior width (default 5%). Existing median and relative endpoint criteria remain in place. These thresholds are explicit diagnostic choices, not universal scientific guarantees.

### Medium: verdict used only the last pair and ignored coverage flips

**Finding:** coarse-to-fine comparisons could change 68% interval coverage classifications while the final pair still passed.

**Response:** the verdict now evaluates every adjacent comparison whose coarse grid is at or above the declared production resolution (default 1601 points). Any 68% coverage-classification flip at production resolution forces a FAIL. A threshold miss without a coverage flip produces REVIEW rather than PASS. The audit also records whether the maximum absolute change across posterior summaries shrinks from one production-resolution comparison to the next.

### Low: boundary mass used a node mask and floating-point-sensitive edge inclusion

**Finding:** a node mask could select only one point on the 801-point grid and report zero mass even when the intended fixed-width edge region contained nonzero posterior probability.

**Response:** boundary mass is now computed from the interpolated posterior CDF at a fixed H0-width from each endpoint. This avoids node-count dependence and exact floating-point equality tests. The diagnostic remains a flag, not a universal release threshold for all models.

### Low: multi-resolution test checked structure rather than numerical behavior

**Response:** regression assertions now check the production-resolution comparisons, width-normalized changes, coverage stability, and decreasing width-normalized endpoint changes. These tests are intended to catch regressions in the convergence contract; they do not substitute for a full campaign on scientific-scale artifacts.

### Medium: stale CI evidence

**Finding:** CI evidence cited in the pre-existing release-gate document was for an older commit, not the reviewed head.

**Response:** CI for the original reviewed snapshot was separately observed as successful on the repository's Python 3.10/3.11/3.12 matrix. Subsequent convergence changes have their own GitHub Actions run. **Do not treat the revised code as verified until that run completes successfully on the latest head.**

## Independent review results preserved as reviewer-reported evidence

The reviewer reported:
- 806 tests passed in 511 seconds on Python 3.13 with NumPy 2.5.3 at the reviewed head; this is not a Python 3.10–3.12 run.
- Nested grids of 801, 1601, 3201 and 6401 points showed maximum summary changes of approximately 9.3e-3, 2.3e-3 and 6.8e-4, with successive ratios around 4.17 and 4.16. The reviewer interpreted this as quadratic discretization convergence for the tested setup.
- An independent posterior reconstruction differed from production posterior means by at most 0.0012 over 15 datasets.
- A separate 1,000-draw calibration run had mean rank 0.507 ± 0.009 and descriptive chi-square p = 0.96.
- A negative control omitting H0-dependent detection normalization produced mean rank 0.561, a descriptive chi-square p-value near zero, and posterior bias around -1.85 in H0.

These are **results reported by the independent reviewer**, not results independently rerun by this response. The posterior-grid convergence tests discretization under the same likelihood formula; they do not establish that the model is physically complete or correctly specified.

## Remaining limitations and release gates

1. Confirm GitHub Actions passes for the latest response-branch head, including the new regression tests.
2. If CI fails, inspect the failing logs and correct the issue before claiming the convergence audit is validated.
3. Keep the PR open and unmerged until verification is complete.
4. Keep synthetic toy-model results clearly separated from real LVK/GWTC analysis.
5. Real-data validation and the external CHIMERA v2.2 benchmark remain unestablished/NO-GO until compatible inputs and an independently verified data pipeline are available.
6. Report coverage by truth stratum, including boundary and interior behavior; do not use pooled coverage as a substitute for stratified calibration.
7. Keep claims about convergence local to the tested model, data-generating assumptions, grid range, and tolerances.

## Current disposition

**Conditional for the synthetic method study; not validated for real data.** The code-level responses to the identified convergence and boundary-mass concerns have been committed. The final release decision remains contingent on a successful CI run for the exact latest head and inspection of any failures.
