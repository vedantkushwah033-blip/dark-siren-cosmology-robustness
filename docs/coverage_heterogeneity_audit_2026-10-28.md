# Independent reconstruction of fixed-truth coverage heterogeneity — 2026-10-28

## Scope and provenance

This report independently reconstructs the 68% and 95% interval-coverage indicators from the archived raw replicate CSV, without importing the simulation runner, likelihood implementation, or the experiment's summary JSON. The source is the `distance_uncertainty_sensitivity_full/replicates.csv` file in GitHub Actions artifact `numerical-diagnostics-py3.12`, run [40630753593](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40630753593), on PR #4 head `54a6a19011b7b870e0d3849a27d93529507630c1` (workflow's tested merge-ref revision recorded in the artifact metadata: `1420ec1b54266588f863992edf809dabec0bd5dc`).

- Raw records: 32,400 (three inference arms × nine fixed truths × four fractional-error levels × 300 replicates).
- Raw CSV SHA-256: `62af121350d62b141f0f93917258e1967cd51b611234e6c4c9839f3088e2a4da`.
- Record key `(true_h0, fractional_error, replicate, arm)` is unique.
- Coverage was reconstructed using the closed-interval rule `lower <= true_h0 <= upper`.
- The 95% intervals below are Wilson binomial intervals for the descriptive coverage proportions. They quantify Monte Carlo uncertainty conditional on these generated records, not model uncertainty.

The audit script is `scripts/audit_full_run_coverage_heterogeneity.py`; its tests are `tests/test_audit_full_run_coverage_heterogeneity.py`. The script deliberately does not label a pooled result as calibrated or astrophysically valid.

## Exact reconstruction of the reported discrepancy

| Inference arm | All truth cells: successes / n | All-truth 68% coverage (95% Wilson CI) | Interior truths 55–85: successes / n | Interior 68% coverage (95% Wilson CI) | Exact boundaries 50 and 90: successes / n |
|---|---:|---:|---:|---:|---:|
| A — fixed realized sigma | 6,156 / 10,800 | 57.00% (56.06–57.93%) | 6,156 / 8,400 | 73.29% (72.33–74.22%) | 0 / 2,400 |
| B — fractional sigma(H0) | 6,184 / 10,800 | 57.26% (56.32–58.19%) | 6,184 / 8,400 | 73.62% (72.67–74.55%) | 0 / 2,400 |
| C — fixed reference sigma | 6,163 / 10,800 | 57.06% (56.13–58.00%) | 6,163 / 8,400 | 73.37% (72.41–74.30%) | 0 / 2,400 |

**The reported 57% versus 73% pattern is reproduced from the actual CI raw artifact.** For each arm, every successful interval belongs to an interior-truth replicate: the 2,400 records at the two exact prior endpoints have zero coverage. Adding those 2,400 boundary cases to the denominator reduces the displayed pooled rate from roughly 73% to roughly 57% without changing the numerator.

These two percentages answer different questions. The all-truth number is an equal-cell mixture of fixed-truth experiments that includes exact prior endpoints; it is not a single fixed-parameter coverage estimate. The interior number is also an average across seven distinct truths and four noise fractions, so it is not evidence of uniform coverage over the interior.

## Why the exact endpoints have zero coverage

The stored intervals are central equal-tail posterior intervals. For a continuous posterior supported on `[50, 90]`, its 16th and 84th percentiles lie strictly inside the support; thus a fixed truth exactly equal to 50 or 90 is outside the interval in every replicate. This is a property of this interval rule at the exact boundary, not evidence that a real dark-siren analysis has zero coverage generally.

The separately supplied analytic demo reports approximately 40% endpoint coverage under a different clipped numerical-interval construction. That number must **not** be transferred to this production artifact: it uses a different interval construction. The production raw records and independent reconstruction above show 0/1,200 coverage per arm at each endpoint. Any future boundary-method comparison must define the interval construction explicitly and test it against an analytic reference before interpreting the result.

## Interior coverage is heterogeneous

Coverage by true H0, averaged over the four fractional-error settings (1,200 replicates per arm and truth), is:

| True H0 | A: fixed realized sigma | B: fractional sigma(H0) | C: fixed reference sigma |
|---:|---:|---:|---:|
| 55 | 66.67% | 67.33% | 68.58% |
| 60 | 70.67% | 70.83% | 70.58% |
| 65 | 70.33% | 70.17% | 69.58% |
| 70 | 71.50% | 72.50% | 71.50% |
| 75 | 78.00% | 78.17% | 78.58% |
| 80 | 79.08% | 79.83% | 80.33% |
| 85 | 76.75% | 76.50% | 74.42% |

The interior mean therefore conceals systematic variation across truth cells. The per-truth figures also average over four different noise fractions, so they should be accompanied by the truth-by-noise table when making any substantive comparison. The full audit JSON preserves both truth-only and truth-by-fraction breakdowns, including 95% coverage.

## Numerical grid-resolution finding is distinct

The same Actions artifact includes a separate deterministic stress test comparing a 1,601-point H0 grid (spacing 0.025) with a 400,001-point reference grid (spacing 0.0001) for deliberately narrow single-host posteriors. It flags four of five tested noise scales as under-resolved; at the narrowest scales the posterior mean/median can look stable while interval endpoints and posterior SD are materially wrong. This is a valid numerical stress result for the simplified model, not proof that the full production experiment's grid is inaccurate. In fact, the archived 32,400-row experiment separately passed its 1,601-versus-3,201-grid check, with coverage unchanged at reported precision. Keep these claims distinct: one is a narrow-posterior stress test; the other is a stability check for the observed simulated records.

## Interpretation and next gate

1. Do not present 57% as an unqualified calibration estimate. State the truth-grid composition and boundary rule, and report boundary and interior strata separately.
2. Do not present 73% as proof of calibration: it averages heterogeneous fixed truths and noise fractions; several truth/fraction cells depart appreciably from 68%.
3. Report cell-level intervals and the full truth × noise breakdown. The independent dataset replicate is the unit for Monte Carlo uncertainty.
4. Keep this one-known-host, low-redshift Gaussian model clearly separated from catalogue-incompleteness/host-mixture claims. This audit does not establish real LVK validity, novelty, or publication readiness.
5. Before a confirmatory claim, prespecify the primary estimand and truth/noise weighting, define boundary handling, and distinguish a numerical-resolution check from scientific model adequacy.


## Follow-up: matched frequentist reference interval

A separate independent calculation now compares the stored Bayesian equal-tail intervals with a Gaussian distance confidence interval inverted through \(d=cz/H_0\), using the same Arm A observations and realized sigma. The interior-truth Bayesian coverage is 73.29% at nominal 68% and 97.12% at nominal 95%; the inverted Gaussian reference gives 67.62% and 95.12%, respectively. At the central truth \(H_0=70\), the Bayesian coverage is 71.50% / 98.17%, versus 66.33% / 95.50% for the reference intervals. This supports the methodological distinction between Bayesian credible coverage and frequentist confidence coverage; it does not prove universal posterior calibration. The reference uses Arm A's oracle-like known sigma and is not an operational correction. Full method, counts, Wilson intervals, assumptions, and caveats: [fixed-truth coverage reference-interval audit](fixed_truth_coverage_reference_interval_audit_2026-10-28.md).
