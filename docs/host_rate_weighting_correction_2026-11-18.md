# Host-rate weighting correction — follow-up record (18 November 2026)

## Finding

A source-level trace found that the density smoke runner supplied the same host/merger-rate function b(z) twice:
1. as parent/catalogue sample-level base weights to the five-arm builder; and
2. as the post-KDE merger_weight multiplier to the density estimator.

The truth density was defined proportional to g(z)b(z), so the prior implementation's large-sample target was proportional to g(z)b(z)^2. This made the old A-E density metrics unsuitable as evidence for recovery of the declared target.

## Correction submitted in PR #59

The runner now omits parent_base_weights and catalogue_base_weights when building the five arms. Consequently:
- A and B use unit sample weights;
- C uses 1/pi_true;
- D uses 1/pi_MAP;
- E uses node-specific 1/pi(theta_k).

The common b(z) factor is applied once by the KDE density estimator after the galaxy-density estimate. The independent raw-input auditor was updated to reconstruct exactly these inclusion-only sample weights, and a regression test asserts that the A and B sample weights saved by the runner are all ones.

## Scope and validation boundary

- Code branch: research/fix-host-weight-double-application-20261118
- Pull request: #59 (draft, unmerged)
- Changes: density smoke runner, independent density smoke auditor, and regression test.
- No main-branch changes and no merge.
- Tests have not yet been executed in this environment.
- No seven-gamma run has been dispatched.
- Existing density artifacts generated before this correction must be treated as superseded for the stated g(z)b(z) target and regenerated before any density conclusion or bandwidth selection.

## Remaining follow-up

1. Obtain green exact-head CI for PR #59.
2. Run the corrected single-replicate density smoke and independently audit its raw artifacts.
3. Run the corrected multi-replicate pilot and independently audit all replicate records.
4. Reassess the 11-bandwidth grid and only then prepare a fresh-seed holdout.
5. Keep H0 propagation blocked until density-stage estimands and results are independently validated.
