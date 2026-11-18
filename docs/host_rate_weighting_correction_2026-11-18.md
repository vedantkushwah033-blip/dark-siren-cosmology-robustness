# Host-rate weighting correction integrated into the quadrature branch

Date: 2026-11-18

## Source finding

The absolute-bandwidth density runner previously applied the common host/merger
factor b(z) twice: in parent/catalogue sample weights and again through
`merger_weight` after KDE. Since the declared target is proportional to
g(z)b(z), that implementation instead approached g(z)b(z)^2. Old outputs
under that convention are superseded and must not be used for density claims
or bandwidth selection.

## Correction now ported to PR #60

The quadrature implementation branch now follows the corrected convention
from companion PR #59:

- A (catalogue-only) and B (full-parent oracle) raw sample weights are ones.
- C, D, and E raw sample weights contain only their respective inverse
  inclusion probabilities.
- The common b(z) factor is applied once after KDE through
  `merger_weight`.
- The independent raw-input auditor reconstructs the same convention rather
  than trusting runner-produced metrics.
- A density-runner regression checks A/B sample weights are exactly one.
- A selected-likelihood regression checks the post-KDE rate convention.
- The curved-logit H0 quadrature regression now uses that same convention.

The likelihood API documentation clarifies that pre-KDE and post-KDE
host-rate weighting are not identical finite-bandwidth estimators. Do not
combine results produced under the two conventions.

## Validation status

This correction has been committed to the PR #60 working branch, not merged
to the default branch. The focused curved-DGP test passed on an earlier
pre-correction commit; the corrected density, likelihood, and H0 regressions
must pass exact-head CI before this branch can be considered validated.
The seven-gamma density experiment remains blocked until exact-head tests pass
and a corrected smoke output is independently audited.

## Scope

Synthetic low-redshift Euclidean methodology only. This correction does not
constitute an observational H0 measurement, public LVK/GWTC analysis, or
evidence of estimator superiority.
