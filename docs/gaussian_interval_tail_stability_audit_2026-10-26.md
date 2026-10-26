# Gaussian interval tail-stability audit — 2026-10-26

## Motivation

The toy selection model computes the probability that a noisy distance falls inside an interval. For standardized Gaussian bounds a < b, the probability is P(a < Z < b) = Phi(b) - Phi(a). Direct CDF subtraction can lose precision when both bounds lie in the positive tail (both CDF values are near one) or when both lie far into the negative tail (both CDF values are tiny and close). The positive-tail survival-function identity was already used in the implementation. A focused review identified that the negative-tail branch still used direct subtraction, despite the module's numerical-robustness role.

## Change

In src/dark_siren_robustness/selection.py, the helper now uses three cases:

- If a >= 0, use Phi(-a) - Phi(-b), avoiding subtraction of values near one.
- If b <= 0, use the log-CDF difference: Phi(b) * [-expm1(log(Phi(a)) - log(Phi(b)))]. This preserves relative precision when the log-CDF values are close, while still returning zero if the true interval mass is below floating-point representability.
- If the interval crosses zero, use the direct CDF difference, where cancellation is not severe.

Non-finite standardized bounds now fail explicitly instead of flowing into an apparently ordinary interval probability. The public detection-efficiency routine continues to fail closed when the resulting efficiency is not representable as a positive floating-point probability; it does not substitute an artificial machine-tiny value.

## Regression coverage and provenance

Added tests cover:

1. A narrow, representable negative-tail interval near -37.5, compared against a log-CDF/expm1 reference.
2. Rejection of non-finite interval bounds.
3. Existing moderate negative-tail agreement with direct CDF differences.
4. Existing positive-tail preservation and explicit fail-closed behavior when probability underflows.

The code and tests are in commits af8cf99b12b9f3bc013e4d1be4a489be330cc77c and 1e6d207b69f9b8b0cbe365ad6ef3a2614eab30ea. The numerical-diagnostics workflow passed on Python 3.10, 3.11, and 3.12 for the test commit:

- [Run 40271904706 — negative-tail regression tests](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40271904706)
- [Run 40271902279 — interval-probability implementation](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40271902279)

## Scope

This is a numerical-stability improvement for a simplified Gaussian selection model. It does not make the selection function a realistic gravitational-wave detector model, validate a galaxy catalogue, or establish astrophysical correctness. The change protects interval-probability calculations in difficult tails; scientific conclusions still require model validation and stress tests under relevant assumptions.
