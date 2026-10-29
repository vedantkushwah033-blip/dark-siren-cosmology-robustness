# Prior-support grid guard audit — 2026-10-29

## Finding

A source-level review of `src/dark_siren_robustness/hardened_inference.py` found that `posterior_h0` checked only that at least two H0 grid points fell inside the declared prior interval. It did not require the computational grid to span the full prior support or to contain its endpoints.

This permitted two silent failure modes:

1. A grid such as [60, 70, 80] with a declared prior [50, 90] could be normalized over only the represented subset, despite the prior claiming support over [50, 90].
2. A grid that spanned [50, 90] but omitted either exact endpoint could silently approximate the discontinuous prior boundary using interior nodes. Since the function returns density only at the supplied grid nodes, downstream summaries cannot recover the omitted endpoint segments exactly from that returned array alone.

These are API-level numerical correctness risks, distinct from the already documented extreme narrow-likelihood grid-resolution stress test.

## Change made

The function now fails explicitly unless:
- the H0 grid spans the full declared prior interval; and
- both prior endpoints occur exactly in the grid.

The error message instructs callers to use an endpoint-inclusive grid. This is intentionally a fail-fast contract rather than silently changing the returned grid shape or interpolating posterior summaries behind the caller's back.

## Regression tests added

Four parameterized cases cover:
- a grid missing the lower prior range;
- a grid missing the upper prior range;
- a spanning grid missing exact prior endpoints; and
- a grid with neither exact endpoint represented.

Existing endpoint-aligned cases in the hardened-inference test suite remain the compatibility baseline. The targeted hardened-inference and prior-support grid-guard tests completed successfully in all three CI environments (Python 3.10, 3.11, and 3.12). The full `tests` workflow also completed successfully on commit `c8a56e4f662dc60b74a564620c958f4801322c28` (run 40900854372). This confirms the repository's configured workflow passed; it does not independently establish scientific validity or real-data performance. See [the workflow run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40900854372).

## Scope and limitations

This change does not establish correctness of every posterior calculation, resolve all discretization error, or validate the toy likelihood physically. It prevents a specific class of silent prior-support truncation/misrepresentation. Callers with off-grid prior endpoints must now construct an endpoint-inclusive grid before calling `posterior_h0`.

The separate selection-conditioned implementation already constructs an endpoint-inclusive normalization grid and returns it explicitly; this patch does not modify that API.

## Provenance

- Code commit: `e03de64df538713ab9411d7be27a9b17a686b5ef`
- Regression-test commit: `160afd5b0488039f419598b21acad4a5d0ca1caf`
- CI workflow update: `c8a56e4f662dc60b74a564620c958f4801322c28`
- Verified CI run: [40900854372](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40900854372), successful on Python 3.10, 3.11, and 3.12
- Branch: `research/audit-posterior-grid-integrity`
