# Prior-support grid guard audit — 2026-10-29

## Finding

A source-level review found inconsistent prior-support contracts across three posterior APIs: `hardened_inference.posterior_h0`, `inference.posterior_h0_low_z`, and `cosmological_inference.posterior_h0_from_distance_grid`. The first checked only that at least two H0 grid points fell inside the prior; the other two likewise permitted incomplete or off-grid support and the latter also accepted non-positive H0 grid values. None consistently required the supplied computational grid to span the full prior and represent both endpoints.

This permitted two silent failure modes:

1. A grid such as [60, 70, 80] with a declared prior [50, 90] could be normalized over only the represented subset, despite the prior claiming support over [50, 90].
2. A grid that spanned [50, 90] but omitted either exact endpoint could silently approximate the discontinuous prior boundary using interior nodes. Since the function returns density only at the supplied grid nodes, downstream summaries cannot recover the omitted endpoint segments exactly from that returned array alone.

These are API-level numerical correctness risks, distinct from the already documented extreme narrow-likelihood grid-resolution stress test.

## Change made

All three posterior APIs now fail explicitly unless:
- the H0 grid is finite, strictly increasing, and physically positive;
- the H0 grid spans the full declared prior interval; and
- both prior endpoints occur exactly in the grid.

The exact-endpoint requirement is deliberately stricter than `dataset_posterior_h0`, which interpolates log likelihood to off-grid endpoints and normalizes on an augmented internal grid. Its returned density remains sampled on the original grid, so callers should follow its documented integration contract rather than assume trapezoidal integration over the returned array is exact when prior endpoints are off-grid.

The error message instructs callers to use an endpoint-inclusive grid. This is intentionally a fail-fast contract rather than silently changing the returned grid shape or interpolating posterior summaries behind the caller's back.

## Regression tests added

Regression tests cover incomplete lower/upper prior ranges, missing exact endpoints, and non-positive H0 grids in both the low-z and cosmology-aware APIs, alongside the original four hardened-inference cases.

Existing endpoint-aligned cases in the hardened-inference test suite remain the compatibility baseline. The targeted hardened-inference and prior-support grid-guard tests completed successfully in all three CI environments (Python 3.10, 3.11, and 3.12). The full `tests` workflow also completed successfully on commit `c8a56e4f662dc60b74a564620c958f4801322c28` (run 40900854372). This confirms the repository's configured workflow passed; it does not independently establish scientific validity or real-data performance. See [the workflow run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40900854372).

## Scope and limitations

This change does not establish correctness of every posterior calculation, resolve all discretization error, or validate the toy likelihood physically. It prevents a specific class of silent prior-support truncation/misrepresentation. Callers with off-grid prior endpoints must now construct an endpoint-inclusive grid before calling `posterior_h0`.

The selection-conditioned implementation already constructs an endpoint-inclusive normalization grid and returns it explicitly; this patch does not modify that API. The low-z and cosmology-aware guards are new commits, and their expanded tests are being checked by CI; no pass is claimed for those additions until the current workflow completes.

## Provenance

- Code commit: `e03de64df538713ab9411d7be27a9b17a686b5ef`
- Regression-test commit: `160afd5b0488039f419598b21acad4a5d0ca1caf`
- CI workflow update: `c8a56e4f662dc60b74a564620c958f4801322c28`
- Verified CI run: [40900854372](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40900854372), successful on Python 3.10, 3.11, and 3.12
- Additional low-z guard commit: `c520199890c15025e91d27c0e136002a973879de`
- Additional cosmology-aware guard commit: `fb378a6f75b6cd5256eff03d989050fc96720429`
- Low-z regression tests: `913e6c3ee339938bc7a2ec06eecdcdb5acaea497`
- Cosmology-aware regression tests: `f8e74da2e0c98bd91d9f71c19011befff3240e18`
- CI workflow update for the additional tests: `d25b34a1ec824e88b4b38c224bfabe4402f6ed61`
- Branch: `research/audit-posterior-grid-integrity`
