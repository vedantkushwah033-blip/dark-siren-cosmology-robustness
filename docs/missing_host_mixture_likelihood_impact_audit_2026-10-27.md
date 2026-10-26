# Impact audit: missing-host mixture likelihood correction

**Audit date:** 2026-10-27  
**Status (updated 2026-10-27):** the host-misspecification audit has been rerun from the corrected implementation and its independent reconstruction passed. That experiment's historical artifact was generated at a source revision that predates the transient regression and matches the corrected rerun to floating-point precision. The other four directly identified hidden-host experiments still require separate regenerated artifacts and audits; do not treat this one successful rerun as clearing them.

## Finding

For a catalogue component with normalized internal host weights and a missing-population component with normalized internal weights, the intended mixture is

\[
L(H_0)=(1-f_{\rm miss})L_{\rm cat}(H_0)+f_{\rm miss}L_{\rm miss}(H_0).
\]

The direct `posterior_h0(..., missing_weight=f, missing_redshifts=...)` implementation previously applied the catalogue mass factor twice. That gave a relative mixture proportional to `(1-f)^2 L_cat + f L_missing`, rather than the intended expression. Commit `bf0d14fceecdae63a9f05a9913ae36f192e5d8c0` removed the extra catalogue factor. Commit `aadd8532087a56925c5fe9aa6edc9de1ffcfd82f` adds a regression test comparing the API output against an explicitly concatenated full population, with a tight numerical tolerance.

This is a model-correctness change, not a cosmetic refactor. Any previously generated posterior or summary that used the affected direct API with a nonzero missing fraction must be treated as provisional until rerun from the corrected revision.

## Directly identified experiment entry points

The following branch files call the corrected API with a nonzero missing fraction and therefore require regenerated artifacts if their published outputs were produced before the correction:

- `scripts/run_hidden_host_assumption_sensitivity.py` — associated report: `docs/hidden_host_assumption_sensitivity_results_2026-10-26.md`.
- `scripts/run_host_misspecification_audit.py` — associated report: `docs/host_misspecification_results_2026-10-26.md`.
- `scripts/run_paired_host_misspecification.py` — associated paired host-misspecification results/protocol.
- `scripts/run_host_misspecification_sensitivity_surface.py` — associated report: `docs/host_misspecification_sensitivity_surface_results_2026-10-26.md`.
- `scripts/run_host_prior_support_factorial.py` — hidden-host × prior-support factorial results and its independent audit.

This list is a confirmed minimum, not a claim that all downstream callers have been exhaustively enumerated. Other scripts or wrappers may call the API indirectly; inspect their call graph before declaring the impact inventory complete. The `run_oracle_missing_host_control.py` script uses a different `missing_host_weighting` argument and is not assumed affected by this particular code change without tracing its implementation.

## Rerun evidence and corrected scope

The host-misspecification stress test was rerun on source revision `0ac70cba1708be82c7b994145f4b0507101afa05` in [workflow run 40416457719](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40416457719). The run passed the experiment, file validation, and independent reconstruction steps for 4,000 records. Its artifact is available at [artifact 12694878568](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40416457719/artifacts/12694878568), ZIP SHA-256 `e56cbeebe4b353eab78bfb533117142fc138f163e0719575264ee416ba910c20`.

A direct comparison against the previously published artifact (ZIP SHA-256 `c627c12bcad2c41db75769f3566dd2a5cc1b1bd3e9988332760ac901462f094c`) found identical scenario-level metrics at the report's displayed precision and only floating-point-rounding differences in per-run posterior summaries (maximum absolute difference below `1.5e-13`). The old report names source revision `5caad04955d735583f77d2c56e455aa28ad43b40`; its linear-space likelihood predates the later log-space rewrite that briefly introduced the extra `log1p(-missing_weight)` factor. Thus this specific historical experiment is reproduced by the corrected run and should not be described as numerically invalidated by that later transient bug.

This finding narrows the earlier blanket warning; it does not clear other reports. The direct-call inventory below remains a minimum set for separate revision/artifact verification. For each other experiment, inspect its recorded source revision and regenerate it from the corrected code before reusing its numerical claims.

## What is not established by this audit

- No numerical claim is made here about how much any historical result will change.
- Existing independent reconstruction audits establish that an artifact matches the implementation and seed that generated it; they do not protect against a later-discovered error in the implementation itself.
- Passing unit tests does not establish that the model's assumptions describe real gravitational-wave events or galaxy catalogues.
- Known-host-only experiments that never invoke the missing-host component are not automatically affected by this fix, but should still be linked to the exact revision used for their results.

## Reproduction commands

Run from the repository root in a clean environment after installing the package and development dependencies (`python -m pip install -e ".[dev]"`). Keep each output directory separate; do not overwrite or mix pre-correction artifacts with regenerated files.

```bash
python scripts/run_hidden_host_assumption_sensitivity.py --replicates 200 --seed 20261032 --bootstrap-resamples 1000 --output-dir results/hidden-host-assumption-sensitivity
python scripts/audit_hidden_host_assumption_sensitivity.py --input-dir results/hidden-host-assumption-sensitivity --replicates 200 --seed 20261032 --bootstrap-resamples 1000

python scripts/run_host_misspecification_audit.py --replicates 1000 --seed 20261027 --output-dir results/host-misspecification-audit
python scripts/audit_host_misspecification_audit.py --input-dir results/host-misspecification-audit --replicates 1000 --seed 20261027

python scripts/run_host_misspecification_sensitivity_surface.py --replicates 200 --seed 20261029 --bootstrap-resamples 1000 --output-dir results/host-misspecification-sensitivity-surface
python scripts/audit_host_misspecification_sensitivity_surface.py --input-dir results/host-misspecification-sensitivity-surface --replicates 200 --seed 20261029 --bootstrap-resamples 1000

python scripts/run_paired_host_misspecification.py --replicates 1000 --seed 20261027 --bootstrap-resamples 4000 --output-dir results/paired-host-misspecification
python scripts/audit_paired_host_misspecification.py --input-dir results/paired-host-misspecification --replicates 1000 --seed 20261027 --bootstrap-resamples 4000

python scripts/run_host_prior_support_factorial.py --replicates 200 --seed 20261030 --bootstrap-resamples 1000 --output-dir results/host-prior-support-factorial
python scripts/audit_host_prior_support_factorial.py --input-dir results/host-prior-support-factorial --replicates 200 --seed 20261030 --bootstrap-resamples 1000
```

The commands above use the runners' declared default design sizes and seeds. A successful process exit is necessary but not sufficient: retain the raw records and summaries, capture the tested Git revision and dependency versions, and verify each independent auditor's success before updating any result narrative. If any audit fails, stop and investigate rather than relaxing tolerances.

## Required release gate

1. Confirm the explicit full-mixture regression test passes on Python 3.10, 3.11, and 3.12.
2. Host-misspecification audit: **completed** on revision `0ac70cba1708be82c7b994145f4b0507101afa05`; its historical artifact matches at reported precision. Rerun each of the other affected experiments from the corrected commit with its declared seed and design.
3. Run the independent auditor against the newly generated raw artifacts; record hashes, source revision, runtime versions, and audit result.
4. Recompute marginal summaries and paired intervals from the regenerated raw records; record any historical artifact that is shown to predate the regression and match the corrected rerun.
5. Replace historical interpretations only after the regenerated results pass audit; preserve old artifacts as superseded for traceability.
6. Do not claim external-method validation or publication readiness until a benchmark with verified input provenance and an independently checked reference implementation has been completed.

## Provenance

- Likelihood correction: [commit bf0d14f](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/commit/bf0d14fceecdae63a9f05a9913ae36f192e5d8c0)
- Explicit mixture regression: [commit aadd853](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/commit/aadd8532087a56925c5fe9aa6edc9de1ffcfd82f)
- Working pull request: [PR #3](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/3)
