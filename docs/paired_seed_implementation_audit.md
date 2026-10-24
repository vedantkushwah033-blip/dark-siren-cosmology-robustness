# Paired-seed implementation audit

**Review date:** 2026-10-24  
**Scope:** `scripts/run_paired_seed_sensitivity.py`, the pinned original simulator at commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`, and the paired-seed protocol.  
**Status:** Static source review; this document does not claim that the paired-seed experiment has been executed.

## Findings

### 1. The pinned simulator commit is retrievable

The original repository contains commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. The relevant function is `src/experiments.py::run_repeated`.

### 2. What the same seed actually shares

For replicate index `r`, `run_repeated` sets `s = seed + r` and uses:

- catalogue seed: `s`;
- host-selection seed: `s + 100000`;
- gravitational-wave measurement seed: `s + 200000`;
- catalogue-mask seed: `s + 300000`.

The catalogue is generated before the scenario mask is applied. Therefore, with the same base seed, replicate count, galaxy count, and true H0, catalogue generation is deterministic and uses the same seed across scenarios. The random-number generators are invoked in separate functions with explicit seeds.

The host-selection and measurement draws are shared in the intended sense only when host-generation settings match. When host-generation weighting changes, the selected host can change, and therefore the injected event distance can change. That is expected for the host-weighting contrast, not a failure of pairing.

### 3. Interpretation of the three prespecified contrasts

- **Random 70% versus complete catalogue:** same parent catalogue and same host-generation model; the mask changes. This is a paired comparison if the generated host and GW measurement draws remain identical, as intended by the seed scheme.
- **Redshift-dependent 70% versus random 70%:** same parent catalogue and same host-generation model; the selection transformation changes. This is a paired comparison of the two specified mechanisms, conditional on this simulator.
- **Uniform versus matched luminosity inference weighting:** both scenarios generate hosts using luminosity weighting, so the host-generation model is held fixed while inference weights change. This contrast is paired under the current seed scheme.

Using identical seeds alone is not a universal proof of pairing. The above assessment depends on the inspected implementation and its deterministic, explicitly seeded calls.

### 4. Known limitations and validation gaps

- The paired-seed script's current unit tests test contrast arithmetic and alignment, but do not execute the original simulator to test shared parent catalogues.
- The output metadata describes the pairing strategy but does not save a parent-catalogue fingerprint or explicit RNG-stream identifiers per replicate.
- The original simulator calls its runs independent mock universes within one scenario. Across scenarios, the same seed reconstructs corresponding universes; this distinction should be stated clearly.
- The simulator skips a replicate if no observed galaxies remain. With 500 galaxies and a 70% target completeness this is unlikely, but output row counts must be checked rather than assumed.
- The paired script calculates a mean and standard error for all metrics. Coverage differences are paired binary outcomes; a final report should add an interval suited to paired binary data rather than relying only on the standard error.
- The current true value is fixed at H0 = 70. This fixed-parameter coverage study is not prior-based simulation-based calibration (SBC).
- These tests concern a synthetic simulator. They do not validate the realism of the galaxy population, the survey selection function, or real LVK data analysis.

## Recommended next steps

1. Run the paired script at a small replicate count (e.g. 20–30) against the pinned simulator.
2. Validate all expected scenario/replicate counts, finite posterior summaries, interval ordering, realized completeness, and deterministic reruns.
3. Add an integration test that confirms identical parent catalogue arrays for matched seeds and scenario configurations, and confirms that the masks implement their intended transformations.
4. Preserve raw per-run results and metadata before increasing the replicate count.
5. For coverage, report scenario-level binomial uncertainty and paired coverage differences with a paired-data interval.
6. Only after those checks pass, run the full prespecified experiment and interpret effect sizes alongside Monte Carlo uncertainty.

## Scope of this audit

This is a source-level implementation audit, not a run report. No numerical result, calibration claim, or scientific conclusion is established by this document.
