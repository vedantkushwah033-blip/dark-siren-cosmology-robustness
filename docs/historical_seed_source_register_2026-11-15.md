# Historical Actions artifact source register — 2026-11-15

Status: **PARTIAL / NOT COMPLETE**. This register records known density-pilot workflow artifacts, not the complete history of every development run.

| Workflow run | Artifact ID | Artifact name | Bytes | Full GitHub SHA-256 | Status |
|---|---:|---|---:|---|---|
| 43619594762 | 14124565039 | absolute-width-density-pilot-floor-010-17466326a4be59ca7698c6c1f34529801e6e5507 | 66,001,638 | caeec4bad484aa7baf1b963c0f914900f96de39db8e3bb1bc1fe0bbb0777e753 | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43619594762 | 14123753124 | absolute-width-density-pilot-no-floor-17466326a4be59ca7698c6c1f34529801e6e5507 | 63,634,331 | 9a3ba516f9f63919d0aa49e57d5636cafcc322aed4aa752e90b9fc578c9d1c4b | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43619594762 | 14123414834 | absolute-width-density-pilot-floor-025-17466326a4be59ca7698c6c1f34529801e6e5507 | 69,496,090 | a78a66f79d8cab8ef0179d55ba07881c302666b8d3cc12add2b93b9d0811cdb7 | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43617052462 | 14123894369 | absolute-width-expanded-grid-no-floor-d94fdd62e514c81c20444196e166d2ebbb341a40 | 97,825,684 | 26f9ccc6d711cff45ea7ac4a130eba48dcde7b243185d2f6c2fc2bb6c9a8af86 | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43617052462 | 14123867382 | absolute-width-expanded-grid-floor-010-d94fdd62e514c81c20444196e166d2ebbb341a40 | 100,214,396 | e83f944b22cc153c1b6354515f42ed25978db2cf7b94b2e66c13f16e04768d1a | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43617052462 | 14123634445 | absolute-width-expanded-grid-floor-025-d94fdd62e514c81c20444196e166d2ebbb341a40 | 103,672,379 | f7eb0576f0b9dc5badafb2a3b3b136bc125106856345df4749e173f94669c627 | Downloaded; SHA-256 verified; 248 JSON/CSV files scanned |
| 43619577181 | 14123904013 | five-arm-kde-bandwidth-pilot | 8,157,332 | c840a912f1cfc7ddcef930d2c92984f748d34859452dc54112cda6df76985d36 | Downloaded; partial smoke/engineering payload inspected |
| 43619577181 | 14123779398 | numerical-diagnostics-py3.12 | 74,359,526 | 3fc0e848d45526cfbd222c6fc2fa9684a1fc4408390dc2c4b5538f85eb715dff | Available; extraction pending |
| 43619577181 | 14123673249 | numerical-diagnostics-py3.11 | 12,279,818 | b8b10c69a19857b9dfd837cd2b96ccf2531ee7b86475e171c8c3431370719206 | Available; extraction pending |
| 43619577181 | 14122914778 | numerical-diagnostics-py3.10 | 12,280,279 | f73e22abd9c6e88218052dd47a6dcb514e0f40454ebaab95a0a018c2c69f0e6c | Available; extraction pending |
| 43635486478 | — | Focused curvature CI run | — | — | No artifacts returned; failed job logs and checked-out commit are separate evidence |

## Partial artifact inspection

Artifact 14123904013 was downloaded and its ZIP directory listing inspected. It contains a parent-density pilot configuration, a selection configuration, absolute-bandwidth smoke inputs, and a two-replicate smoke payload. The seed values observed in its JSON/CSV configuration and ledger files include:

- parent-density pilot seed: 20261114;
- selection split seed: 20261115;
- absolute-bandwidth smoke seed pair: 20261114 / 20261115;
- two-replicate smoke population/validation pairs: 20261120 / 20261121 and 20261130 / 20261131.

The smoke replicate config also contains `master_seed=20261120` and `seed_stride=10`. The stride is a design parameter, not a random seed value; the inventory parser excludes `seed_stride`.

This is not a complete historical inventory: only one of the listed artifacts has been downloaded and inspected, and additional historical workflow runs/artifacts may exist. The candidate training/holdout schedule therefore remains **unverified**. Do not treat the absence of overlap in this partial sample as evidence of global disjointness.

## Required next source-register work

1. Enumerate all relevant workflow runs, including older floor-family, absolute-width, curvature, smoke, and manual pilot runs.
2. Record every artifact ID, full GitHub digest, run ID, head SHA, expiration state, and extraction result.
3. Download and verify each artifact against its GitHub digest before extraction.
4. Preserve the original ZIP and extraction log; scan only a copy.
5. Reconcile child-level seeds with parent configs and ledgers; distinguish actual seed values from master-seed derivation parameters such as stride.
6. Record expired, missing, or non-enumerable artifacts as unresolved blockers, never as empty seed sets.

**Gate decision: NOT VERIFIED. No training run and no holdout run is authorized by this partial register.**


## Verified artifact batch update — 2026-11-15

The six ZIP artifacts from runs 43619594762 and 43617052462 have now been downloaded. Their local SHA-256 digests match the GitHub artifact metadata exactly:

- Run 43619594762: all three artifacts verified; each contains 248 supported JSON/CSV files, 81 seed-field occurrences, 40 distinct integer seed values, and zero parse errors.
- Run 43617052462: all three artifacts verified; each contains 248 supported JSON/CSV files, 81 seed-field occurrences, 40 distinct integer seed values, and zero parse errors.
- All six artifacts contain the same 40 distinct values: the 20 population/validation seed pairs beginning with 20261215/20261216 and continuing in increments of 10. The expanded-grid artifacts reuse the same paired replicate seeds as the five-width artifacts; they are not independent seed draws.
- Across these six archives, the candidate training schedule and candidate holdout schedule have no numerical seed-value overlap with the extracted historical seed set.
- Duplicate references across floor arms, ledgers, and child configurations are provenance repetition, not extra independent seeds.

This establishes an internally consistent inventory for these six artifact archives only. It does **not** establish global historical completeness: the smoke/engineering artifact, numerical-diagnostics artifacts, curvature workflow artifacts, earlier smoke runs, and any other relevant historical Actions runs remain to be reviewed. The candidate schedule remains non-authorized until the source register is complete and every relevant artifact is accounted for.
