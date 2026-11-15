# Historical Actions artifact source register — 2026-11-15

Status: **PARTIAL / NOT COMPLETE**. This register records artifacts that have been identified from the known density-pilot workflow runs. It is not a complete history of every development run in the repository.

| Workflow run | Artifact ID | Artifact name | Size | GitHub digest | Current status |
|---|---:|---|---:|---|---|
| 43619594762 | 14124565039 | absolute-width-density-pilot-floor-010-17466326a4be59ca7698c6c1f34529801e6e5507 | 66,001,638 B | caee...7753 | Available; extraction pending |
| 43619594762 | 14123753124 | absolute-width-density-pilot-no-floor-17466326a4be59ca7698c6c1f34529801e6e5507 | 63,634,331 B | 9a3b...1c4b | Available; extraction pending |
| 43619594762 | 14123414834 | absolute-width-density-pilot-floor-025-17466326a4be59ca7698c6c1f34529801e6e5507 | 69,496,090 B | a78a...1db7 | Available; extraction pending |
| 43617052462 | 14123894369 | absolute-width-expanded-grid-no-floor-d94fdd62e514c81c20444196e166d2ebbb341a40 | 97,825,684 B | 26f9...af86 | Available; extraction pending |
| 43617052462 | 14123867382 | absolute-width-expanded-grid-floor-010-d94fdd62e514c81c20444196e166d2ebbb341a40 | 100,214,396 B | e83f...6d1a | Available; extraction pending |
| 43617052462 | 14123634445 | absolute-width-expanded-grid-floor-025-d94fdd62e514c81c20444196e166d2ebbb341a40 | 103,672,379 B | f7eb...c627 | Available; extraction pending |
| 43619577181 | 14123904013 | five-arm-kde-bandwidth-pilot | 8,157,332 B | c840...d36 | Downloaded; partial smoke/engineering payload inspected |
| 43619577181 | 14123779398 | numerical-diagnostics-py3.12 | 74,359,526 B | 3fc0...5dff | Available; extraction pending |
| 43619577181 | 14123673249 | numerical-diagnostics-py3.11 | 12,279,818 B | b8b1...9206 | Available; extraction pending |
| 43619577181 | 14122914778 | numerical-diagnostics-py3.10 | 12,280,279 B | f73e...0e6c | Available; extraction pending |
| 43635486478 | — | Focused curvature CI run | — | — | No workflow artifacts returned; failed job logs and checked-out commit remain separate evidence |

The digest column abbreviates GitHub's artifact SHA-256 value for readability. The full digest is available from the workflow-artifact metadata endpoint and should be recorded in the machine-readable extraction manifest before any artifact is accepted.

## Partial artifact inspection

Artifact 14123904013 was downloaded and its ZIP directory listing inspected. It contains a parent-density pilot configuration, a selection configuration, absolute-bandwidth smoke inputs, and a two-replicate smoke payload. Seed-like values observed include:

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
