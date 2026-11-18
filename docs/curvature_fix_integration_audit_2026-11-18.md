# Curvature-branch integration audit for host-rate weighting correction

Date: 2026-11-18

Status: **SOURCE-LEVEL AUDIT COMPLETE; CI AND EXECUTABLE ARTIFACT VALIDATION NOT COMPLETE**

## Purpose

Follow up the host-rate double-application finding in PR #59 and determine whether the corrected implementation is actually present in the branch that contains the seven-gamma curvature workflow.

## Exact repository evidence inspected

- Correction branch: `research/fix-host-weight-double-application-20261118`.
- PR: [#59](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/59).
- Correction head observed: `4902016d817de66be450d02230ce06acab86b738`.
- Base branch: `research/curvature-known-probability-fix-20261114`, whose observed base SHA was `95e8164bd1c17756e5416566ec3a28475b15f1b9`.
- GitHub compare reports the correction head is five commits ahead of the base and changes four files: the smoke runner, its independent auditor, its regression tests, and a correction note.
- The seven-gamma workflow file was fetched from the base branch. It runs `scripts/run_absolute_bandwidth_density_replicates.py`, which invokes the single-replicate smoke runner, followed by the multi-replicate auditor.

## Source-level result

### Corrected branch

The corrected smoke runner still computes and archives `base_p` and `base_c` for audit context, but no longer supplies them as `parent_base_weights` or `catalogue_base_weights` to `build_five_arm_weights`.

Consequently, the five-arm builder receives its default equal base weights:
- A: unit catalogue weights;
- B: unit parent weights;
- C: (1/\pi_{true});
- D: (1/\pi_{MAP});
- E: node-specific (1/\pi(\theta_k)).

The runner then passes the common (b(z)) factor once through `estimate_host_density(..., merger_weight=merger)`.

The corrected independent auditor reconstructs A and B as unit vectors and C/D/E from the corresponding inverse inclusion probabilities. The correction branch also contains a regression test that runs the smoke and auditor and asserts A/B weights are all ones.

The multi-replicate auditor calls the independent single-replicate auditor for every successful child, so the corrected child-level convention will be checked once the corrected single-replicate runner and auditor are integrated into that branch.

### Curvature base branch

The currently inspected curvature base branch still passes `parent_base_weights=base_p` and `catalogue_base_weights=base_c` into the five-arm builder. Its auditor still reconstructs A/B as those base-weight vectors. Its test file does not contain the new single-application regression test.

Therefore, **the seven-gamma workflow on the curvature base branch is still wired to the pre-correction implementation**. Merely having PR #59 open does not change the code checked out by that workflow.

## CI status

For correction head `4902016d817de66be450d02230ce06acab86b738`, GitHub reports workflow run `44360570755` (run number 5553, workflow `tests`) as **completed / cancelled**. The jobs and artifacts endpoints returned empty arrays. The combined-status endpoint returned no status entries.

This is not a pass, a fail, or evidence that the tests executed. No corrected smoke artifact is available from that run.

A local clone was also attempted, but the runtime could not resolve `github.com`; consequently no repository tests were executed locally in this session.

## Independent deterministic calculation

Using the declared parent intensity
`g(z)=1+2.2\exp[-(z-0.105)^2/(2(0.025)^2)]+1.6\exp[-(z-0.275)^2/(2(0.040)^2)]`
and
`b(z)=1+1.5(z-0.01)/(0.40-0.01)`
on 200,001 equally spaced points over [0.01, 0.40], I independently normalized the two target densities:

- Intended: proportional to (g(z)b(z)).
- Duplicated: proportional to (g(z)b(z)^2).

Results:
- L1 density distance: 0.1846197162.
- Total-variation distance: 0.0923098581.
- Mean redshift, intended target: 0.2251874658.
- Mean redshift, duplicated target: 0.2469060725.
- Mean-redshift shift: +0.0217186068.

These are deterministic integrations of the declared functions, not Monte Carlo results.

## Decision

**NO-GO for dispatching the seven-gamma workflow at the currently inspected curvature-base revision.**

Required sequence:
1. Obtain an executable exact-head CI result for the correction.
2. Run the corrected one-replicate smoke on the correction head and independently audit the raw artifact.
3. Integrate the validated correction into the curvature research branch without changing main.
4. Verify the multi-replicate runner/auditor uses the corrected child implementation and run a small multi-replicate smoke with independent outer audit.
5. Only then consider the seven-gamma density-only workflow.

No code was changed in this audit, no workflow was dispatched, no PR was merged, and no H0 inference was run.
