# Consolidated Conversation Checkpoint Index — 2026-11-26

This directory preserves the latest dark-siren project status and agent instructions from the ChatGPT conversation on 2026-11-26.

## Documents

1. [Latest v2 project status and evidence register](latest_v2_project_status_and_evidence_register.md)
   - Latest user-reported /workspace/v2 state.
   - Pilot design and 1,500-cell results.
   - IPW stability and misspecification findings.
   - Reported 475-test suite, v2 tests, scripts and failure resolutions.
   - D3 production-integration report claims and FPC discrepancy.
   - Oracle undercoverage, selection, quadrature, publication-readiness and scientific limitations.
   - Evidence labels and no-go decision.

2. [Agent work orders and coordination](agent_work_orders_and_coordination.md)
   - OpenHands independent audit.
   - OpenCode mathematical diagnosis and roadmap.
   - Consensus literature and novelty review.
   - Coordination rules and campaign restrictions.

3. Earlier conversation checkpoint:
   - `research/project_completion/conversation_checkpoint_2026-11-26.md`
   - Contains the prior D3/Arm E mathematical reconciliation, DGP/estimands, finite-frame versus superpopulation variance derivation, literature leads, and earlier implementation checkpoint.

## What this checkpoint does and does not save

Saved remotely on branch `research/conversation-checkpoint-20261126`:
- This documentation index.
- The latest project status/evidence register.
- The three-agent work orders and coordination protocol.
- The earlier conversation checkpoint already on this branch.

Not saved by these documentation commits:
- The local `/workspace/v2` source tree.
- Uncommitted source changes or local Git objects.
- The 1,500-cell pilot raw outputs, JSON artifacts, logs, or ledger.
- The 230 MB CI artifact.
- Local virtual environment or dependency state.
- Any unpushed OpenCode/OpenHands worktree content.

Remote verification in the session:
- Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`.
- Default branch: `main`.
- The reported v2 branch was not returned by GitHub branch search.
- The reported HEAD SHA `e06b88eb87fc4851105f34f9680de2f9deece7e6` exists remotely, but the inspected commit is a revert deleting `scripts/run_degree_matched_five_arm_quadrature_sensitivity.py`; it does not prove the reported local v2 work is present.
- The documentation branch is separate; main was not modified by these checkpoint writes.

## Current decision

**NO-GO for confirmatory campaign.**

Resolve the validation-frame/FPC denominator discrepancy and reproduce/diagnose oracle undercoverage before additional campaign-scale simulations. Do not claim posterior calibration, estimator superiority, novelty, observational inference, or publication readiness.

## Evidence convention

All latest v2 test, pilot, D3 integration, coverage, selection, and quadrature values are marked USER-REPORTED in the status register until the actual local source and artifacts are uploaded/pushed or otherwise made accessible and independently inspected.
