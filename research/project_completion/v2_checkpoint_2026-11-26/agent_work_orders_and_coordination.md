# Agent Work Orders and Coordination — 2026-11-26

This document preserves the current agent assignments drafted in the ChatGPT conversation. These are instructions for the user to provide to external agents; ChatGPT has not launched OpenHands, OpenCode, or Consensus.

## Coordination protocol

1. OpenHands performs the independent repository/scientific audit first.
2. OpenCode does not edit the shared worktree while OpenHands is auditing it.
3. After the OpenHands report is available, OpenCode reads it, independently checks critical conclusions, and produces the mathematical diagnosis/implementation roadmap.
4. Consensus performs an independent literature and novelty review; it does not validate code.
5. Compare disagreements across all three reports against source code, equations, raw artifacts, and reproducible commands.
6. No implementation integration or new campaign until an explicit scientific design approval.
7. Do not let two agents edit the same worktree simultaneously.

## Work order A — OpenHands: independent scientific and repository audit

### Starting state (reported, verify before use)
- Worktree: /workspace/v2
- Branch: research/v2-finite-frame-arm-e-2026-11-25
- HEAD: e06b88eb87fc4851105f34f9680de2f9deece7e6

### Mission
First verify the repository, branch, commit and working-tree status. Do not assume the reported identifiers are correct.

**Priority 1: D3 sampling-frame/FPC discrepancy**
- Inspect frozen config, parent generation, validation-frame generation, validation-sample selection, validation exclusion, D3 covariance implementation, pilot artifacts and ledger.
- Determine whether actual SRSWOR frame size is 1,500 or 4,000 (or whether those numbers refer to different populations/experiments).
- Derive the FPC denominator from actual design. Do not change the frozen design to match a report.
- Check that D3 covariance and downstream propagation use the same design.

**Priority 2: production integration**
Trace executable path:
beta estimate -> inclusion probabilities -> normalized IPW weights -> weighted redshift density -> KDE -> event likelihood -> selection normalization -> H0 posterior.

For every stage record source file/function, inputs/outputs, whether uncertainty is propagated or merely computed, whether pilot executes it, and relevant tests/artifacts. Distinguish implementation, unit tests, numerical validation, independent reconstruction and scientific validation.

**Priority 3: oracle undercoverage**
- Reproduce the reported true-pi/oracle nominal 68% coverage of roughly 0.60–0.65 from saved artifacts.
- Investigate interval construction, likelihood, prior/boundaries, event-generation model, host density/KDE, selection normalization, quadrature, Monte Carlo precision, and whether truths/data are generated consistently with the inferential model.
- Separate demonstrated causes from hypotheses. Do not blame IPW by default.

**Priority 4: provenance**
- Verify 1,500 pilot cells against config and artifacts.
- Verify ledger completeness and status classification.
- Reproduce the reported 475-test run if practical; inspect the two ignored collection-error modules and verify their provenance.
- Verify the three-seed Route A/Route B comparison and 3.5% / 2.8% delta-method-vs-MC agreements.
- Identify claims supported only by reports rather than reproducible artifacts.

### Required OpenHands deliverables
- repository_state.md
- fpc_sampling_design_audit.md
- d3_production_trace.md
- oracle_undercoverage_audit.md
- pilot_and_test_provenance.md
- claim_evidence_matrix.md
- discrepancy_register.md
- go_no_go_recommendation.md

Final recommendation must be GO / CONDITIONAL GO / NO-GO and list conditions.

### OpenHands restrictions
- No main modifications.
- No reset, clean, delete, overwrite, commit, push, or merge.
- No changes to preregistered configurations.
- No confirmatory simulations or long campaign.
- Preserve existing work and artifacts.
- Report exact commands and actual outputs.
- If artifact/repository inaccessible, report it rather than guessing.

## Work order B — OpenCode: mathematical diagnosis and implementation roadmap

### Starting state (reported, verify before use)
- Worktree: /workspace/v2
- Branch: research/v2-finite-frame-arm-e-2026-11-25
- HEAD: e06b88eb87fc4851105f34f9680de2f9deece7e6

Do not edit the shared worktree while OpenHands is auditing. After the audit is delivered, read it and independently verify critical conclusions.

### Mission A: D3 sampling design
Derive the actual fixed-frame SRSWOR variance:
Var_d(beta_hat | F_M) approximately A_M^(-1) [(1-n/M) S_psi,M/n] A_M^(-T).

Establish what F_M is, whether validation records are sampled from the parent or a separate frame, whether validation exclusion changes only the analysis catalogue, whether finite-frame score covariance is correct, and whether the target is conditional design variance or unconditional superpopulation variance. Do not conflate these targets.

### Mission B: oracle undercoverage fault tree
Investigate:
1. Likelihood misspecification.
2. Posterior interval construction/asymmetry.
3. Prior effects and H0 boundaries.
4. Event-generation assumptions.
5. Host-population density and KDE approximation.
6. GW detection selection and normalization.
7. Numerical integration/quadrature.
8. Monte Carlo uncertainty and coverage-estimator precision.
9. Differences between DGP and inference model.

For each candidate, state mechanism, required code/artifact, evidence for/against, and smallest discriminating diagnostic. No unsupported causal conclusion.

### Mission C: calibration validation strategy
Design SBC correctly, but do not implement until DGP and inferential model are reconciled:
- Draw H0 and nuisance parameters from the actual prior.
- Generate validation data, catalogue inclusion, GW events and detection consistently conditional on truths.
- Compute posterior ranks.
- Handle discrete/numerically approximated posterior CDFs.
- Set replicate count and report MC uncertainty.
- Distinguish SBC from frequentist coverage at fixed true H0.

### Mission D: D3 and Arm E
Keep separate unless a valid decomposition is derived:
- D3 = fixed-frame frequentist design uncertainty.
- Arm E = Bayesian uncertainty conditional on validation observations/model/prior.
- Do not add covariance matrices or inject D3 into Arm E by default.
- Verify shared beta across all events and beta-dependent selection normalization inside the integral.

### Mission E: prioritized roadmap
After OpenHands audit, prioritize:
1. FPC and sampling-frame reconciliation.
2. Oracle undercoverage diagnosis.
3. SBC only after DGP/inference reconciliation.
4. Detection-limited selection validation.
5. Arm E R7 prior sensitivity and remaining certification.
6. AIPW only after calibration diagnosis.
7. Separately preregister stress designs for separation/extreme weights.

For each change specify files/functions, rationale, tests, acceptance criteria, computational cost, and whether preregistration changes.

### Required OpenCode deliverables
- research/next_phase_coverage_audit/mathematical_fault_tree.md
- research/next_phase_coverage_audit/sbc_design.md
- research/next_phase_coverage_audit/d3_arm_e_reconciliation.md
- research/next_phase_coverage_audit/implementation_roadmap.md
- research/next_phase_coverage_audit/open_questions.md

No paper draft or novelty claim. End with the single most important scientific uncertainty to resolve before any campaign.

### OpenCode restrictions
- No main modifications or destructive Git operations.
- No reset, clean, delete, overwrite, commit, push, or merge.
- No confirmatory/large simulation campaign.
- No silent frozen-config changes.
- Do not weaken tests.
- Do not claim unrun tests passed.
- Record exact branch, commit, modified files, commands and results.

## Work order C — Consensus: literature and novelty audit

Conduct independent critical literature review. Do not claim code validation.

### Core question
The synthetic pilot reports that true inclusion probabilities still yield only about 60–65% coverage for nominal 68% posterior intervals. Estimated inclusion probabilities add residual bias under logistic misspecification, but cannot explain the oracle's coverage deficit by themselves.

### Review topics
1. Bayesian credible-set frequentist coverage, posterior calibration under misspecification, SBC, posterior contraction and hierarchical population inference.
2. Dark-siren catalogue incompleteness and reconstruction: Gair et al.; Borghi et al. 2026; Hanselman et al. 2025; Tagliazucchi et al. 2026; Leyde et al.; other direct work.
3. Candidate undercoverage mechanisms: likelihood misspecification, interval construction, prior/boundaries, finite-event uncertainty, selection normalization, KDE/density reconstruction, quadrature, DGP/inference mismatch.
4. SBC methodology: prior-drawn truths, conditional data generation, posterior ranks, numerical CDFs, rank histograms, Monte Carlo precision, distinction from fixed-truth frequentist coverage.
5. Novelty: whether prior work combines finite validation samples, misspecified inclusion models, shared inclusion nuisance across GW events, parameter-dependent GW selection normalization, downstream H0 calibration/coverage, and pre-inference diagnostics.

### Required outputs
- oracle_undercoverage_literature.md
- dark_siren_calibration_prior_art.md
- sbc_methodology_recommendations.md
- novelty_risk_assessment.md
- annotated_bibliography.md

For important claims, provide full citation, DOI/arXiv/publisher link, exact section/equation/result, and whether full text or abstract only was inspected. Use primary papers, include contradictory work, and identify prior art that weakens novelty.

### Consensus restrictions
- No invented papers/equations/results/citations.
- Do not equate SBC with frequentist coverage.
- Oracle undercoverage does not prove a specific likelihood defect.
- No “first” claim based on search absence.
- Mark inaccessible full texts as not fully verified.

Conclude with three plausible explanations for oracle undercoverage, smallest discriminating experiments, and strongest defensible research question after prior-art review.

## Current campaign gate
NO-GO for confirmatory/large simulation campaign until:
- the actual validation frame and FPC denominator are resolved;
- oracle undercoverage is reproduced and diagnosed to a defensible extent;
- the generative and inferential models are reconciled;
- independent audit findings are reviewed.

The three agents have complementary responsibilities; none of their reports alone authorizes production integration or a confirmatory campaign.
