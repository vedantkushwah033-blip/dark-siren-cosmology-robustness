# External AI collaboration protocol and task allocation — 29 October 2026

## Purpose

Coordinate independent assistance from Perplexity, Claude and Gemini on the Dark-Siren Cosmology Robustness project without duplicating work, conflicting edits, or inflating claims.

Repository: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness

Current scope: a controlled synthetic study of catalogue incompleteness, host-density estimation, selection-function misspecification, numerical integration and posterior calibration. It is not a real LVK/GWTC H0 measurement.

## Rules for all collaborators

1. Start by inspecting the current repository state, relevant branches, open PRs, README, manuscript, audit reports, raw artifacts, tests and workflows. Do not assume earlier status summaries are current.
2. Do not merge PRs, retarget PRs, force-push, rewrite history, delete data or change protected/default branches.
3. Do not commit to the same working branch as another collaborator concurrently. Work in separate branches or return a patch/diff and exact file/line recommendations.
4. Do not fabricate tool execution, test outcomes, papers, DOIs, data access, artifact hashes or review status. Clearly separate verified facts, inferred risks, and proposed work.
5. Do not run huge simulation campaigns before inspecting the existing protocol, computing cost, estimands, seed allocation, power/Monte Carlo precision and available artifacts.
6. Every new result must be reproducible: code/configuration, seeds, versions, raw replicate records, immutable artifact/hash, independent audit, uncertainty intervals and documented failure cases.
7. Treat real-data/CHIMERA v2.2 validation as NO-GO until compatible version-pinned inputs and a complete auditable run are actually available.
8. The novelty boundary is narrow. Catalogue incompleteness, GP reconstruction, clustering-based completion and selection-aware inference are established areas. No “first ever,” Hubble-tension resolution, universal weighting superiority, real H0 measurement or observational validation claims.
9. Prefer minimal, high-value changes. Report blockers explicitly. Do not merely produce a plan if a concrete read-only audit or isolated patch can be completed safely.
10. End with changed files/commit or exact patch, commands actually run, test outputs, unresolved risks and the single highest-value next action.

## Task allocation

### Perplexity — literature and external-evidence lead
- Task P1: source-backed novelty and assumptions matrix. Search primary papers, publisher pages, supplements, official code/docs and data documentation; compare methods by data regime, selection likelihood, galaxy-field treatment, assumptions and evaluation. Verify bibliographic details; do not infer novelty from abstracts alone.
- Task P2: external benchmark feasibility dossier. Inspect official CHIMERA v2.2 release/docs/tutorials and required data provenance. Produce an exact dependency/input checklist, compatible version/commit, minimal reproducible command sequence and explicit GO/NO-GO decision. Do not silently substitute old versions or claim execution without logs.

### Claude — code and software-quality lead
- Task C1: read-only adversarial code audit of likelihood conditioning, selection normalization, prior support/boundaries, numerical integration, KDE bandwidth convention, seed reproducibility, error handling and scientific claim-to-code traceability. Give file/line references and minimal fixes; independently verify with tests where possible.
- Task C2: propose and, if tool access permits, implement one isolated, high-value software improvement on a separate branch: e.g. a release-manifest validator that maps manuscript tables to immutable artifacts and verifies hashes/schema/replicate counts. Include tests for corrupted/missing artifacts and avoid rewriting scientific outputs.

### Gemini — independent statistical design and red-team lead
- Task G1: independently reconstruct the estimands and audit the statistical interpretation: paired contrasts, Monte Carlo uncertainty, Wilson coverage intervals, bootstrap design, overlapping seeds, boundary mass and negative controls. Identify any pseudo-replication, conditioning errors or claims unsupported by the design.
- Task G2: design one bounded next-stage synthetic experiment that adds meaningful insight without scope creep. State the scientific question, data-generating process, estimand, factorial cells, replicate count justified by precision/power, seed scheme, preregistered pass/fail gates, compute estimate, artifact format and stopping rule. Do not run it until the design is reviewed against existing experiments; avoid repeating existing cells.

## Integration policy

Do not blindly combine outputs. The project owner should reconcile disagreements using code, raw artifacts and primary sources. Prioritize:
1. critical statistical/code defects;
2. manuscript-artifact mismatch;
3. reproducibility and release gates;
4. external benchmark readiness;
5. only then a new experiment.

A collaborator's recommendation is not validation until independently checked. Keep all external collaborators' work in separate branches/patches until reviewed.

## Required final report template

- Exact repository/branch/commit reviewed:
- Task completed:
- Files changed (or no changes):
- Evidence inspected:
- Commands actually run and exact results:
- Findings ranked Critical/High/Medium/Low:
- Reproducibility and uncertainty limitations:
- Conflicts with other collaborator work:
- Remaining blockers:
- Recommended next action:
