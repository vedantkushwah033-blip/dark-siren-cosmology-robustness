# Multi-frame benchmark implementation contract (v3 architecture)

**Status:** engineering specification; not an executed experiment or preregistration amendment.
**Scope:** linear finite-population arithmetic-mean benchmark under SRSWOR.
**Release rule:** do not describe a run as confirmatory until the frozen design, runner, auditor, and manuscript all implement this contract and all release gates pass.

## 1. Why v2 cannot simply be scaled up

The current v2 pilot schema represents one generated finite frame per scenario. Increasing `scenarios` repeats distribution labels cyclically; it does not create a proper independently replicated frame hierarchy within each scenario. Therefore, raising the scenario count or repetitions is not a valid substitute for a multi-frame design.

A genuine multi-frame implementation must add an explicit `frame_id` dimension and independently seed each frame within each scenario. The auditor must reconstruct that hierarchy and reject missing, duplicated, or relabelled frames.

## 2. Locked hierarchy and identifiers

The confirmatory key is:

`(scenario_id, frame_id, sample_size, outer_replicate, method)`

The frame table key is:

`(scenario_id, frame_id, population_unit)`

The inner-attempt table (or canonical serialized inner summary attached to each studentized outer attempt) is keyed by:

`(scenario_id, frame_id, sample_size, outer_replicate, method, inner_attempt)`

Identifiers are zero-based integers or stable strings with one documented canonical representation. No key may be inferred from row order.

## 3. Seed lineage

Use a root `numpy.random.SeedSequence(root_seed)` with explicitly documented child spawning:

1. Spawn one child per scenario in frozen scenario order.
2. For each scenario child, spawn one child per frame in increasing `frame_id`.
3. Use each frame child only to generate that frame's finite population.
4. Derive outer sample seeds deterministically from the frame child and the prespecified sample-size / replicate ordering, or use a documented independent child tree.
5. Derive method-specific bootstrap seeds from a collision-resistant tuple of root seed, scenario, frame, sample size, outer replicate, and method-specific domain-separation constant.

The exact derivation algorithm must be shared in a small pure function and independently reimplemented in the auditor. Record the algorithm identifier and NumPy version in the manifest. Never rely on incidental iteration order or global RNG state.

## 4. Paired outer samples

For each frame, sample size, and outer replicate, generate one SRSWOR index vector and reuse that exact sample for every method. This is essential for paired method contrasts.

Store the sample indices once in an outer-sample table, not redundantly in each method row. The auditor reconstructs indices from the seed and checks uniqueness, range, sample size, and equality across all method records linked to that outer sample.

Method-specific random streams must be independent of method execution order. Reordering the method list must not change any method's output.

## 5. Scenario manifest

Freeze an explicit ordered list of scenarios, each containing:
- scenario ID and human-readable name;
- distribution family and all parameters;
- population size;
- number of independent frames;
- sample sizes / sampling fractions;
- outer repetitions per frame and cell;
- nominal interval level;
- bootstrap replicate count;
- methods and primary reference method;
- frame and outer seed derivation version;
- primary estimand, contrast sign convention, and multiplicity family.

Do not infer scenario definitions from cyclic indexing. Reject unknown or omitted scenario IDs. If finite-frame shape strata are added, define their construction before outcomes and retain the unstratified generating-regime analysis as the primary target.

## 6. Required summaries

Retain frame-level summaries before aggregation. For every scenario × sample-size × method cell, report:
- number of generated and valid frames;
- per-frame attempted, successful, failed, covered and successful-and-covered counts;
- conditional coverage (undefined if no successes);
- interval-success probability;
- unconditional successful-and-covering probability;
- interval-width distribution among successful intervals;
- canonical failure categories;
- paired frame-level method contrasts;
- equal-frame-weighted mean, standard deviation, median, quantiles and range.

Do not pool all outer replicates into a single binomial denominator for across-frame inference. The frame is the independent unit for the superpopulation generalization target.

## 7. Primary contrast and uncertainty

Before confirmatory execution, select exactly one primary contrast and one primary estimand. A defensible default is the equal-frame-weighted mean paired difference in unconditional successful-and-covering probability between studentized pseudo-population and Wald SRSWOR intervals, with sign defined as studentized minus Wald.

Report the paired frame-level contrast distribution and a prespecified frame-level confidence interval. If using a frame bootstrap, resample whole frames within scenario; never resample individual outer attempts as though they were independent frames. Validate interval coverage of the frame-level procedure in a separate design simulation where feasible.

Conditional coverage, interval success, widths, RMSE, and scenario interactions remain secondary unless explicitly promoted before the design lock. Apply the preregistered multiplicity procedure to the declared secondary family.

This default is a recommendation for the protocol, not an assertion that the estimand is scientifically optimal. The research team must justify it before lock.

## 8. Independent auditor requirements

The auditor must not import the runner's private calculations for scientific reconstruction. It must independently:
- regenerate every finite frame from its seed and scenario specification;
- verify frame hashes, finite-frame means and population diagnostics;
- reconstruct every outer SRSWOR sample and ensure method pairing;
- independently recompute all interval endpoints and failures;
- replay inner pseudo-population attempts and failure categories;
- verify complete expected key sets and reject duplicates / omissions;
- reconstruct frame-level and across-frame summaries;
- reproduce the primary paired contrast and its uncertainty;
- check the frozen manifest against actual artifact dimensions and code-version metadata.

Add adversarial tests that alter one frame seed, one frame value, one outer sample index, one method endpoint, one inner failure code, one frame-level summary, and one primary contrast while refreshing file checksums. Each must be rejected for a scientific reason.

## 9. Staged implementation and stop/go gates

**Stage A — schema and seed primitives:** implement pure seed derivation, manifest validation, and v3 schemas; unit-test collision resistance, deterministic replay, and method-order invariance.

**Stage B — small multi-frame pilot:** run at least a few frames per scenario with tiny repetitions; audit every level; label `PILOT_NOT_CONFIRMATORY`. This tests engineering, not scientific performance.

**Stage C — planning simulation:** use separate planning frames to estimate frame-level variance and computational cost. Evaluate planner sensitivity and select frame count / within-frame repetitions. Do not reuse planning frames in confirmatory inference.

**Stage D — design lock:** commit immutable manifest, analysis plan, seeds, primary contrast, multiplicity family, and deviation rules. Record a hash of the locked design.

**Stage E — confirmatory run:** execute only the frozen manifest. No optional stopping or post-hoc frame exclusion.

**Stage F — independent release audit:** audit the complete artifacts, reconstruct results, and reconcile every manuscript number against machine-readable output.

If any stage fails, stop progression to the next stage; fix the defect and rerun the affected stage. A passing v2 pilot does not waive any v3 gate.

## 10. Explicit non-claims

This architecture is not a new bootstrap theorem, does not guarantee universal coverage, and does not establish superiority over Wald intervals. It does not transfer to nonlinear dark-siren posterior inference or constitute a cosmological measurement. Its contribution, if supported by results, will be a carefully controlled and auditable empirical comparison under the explicitly defined finite-population designs.
