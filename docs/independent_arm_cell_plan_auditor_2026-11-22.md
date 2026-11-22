# Independent aggregate cell-plan auditor

Added 2026-11-22.

`scripts/audit_arm_cell_plan.py` is a read-only, standard-library auditor. It
does not import the checkpoint store, lifecycle ledger, coordinator, or cell
evaluator. It independently checks the canonical immutable-plan hash, exact
plan/ledger key coverage, planned provenance, four-state reason/audit rules,
timestamp ordering, and each SUCCESS/NON_EVALUABLE row's referenced audit
report and output directory.

The referenced cell report must say PASS for the exact cell key. The auditor
distinguishes an incomplete plan from a completed plan containing FAILED or
NON_EVALUABLE cells. It never edits or repairs a checkpoint.

This is an aggregate structural audit, not a replacement for
`audit_arm_cell_from_bundle.py`: numerical reconstruction remains the job of
the per-cell independent auditor. It also does not define scientific
NON_EVALUABLE thresholds; those must be predeclared and justified separately.
