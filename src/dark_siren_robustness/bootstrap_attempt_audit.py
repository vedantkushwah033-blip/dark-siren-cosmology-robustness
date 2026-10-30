"""Canonical accounting helpers for inner bootstrap attempts.

This module validates and serializes attempt-level diagnostics returned by
studentized finite-population bootstrap routines. It does not independently
audit the draws; the artifact auditor must replay those from their seed.
"""
from __future__ import annotations

import json
from collections import Counter
from typing import Any


def canonical_inner_attempt_record(result: dict[str, Any]) -> dict[str, object]:
    """Return validated, JSON-stable inner-attempt accounting.

    The counts are checked against both the result counters and the per-draw
    success/failure arrays. Failure-code keys are sorted for deterministic
    serialization. A successful attempt must have an empty failure code and
    a failed attempt must have a non-empty failure code.
    """
    required = ("n_attempted", "n_successful", "n_failed", "success_mask", "failure_codes")
    missing = [key for key in required if key not in result]
    if missing:
        raise ValueError(f"bootstrap result missing attempt fields: {missing}")

    attempted = result["n_attempted"]
    successful = result["n_successful"]
    failed = result["n_failed"]
    for name, value in (("n_attempted", attempted), ("n_successful", successful), ("n_failed", failed)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a non-negative integer")
    if attempted != successful + failed:
        raise ValueError("inner attempt counts do not reconcile")

    success_mask = list(result["success_mask"])
    failure_codes = [str(code) for code in result["failure_codes"]]
    if len(success_mask) != attempted or len(failure_codes) != attempted:
        raise ValueError("inner per-attempt arrays differ from n_attempted")

    observed_success = sum(bool(value) for value in success_mask)
    if observed_success != successful:
        raise ValueError("success_mask does not match n_successful")

    failures = []
    for ok, code in zip(success_mask, failure_codes):
        if bool(ok) and code:
            raise ValueError("successful inner attempt has a failure code")
        if not bool(ok):
            if not code:
                raise ValueError("failed inner attempt lacks a failure code")
            failures.append(code)

    counts = dict(sorted(Counter(failures).items()))
    if sum(counts.values()) != failed:
        raise ValueError("failure-code histogram does not match n_failed")

    return {
        "inner_attempts": attempted,
        "inner_successes": successful,
        "inner_failures": failed,
        "inner_failure_code_counts": json.dumps(counts, sort_keys=True, separators=(",", ":")),
    }
