"""Disjoint finite-validation split for synthetic catalogue-inclusion studies.

Validation objects are sampled uniformly from the full parent frame before
conditioning on the retention label. The returned analysis catalogue excludes
every validation object, preventing validation-label reuse in the primary arm.
Hidden true inclusion probabilities are deliberately not returned.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ValidationSplit:
    validation: pd.DataFrame
    analysis_catalogue: pd.DataFrame
    validation_object_ids: np.ndarray
    analysis_object_ids: np.ndarray
    n_parent: int
    n_validation: int


def split_disjoint_validation(
    parent_frame: pd.DataFrame,
    *,
    n_validation: int,
    seed: int,
) -> ValidationSplit:
    """Draw a simple random validation sample from all parent objects.

    Required input columns are object_id, redshift, and retained. Any other
    columns (including true inclusion probabilities) are ignored. Validation
    selection is independent of retained outcomes by construction.
    """
    required = {"object_id", "redshift", "retained"}
    missing = required.difference(parent_frame.columns)
    if missing:
        raise ValueError(f"parent_frame is missing required columns: {sorted(missing)}")
    if parent_frame.empty or parent_frame["object_id"].isna().any():
        raise ValueError("parent_frame must be nonempty with non-null object IDs")
    if parent_frame["object_id"].duplicated().any():
        raise ValueError("object_id values must be unique")
    z = parent_frame["redshift"].to_numpy(dtype=float)
    y = parent_frame["retained"].to_numpy()
    if not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("redshift must be finite and positive")
    if np.any((y != 0) & (y != 1)):
        raise ValueError("retained must contain only binary 0/1 labels")
    if not isinstance(n_validation, int) or isinstance(n_validation, bool):
        raise ValueError("n_validation must be an integer")
    if not 5 <= n_validation < len(parent_frame):
        raise ValueError("n_validation must be at least 5 and smaller than parent size")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    rng = np.random.default_rng(seed)
    positions = rng.choice(len(parent_frame), size=n_validation, replace=False)
    validation_mask = np.zeros(len(parent_frame), dtype=bool)
    validation_mask[positions] = True
    base = parent_frame.loc[:, ["object_id", "redshift", "retained"]].reset_index(drop=True)
    validation = base.loc[validation_mask].reset_index(drop=True)
    analysis = base.loc[(~validation_mask) & (base["retained"].to_numpy() == 1),
                        ["object_id", "redshift"]].reset_index(drop=True)
    validation_ids = validation["object_id"].to_numpy(copy=True)
    analysis_ids = analysis["object_id"].to_numpy(copy=True)
    if np.intersect1d(validation_ids, analysis_ids).size:
        raise ArithmeticError("Validation and analysis object IDs overlap")
    return ValidationSplit(
        validation=validation,
        analysis_catalogue=analysis,
        validation_object_ids=validation_ids,
        analysis_object_ids=analysis_ids,
        n_parent=int(len(parent_frame)),
        n_validation=int(n_validation),
    )
