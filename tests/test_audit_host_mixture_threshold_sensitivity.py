import csv, json
from pathlib import Path
import pytest
from scripts.audit_host_mixture_threshold_sensitivity import audit


def test_audit_rejects_missing_artifact(tmp_path):
    with pytest.raises(FileNotFoundError):
        audit(tmp_path, grid_sensitivity=False)
