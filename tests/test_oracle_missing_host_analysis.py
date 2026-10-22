import importlib.util
from pathlib import Path


def load_script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_oracle_analysis_declares_four_prespecified_pairs():
    module = load_script("analyze_oracle_missing_host.py")
    assert len(module.PAIRS) == 4
    assert all(len(pair) == 3 for pair in module.PAIRS)
    assert all("oracle" in pair[0] or "observed_only" in pair[0] or "inference_weights" in pair[0]
               for pair in module.PAIRS)


def test_oracle_analysis_uses_same_validated_bootstrap_helper():
    module = load_script("analyze_oracle_missing_host.py")
    assert callable(module.analysis.bootstrap_mean_difference)
    result = module.analysis.bootstrap_mean_difference([0.0, 0.0], [1.0, 1.0],
                                                        __import__("random").Random(7), 200)
    assert result["difference_B_minus_A"] == 1.0
