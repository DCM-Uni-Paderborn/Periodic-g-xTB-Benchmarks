import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "campaigns" / "exchange-acceleration-20260716"


def load_json(name: str):
    return json.loads((CAMPAIGN / name).read_text())


def test_validation_matrix_is_self_consistent():
    matrix = load_json("validation_matrix.json")
    allowed = set(matrix["status_vocabulary"])
    modules = {entry["id"]: entry for entry in matrix["modules"]}

    assert matrix["oracle"] == "explicit-expanded-full-mesh"
    assert len(modules) == len(matrix["modules"])
    assert all(entry["status"] in allowed for entry in modules.values())

    for case in matrix["cases"]:
        assert case["status"] in allowed
        assert case["modules"]
        assert set(case["modules"]) <= set(modules)


def test_cache_runtime_archive_matches_provenance_counts():
    provenance = load_json("provenance/cache_runtime_r3.json")
    runs = CAMPAIGN / "raw" / "cache_runtime_r3" / "runs"
    run_dirs = sorted(path for path in runs.iterdir() if path.is_dir())
    returncodes = [(path / "returncode.txt").read_text().strip() for path in run_dirs]

    assert len(run_dirs) == provenance["launcher_count"] == 19
    assert returncodes.count("0") == provenance["returncodes"]["zero"] == 18
    assert len(returncodes) - returncodes.count("0") == provenance["returncodes"]["nonzero"] == 1
    assert sum((path / "SHA256SUMS.initial").is_file() for path in run_dirs) == 19
    assert sum((path / "SHA256SUMS.final").is_file() for path in run_dirs) == 19


def test_paper_evidence_has_reproducible_derivation():
    assert (CAMPAIGN / "scripts" / "compare_cache_runtime.py").is_file()
    comparison = CAMPAIGN / "derived" / "cache_runtime_r3_comparison.txt"
    text = comparison.read_text()
    assert "BASELINE_PHYSICAL_COMPONENTS" in text
    assert "SYMMETRY_PHYSICAL_COMPONENTS" in text
    assert "BASELINE_SYMMETRY_METADATA" in text
