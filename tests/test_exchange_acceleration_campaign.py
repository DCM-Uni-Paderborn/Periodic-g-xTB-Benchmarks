import json
import hashlib
import subprocess
import sys
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


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def test_cp2k_block_helper_status_is_scoped_conservatively():
    matrix = load_json("validation_matrix.json")
    modules = {entry["id"]: entry for entry in matrix["modules"]}

    assert modules["cp2k_block_expansion_foldback_helpers"]["status"] == "passed"
    assert modules["streamed_symmetry_stars"]["status"] == "implementation_in_progress"
    assert "full-array oracle" in modules["cp2k_block_expansion_foldback_helpers"][
        "current_evidence"
    ]


def test_cp2k_block_helper_raw_archive_and_manifests():
    provenance = load_json("provenance/cp2k_block_helpers.json")
    raw = CAMPAIGN / "raw" / "cp2k_block_helpers"
    run_dirs = sorted(path for root in raw.iterdir() for path in root.iterdir() if path.is_dir())

    assert len(run_dirs) == 10
    assert sum((path / "returncode.txt").read_text().strip() == "0" for path in run_dirs) == 10
    assert sum((path / "SHA256SUMS.initial").is_file() for path in run_dirs) == 10
    assert sum((path / "SHA256SUMS.final").is_file() for path in run_dirs) == 10
    assert provenance["runtime"]["zero_returncodes"] == 10

    for run in run_dirs:
        output = (run / "cp2k.out").read_text(errors="replace")
        assert "PROGRAM ENDED AT" in output
        final_manifest = (run / "SHA256SUMS.final").read_text().splitlines()
        for line in final_manifest:
            expected, relative = line.split(maxsplit=1)
            assert digest(run / relative) == expected

        initial_manifest = (run / "SHA256SUMS.initial").read_text().splitlines()
        input_hash, input_path = initial_manifest[0].split(maxsplit=1)
        binary_hash, binary_path = initial_manifest[1].split(maxsplit=1)
        assert Path(input_path).name == "input.inp"
        assert digest(run / "input.inp") == input_hash
        expected_binary = (
            provenance["cp2k"]["baseline_binary_sha256"]
            if "baseline" in run.parent.name
            else provenance["cp2k"]["current_binary_sha256"]
        )
        assert Path(binary_path).name == "cp2k.psmp"
        assert binary_hash == expected_binary


def test_cp2k_block_helper_comparison_is_reproducible():
    script = CAMPAIGN / "scripts" / "compare_cp2k_block_helpers.py"
    derived_json = CAMPAIGN / "derived" / "cp2k_block_helpers_comparison.json"
    derived_text = CAMPAIGN / "derived" / "cp2k_block_helpers_comparison.txt"

    generated_json = subprocess.run(
        [sys.executable, str(script), "--format", "json"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    generated_text = subprocess.run(
        [sys.executable, str(script), "--format", "text"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    assert generated_json == derived_json.read_text()
    assert generated_text == derived_text.read_text()
    result = json.loads(generated_json)
    assert result["all_passed"] is True
    assert result["maxima"]["energy_sequence_hartree_max_abs"] < 1e-13
    assert result["maxima"]["atomic_force_hartree_per_bohr_max_abs"] < 1e-14
    assert "production streamed exchange module is not qualified" in result["scope_note"]


def test_cp2k_block_helper_provenance_is_self_verifying():
    provenance_dir = CAMPAIGN / "provenance" / "cp2k_block_helpers"
    manifest = provenance_dir / "SHA256SUMS"
    entries = manifest.read_text().splitlines()

    assert len(entries) == 90
    for line in entries:
        expected, relative = line.split(maxsplit=1)
        assert digest(CAMPAIGN / relative) == expected

    patch = (provenance_dir / "cp2k_worktree.patch").read_text()
    assert "tb_expand_gxtb_kpoint_block" in patch
    assert "tb_accumulate_gxtb_kpoint_response" in patch
    assert "tb_validate_gxtb_kpoint_block_helpers" in patch
    assert "Blockwise CP2K symmetry fold disagrees with the full-array oracle" in patch

    first_build = (provenance_dir / "terok_build_and_launch" / "build_helpers.log").read_text()
    retry = (provenance_dir / "terok_build_and_launch" / "build_helpers_retry.log").read_text()
    assert "Error: Line truncated" in first_build
    assert "[ 99%] Built target cp2k-bin" in retry
