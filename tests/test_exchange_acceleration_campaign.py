import json
import hashlib
import subprocess
import sys
import unittest
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

    helpers = modules["cp2k_block_expansion_foldback_helpers"]
    assert helpers["status"] == "passed"
    assert modules["cp2k_streamed_star_consumer"]["status"] == "passed"
    assert {
        "physical_overlap_expansion",
        "weighted_real_adjoint",
        "time_reversal",
        "K290",
        "SPGLIB",
        "energy",
        "forces",
        "stress",
    } <= set(helpers["required_observables"])
    assert "full-array oracle" in helpers["current_evidence"]


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


def test_provider_forward_status_is_scoped_conservatively():
    matrix = load_json("validation_matrix.json")
    modules = {entry["id"]: entry for entry in matrix["modules"]}

    assert modules["regular_grid_bvk_cache"]["status"] == "passed"
    assert modules["provider_matrix_lean_forward_stream"]["status"] == "passed"
    assert modules["provider_bounded_ao_image_batching"]["status"] == "passed"
    assert modules["provider_streamed_reverse"]["status"] == "passed"
    assert modules["cp2k_streamed_star_consumer"]["status"] == "passed"
    assert "does not retain full-k density or overlap input arrays" in modules[
        "provider_matrix_lean_forward_stream"
    ]["scope"]
    assert "does not imply k-independent total CP2K memory" in modules[
        "provider_bounded_ao_image_batching"
    ]["scope"]

    assert {
        "energy",
        "fock",
        "overlap_adjoint",
        "forces",
        "stress",
        "cold_warm_identity",
        "invalidation",
    } <= set(modules["regular_grid_bvk_cache"]["required_observables"])
    assert {
        "energy",
        "shell_potential",
        "fock",
        "no_retained_full_k_space_density_overlap_arrays",
        "push_order_permutation",
        "twist_and_physical_kpoint_permutation",
        "negative_recovery_paths",
        "9x9x1_dense_oracle_identity",
    } <= set(modules["provider_matrix_lean_forward_stream"]["required_observables"])
    evidence = modules["provider_matrix_lean_forward_stream"]["current_evidence"]
    assert "not a bounded-memory implementation" in evidence
    assert "amat_r, cmat_r and vmat_r" in evidence
    assert "two dense Nk x Nk phase tables" in evidence

    assert {
        "bounded_peak_memory",
        "r_image_batching",
        "batched_or_on_demand_phases",
        "energy",
        "shell_potential",
        "fock",
        "push_order_permutation",
        "9x9x1_dense_oracle_identity",
    } <= set(modules["provider_bounded_ao_image_batching"]["required_observables"])
    assert {
        "overlap_adjoint",
        "forces",
        "stress",
        "reverse_pull_permutation",
        "no_retained_full_k_space_density_overlap_arrays",
        "state_machine_negative_tests",
    } <= set(modules["provider_streamed_reverse"]["required_observables"])
    assert {
        "irreducible_to_star_block_expansion",
        "provider_stream_push_pull",
        "weighted_real_adjoint_foldback",
        "energy",
        "fock",
        "overlap_adjoint",
        "forces",
        "stress",
        "push_pull_permutation",
        "state_machine_negative_tests",
    } <= set(modules["cp2k_streamed_star_consumer"]["required_observables"])


def test_provider_forward_focused_raw_record():
    raw = CAMPAIGN / "raw" / "save_tblite_provider_forward"
    stdout = (raw / "focused_exchange.stdout").read_text(errors="replace")
    stderr = (raw / "focused_exchange.stderr").read_text(errors="replace")

    assert (raw / "returncode.txt").read_text().strip() == "0"
    assert stdout == ""
    assert "bvk_exchange_supercell [PASSED]" in stderr
    assert "[FAILED]" not in stderr
    assert stderr.count("Fortran runtime warning: An array temporary was created") == 42


def test_provider_forward_summary_is_reproducible():
    script = CAMPAIGN / "scripts" / "summarize_save_tblite_provider_forward.py"
    derived_json = CAMPAIGN / "derived" / "save_tblite_provider_forward_summary.json"
    derived_text = CAMPAIGN / "derived" / "save_tblite_provider_forward_summary.txt"

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
    assert result["all_scoped_gates_passed"] is True
    assert result["provider_cache_planner"]["status"] == "passed"
    assert result["matrix_lean_forward_stream"]["status"] == "passed"
    assert (
        result["matrix_lean_forward_stream"]["k_space_input_storage_query"]
        ["reduced_no_retained_full_k_space_density_overlap_assertion_count"]
        == 6
    )
    assert result["matrix_lean_forward_stream"]["order_twist_and_large_mesh"][
        "mesh"
    ] == [9, 9, 1]
    assert result["matrix_lean_forward_stream"]["k_space_input_storage_query"][
        "query_scope"
    ] == "stream%density and stream%overlap allocations only"
    assert "stream%amat_r" in result["matrix_lean_forward_stream"][
        "k_space_input_storage_query"
    ]["not_measured"]
    assert "cache%bvk_phase_inverse" in result["matrix_lean_forward_stream"][
        "k_space_input_storage_query"
    ]["not_measured"]
    assert (
        result["not_qualified"]["true_bounded_memory_r_image_batching"]["status"]
        == "implementation_in_progress"
    )
    assert (
        result["not_qualified"]["reduced_memory_reverse_stream"]["status"]
        == "implementation_in_progress"
    )
    assert (
        result["not_qualified"]["cp2k_consumer_integration"]["status"]
        == "implementation_in_progress"
    )


def test_provider_forward_exact_source_and_patch_provenance():
    provenance = load_json("provenance/save_tblite_provider_forward.json")
    directory = CAMPAIGN / "provenance" / "save_tblite_provider_forward"
    snapshot = directory / "source_snapshot"

    assert digest(directory / "save_tblite_tested.patch") == provenance[
        "save_tblite"
    ]["worktree_patch_sha256"]
    for relative, expected in provenance["save_tblite"]["source_file_sha256"].items():
        assert digest(snapshot / relative) == expected
    assert digest(directory / "CMakeCache.txt") == provenance["build"][
        "cmake_cache_sha256"
    ]

    patch = (directory / "save_tblite_tested.patch").read_text()
    assert "cp2k_exchange_stream_has_full_mesh_storage" in patch
    assert "get_KFock_stream_apply" in patch
    assert "integer, parameter :: nmesh_large(3) = [9, 9, 1]" in patch
    assert "g-xTB exchange stream reverse requires oracle mode" in patch
    qualification = provenance["qualification"]
    assert qualification["matrix_lean_forward_stream"] == "passed"
    assert (
        qualification["true_bounded_memory_r_image_batching"]
        == "implementation_in_progress"
    )
    query = provenance["runtime"]["k_space_input_storage_query"]
    assert query["implementation_scope"].endswith("only")
    assert "total process memory" in query["not_measured"]


def test_earlier_terok_record_is_preserved_but_not_qualification_basis():
    raw = CAMPAIGN / "raw" / "save_tblite_provider_cache_terok_earlier"
    log = (raw / "CTest-LastTest.log").read_text(errors="replace")
    summary = load_json("derived/save_tblite_provider_forward_summary.json")

    assert log.count("[PASSED]") == 30
    assert "Test Passed." in log
    assert summary["historical_terok_record"]["preserved"] is True
    assert summary["historical_terok_record"]["qualification_basis"] is False


def test_provider_forward_archive_manifest_is_self_verifying():
    manifest = (
        CAMPAIGN
        / "provenance"
        / "save_tblite_provider_forward"
        / "SHA256SUMS"
    )
    entries = manifest.read_text().splitlines()

    assert len(entries) == 21
    for line in entries:
        expected, relative = line.split(maxsplit=1)
        assert digest(CAMPAIGN / relative) == expected


def test_acceleration_catalog_matches_the_existing_curated_archives():
    matrix = load_json("validation_matrix.json")
    modules = {entry["id"]: entry for entry in matrix["modules"]}
    cases = {entry["id"]: entry for entry in matrix["cases"]}

    expected_modules = {
        "regular_mesh_fft": "validation/accelerated_exchange/mixed_radix_fft_20260717",
        "replicated_importer_disjoint_pull": (
            "validation/accelerated_exchange/cp2k_distributed_images_20260717"
        ),
        "distributed_nonlinear_exchange_kernel": (
            "validation/accelerated_exchange/cp2k_distributed_images_20260717"
        ),
        "streamed_symmetry_covariance_check": (
            "validation/accelerated_exchange/cp2k_streamed_star_memory"
        ),
    }
    for module_id, evidence in expected_modules.items():
        entry = modules[module_id]
        assert entry["status"] == "passed"
        assert entry["evidence"] == evidence
        assert (ROOT / evidence).is_dir()

    assert cases["regular_mesh_fft_oracle_matrix"]["status"] == "passed"
    assert cases["regular_mesh_fft_oracle_matrix"]["maximum_residuals"][
        "energy_sequence_hartree"
    ] < 1.0e-12
    assert cases["streamed_symmetry_covariance_matrix"]["runs"] == 48
    assert cases["distributed_nonlinear_oracle_matrix"]["pairs"] == 30
    assert cases["distributed_nonlinear_oracle_matrix"]["faults"] == 6
    assert cases["distributed_kernel_scaling"]["status"] == "pending"


def load_tests(loader, standard_tests, pattern):
    """Expose the assertion-style campaign checks to unittest discovery."""
    del loader, standard_tests, pattern
    suite = unittest.TestSuite()
    for name, check in sorted(globals().items()):
        if name.startswith("test_") and callable(check):
            suite.addTest(unittest.FunctionTestCase(check, description=name))
    return suite


if __name__ == "__main__":
    unittest.main()
