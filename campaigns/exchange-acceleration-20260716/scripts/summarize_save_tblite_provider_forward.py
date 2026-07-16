#!/usr/bin/env python3
"""Summarize the exact save_tblite cache and reduced-forward qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parent.parent
RAW = CAMPAIGN / "raw" / "save_tblite_provider_forward"
HISTORICAL = CAMPAIGN / "raw" / "save_tblite_provider_cache_terok_earlier"
PROVENANCE = CAMPAIGN / "provenance" / "save_tblite_provider_forward"
SNAPSHOT = PROVENANCE / "source_snapshot"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_all(text: str, markers: tuple[str, ...]) -> dict[str, bool]:
    return {marker: marker in text for marker in markers}


def build_result() -> dict[str, object]:
    stdout = (RAW / "focused_exchange.stdout").read_text(errors="replace")
    stderr = (RAW / "focused_exchange.stderr").read_text(errors="replace")
    output = stdout + stderr
    returncode = int((RAW / "returncode.txt").read_text().strip())
    test_source_path = SNAPSHOT / "test" / "unit" / "test_exchange.f90"
    compat_path = SNAPSHOT / "src" / "tblite" / "cp2k_compat.f90"
    fock_path = SNAPSHOT / "src" / "tblite" / "exchange" / "fock.f90"
    test_source = test_source_path.read_text()
    compat_source = compat_path.read_text()
    fock_source = fock_path.read_text()

    runtime_passed = (
        returncode == 0
        and "bvk_exchange_supercell [PASSED]" in output
        and "[FAILED]" not in output
    )
    cache_markers = require_all(
        test_source + compat_source,
        (
            "integer, parameter :: nmesh_large(3) = [9, 9, 1]",
            "g-xTB k-point mesh contains a duplicate grid point",
            "g-xTB exchange requires a complete regular k mesh",
            "bvk_plan_builds",
            "bvk_input_to_grid",
            "bvk_grid_to_input",
            "bvk_matches",
        ),
    )
    forward_markers = require_all(
        test_source + compat_source + fock_source,
        (
            "logical function cp2k_exchange_stream_has_full_mesh_storage",
            "selected_mode = cp2k_exchange_stream_reduced",
            "duplicate g-xTB exchange stream block",
            "g-xTB exchange stream has missing blocks",
            "g-xTB exchange stream state was invalidated before apply",
            "stream_permuted",
            "get_KFock_stream_apply",
            "fock_stream-fock_cp2k",
            "fock_large-fock_large_oracle",
        ),
    )
    storage = {
        "reduced_no_retained_full_k_space_density_overlap_assertion_count": test_source.count(
            "if (cp2k_exchange_stream_has_full_mesh_storage("
        ),
        "oracle_retained_full_k_space_density_overlap_assertion_count": test_source.count(
            "if (.not.cp2k_exchange_stream_has_full_mesh_storage("
        ),
        "query_implementation_present": (
            "has_storage = allocated(stream%density) .or. allocated(stream%overlap)"
            in compat_source
        ),
        "query_scope": "stream%density and stream%overlap allocations only",
        "not_measured": [
            "stream%amat_r",
            "stream%cmat_r",
            "stream%vmat_r",
            "cache%bvk_phase_forward",
            "cache%bvk_phase_inverse",
            "total_process memory",
        ],
    }
    storage["passed"] = (
        runtime_passed
        and storage[
            "reduced_no_retained_full_k_space_density_overlap_assertion_count"
        ]
        >= 5
        and storage[
            "oracle_retained_full_k_space_density_overlap_assertion_count"
        ]
        >= 1
        and storage["query_implementation_present"]
    )
    negative_recovery = {
        "duplicate_push": "Duplicate exchange stream block was accepted" in test_source,
        "missing_push_then_complete": "Incomplete exchange stream was applied" in test_source,
        "stale_onsite_then_restore": "Invalidated exchange stream was applied" in test_source,
        "incomplete_grid_then_valid_retry": (
            "Incomplete whole-mesh exchange grid was accepted" in test_source
        ),
        "duplicate_grid_then_valid_retry": (
            "Duplicate whole-mesh exchange grid point was accepted" in test_source
        ),
    }
    negative_recovery["passed"] = runtime_passed and all(negative_recovery.values())
    order_and_mesh = {
        "arbitrary_push_order": "ik = korder(j)" in test_source,
        "common_twist": "A common twist and a permutation" in test_source,
        "physical_kpoint_permutation": "stream_permuted" in test_source,
        "mesh": [9, 9, 1],
        "mesh_points": 81,
        "large_mesh_reduced_forward_vs_dense_oracle": (
            "fock_large-fock_large_oracle" in test_source
            and "cp2k_exchange_stream_begin(stream_large" in test_source
        ),
    }
    order_and_mesh["passed"] = runtime_passed and all(
        value for key, value in order_and_mesh.items() if key not in {"mesh", "mesh_points"}
    )
    current_passed = (
        runtime_passed
        and all(cache_markers.values())
        and all(forward_markers.values())
        and storage["passed"]
        and negative_recovery["passed"]
        and order_and_mesh["passed"]
    )

    historical_log = (HISTORICAL / "CTest-LastTest.log").read_text(errors="replace")
    historical_passed = re.findall(r"\.\.\.\s+(\S+) \[PASSED\]", historical_log)
    return {
        "schema_version": 1,
        "component": "save-tblite-provider-cache-and-matrix-lean-forward-stream",
        "oracle": "unchanged-complete-mesh-and-dense-transform-provider-path",
        "focused_runtime": {
            "command": (RAW / "command.txt").read_text().strip(),
            "returncode": returncode,
            "test": "bvk_exchange_supercell",
            "passed": runtime_passed,
            "array_temporary_warning_count": output.count(
                "Fortran runtime warning: An array temporary was created"
            ),
            "stdout_sha256": sha256(RAW / "focused_exchange.stdout"),
            "stderr_sha256": sha256(RAW / "focused_exchange.stderr"),
        },
        "provider_cache_planner": {
            "status": "passed" if current_passed else "failed",
            "markers": cache_markers,
        },
        "matrix_lean_forward_stream": {
            "status": "passed" if current_passed else "failed",
            "observables": ["energy", "shell_potential", "fock"],
            "classification": (
                "no retained full-k-space density/overlap arrays; not a "
                "bounded-memory implementation"
            ),
            "markers": forward_markers,
            "k_space_input_storage_query": storage,
            "negative_and_recovery_tests": negative_recovery,
            "order_twist_and_large_mesh": order_and_mesh,
        },
        "not_qualified": {
            "true_bounded_memory_r_image_batching": {
                "status": "implementation_in_progress",
                "evidence": (
                    "Reduced forward mode still retains three nao x nao x Nk x nspin "
                    "BvK-image tensors (amat_r, cmat_r, vmat_r) and the cache retains "
                    "two dense Nk x Nk phase tables. R/image batching and batched or "
                    "on-demand phases are not implemented."
                ),
            },
            "reduced_memory_reverse_stream": {
                "status": "implementation_in_progress",
                "evidence": (
                    "The tested reduced stream intentionally rejects reverse_apply; "
                    "only the retained full-mesh oracle supports the reverse API."
                ),
                "diagnostic_present": (
                    "g-xTB exchange stream reverse requires oracle mode" in test_source
                ),
            },
            "cp2k_consumer_integration": {
                "status": "implementation_in_progress",
                "evidence": (
                    "This is a provider-unit qualification; no production CP2K "
                    "consumer invoked the stream transaction."
                ),
            },
        },
        "historical_terok_record": {
            "preserved": True,
            "qualification_basis": False,
            "passed_subtests": len(historical_passed),
            "reported_total": 30,
            "reason": (
                "The executable at the recorded path was replaced while the test was "
                "running, so it is retained as raw history but not used for byte-exact "
                "source qualification."
            ),
        },
        "all_scoped_gates_passed": current_passed,
        "scope_note": (
            "Only the provider cache/planner and matrix-lean forward stream are "
            "qualified. The storage query proves only that full-k-space density and "
            "overlap input arrays are not retained. True bounded-memory R/image "
            "batching, reduced-memory reverse, and actual CP2K integration remain in "
            "progress."
        ),
    }


def render_text(result: dict[str, object]) -> str:
    runtime = result["focused_runtime"]
    forward = result["matrix_lean_forward_stream"]
    storage = forward["k_space_input_storage_query"]
    large = forward["order_twist_and_large_mesh"]
    lines = [
        "SAVE_TBLITE_PROVIDER_MATRIX_LEAN_FORWARD_QUALIFICATION",
        f"test\t{runtime['test']}",
        f"returncode\t{runtime['returncode']}",
        f"provider_cache_planner\t{result['provider_cache_planner']['status']}",
        f"matrix_lean_forward_stream\t{forward['status']}",
        "storage_query_scope\tfull-k-space density/overlap arrays only",
        f"reduced_no_retained_kspace_input_queries\t{storage['reduced_no_retained_full_k_space_density_overlap_assertion_count']}",
        f"oracle_retained_kspace_input_queries\t{storage['oracle_retained_full_k_space_density_overlap_assertion_count']}",
        f"large_mesh\t{'x'.join(str(value) for value in large['mesh'])}",
        f"negative_recovery\t{str(forward['negative_and_recovery_tests']['passed']).lower()}",
        f"twist_permutation\t{str(large['passed']).lower()}",
        "true_bounded_memory_r_image_batching\timplementation_in_progress",
        "reduced_memory_reverse_stream\timplementation_in_progress",
        "cp2k_consumer_integration\timplementation_in_progress",
        f"ALL_SCOPED_GATES_PASSED\t{str(result['all_scoped_gates_passed']).lower()}",
        f"SCOPE\t{result['scope_note']}",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "text"), default="json")
    args = parser.parse_args()
    result = build_result()
    if args.format == "json":
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(render_text(result), end="")


if __name__ == "__main__":
    main()
