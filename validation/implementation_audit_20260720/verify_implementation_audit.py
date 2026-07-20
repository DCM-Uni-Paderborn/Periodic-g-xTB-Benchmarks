#!/usr/bin/env python3
"""Aggregate the exact Part-I periodic g-xTB validation gates."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPECTED_PHASES = (
    "II", "III", "IV", "VI", "VII", "VIII",
    "IX", "XI", "XIII", "XIV", "XV", "XVII",
)
ADAPTIVE_TOLERANCE_KJ_MOL_PER_WATER = 0.10
MAXIMUM_ADAPTIVE_MESH = 8
PENDING_DIAGNOSTIC_ENDPOINTS = ()
QUALIFIED_CP2K_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
SEIDLER_PACKAGE = Path("DMC-ICE13/reproduction/seidler_dmc13_recalculation")


def load_json(relative: str) -> dict:
    with (ROOT / relative).open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def status_passes(payload: dict) -> bool:
    return str(payload.get("status", "")).lower() in {"pass", "passed"}


def main() -> None:
    gates = {
        "binary_provider_identity": "validation/binary_provider_identity_20260720/verification.json",
        "qualified_build_head_delta": "validation/qualified_build_head_delta_20260720/verification.json",
        "macdonald_gamma_bvk_mesh_equivalence": "validation/macdonald_bvk_mesh_equivalence_20260720/verification.json",
        "native_cli_full_k111_to_k444_parity": "validation/native_cli_full_parity_20260720/verification.json",
        "native_cli_inprocess_derivative_parity": "validation/native_cli_inprocess_derivatives_20260720/verification.json",
        "native_restart_equivalence": "validation/gxtb_restart_equivalence_20260720/verification.json",
        "restart_operational_safety": "validation/restart_operational_safety_20260720/verification.json",
        "seidler_package_tracking": "validation/seidler_package_tracking_20260720/verification.json",
        "relative_energy_postprocessing": "validation/relative_energy_postprocessing_20260720/verification.json",
        "native_cli_tight_k222": "validation/tight_parity_k222_20260720/verification.json",
        "accuracy_equivalence": "validation/accuracy_equivalence_20260720/verification.json",
        "native_gamma_bvk_oracle": "validation/explicit_cp2k_gamma_supercell_oracle_20260719/verification.reproduced.json",
        "native_gamma_bvk_oracle_ice_vii": "validation/explicit_cp2k_gamma_supercell_oracle_20260719/verification-vii.json",
        "geometry_equivalence": "validation/geometry_equivalence_20260720/verification.json",
        "mstore_accuracy_equivalence": "validation/mstore_accuracy_equivalence_20260720/verification.json",
        "mstore_pbc_component_ablation": "validation/mstore_pbc_component_ablation_20260720/verification.json",
        "mstore_inorganic_k444_partial": "DMC-ICE13/reproduction/seidler_dmc13_recalculation/evidence/mstore_inorganic_k444_partial/verification.json",
        "wigner_seitz_self_image_attribution": "validation/wigner_seitz_self_image_attribution_20260720/verification.json",
        "second_order_mic_attribution": "validation/second_order_mic_attribution_20260720/verification.json",
        "dmc13_discrepancy_attribution": "validation/dmc13_discrepancy_attribution_20260720/verification.json",
        "three_route_k333_closure": "validation/three_route_k333_closure_20260719/summary.json",
        "lowk_derivatives_and_partial_pbc": "validation/gxtb_final_lowk_derivatives_20260719/verification.reproduced.json",
        "energy_component_ablation": "validation/dmc13_k222_viii_component_ablation_20260719/verification.json",
        "derivative_component_ablation": "validation/dmc13_k222_xvii_derivative_component_ablation_20260719/verification.json",
        "dmc13_xvii_full_derivatives": "validation/dmc13_xvii_full_derivatives_20260718/verification.json",
        "dmc13_xvii_same_binary_derivatives": "validation/dmc13_xvii_full_derivatives_20260718/current_binary_requalification/verification.json",
        "provider_component_attribution": "validation/provider_component_attribution_20260719/verification.json",
        "h0_anisotropy_attribution": "validation/pbc_h0_anisotropy_attribution_20260719/verification.json",
        "cecl3_tolerance_recheck": "validation/cecl3_tolerance_recheck_20260720/verification.reproduced.json",
    }
    gate_results = {}
    for name, relative in gates.items():
        path = ROOT / relative
        payload = load_json(relative)
        gate_results[name] = {
            "status": payload.get("status"),
            "passed": status_passes(payload),
            "file": relative,
            "sha256": sha256(path),
        }

    source_relative = "validation/save_tblite_periodic_source_tests_20260719/results.json"
    source = load_json(source_relative)
    required_source_tests = [
        "current_h0_diamond",
        "current_h0_supercell",
        "current_h0_gradient",
        "current_wignerseitz",
        "current_exchange",
        "current_acp",
        "current_coulomb_charge",
        "current_coulomb_multipole",
        "current_dispersion",
        "current_repulsion",
        "current_gxtb",
    ]
    periodic_source_pass = all(
        source["tests"][name]["returncode"] == 0
        and source["tests"][name]["failed_count"] == 0
        for name in required_source_tests
    )
    cecl3 = source["cecl3_nonperiodic_finite_difference"]
    inherited_cecl3_match = (
        cecl3["identical_difference_components"]
        and math.isclose(
            cecl3["current_max_abs_hartree_per_bohr"],
            cecl3["pbc_max_abs_hartree_per_bohr"],
            rel_tol=0.0,
            abs_tol=1.0e-18,
        )
    )
    source_gate_pass = periodic_source_pass and inherited_cecl3_match
    gate_results["periodic_source_tests"] = {
        "status": "PASS" if source_gate_pass else "FAIL",
        "passed": source_gate_pass,
        "file": source_relative,
        "sha256": sha256(ROOT / source_relative),
        "required_test_groups": required_source_tests,
        "inherited_nonperiodic_cecl3_residual_identical_to_author_pbc": inherited_cecl3_match,
    }

    adaptive_relative = "DMC-ICE13/data/dmc_ice13_gxtb_current_adaptive_set.csv"
    statistics_relative = "DMC-ICE13/data/dmc_ice13_gxtb_current_adaptive_statistics.csv"
    progress_relative = "DMC-ICE13/data/dmc_ice13_gxtb_phasewise_progress.csv"
    native_absolute_relative = (
        "DMC-ICE13/reproduction/seidler_dmc13_recalculation/"
        "tables/cp2k_native_absolute_energies_by_mesh.csv"
    )
    with (ROOT / adaptive_relative).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (ROOT / statistics_relative).open(newline="", encoding="utf-8") as handle:
        statistics = next(csv.DictReader(handle))
    with (ROOT / progress_relative).open(newline="", encoding="utf-8") as handle:
        progress_rows = list(csv.DictReader(handle))
    with (ROOT / native_absolute_relative).open(newline="", encoding="utf-8") as handle:
        native_absolute_rows = list(csv.DictReader(handle))

    native_endpoint_checks = []
    for row in native_absolute_rows:
        output = ROOT / SEIDLER_PACKAGE / row["raw_output"]
        status_file = output.parent / "exit_status"
        output_text = (
            output.read_text(encoding="utf-8", errors="replace")
            if output.is_file()
            else ""
        )
        passed = (
            row.get("qualification") == "PASS"
            and row.get("cp2k_binary_sha256") == QUALIFIED_CP2K_SHA256
            and row.get("exit_status") == "0"
            and row.get("normal_termination_qualification") == "PASS"
            and status_file.is_file()
            and status_file.read_text(encoding="utf-8", errors="replace").strip() == "0"
            and "PROGRAM ENDED AT" in output_text
            and output.is_file()
            and sha256(output) == row.get("output_sha256")
        )
        native_endpoint_checks.append({
            "mesh_n": int(row["mesh_n"]),
            "phase": row["phase"],
            "passed": passed,
            "raw_output": row["raw_output"],
        })
    native_endpoint_qualification_pass = bool(native_endpoint_checks) and all(
        item["passed"] for item in native_endpoint_checks
    )
    gate_results["native_endpoint_build_and_termination"] = {
        "status": "PASS" if native_endpoint_qualification_pass else "FAIL",
        "passed": native_endpoint_qualification_pass,
        "file": native_absolute_relative,
        "sha256": sha256(ROOT / native_absolute_relative),
        "qualified_cp2k_sha256": QUALIFIED_CP2K_SHA256,
        "endpoint_count": len(native_endpoint_checks),
        "failed_endpoints": [
            item for item in native_endpoint_checks if not item["passed"]
        ],
        "requirements": [
            "exact qualified CP2K binary SHA-256",
            "archived process exit status 0",
            "PROGRAM ENDED AT marker",
            "raw output SHA-256 matches table",
        ],
    }

    phases = [row["phase"] for row in rows]
    phase_matrix_valid = (
        len(phases) == len(EXPECTED_PHASES)
        and len(phases) == len(set(phases))
        and set(phases) == set(EXPECTED_PHASES)
    )
    convergence_values_valid = all(
        row["phase_converged"].strip().lower() in {"true", "false"}
        for row in rows
    )
    native_by_key = {
        (int(row["mesh_n"]), row["phase"]): row
        for row in native_absolute_rows
    }
    native_keys_unique = len(native_by_key) == len(native_absolute_rows)
    same_mesh_phase_ih_pairs_valid = native_keys_unique and all(
        (int(row["mesh_n"]), row["phase"]) in native_by_key
        and (int(row["mesh_n"]), "Ih") in native_by_key
        and native_by_key[(int(row["mesh_n"]), row["phase"])]["qualification"]
        == "PASS"
        and native_by_key[(int(row["mesh_n"]), "Ih")]["qualification"] == "PASS"
        and native_by_key[(int(row["mesh_n"]), row["phase"])][
            "cp2k_binary_sha256"
        ]
        == QUALIFIED_CP2K_SHA256
        and native_by_key[(int(row["mesh_n"]), "Ih")]["cp2k_binary_sha256"]
        == QUALIFIED_CP2K_SHA256
        and native_by_key[(int(row["mesh_n"]), row["phase"])]["exit_status"] == "0"
        and native_by_key[(int(row["mesh_n"]), "Ih")]["exit_status"] == "0"
        for row in rows
    )
    convergence_criterion_valid = convergence_values_valid
    for row in rows:
        converged = row["phase_converged"].strip().lower() == "true"
        delta_text = row["absolute_adjacent_delta_kj_mol_per_water"].strip()
        delta = None if not delta_text else float(delta_text)
        if converged:
            convergence_criterion_valid = (
                convergence_criterion_valid
                and delta is not None
                and delta <= ADAPTIVE_TOLERANCE_KJ_MOL_PER_WATER
            )
        elif delta is not None:
            convergence_criterion_valid = (
                convergence_criterion_valid
                and delta > ADAPTIVE_TOLERANCE_KJ_MOL_PER_WATER
            )

    errors = [float(row["error_kj_mol_per_water"]) for row in rows]
    abs_errors = [abs(value) for value in errors]
    calculated = {
        "phase_count": len(rows),
        "converged_phase_count": sum(
            row["phase_converged"].strip().lower() == "true" for row in rows
        ),
        "largest_mesh_n": max(int(row["mesh_n"]) for row in rows),
        "me_kj_mol_per_water": sum(errors) / len(errors),
        "mae_kj_mol_per_water": sum(abs_errors) / len(abs_errors),
        "rmse_kj_mol_per_water": math.sqrt(
            sum(value * value for value in errors) / len(errors)
        ),
        "maxae_kj_mol_per_water": max(abs_errors),
    }
    declared_final_text = statistics["final_result"].strip().lower()
    declared_final_valid = declared_final_text in {"true", "false"}
    declared_final = declared_final_text == "true"
    final_state_consistent = (
        declared_final_valid
        and declared_final == (calculated["converged_phase_count"] == len(EXPECTED_PHASES))
    )
    adaptive_match = (
        phase_matrix_valid
        and same_mesh_phase_ih_pairs_valid
        and convergence_criterion_valid
        and final_state_consistent
        and calculated["phase_count"] == int(statistics["phase_count"])
        and calculated["converged_phase_count"] == int(statistics["converged_phase_count"])
        and calculated["largest_mesh_n"] == int(statistics["largest_mesh_n"])
        and all(
            math.isclose(
                calculated[key],
                float(statistics[key]),
                rel_tol=0.0,
                # The machine-readable statistics are written with twelve
                # decimal places.  Allow one unit in the final written place
                # when recomputing them from the likewise rounded phase rows.
                abs_tol=1.0e-12,
            )
            for key in (
                "me_kj_mol_per_water",
                "mae_kj_mol_per_water",
                "rmse_kj_mol_per_water",
                "maxae_kj_mol_per_water",
            )
        )
    )
    gate_results["adaptive_statistics_internal_consistency"] = {
        "status": "PASS" if adaptive_match else "FAIL",
        "passed": adaptive_match,
        "files": [adaptive_relative, statistics_relative, native_absolute_relative],
        "sha256": {
            adaptive_relative: sha256(ROOT / adaptive_relative),
            statistics_relative: sha256(ROOT / statistics_relative),
            native_absolute_relative: sha256(ROOT / native_absolute_relative),
        },
        "phase_matrix_valid": phase_matrix_valid,
        "native_endpoint_keys_unique": native_keys_unique,
        "same_mesh_phase_and_ih_exact_build_pairs_valid": same_mesh_phase_ih_pairs_valid,
        "convergence_criterion_valid": convergence_criterion_valid,
        "adaptive_tolerance_kj_mol_per_water": ADAPTIVE_TOLERANCE_KJ_MOL_PER_WATER,
        "final_state_consistent": final_state_consistent,
    }

    expected_progress_limits = [6, 7, 8]
    progress_limits = [int(row["mesh_limit_n"]) for row in progress_rows]
    progress_pass_counts = [
        int(row["converged_phase_count"]) for row in progress_rows
    ]
    progress_match = (
        progress_limits == expected_progress_limits
        and len(progress_rows) == len(expected_progress_limits)
        and all(int(row["phase_count"]) == len(EXPECTED_PHASES) for row in progress_rows)
        and all(row["qualification"] == "PASS" for row in progress_rows)
        and progress_pass_counts == sorted(progress_pass_counts)
        and progress_pass_counts[-1] == calculated["converged_phase_count"]
        and all(
            math.isfinite(float(row[field]))
            for row in progress_rows
            for field in (
                "me_kj_mol_per_water",
                "mae_kj_mol_per_water",
                "rmse_kj_mol_per_water",
                "maxae_kj_mol_per_water",
            )
        )
        and math.isclose(
            float(progress_rows[-1]["me_kj_mol_per_water"]),
            calculated["me_kj_mol_per_water"],
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            float(progress_rows[-1]["mae_kj_mol_per_water"]),
            calculated["mae_kj_mol_per_water"],
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            float(progress_rows[-1]["rmse_kj_mol_per_water"]),
            calculated["rmse_kj_mol_per_water"],
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            float(progress_rows[-1]["maxae_kj_mol_per_water"]),
            calculated["maxae_kj_mol_per_water"],
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    )
    gate_results["phasewise_progress_internal_consistency"] = {
        "status": "PASS" if progress_match else "FAIL",
        "passed": progress_match,
        "file": progress_relative,
        "sha256": sha256(ROOT / progress_relative),
        "mesh_limits": progress_limits,
        "converged_phase_counts": progress_pass_counts,
        "latest_matches_current_adaptive_statistics": progress_match,
    }

    ih_meshes = {
        int(row["mesh_n"])
        for row in native_absolute_rows
        if row["phase"] == "Ih" and row["qualification"] == "PASS"
    }
    pending_science_endpoints = []
    capped_unresolved_phases = []
    for row in rows:
        if row["phase_converged"].strip().lower() == "true":
            continue
        selected_mesh = int(row["mesh_n"])
        delta_text = row["absolute_adjacent_delta_kj_mol_per_water"].strip()
        required_mesh = selected_mesh + 1 if delta_text or selected_mesh == 1 else selected_mesh - 1
        if required_mesh > MAXIMUM_ADAPTIVE_MESH:
            capped_unresolved_phases.append(row["phase"])
            continue
        endpoint = f"ice {row['phase']} {required_mesh}x{required_mesh}x{required_mesh}"
        if required_mesh not in ih_meshes:
            endpoint += (
                f" after same-build ice Ih "
                f"{required_mesh}x{required_mesh}x{required_mesh}"
            )
        pending_science_endpoints.append(endpoint)

    all_completed_gates_pass = all(item["passed"] for item in gate_results.values())
    output = {
        "schema": "periodic-gxtb-part-i-implementation-audit-v2",
        "status": "PASS" if all_completed_gates_pass else "FAIL",
        "completed_gate_count": len(gate_results),
        "completed_gates": gate_results,
        "current_adaptive_set": calculated,
        "current_adaptive_set_is_final": declared_final,
        "pending_science_endpoints": pending_science_endpoints,
        "maximum_adaptive_mesh": MAXIMUM_ADAPTIVE_MESH,
        "capped_unresolved_phases": capped_unresolved_phases,
        "pending_diagnostic_endpoints": list(PENDING_DIAGNOSTIC_ENDPOINTS),
        "interpretation": (
            "All completed exact implementation gates pass and the DMC-ICE13 adaptive "
            "statistic is final."
            if (
                declared_final
                and not pending_science_endpoints
                and not PENDING_DIAGNOSTIC_ENDPOINTS
            )
            else "All completed exact implementation gates pass. The DMC-ICE13 adaptive "
            "statistic remains provisional while listed sub-cap endpoints are pending or "
            "a phase remains unresolved at the declared mesh cap; the separately listed "
            "diagnostic endpoint is not yet a completed gate."
        ),
    }
    output_path = HERE / "verification.json"
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))
    if not all_completed_gates_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
