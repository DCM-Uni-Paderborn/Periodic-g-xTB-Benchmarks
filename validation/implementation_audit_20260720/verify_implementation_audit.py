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
        "geometry_equivalence": "validation/geometry_equivalence_20260720/verification.json",
        "mstore_accuracy_equivalence": "validation/mstore_accuracy_equivalence_20260720/verification.json",
        "dmc13_discrepancy_attribution": "validation/dmc13_discrepancy_attribution_20260720/verification.json",
        "three_route_k333_closure": "validation/three_route_k333_closure_20260719/summary.json",
        "lowk_derivatives_and_partial_pbc": "validation/gxtb_final_lowk_derivatives_20260719/verification.reproduced.json",
        "energy_component_ablation": "validation/dmc13_k222_viii_component_ablation_20260719/verification.json",
        "derivative_component_ablation": "validation/dmc13_k222_xvii_derivative_component_ablation_20260719/verification.json",
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
    native_absolute_relative = (
        "DMC-ICE13/reproduction/seidler_dmc13_recalculation/"
        "tables/cp2k_native_absolute_energies_by_mesh.csv"
    )
    with (ROOT / adaptive_relative).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (ROOT / statistics_relative).open(newline="", encoding="utf-8") as handle:
        statistics = next(csv.DictReader(handle))
    with (ROOT / native_absolute_relative).open(newline="", encoding="utf-8") as handle:
        native_absolute_rows = list(csv.DictReader(handle))

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
        "convergence_criterion_valid": convergence_criterion_valid,
        "adaptive_tolerance_kj_mol_per_water": ADAPTIVE_TOLERANCE_KJ_MOL_PER_WATER,
        "final_state_consistent": final_state_consistent,
    }

    ih_meshes = {
        int(row["mesh_n"])
        for row in native_absolute_rows
        if row["phase"] == "Ih" and row["qualification"] == "PASS"
    }
    pending_science_endpoints = []
    for row in rows:
        if row["phase_converged"].strip().lower() == "true":
            continue
        selected_mesh = int(row["mesh_n"])
        delta_text = row["absolute_adjacent_delta_kj_mol_per_water"].strip()
        required_mesh = selected_mesh + 1 if delta_text or selected_mesh == 1 else selected_mesh - 1
        endpoint = f"ice {row['phase']} {required_mesh}x{required_mesh}x{required_mesh}"
        if required_mesh not in ih_meshes:
            endpoint += (
                f" after same-build ice Ih "
                f"{required_mesh}x{required_mesh}x{required_mesh}"
            )
        pending_science_endpoints.append(endpoint)

    all_completed_gates_pass = all(item["passed"] for item in gate_results.values())
    output = {
        "schema": "periodic-gxtb-part-i-implementation-audit-v1",
        "status": "PASS" if all_completed_gates_pass else "FAIL",
        "completed_gate_count": len(gate_results),
        "completed_gates": gate_results,
        "current_adaptive_set": calculated,
        "current_adaptive_set_is_final": declared_final,
        "pending_science_endpoints": pending_science_endpoints,
        "pending_diagnostic_endpoint": None,
        "interpretation": (
            "All completed exact implementation gates pass and the DMC-ICE13 adaptive "
            "statistic is final."
            if declared_final and not pending_science_endpoints
            else "All completed exact implementation gates pass. The DMC-ICE13 adaptive "
            "statistic remains provisional until the listed phase-local endpoints finish."
        ),
    }
    output_path = HERE / "verification.json"
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))
    if not all_completed_gates_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
