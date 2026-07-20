#!/usr/bin/env python3
"""Independently quantify the DMC-ICE13 discrepancy hierarchy."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation/tables"

NATIVE = "CP2K-native pbc provider"
CURRENT = "current pbc CLI"
MSTORE = "historical mstore-inorganic CLI"


def load_json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def statistics(rows: list[dict[str, str]]) -> dict[str, float | int]:
    errors = [float(row["error_kj_mol_per_H2O"]) for row in rows]
    return {
        "phase_count": len(errors),
        "me_kj_mol_per_H2O": sum(errors) / len(errors),
        "mae_kj_mol_per_H2O": sum(abs(value) for value in errors) / len(errors),
        "rmse_kj_mol_per_H2O": math.sqrt(
            sum(value * value for value in errors) / len(errors)
        ),
        "maxae_kj_mol_per_H2O": max(abs(value) for value in errors),
    }


def main() -> None:
    parity_path = TABLES / "pbc_cli_vs_cp2k_native_absolute_parity.csv"
    all_branch_path = TABLES / "all_branch_relative_energy_comparison.csv"
    branch_statistics_path = TABLES / "branch_comparison_statistics.csv"
    author_path = TABLES / "author_pbc_relative_energies.csv"
    mstore_path = TABLES / "mstore_vs_pbc_relative_differences.csv"
    adaptive_path = ROOT / "DMC-ICE13/data/dmc_ice13_gxtb_current_adaptive_set.csv"
    adaptive_statistics_path = (
        ROOT / "DMC-ICE13/data/dmc_ice13_gxtb_current_adaptive_statistics.csv"
    )

    parity_rows = load_csv(parity_path)
    all_branch_rows = load_csv(all_branch_path)
    archived_statistics = load_csv(branch_statistics_path)
    author_rows = load_csv(author_path)
    mstore_rows = load_csv(mstore_path)
    adaptive_rows = load_csv(adaptive_path)
    adaptive_declared = load_csv(adaptive_statistics_path)[0]

    parity_max = max(
        (abs(float(row["native_minus_cli_Ha_per_primitive"])), row)
        for row in parity_rows
    )
    by_route = {
        (row["method"], int(row["mesh_n"]), row["phase"]): row
        for row in all_branch_rows
    }
    relative_differences = []
    for mesh in range(1, 5):
        for phase in sorted(
            row["phase"]
            for row in all_branch_rows
            if row["method"] == NATIVE and int(row["mesh_n"]) == mesh
        ):
            native = float(
                by_route[(NATIVE, mesh, phase)]["relative_energy_kj_mol_per_H2O"]
            )
            current = float(
                by_route[(CURRENT, mesh, phase)]["relative_energy_kj_mol_per_H2O"]
            )
            relative_differences.append((abs(native - current), mesh, phase))
    relative_max = max(relative_differences)

    reconstructed_statistics = {}
    statistics_match = True
    archived_index = {
        (row["method"], int(row["mesh_n"])): row for row in archived_statistics
    }
    for method in (NATIVE, CURRENT, MSTORE):
        meshes = sorted(
            mesh for archived_method, mesh in archived_index if archived_method == method
        )
        for mesh in meshes:
            subset = [
                row
                for row in all_branch_rows
                if row["method"] == method and int(row["mesh_n"]) == mesh
            ]
            calculated = statistics(subset)
            reconstructed_statistics[f"{method}|{mesh}"] = calculated
            archived = archived_index[(method, mesh)]
            statistics_match = statistics_match and all(
                math.isclose(
                    float(calculated[key]),
                    float(archived[key]),
                    rel_tol=0.0,
                    abs_tol=1.0e-9,
                )
                for key in (
                    "phase_count",
                    "me_kj_mol_per_H2O",
                    "mae_kj_mol_per_H2O",
                    "rmse_kj_mol_per_H2O",
                    "maxae_kj_mol_per_H2O",
                )
            )

    author_max_by_mesh = {}
    for mesh in (2, 3):
        values = [
            (
                abs(float(row["author_minus_current_kj_mol_per_water"])),
                row["phase"],
                float(row["author_minus_current_kj_mol_per_water"]),
            )
            for row in author_rows
            if int(row["mesh_n"]) == mesh
        ]
        magnitude, phase, signed = max(values)
        author_max_by_mesh[str(mesh)] = {
            "maximum_absolute_shift_kj_mol_per_H2O": magnitude,
            "phase": phase,
            "signed_shift_kj_mol_per_H2O": signed,
        }

    mstore_max_by_mesh = {}
    mstore_mae_by_mesh = {}
    current_mae_by_mesh = {}
    for mesh in (2, 3):
        subset = [row for row in mstore_rows if int(row["mesh_n"]) == mesh]
        magnitude, phase, signed = max(
            (
                abs(float(row["mstore_minus_pbc_kj_mol_per_H2O"])),
                row["phase"],
                float(row["mstore_minus_pbc_kj_mol_per_H2O"]),
            )
            for row in subset
        )
        mstore_max_by_mesh[str(mesh)] = {
            "maximum_absolute_shift_kj_mol_per_H2O": magnitude,
            "phase": phase,
            "signed_shift_kj_mol_per_H2O": signed,
        }
        mstore_mae_by_mesh[str(mesh)] = sum(
            abs(float(row["mstore_absolute_error_kj_mol_per_H2O"])) for row in subset
        ) / len(subset)
        current_mae_by_mesh[str(mesh)] = sum(
            abs(float(row["pbc_absolute_error_kj_mol_per_H2O"])) for row in subset
        ) / len(subset)

    adaptive_errors = [float(row["error_kj_mol_per_water"]) for row in adaptive_rows]
    adaptive_mae = sum(abs(value) for value in adaptive_errors) / len(adaptive_errors)
    adaptive_converged = sum(
        row["phase_converged"].strip().lower() == "true" for row in adaptive_rows
    )

    full_parity = load_json("validation/native_cli_full_parity_20260720/verification.json")
    provider = load_json("validation/provider_component_attribution_20260719/verification.json")
    h0 = load_json("validation/pbc_h0_anisotropy_attribution_20260719/verification.json")
    mstore_accuracy = load_json("validation/mstore_accuracy_equivalence_20260720/verification.json")

    checks = {
        "complete_cli_native_matrix": len(parity_rows) == 52
        and {int(row["mesh_n"]) for row in parity_rows} == {1, 2, 3, 4}
        and all(row["cli_accuracy_qualification"] == "PASS" for row in parity_rows),
        "cli_native_absolute_max_reproduces": math.isclose(
            parity_max[0],
            float(full_parity["maximum_absolute_difference_Ha_per_primitive"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "cli_native_relative_max_reproduces": math.isclose(
            relative_max[0],
            float(full_parity["maximum_relative_difference_kJ_mol_per_H2O"]),
            rel_tol=0.0,
            # The all-route table is written with twelve decimal places, so
            # differencing two rows can accumulate two final-place units.
            abs_tol=2.0e-11,
        ),
        "all_branch_statistics_reproduce": statistics_match,
        "same_accuracy_branch_comparison": all(
            row["same_cli_accuracy"].strip().lower() == "true" for row in mstore_rows
        ),
        "mstore_accuracy_sensitivity_is_numerical": float(
            mstore_accuracy["maximum_accuracy_sensitivity_hartree_supercell"]
        ) < 1.0e-10,
        "pbc_revision_shift_exceeds_interface_residual": author_max_by_mesh["3"][
            "maximum_absolute_shift_kj_mol_per_H2O"
        ]
        > 1000.0 * relative_max[0],
        "mstore_shift_exceeds_pbc_revision_shift": mstore_max_by_mesh["3"][
            "maximum_absolute_shift_kj_mol_per_H2O"
        ]
        > 100.0
        * author_max_by_mesh["3"]["maximum_absolute_shift_kj_mol_per_H2O"],
        "provider_gap_is_h0_attributed": float(h0["fraction_of_provider_gap_accounted_for"])
        > 0.999999,
        "exchange_ablation_collapses_provider_gap": float(
            provider["gap_reduction_when_exchange_disabled_percent"]
        )
        > 98.0,
        "adaptive_statistics_reproduce": len(adaptive_rows) == 12
        and adaptive_converged == int(adaptive_declared["converged_phase_count"])
        and math.isclose(
            adaptive_mae,
            float(adaptive_declared["mae_kj_mol_per_water"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "adaptive_result_remains_provisional": adaptive_declared["final_result"].strip().lower()
        == "false",
    }
    passed = all(checks.values())
    payload = {
        "schema": "periodic-gxtb-dmc13-discrepancy-attribution-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "cp2k_native_vs_current_pbc": {
            "maximum_absolute_difference_Ha_per_primitive": parity_max[0],
            "maximum_absolute_difference_mesh": int(parity_max[1]["mesh_n"]),
            "maximum_absolute_difference_phase": parity_max[1]["phase"],
            "maximum_relative_difference_kj_mol_per_H2O": relative_max[0],
            "maximum_relative_difference_mesh": relative_max[1],
            "maximum_relative_difference_phase": relative_max[2],
        },
        "author_pbc_vs_current_pbc": author_max_by_mesh,
        "mstore_inorganic_vs_current_pbc": {
            "maximum_shift_by_mesh": mstore_max_by_mesh,
            "mstore_sparse_mesh_mae_kj_mol_per_H2O": mstore_mae_by_mesh,
            "current_pbc_sparse_mesh_mae_kj_mol_per_H2O": current_mae_by_mesh,
            "maximum_accuracy_sensitivity_Ha_supercell": mstore_accuracy[
                "maximum_accuracy_sensitivity_hartree_supercell"
            ],
        },
        "provider_path_attribution": {
            "fraction_of_provider_gap_accounted_for_by_h0": h0[
                "fraction_of_provider_gap_accounted_for"
            ],
            "gap_reduction_when_exchange_disabled_percent": provider[
                "gap_reduction_when_exchange_disabled_percent"
            ],
        },
        "current_adaptive_dmc13": {
            "mae_kj_mol_per_H2O": adaptive_mae,
            "converged_phase_count": adaptive_converged,
            "phase_count": len(adaptive_rows),
            "final": False,
        },
        "source_files": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                parity_path,
                all_branch_path,
                branch_statistics_path,
                author_path,
                mstore_path,
                adaptive_path,
                adaptive_statistics_path,
            )
        },
        "classification": (
            "The complete current-provider CLI/native residual is numerical and far "
            "too small to explain a kJ/mol-scale DMC difference. The pbc source "
            "snapshots are measurably distinct, while mstore-inorganic is a much "
            "larger provider-model change dominated in the tested ablation by the "
            "periodic H0/exchange path. The lower previously quoted author result "
            "still requires exact source and executable provenance; the ongoing "
            "dense-mesh adaptive pbc result remains provisional."
        ),
    }
    output = HERE / "verification.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
