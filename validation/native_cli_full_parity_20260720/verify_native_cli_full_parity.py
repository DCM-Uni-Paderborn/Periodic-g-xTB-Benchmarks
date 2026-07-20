#!/usr/bin/env python3
"""Verify complete direct save_tblite CLI versus CP2K-native energy parity."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path


getcontext().prec = 50

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLE = (
    ROOT
    / "DMC-ICE13/reproduction/seidler_dmc13_recalculation/tables"
    / "pbc_cli_vs_cp2k_native_absolute_parity.csv"
)

EXPECTED_MESHES = (1, 2, 3, 4)
EXPECTED_PHASES = (
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
)
HARTREE_TO_KJMOL = Decimal("2625.4996394798254")
ABSOLUTE_TOLERANCE_HARTREE_PER_PRIMITIVE = Decimal("2e-7")
RELATIVE_TOLERANCE_KJMOL_PER_WATER = Decimal("5e-5")
REPORTED_DIFFERENCE_RESIDUAL_TOLERANCE_HARTREE = Decimal("5e-15")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with TABLE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    expected_keys = {
        (mesh, phase) for mesh in EXPECTED_MESHES for phase in EXPECTED_PHASES
    }
    observed_keys = {(int(row["mesh_n"]), row["phase"]) for row in rows}
    duplicate_free = len(observed_keys) == len(rows)

    water_counts: dict[str, int] = {}
    hashes_valid = True
    accuracy_valid = True
    qualifications_valid = True
    reported_difference_residuals: list[Decimal] = []
    absolute_differences: list[tuple[Decimal, int, str]] = []
    indexed: dict[tuple[int, str], dict[str, Decimal | int]] = {}

    for row in rows:
        mesh = int(row["mesh_n"])
        phase = row["phase"]
        waters = int(row["water_molecules_primitive"])
        previous_waters = water_counts.setdefault(phase, waters)
        if previous_waters != waters:
            raise ValueError(f"Inconsistent primitive water count for {phase}")

        native = Decimal(row["cp2k_native_energy_Ha_per_primitive"])
        cli = Decimal(row["pbc_cli_energy_Ha_per_primitive"])
        reported_difference = Decimal(row["native_minus_cli_Ha_per_primitive"])
        calculated_difference = native - cli
        residual = abs(calculated_difference - reported_difference)
        reported_difference_residuals.append(residual)
        absolute_differences.append((abs(calculated_difference), mesh, phase))
        indexed[(mesh, phase)] = {
            "native": native,
            "cli": cli,
            "waters": waters,
        }

        hashes_valid = hashes_valid and all(
            SHA256_RE.fullmatch(row[column]) is not None
            for column in (
                "cp2k_output_sha256",
                "pbc_cli_json_sha256",
                "poscar_sha256",
            )
        )
        accuracy_valid = accuracy_valid and Decimal(row["cli_accuracy"]) == Decimal("0.1")
        qualifications_valid = (
            qualifications_valid and row["cli_accuracy_qualification"] == "PASS"
        )

    relative_differences: list[tuple[Decimal, int, str]] = []
    per_mesh: dict[str, dict[str, str | int]] = {}
    for mesh in EXPECTED_MESHES:
        ih = indexed[(mesh, "Ih")]
        mesh_relative: list[tuple[Decimal, int, str]] = []
        for phase in EXPECTED_PHASES:
            point = indexed[(mesh, phase)]
            native_relative = (
                Decimal(point["native"]) / Decimal(point["waters"])
                - Decimal(ih["native"]) / Decimal(ih["waters"])
            ) * HARTREE_TO_KJMOL
            cli_relative = (
                Decimal(point["cli"]) / Decimal(point["waters"])
                - Decimal(ih["cli"]) / Decimal(ih["waters"])
            ) * HARTREE_TO_KJMOL
            item = (abs(native_relative - cli_relative), mesh, phase)
            relative_differences.append(item)
            mesh_relative.append(item)

        mesh_absolute = max(item for item in absolute_differences if item[1] == mesh)
        mesh_relative_max = max(mesh_relative)
        per_mesh[str(mesh)] = {
            "point_count": len(EXPECTED_PHASES),
            "maximum_absolute_difference_Ha_per_primitive": str(mesh_absolute[0]),
            "maximum_absolute_difference_phase": mesh_absolute[2],
            "maximum_relative_difference_kJ_mol_per_H2O": str(mesh_relative_max[0]),
            "maximum_relative_difference_phase": mesh_relative_max[2],
        }

    maximum_absolute = max(absolute_differences)
    maximum_relative = max(relative_differences)
    maximum_reported_residual = max(reported_difference_residuals, default=Decimal("0"))

    checks = {
        "exactly_52_unique_required_points": (
            duplicate_free and observed_keys == expected_keys and len(rows) == 52
        ),
        "all_three_provenance_hashes_are_sha256": hashes_valid,
        "all_cli_calculations_use_accuracy_0.1": accuracy_valid,
        "all_cli_accuracy_qualifications_pass": qualifications_valid,
        "reported_absolute_differences_reproduce": (
            maximum_reported_residual
            <= REPORTED_DIFFERENCE_RESIDUAL_TOLERANCE_HARTREE
        ),
        "all_absolute_energy_differences_within_tolerance": (
            maximum_absolute[0] <= ABSOLUTE_TOLERANCE_HARTREE_PER_PRIMITIVE
        ),
        "all_same_mesh_ih_referenced_differences_within_tolerance": (
            maximum_relative[0] <= RELATIVE_TOLERANCE_KJMOL_PER_WATER
        ),
    }
    passed = all(checks.values())
    output = {
        "schema": "periodic-gxtb-native-cli-full-parity-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "source_table": str(TABLE.relative_to(ROOT)),
        "source_table_sha256": sha256(TABLE),
        "point_count": len(rows),
        "meshes": list(EXPECTED_MESHES),
        "phase_count_per_mesh": len(EXPECTED_PHASES),
        "absolute_tolerance_Ha_per_primitive": str(
            ABSOLUTE_TOLERANCE_HARTREE_PER_PRIMITIVE
        ),
        "relative_tolerance_kJ_mol_per_H2O": str(
            RELATIVE_TOLERANCE_KJMOL_PER_WATER
        ),
        "maximum_absolute_difference_Ha_per_primitive": str(maximum_absolute[0]),
        "maximum_absolute_difference_mesh": maximum_absolute[1],
        "maximum_absolute_difference_phase": maximum_absolute[2],
        "maximum_relative_difference_kJ_mol_per_H2O": str(maximum_relative[0]),
        "maximum_relative_difference_mesh": maximum_relative[1],
        "maximum_relative_difference_phase": maximum_relative[2],
        "maximum_reported_difference_reproduction_residual_Ha": str(
            maximum_reported_residual
        ),
        "per_mesh": per_mesh,
        "interpretation": (
            "The complete same-provider direct save_tblite CLI and CP2K-native "
            "matrix agrees within the declared absolute and independently "
            "reconstructed same-mesh relative-energy tolerances."
        ),
    }
    rendered = json.dumps(output, indent=2, sort_keys=True) + "\n"
    (HERE / "verification.json").write_text(rendered, encoding="utf-8")
    (HERE / "verification.stdout").write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
