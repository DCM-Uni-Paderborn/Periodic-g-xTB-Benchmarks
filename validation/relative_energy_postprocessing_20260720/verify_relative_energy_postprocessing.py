#!/usr/bin/env python3
"""Rebuild DMC-ICE13 relative energies from raw outputs with Decimal."""

from __future__ import annotations

import csv
from decimal import Decimal, getcontext
import hashlib
import json
import re
from pathlib import Path


getcontext().prec = 50
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
TABLES = PACKAGE / "tables"
RAW = PACKAGE / "raw"
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
NONREFERENCE = PHASES[1:]
HARTREE_TO_KJMOL = Decimal("2625.4996394798254")
TABLE_TOLERANCE_KJMOL = Decimal("1e-9")
ENERGY_PATTERN = re.compile(
    r"ENERGY\| Total FORCE_EVAL .*?([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s*$"
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cp2k_energy(path: Path) -> Decimal:
    value: Decimal | None = None
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ENERGY_PATTERN.search(line)
        if match:
            value = Decimal(match.group(1))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if value is None or not ended:
        raise RuntimeError(f"incomplete raw CP2K output: {path}")
    return value


def primitive_counts(phase: str) -> tuple[int, int]:
    path = PACKAGE / "structures/primitive" / phase / "POSCAR"
    lines = path.read_text(encoding="utf-8").splitlines()
    symbols = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    by_symbol = dict(zip(symbols, counts, strict=True))
    waters = by_symbol.get("O", 0)
    atoms = sum(counts)
    if waters <= 0 or by_symbol.get("H", 0) != 2 * waters or atoms != 3 * waters:
        raise AssertionError(f"not a pure H2O primitive structure: {phase}")
    return atoms, waters


def poscar_atom_count(path: Path) -> int:
    lines = path.read_text(encoding="utf-8").splitlines()
    return sum(int(value) for value in lines[6].split())


def relative(energies: dict[str, Decimal], waters: dict[str, int]) -> dict[str, Decimal]:
    ih = energies["Ih"] / Decimal(waters["Ih"])
    return {
        phase: (energies[phase] / Decimal(waters[phase]) - ih) * HARTREE_TO_KJMOL
        for phase in NONREFERENCE
    }


def statistics(values: list[Decimal]) -> dict[str, Decimal]:
    count = Decimal(len(values))
    mean = sum(values) / count
    mae = sum(abs(value) for value in values) / count
    rmse = (sum(value * value for value in values) / count).sqrt()
    return {
        "me_kj_mol_per_H2O": mean,
        "mae_kj_mol_per_H2O": mae,
        "rmse_kj_mol_per_H2O": rmse,
        "maxae_kj_mol_per_H2O": max(abs(value) for value in values),
    }


def main() -> None:
    primitive = {phase: primitive_counts(phase) for phase in PHASES}
    waters = {phase: value[1] for phase, value in primitive.items()}
    references = {
        row["phase"]: Decimal(row["reference_relative_energy_kJmol_per_H2O"])
        for row in rows(TABLES / "dmc_reference_relative_energies.csv")
    }
    if set(references) != set(PHASES):
        raise AssertionError("incomplete DMC reference table")

    generated_relative = {
        (int(row["mesh_n"]), row["phase"]): row
        for row in rows(TABLES / "cp2k_native_relative_energies_by_mesh.csv")
        if row["qualification"] == "PASS"
    }
    generated_statistics = {
        (row["method"], int(row["mesh_n"])): row
        for row in rows(TABLES / "branch_comparison_statistics.csv")
    }

    native_relative: dict[tuple[int, str], Decimal] = {}
    native_table_residuals: list[Decimal] = []
    native_statistics_residuals: list[Decimal] = []
    for mesh in range(1, 6):
        energies = {
            phase: cp2k_energy(
                RAW / "cp2k_native" / f"k{mesh}{mesh}{mesh}-reduced" / phase / "cp2k.out"
            )
            for phase in PHASES
        }
        rel = relative(energies, waters)
        errors: list[Decimal] = []
        for phase, value in rel.items():
            native_relative[(mesh, phase)] = value
            table = generated_relative[(mesh, phase)]
            native_table_residuals.extend(
                [
                    abs(value - Decimal(table["relative_energy_kj_mol_per_H2O"])),
                    abs((value - references[phase]) - Decimal(table["error_kj_mol_per_H2O"])),
                    abs(abs(value - references[phase]) - Decimal(table["absolute_error_kj_mol_per_H2O"])),
                ]
            )
            errors.append(value - references[phase])
        calculated = statistics(errors)
        table_stats = generated_statistics[("CP2K-native pbc provider", mesh)]
        for key, value in calculated.items():
            native_statistics_residuals.append(abs(value - Decimal(table_stats[key])))

    all_branch = {
        (row["method"], int(row["mesh_n"]), row["phase"]): row
        for row in rows(TABLES / "all_branch_relative_energy_comparison.csv")
        if row["qualification"] == "PASS"
    }
    complete_cli_meshes: list[int] = []
    incomplete_cli_meshes: list[int] = []
    cli_table_residuals: list[Decimal] = []
    cli_native_relative_differences: list[Decimal] = []
    cli_statistics_residuals: list[Decimal] = []
    for mesh in range(1, 5):
        runs = [RAW / "current_pbc_cli" / f"cli-k{mesh}{mesh}{mesh}" / phase for phase in PHASES]
        complete = all(
            (run / "tblite.json").is_file()
            and (run / "exit_status").is_file()
            and (run / "exit_status").read_text(encoding="utf-8").strip() == "0"
            for run in runs
        )
        if not complete:
            incomplete_cli_meshes.append(mesh)
            continue
        complete_cli_meshes.append(mesh)
        energies: dict[str, Decimal] = {}
        for phase, run in zip(PHASES, runs, strict=True):
            payload = json.loads(
                (run / "tblite.json").read_text(encoding="utf-8"),
                parse_float=Decimal,
            )
            poscar = run / "POSCAR"
            expected_atoms = primitive[phase][0] * mesh**3
            if poscar_atom_count(poscar) != expected_atoms:
                raise AssertionError(f"wrong BvK atom count: mesh={mesh} phase={phase}")
            energies[phase] = payload["energy"] / Decimal(mesh**3)
        rel = relative(energies, waters)
        errors: list[Decimal] = []
        for phase, value in rel.items():
            table = all_branch[("current pbc CLI", mesh, phase)]
            cli_table_residuals.extend(
                [
                    abs(value - Decimal(table["relative_energy_kj_mol_per_H2O"])),
                    abs((value - references[phase]) - Decimal(table["error_kj_mol_per_H2O"])),
                    abs(abs(value - references[phase]) - Decimal(table["absolute_error_kj_mol_per_H2O"])),
                ]
            )
            cli_native_relative_differences.append(abs(value - native_relative[(mesh, phase)]))
            errors.append(value - references[phase])
        calculated = statistics(errors)
        table_stats = generated_statistics[("current pbc CLI", mesh)]
        for key, value in calculated.items():
            cli_statistics_residuals.append(abs(value - Decimal(table_stats[key])))

    checks = {
        "all_primitive_water_counts_valid": len(primitive) == len(PHASES),
        "native_gamma_through_k555_complete": len(native_relative) == 5 * len(NONREFERENCE),
        "current_cli_complete_through_k444": all(mesh in complete_cli_meshes for mesh in (1, 2, 3, 4)),
        "no_incomplete_cli_mesh": not incomplete_cli_meshes,
        "native_relative_tables_reproduce": max(native_table_residuals) <= TABLE_TOLERANCE_KJMOL,
        "native_statistics_reproduce": max(native_statistics_residuals) <= TABLE_TOLERANCE_KJMOL,
        "cli_relative_tables_reproduce": max(cli_table_residuals) <= TABLE_TOLERANCE_KJMOL,
        "cli_statistics_reproduce": max(cli_statistics_residuals) <= TABLE_TOLERANCE_KJMOL,
    }
    passed = all(checks.values())
    output = {
        "schema": "periodic-gxtb-relative-energy-postprocessing-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "hartree_to_kjmol": str(HARTREE_TO_KJMOL),
        "complete_current_cli_meshes": complete_cli_meshes,
        "incomplete_current_cli_meshes": incomplete_cli_meshes,
        "maximum_native_table_reproduction_residual_kj_mol_per_H2O": str(max(native_table_residuals)),
        "maximum_native_statistics_reproduction_residual_kj_mol_per_H2O": str(max(native_statistics_residuals)),
        "maximum_cli_table_reproduction_residual_kj_mol_per_H2O": str(max(cli_table_residuals)),
        "maximum_cli_statistics_reproduction_residual_kj_mol_per_H2O": str(max(cli_statistics_residuals)),
        "maximum_current_cli_minus_native_relative_energy_kj_mol_per_H2O": str(max(cli_native_relative_differences)),
        "source_files": {
            "dmc_references_sha256": sha256(TABLES / "dmc_reference_relative_energies.csv"),
            "generated_native_relative_sha256": sha256(TABLES / "cp2k_native_relative_energies_by_mesh.csv"),
            "generated_statistics_sha256": sha256(TABLES / "branch_comparison_statistics.csv"),
        },
    }
    (HERE / "verification.json").write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
