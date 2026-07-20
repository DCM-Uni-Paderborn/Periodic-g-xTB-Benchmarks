#!/usr/bin/env python3
"""Close the author-pbc/current-CLI/CP2K-native 3^3 energy triangle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
PACKAGE = REPOSITORY / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
AUTHOR_ABSOLUTE = PACKAGE / "tables/author_pbc_absolute_energies.csv"
CURRENT_ABSOLUTE = PACKAGE / "tables/current_absolute_energies_by_mesh.csv"
NATIVE = PACKAGE / "raw/cp2k_native/k333-reduced"
CURRENT_CLI = PACKAGE / "raw/current_pbc_cli/cli-k333"
STRUCTURES = PACKAGE / "structures/primitive"
REFERENCES = PACKAGE / "tables/dmc_reference_relative_energies.csv"
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
NONREFERENCE = PHASES[1:]
HARTREE_TO_KJMOL = 2625.4996394799
QUALIFIED_CP2K_SHA256 = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
QUALIFIED_PBC_CLI_SHA256 = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)


def recorded_digest(path: Path) -> str:
    fields = path.read_text(encoding="utf-8").split()
    if not fields or not re.fullmatch(r"[0-9a-fA-F]{64}", fields[0]):
        raise AssertionError(f"invalid SHA-256 sidecar: {path}")
    return fields[0].lower()


def cp2k_energy(run: Path) -> float:
    path = run / "cp2k.out"
    if (run / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise AssertionError(f"nonzero native exit status: {run}")
    if recorded_digest(run / "binary.sha256") != QUALIFIED_CP2K_SHA256:
        raise AssertionError(f"unqualified native executable: {run}")
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values or not math.isfinite(values[-1]):
        raise AssertionError(f"incomplete native output: {path}")
    return values[-1]


def current_cli_energy(run: Path) -> float:
    if (run / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise AssertionError(f"nonzero CLI exit status: {run}")
    if recorded_digest(run / "binary.sha256") != QUALIFIED_PBC_CLI_SHA256:
        raise AssertionError(f"unqualified CLI executable: {run}")
    poscar = run / "POSCAR"
    if recorded_digest(run / "input.sha256") != hashlib.sha256(poscar.read_bytes()).hexdigest():
        raise AssertionError(f"CLI input hash mismatch: {run}")
    result = json.loads((run / "tblite.json").read_text(encoding="utf-8"))
    energy = float(result["energy"])
    if not math.isfinite(energy):
        raise AssertionError(f"non-finite CLI energy: {run}")
    return energy / 27.0


def primitive_water_count(poscar: Path) -> int:
    lines = [line.strip() for line in poscar.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 7:
        raise AssertionError(f"incomplete POSCAR: {poscar}")
    total_atoms = sum(int(value) for value in lines[6].split())
    divisor = 3
    if total_atoms % divisor:
        raise AssertionError(f"nonintegral primitive water count: {poscar}")
    return total_atoms // divisor


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with AUTHOR_ABSOLUTE.open(newline="", encoding="utf-8") as handle:
        author_rows = {
            row["phase"]: row for row in csv.DictReader(handle)
            if int(row["mesh_n"]) == 3
        }
    with CURRENT_ABSOLUTE.open(newline="", encoding="utf-8") as handle:
        current_rows = {
            row["phase"]: row for row in csv.DictReader(handle)
            if int(row["mesh_n"]) == 3
        }
    with REFERENCES.open(newline="", encoding="utf-8") as handle:
        references = {
            row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
            for row in csv.DictReader(handle)
        }
    if (
        set(author_rows) != set(PHASES)
        or set(current_rows) != set(PHASES)
        or set(references) != set(PHASES)
    ):
        raise AssertionError("phase coverage differs across the three-route inputs")

    absolute_rows: list[dict[str, object]] = []
    energies: dict[str, dict[str, float]] = {}
    waters: dict[str, int] = {}
    for phase in PHASES:
        author = float(author_rows[phase]["author_pbc_Ha_per_primitive"])
        current = current_cli_energy(CURRENT_CLI / phase)
        recorded_current = float(current_rows[phase]["save_tblite_cli_energy_Ha_per_primitive"])
        if abs(current - recorded_current) > 5.0e-13:
            raise AssertionError(f"current CLI table/raw mismatch for {phase}")
        native = cp2k_energy(NATIVE / phase)
        if abs(native - current) > 2.0e-7:
            raise AssertionError(
                f"current CLI/native mismatch for {phase}: {native-current:+.6e} Ha"
            )
        waters[phase] = primitive_water_count(STRUCTURES / phase / "POSCAR")
        energies[phase] = {"author": author, "current": current, "native": native}
        absolute_rows.append(
            {
                "phase": phase,
                "water_count_primitive": waters[phase],
                "author_pbc_Ha": f"{author:.15f}",
                "current_cli_Ha": f"{current:.15f}",
                "cp2k_native_Ha": f"{native:.15f}",
                "native_minus_current_Ha": f"{native-current:+.15e}",
                "author_minus_current_Ha": f"{author-current:+.15e}",
                "author_minus_native_Ha": f"{author-native:+.15e}",
            }
        )

    relative_rows: list[dict[str, object]] = []
    errors: dict[str, list[float]] = {route: [] for route in ("author", "current", "native")}
    relative_values: dict[str, dict[str, float]] = {route: {} for route in errors}
    for phase in NONREFERENCE:
        row: dict[str, object] = {
            "phase": phase,
            "reference_kJ_mol_per_H2O": f"{references[phase]:.12f}",
        }
        for route in errors:
            relative = (
                energies[phase][route] / waters[phase]
                - energies["Ih"][route] / waters["Ih"]
            ) * HARTREE_TO_KJMOL
            error = relative - references[phase]
            relative_values[route][phase] = relative
            errors[route].append(error)
            row[f"{route}_relative_kJ_mol_per_H2O"] = f"{relative:.12f}"
            row[f"{route}_error_kJ_mol_per_H2O"] = f"{error:+.12f}"
        row["native_minus_current_relative_kJ_mol_per_H2O"] = (
            f"{relative_values['native'][phase]-relative_values['current'][phase]:+.12e}"
        )
        row["author_minus_current_relative_kJ_mol_per_H2O"] = (
            f"{relative_values['author'][phase]-relative_values['current'][phase]:+.12f}"
        )
        relative_rows.append(row)

    maximum_native_current_absolute = max(
        abs(float(row["native_minus_current_Ha"])) for row in absolute_rows
    )
    maximum_native_current_relative = max(
        abs(float(row["native_minus_current_relative_kJ_mol_per_H2O"]))
        for row in relative_rows
    )
    summary = {
        "mesh": "3x3x3",
        "phase_count_including_Ih": len(PHASES),
        "relative_phase_count": len(NONREFERENCE),
        "status": "PASS",
        "routes": {
            route: {
                "mae_kJ_mol_per_H2O": sum(map(abs, values)) / len(values),
                "rmse_kJ_mol_per_H2O": math.sqrt(
                    sum(value * value for value in values) / len(values)
                ),
            }
            for route, values in errors.items()
        },
        "maximum_absolute_native_minus_current_Ha": maximum_native_current_absolute,
        "maximum_absolute_native_minus_current_relative_kJ_mol_per_H2O": maximum_native_current_relative,
        "classification": (
            "current CLI and CP2K-native are energetically identical within the "
            "qualified numerical tolerance; the author-pbc/current shift is a "
            "provider-model revision, not a CP2K k-point integration error"
        ),
    }
    if maximum_native_current_relative > 5.0e-5:
        raise AssertionError("Ih-referenced CLI/native difference exceeds tolerance")

    write_csv(
        HERE / "absolute_energies.csv",
        tuple(absolute_rows[0]),
        absolute_rows,
    )
    write_csv(
        HERE / "relative_energies_and_errors.csv",
        tuple(relative_rows[0]),
        relative_rows,
    )
    (HERE / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
