#!/usr/bin/env python3
"""Close the author-pbc/current-CLI/CP2K-native 3^3 energy triangle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
PACKAGE = REPOSITORY / "DMC-ICE13/reproduction/save_tblite_direct_dmc13"
PROVIDER = PACKAGE / "validation/provider_revision_bvk_ab_20260718"
ABSOLUTE = PROVIDER / "full_k333_absolute_energy_comparison.csv"
NATIVE = PACKAGE / "results/current_cp2k_native/k333"
STRUCTURES = PACKAGE / "structures/k333"
REFERENCES = PACKAGE / "tables/dmc_reference_relative_energies.csv"
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
NONREFERENCE = PHASES[1:]
HARTREE_TO_KJMOL = 2625.4996394799
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)


def run_gate(path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, str(path)],
        cwd=path.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise AssertionError(
            f"prerequisite gate failed: {path}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def cp2k_energy(path: Path) -> float:
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


def primitive_water_count(poscar: Path) -> int:
    lines = [line.strip() for line in poscar.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 7:
        raise AssertionError(f"incomplete POSCAR: {poscar}")
    total_atoms = sum(int(value) for value in lines[6].split())
    divisor = 3 * 27
    if total_atoms % divisor:
        raise AssertionError(f"nonintegral primitive water count: {poscar}")
    return total_atoms // divisor


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_manifest() -> None:
    for line in (HERE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split(maxsplit=1)
        path = HERE / relative.strip()
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"SHA-256 mismatch: {relative.strip()}")


def main() -> None:
    run_gate(PROVIDER / "verify_comparison.py")
    run_gate(PACKAGE / "tools/verify_absolute_energy_parity.py")

    with ABSOLUTE.open(newline="", encoding="utf-8") as handle:
        provider_rows = {row["phase"]: row for row in csv.DictReader(handle)}
    with REFERENCES.open(newline="", encoding="utf-8") as handle:
        references = {
            row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
            for row in csv.DictReader(handle)
        }
    if set(provider_rows) != set(PHASES) or set(references) != set(PHASES):
        raise AssertionError("phase coverage differs across the three-route inputs")

    absolute_rows: list[dict[str, object]] = []
    energies: dict[str, dict[str, float]] = {}
    waters: dict[str, int] = {}
    for phase in PHASES:
        provider = provider_rows[phase]
        current = float(provider["current_Ha_per_primitive"])
        author = float(provider["author_Ha_per_primitive"])
        native = cp2k_energy(NATIVE / phase / "cp2k.out")
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
    verify_manifest()
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
