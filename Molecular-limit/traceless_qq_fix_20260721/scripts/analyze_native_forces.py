#!/usr/bin/env python3
"""Build the CP2K-native H2O molecular-limit energy/force/stress table."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NATIVE = ROOT / "native_forces_20260721"
STRESS = ROOT / "limit_fix_requested_geometry_20260721" / "stress" / "raw"
OUTPUT = NATIVE / "analysis"
BOXES = (8, 10, 12, 15, 20, 30, 40, 50, 60, 80, 100)
EH_TO_KJMOL = 2625.4996394799
BOHR_PER_ANGSTROM = 1.88972613288564
AU_PRESSURE_TO_GPA = 2.94210107994716e4


def read_provenance(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def one(pattern: str, text: str, flags: int = 0) -> re.Match[str]:
    match = re.search(pattern, text, flags)
    if match is None:
        raise RuntimeError(f"Pattern not found: {pattern}")
    return match


def force_components(text: str) -> list[float]:
    blocks = re.findall(
        r"\s*FORCES\| Atomic forces \[hartree/bohr\]\s*\n"
        r"\s*FORCES\|\s+Atom\s+x\s+y\s+z\s+\|f\|\s*\n"
        r"((?:\s*FORCES\|\s+\d+[^\n]*\n)+)",
        text,
    )
    if not blocks:
        raise RuntimeError("Missing CP2K atomic-force block")
    result: list[float] = []
    for line in blocks[-1].splitlines():
        values = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", line)
        if len(values) != 4:
            raise RuntimeError(f"Malformed CP2K force row: {line}")
        result.extend(float(value) for value in values[:3])
    if len(result) != 9:
        raise RuntimeError(f"Expected 9 force components, found {len(result)}")
    return result


def matrix_after(label: str, text: str) -> list[float]:
    matches = re.findall(
        rf"{re.escape(label)}\s*\n((?:\s*DEBUG\|[^\n]*\n){{3}})", text
    )
    if not matches:
        raise RuntimeError(f"Missing matrix: {label}")
    values = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", matches[-1])
    if len(values) != 9:
        raise RuntimeError(f"Expected 9 matrix entries, found {len(values)}")
    return [float(value) for value in values]


def validate_case(directory: Path, cp2k_hash: str) -> str:
    output = (directory / "cp2k.out").read_text(errors="replace")
    if "PROGRAM ENDED AT" not in output:
        raise RuntimeError(f"No normal termination: {directory.name}")
    if (directory / "exit_status").read_text().strip() != "0":
        raise RuntimeError(f"Nonzero exit status: {directory.name}")
    if (directory / "qualification_status").read_text().strip() != "qualified":
        raise RuntimeError(f"Not qualified: {directory.name}")
    if (directory / "binary.sha256").read_text().split()[0] != cp2k_hash:
        raise RuntimeError(f"Wrong CP2K hash: {directory.name}")
    proof = (directory / "affinity_preexec.txt").read_text()
    if "expected_cpu=42 allowed=42" not in proof:
        raise RuntimeError(f"Wrong affinity proof: {directory.name}")
    return output


def load_native(case: str, box: int | None, cp2k_hash: str) -> dict[str, object]:
    output = validate_case(NATIVE / "raw" / case, cp2k_hash)
    energy = float(
        one(
            r"ENERGY\| Total FORCE_EVAL .*?([-+]?\d+\.\d+)\s*$",
            output,
            re.MULTILINE,
        ).group(1)
    )
    return {
        "case": case,
        "periodicity": "0D" if box is None else "3D",
        "L_angstrom": box,
        "cp2k_energy_Eh": energy,
        "forces_Eh_per_bohr": force_components(output),
    }


def load_stress(box: int, cp2k_hash: str) -> tuple[float, float]:
    case = f"H2O_stress_L{box:02d}"
    output = validate_case(STRESS / case, cp2k_hash)
    analytical = matrix_after("DEBUG| Analytical pv_virial [a.u.]", output)
    numerical = matrix_after("DEBUG| Numerical pv_virial [a.u.]", output)
    volume_bohr3 = (box * BOHR_PER_ANGSTROM) ** 3
    stress_gpa = [
        value / volume_bohr3 * AU_PRESSURE_TO_GPA for value in analytical
    ]
    return max(abs(value) for value in stress_gpa), max(
        abs(a - n) for a, n in zip(analytical, numerical, strict=True)
    )


def main() -> None:
    provenance = read_provenance(NATIVE / "provenance.env")
    cp2k_hash = provenance["cp2k_sha256"]
    if provenance.get("reference_cli") != "disabled":
        raise RuntimeError("The CP2K-native series unexpectedly enabled REFERENCE_CLI")

    rows = [load_native("H2O_0D", None, cp2k_hash)]
    rows.extend(load_native(f"H2O_L{box:02d}", box, cp2k_hash) for box in BOXES)
    reference = rows[0]
    reference_energy = float(reference["cp2k_energy_Eh"])
    reference_forces = list(reference["forces_Eh_per_bohr"])
    for row in rows:
        delta_energy = float(row["cp2k_energy_Eh"]) - reference_energy
        delta_force = [
            current - ref
            for current, ref in zip(
                list(row["forces_Eh_per_bohr"]), reference_forces, strict=True
            )
        ]
        max_force_au = max(abs(value) for value in delta_force)
        row["signed_delta_E_kJ_per_mol"] = delta_energy * EH_TO_KJMOL
        row["abs_delta_E_kJ_per_mol"] = abs(delta_energy * EH_TO_KJMOL)
        row["max_component_delta_F_Eh_per_bohr"] = max_force_au
        row["max_component_delta_F_Eh_per_angstrom"] = (
            max_force_au * BOHR_PER_ANGSTROM
        )
        if row["L_angstrom"] is None:
            row["max_abs_stress_GPa"] = None
            row["max_abs_analytical_minus_numerical_virial_Eh"] = None
        else:
            stress, virial_difference = load_stress(
                int(row["L_angstrom"]), cp2k_hash
            )
            row["max_abs_stress_GPa"] = stress
            row["max_abs_analytical_minus_numerical_virial_Eh"] = virial_difference
        row["cp2k_sha256"] = cp2k_hash
        row["normal_termination"] = True
        row["qualified"] = True
        del row["forces_Eh_per_bohr"]

    OUTPUT.mkdir(parents=True, exist_ok=True)
    table = OUTPUT / "h2o_molecular_limit_cp2k_native_0d_8_100.csv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    periodic = rows[1:]
    l100 = rows[-1]
    summary = {
        "all_normally_terminated_qualified_same_build": True,
        "boxes_angstrom": list(BOXES),
        "cp2k_sha256": cp2k_hash,
        "libcp2k_sha256": provenance["libcp2k_sha256"],
        "libtblite_sha256": provenance["libtblite_sha256"],
        "save_tblite_source_commit": provenance["save_tblite_source_commit"],
        "geometry_unit": "bohr",
        "force_native_unit": "Eh/bohr",
        "force_display_unit": "Eh/bohr",
        "stress_display_unit": "GPa",
        "molecular_0D_cp2k_energy_Eh": reference_energy,
        "L100_cp2k_energy_Eh": float(l100["cp2k_energy_Eh"]),
        "L100_minus_0D_energy_Eh": (
            float(l100["cp2k_energy_Eh"]) - reference_energy
        ),
        "L100_minus_0D_energy_kJ_per_mol": float(
            l100["signed_delta_E_kJ_per_mol"]
        ),
        "L100_max_component_force_difference_Eh_per_bohr": float(
            l100["max_component_delta_F_Eh_per_bohr"]
        ),
        "L100_max_abs_stress_GPa": float(l100["max_abs_stress_GPa"]),
        "max_abs_analytical_minus_numerical_virial_Eh": max(
            float(row["max_abs_analytical_minus_numerical_virial_Eh"])
            for row in periodic
        ),
        "reference_cli": "disabled",
    }
    (OUTPUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
