#!/usr/bin/env python3
"""Rebuild the qualified 8--50 A energy/force/stress series for Fig. 2."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "current_build_20260721"
BOXES = (8, 10, 12, 15, 20, 30, 40, 50)
EXPECTED_CP2K_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
BOHR_PER_ANGSTROM = 1.88972613288564
AU_PRESSURE_TO_MPA = 2.94210107994716e7


def matrix_after(label: str, text: str) -> list[float]:
    matches = re.findall(
        rf"{re.escape(label)}\s*\n((?:\s*DEBUG\|[^\n]*\n){{3}})", text
    )
    if not matches:
        raise RuntimeError(f"Missing matrix: {label}")
    values = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", matches[-1])
    if len(values) != 9:
        raise RuntimeError(f"Expected 9 values for {label}, found {len(values)}")
    return [float(value) for value in values]


def main() -> None:
    with (BUILD / "results_current.csv").open(newline="") as handle:
        energy_force_rows = {
            int(float(row["L_angstrom"])): row
            for row in csv.DictReader(handle)
            if row["group"] == "main" and row["L_angstrom"]
        }

    rows: list[dict[str, object]] = []
    for box in BOXES:
        case = f"H2O_stress_L{box:02d}"
        directory = BUILD / "raw_stress" / case
        output = (directory / "cp2k.out").read_text()
        if "PROGRAM ENDED AT" not in output:
            raise RuntimeError(f"Run did not terminate normally: {case}")
        if (directory / "qualification_status").read_text().strip() != "qualified":
            raise RuntimeError(f"Run is not qualified: {case}")
        binary_hash = (directory / "binary.sha256").read_text().split()[0]
        if binary_hash != EXPECTED_CP2K_SHA256:
            raise RuntimeError(f"Wrong CP2K binary for {case}: {binary_hash}")

        analytical = matrix_after("DEBUG| Analytical pv_virial [a.u.]", output)
        numerical = matrix_after("DEBUG| Numerical pv_virial [a.u.]", output)
        difference = [a - n for a, n in zip(analytical, numerical, strict=True)]
        volume_bohr3 = (box * BOHR_PER_ANGSTROM) ** 3
        stress_components_mpa = [
            value / volume_bohr3 * AU_PRESSURE_TO_MPA for value in analytical
        ]
        base = energy_force_rows[box]
        rows.append(
            {
                "L_angstrom": box,
                "abs_delta_E_kJ_per_mol": abs(
                    float(base["delta_E_vs_matching_0D_kJ_per_mol"])
                ),
                "signed_delta_E_kJ_per_mol": float(
                    base["delta_E_vs_matching_0D_kJ_per_mol"]
                ),
                "max_component_delta_F_eV_per_angstrom": float(
                    base["max_component_force_delta_vs_matching_0D_eV_per_angstrom"]
                ),
                "max_abs_stress_MPa": max(abs(value) for value in stress_components_mpa),
                "max_abs_analytical_virial_Eh": max(abs(value) for value in analytical),
                "max_abs_analytical_minus_numerical_virial_Eh": max(
                    abs(value) for value in difference
                ),
                "cp2k_sha256": binary_hash,
                "normal_termination": True,
                "qualified": True,
            }
        )

    with (ROOT / "results_energy_force_stress_8_50.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "boxes_angstrom": list(BOXES),
        "case_count": len(rows),
        "all_normally_terminated_qualified_same_build": True,
        "cp2k_sha256": EXPECTED_CP2K_SHA256,
        "max_abs_analytical_minus_numerical_virial_Eh": max(
            float(row["max_abs_analytical_minus_numerical_virial_Eh"])
            for row in rows
        ),
        "stress_metric": "maximum absolute analytical stress-tensor component",
        "stress_unit": "MPa",
    }
    (BUILD / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
