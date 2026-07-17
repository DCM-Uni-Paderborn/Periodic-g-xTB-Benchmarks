#!/usr/bin/env python3
"""Summarize the frozen H2O molecular-limit benchmark artifacts."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent / "raw"
OUT = Path(__file__).resolve().parent / "results.csv"
EH_TO_KJMOL = 2625.4996394799
EH_PER_BOHR_TO_EV_PER_ANG = 51.4220674763

CASES = [
    ("H2O_0D", None),
    ("H2O_L08", 8.0),
    ("H2O_L10", 10.0),
    ("H2O_L12", 12.0),
    ("H2O_L15", 15.0),
    ("H2O_L20", 20.0),
    ("H2O_L30", 30.0),
    ("H2O_L40", 40.0),
    ("H2O_L50", 50.0),
    ("H2O_L60", 60.0),
    ("H2O_L80", 80.0),
    ("H2O_L100", 100.0),
    ("H2O_L150", 150.0),
    ("H2O_L200", 200.0),
]


def one(pattern: str, text: str, flags: int = 0) -> re.Match[str]:
    match = re.search(pattern, text, flags)
    if match is None:
        raise RuntimeError(f"Pattern not found: {pattern}")
    return match


def load_case(name: str, box: float | None) -> dict[str, object]:
    directory = ROOT / name
    output = (directory / "cp2k.out").read_text()
    native = json.loads(next(directory.glob("*_cli.json")).read_text())

    cp2k_energy = float(
        one(r"ENERGY\| Total FORCE_EVAL .*?([-+]?\d+\.\d+)\s*$", output, re.MULTILINE).group(1)
    )
    parity = one(
        r"Energy CP2K/CLI/absdiff:\s+(\S+)\s+(\S+)\s+(\S+)", output
    )
    gradient_parity = one(r"Gradient diff sum/max:\s+(\S+)\s+(\S+)", output)
    steps = int(one(r"SCF run converged in\s+(\d+)\s+steps", output).group(1))

    return {
        "case": name,
        "periodicity": "0D" if box is None else "3D",
        "L_angstrom": box,
        "cp2k_energy_Eh": cp2k_energy,
        "cli_energy_Eh": float(native["energy"]),
        "cp2k_cli_absdiff_Eh": float(parity.group(3)),
        "cp2k_cli_gradient_sumdiff_Eh_per_bohr": float(gradient_parity.group(1)),
        "cp2k_cli_gradient_maxdiff_Eh_per_bohr": float(gradient_parity.group(2)),
        "scf_steps_cp2k": steps,
        "gradient": [float(value) for value in native["gradient"]],
    }


rows = [load_case(name, box) for name, box in CASES]
reference = rows[0]
ref_energy = float(reference["cp2k_energy_Eh"])
ref_cli_energy = float(reference["cli_energy_Eh"])
ref_gradient = list(reference["gradient"])

for row in rows:
    energy_delta = float(row["cp2k_energy_Eh"]) - ref_energy
    cli_energy_delta = float(row["cli_energy_Eh"]) - ref_cli_energy
    gradient_delta = [
        current - ref for current, ref in zip(row["gradient"], ref_gradient, strict=True)
    ]
    row["delta_E_vs_0D_Eh"] = energy_delta
    row["delta_E_vs_0D_kJ_per_mol"] = energy_delta * EH_TO_KJMOL
    row["delta_cli_E_vs_0D_Eh"] = cli_energy_delta
    row["max_component_gradient_delta_vs_0D_Eh_per_bohr"] = max(
        abs(value) for value in gradient_delta
    )
    row["max_component_force_delta_vs_0D_eV_per_angstrom"] = (
        row["max_component_gradient_delta_vs_0D_Eh_per_bohr"]
        * EH_PER_BOHR_TO_EV_PER_ANG
    )
    row["rms_component_gradient_delta_vs_0D_Eh_per_bohr"] = math.sqrt(
        sum(value * value for value in gradient_delta) / len(gradient_delta)
    )
    del row["gradient"]

fieldnames = list(rows[0])
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(
    "case periodicity L/A E_CP2K/Eh E_CLI/Eh dE(3D-0D)/kJmol "
    "CP2K-CLI/Eh max|dF(3D-0D)|/eVAng max|dgrad(CP2K-CLI)|/Eha0"
)
for row in rows:
    box = "--" if row["L_angstrom"] is None else f"{row['L_angstrom']:g}"
    print(
        f"{row['case']:9s} {row['periodicity']:>2s} {box:>4s} "
        f"{row['cp2k_energy_Eh']:.15f} {row['cli_energy_Eh']:.15f} "
        f"{row['delta_E_vs_0D_kJ_per_mol']:+.9f} "
        f"{row['cp2k_cli_absdiff_Eh']:.3e} "
        f"{row['max_component_force_delta_vs_0D_eV_per_angstrom']:.3e} "
        f"{row['cp2k_cli_gradient_maxdiff_Eh_per_bohr']:.3e}"
    )

