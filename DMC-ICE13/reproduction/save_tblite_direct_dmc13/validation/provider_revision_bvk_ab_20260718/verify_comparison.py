#!/usr/bin/env python3
"""Recompute the archived direct-provider BvK comparison from raw JSON."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import statistics


ROOT = Path(__file__).resolve().parent
HARTREE_TO_KJMOL = 2625.499638
N_WATER = {
    "Ih": 12,
    "II": 12,
    "III": 12,
    "IV": 16,
    "VI": 10,
    "VII": 12,
    "VIII": 8,
    "IX": 12,
    "XI": 8,
    "XIII": 28,
    "XIV": 12,
    "XV": 10,
    "XVII": 6,
}


def energy(provider: str, mesh: int, phase: str) -> float:
    path = ROOT / provider / f"k{mesh}{mesh}{mesh}" / phase / "result.json"
    return float(json.loads(path.read_text())["energy"]) / mesh**3


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if abs(actual - expected) > tolerance:
        raise SystemExit(
            f"{label}: actual={actual:.15g} expected={expected:.15g} "
            f"difference={actual - expected:+.3e}"
        )


with (ROOT / "absolute_energy_comparison.csv").open(newline="") as handle:
    for row in csv.DictReader(handle):
        mesh = int(row["mesh_n"])
        phase = row["phase"]
        current = energy("current_save_tblite_cli", mesh, phase)
        seidler = energy("seidler_pbc_cli", mesh, phase)
        close(current, float(row["current_save_tblite_Ha_per_primitive"]), 5.0e-13,
              f"current {mesh} {phase}")
        close(seidler, float(row["seidler_pbc_Ha_per_primitive"]), 5.0e-13,
              f"seidler {mesh} {phase}")
        close(seidler - current,
              float(row["seidler_minus_current_Ha_per_primitive"]), 5.0e-14,
              f"difference {mesh} {phase}")

with (ROOT / "relative_energy_comparison.csv").open(newline="") as handle:
    for row in csv.DictReader(handle):
        mesh = int(row["mesh_n"])
        phase = row["phase"]
        values = {}
        for provider in ("current_save_tblite_cli", "seidler_pbc_cli"):
            values[provider] = (
                energy(provider, mesh, phase) / N_WATER[phase]
                - energy(provider, mesh, "Ih") / N_WATER["Ih"]
            ) * HARTREE_TO_KJMOL
        close(values["current_save_tblite_cli"],
              float(row["current_save_tblite_kj_mol_per_water"]), 5.0e-10,
              f"current relative {mesh} {phase}")
        close(values["seidler_pbc_cli"],
              float(row["seidler_pbc_kj_mol_per_water"]), 5.0e-10,
              f"seidler relative {mesh} {phase}")
        close(values["seidler_pbc_cli"] - values["current_save_tblite_cli"],
              float(row["seidler_minus_current_kj_mol_per_water"]), 5.0e-10,
              f"relative difference {mesh} {phase}")

current_errors = []
author_errors = []
relative_shifts = []
with (ROOT / "full_k222_relative_comparison.csv").open(newline="") as handle:
    for row in csv.DictReader(handle):
        phase = row["phase"]
        values = {}
        for provider in ("current_save_tblite_cli", "seidler_pbc_cli_linux"):
            values[provider] = (
                energy(provider, 2, phase) / N_WATER[phase]
                - energy(provider, 2, "Ih") / N_WATER["Ih"]
            ) * HARTREE_TO_KJMOL
        current = values["current_save_tblite_cli"]
        author = values["seidler_pbc_cli_linux"]
        reference = float(row["dmc_reference_kj_mol_per_water"])
        close(current, float(row["current_save_tblite_kj_mol_per_water"]),
              5.0e-10, f"full current relative 2 {phase}")
        close(author, float(row["author_pbc_kj_mol_per_water"]), 5.0e-10,
              f"full author relative 2 {phase}")
        close(author - current,
              float(row["author_minus_current_kj_mol_per_water"]), 5.0e-10,
              f"full relative difference 2 {phase}")
        close(current - reference,
              float(row["current_error_kj_mol_per_water"]), 5.0e-10,
              f"full current error 2 {phase}")
        close(author - reference,
              float(row["author_error_kj_mol_per_water"]), 5.0e-10,
              f"full author error 2 {phase}")
        current_errors.append(abs(current - reference))
        author_errors.append(abs(author - reference))
        relative_shifts.append((abs(author - current), phase))

summary = json.loads((ROOT / "full_k222_summary.json").read_text())
current_mae = statistics.mean(current_errors)
author_mae = statistics.mean(author_errors)
maximum_shift, maximum_phase = max(relative_shifts)
close(current_mae, summary["current_save_tblite_mae_kj_mol_per_water"],
      5.0e-10, "full current k222 MAE")
close(author_mae, summary["author_pbc_mae_kj_mol_per_water"],
      5.0e-10, "full author k222 MAE")
close(author_mae - current_mae,
      summary["author_minus_current_mae_kj_mol_per_water"], 5.0e-10,
      "full k222 MAE shift")
close(maximum_shift,
      summary["maximum_absolute_relative_energy_shift_kj_mol_per_water"],
      5.0e-10, "full k222 maximum relative-energy shift")
if maximum_phase != summary["maximum_shift_phase"]:
    raise SystemExit(
        f"maximum shift phase: actual={maximum_phase} "
        f"expected={summary['maximum_shift_phase']}"
    )

for phase in ("Ih", "VII", "XVII"):
    linux = energy("seidler_pbc_cli_linux", 2, phase) * 8
    macos = energy("seidler_pbc_cli", 2, phase) * 8
    close(linux, macos, 2.0e-11, f"Linux/macOS author pbc 2 {phase}")

print("provider-revision BvK comparison: all archived values verified")
