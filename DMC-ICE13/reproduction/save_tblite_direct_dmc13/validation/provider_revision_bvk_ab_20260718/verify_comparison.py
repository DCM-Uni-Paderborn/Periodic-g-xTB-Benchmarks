#!/usr/bin/env python3
"""Recompute the archived direct-provider BvK comparison from raw JSON."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
HARTREE_TO_KJMOL = 2625.499638
N_WATER = {"Ih": 12, "VII": 12, "XVII": 6}


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

print("provider-revision BvK comparison: all archived values verified")
