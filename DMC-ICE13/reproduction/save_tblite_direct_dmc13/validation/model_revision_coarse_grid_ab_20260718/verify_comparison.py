#!/usr/bin/env python3
"""Recompute the coarse-grid g-xTB model/revision diagnostic."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
    "XIII", "XIV", "XV", "XVII",
)
N_WATER = {
    "Ih": 12, "II": 12, "III": 12, "IV": 16, "VI": 10, "VII": 12,
    "VIII": 8, "IX": 12, "XI": 8, "XIII": 28, "XIV": 12, "XV": 10,
    "XVII": 6,
}
REFERENCE = {
    "II": 0.31, "III": 1.25, "IV": 3.83, "VI": 1.78, "VII": 4.99,
    "VIII": 4.23, "IX": 0.60, "XI": 0.16, "XIII": 2.12, "XIV": 1.70,
    "XV": 1.74, "XVII": 1.75,
}
HARTREE_TO_KJMOL = 2625.499638
EXPECTED_MAE = {
    (1, "current"): 163.834465930395,
    (1, "authors_exchange"): 130.899820086448,
    (1, "gxtb_v201"): 41.390727637315,
    (1, "dcm_main"): 162.626025116661,
    (2, "current"): 88.681375524804,
    (2, "authors_exchange"): 48.710763687881,
    (2, "gxtb_v201"): 37.201256501330,
    (2, "dcm_main"): 89.396043594685,
}
EXPECTED_EXECUTABLE_HASH = {
    "authors_exchange": "324c2c1e4968eab579fae1bd8571a467d62a8eaf372f2b88906bb0d9f7ba7549",
    "gxtb_v201": "c87471101170b506dae7f54700d5724aad9ce3dc5923e48d5317a4fd8f6cac60",
    "dcm_main": "2af03fdc70875df823038e49319f69751ae4a94dada58ce2960d09d358884bf0",
}


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(
            f"{label}: actual={actual:.15g}, expected={expected:.15g}, "
            f"tolerance={tolerance:.3g}"
        )


def result_path(provider: str, mesh: int, phase: str) -> Path:
    if provider == "current":
        return (
            PACKAGE / "results/current_save_tblite_cli"
            / f"k{mesh}{mesh}{mesh}" / phase / "result.json"
        )
    return ROOT / "raw" / provider / f"k{mesh}{mesh}{mesh}" / phase / "result.json"


def energy(provider: str, mesh: int, phase: str) -> float:
    path = result_path(provider, mesh, phase)
    if not path.is_file():
        raise AssertionError(f"missing result: {path}")
    value = float(json.loads(path.read_text())["energy"]) / mesh**3
    if not math.isfinite(value):
        raise AssertionError(f"non-finite energy: {path}")
    return value


def relative(provider: str, mesh: int, phase: str) -> float:
    return (
        energy(provider, mesh, phase) / N_WATER[phase]
        - energy(provider, mesh, "Ih") / N_WATER["Ih"]
    ) * HARTREE_TO_KJMOL


def executable_hash(provider: str) -> str:
    line = (ROOT / "raw" / provider / "executable.sha256").read_text().splitlines()[0]
    return line.split()[0]


providers = ("current", "authors_exchange", "gxtb_v201", "dcm_main")
computed: dict[tuple[int, str], float] = {}
for mesh in (1, 2):
    for provider in providers:
        errors = [
            abs(relative(provider, mesh, phase) - REFERENCE[phase])
            for phase in PHASES[1:]
        ]
        computed[(mesh, provider)] = sum(errors) / len(errors)
        close(
            computed[(mesh, provider)],
            EXPECTED_MAE[(mesh, provider)],
            5.0e-10,
            f"MAE mesh={mesh} provider={provider}",
        )

for provider, expected_hash in EXPECTED_EXECUTABLE_HASH.items():
    actual_hash = executable_hash(provider)
    if actual_hash != expected_hash:
        raise AssertionError(
            f"executable hash {provider}: actual={actual_hash}, expected={expected_hash}"
        )

with (ROOT / "coarse_mae_summary.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
if len(rows) != len(EXPECTED_MAE):
    raise AssertionError("coarse MAE table has the wrong row count")
for row in rows:
    key = (int(row["mesh"]), row["provider"])
    close(
        float(row["mae_kj_mol_per_water"]),
        computed[key],
        5.0e-10,
        f"summary row {key}",
    )

with (ROOT / "relative_energies_k222.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
if [row["phase"] for row in rows] != list(PHASES[1:]):
    raise AssertionError("relative-energy phase order is incomplete")
columns = {
    "current": "current_kj_mol_per_water",
    "authors_exchange": "authors_exchange_kj_mol_per_water",
    "gxtb_v201": "gxtb_v201_kj_mol_per_water",
    "dcm_main": "dcm_main_kj_mol_per_water",
}
for row in rows:
    phase = row["phase"]
    close(float(row["dmc_reference_kj_mol_per_water"]), REFERENCE[phase],
          5.0e-13, f"DMC reference {phase}")
    for provider, column in columns.items():
        close(float(row[column]), relative(provider, 2, phase),
              7.0e-10, f"k222 relative energy {provider}/{phase}")

platform_deltas = []
for phase in ("Ih", "VII", "XVII"):
    macos = energy("authors_exchange", 2, phase) * 8
    linux_path = (
        ROOT / "raw/authors_exchange_linux_k222" / phase / "result.json"
    )
    linux = float(json.loads(linux_path.read_text())["energy"])
    platform_deltas.append(linux - macos)
max_platform_delta = max(map(abs, platform_deltas))
if max_platform_delta > 2.0e-11:
    raise AssertionError(
        f"authors_exchange Linux/macOS mismatch: {max_platform_delta:.3e} Ha"
    )

print("mesh provider MAE_kJ_mol_per_water")
for key in sorted(computed):
    print(key[0], key[1], f"{computed[key]:.12f}")
print(f"authors_exchange Linux/macOS max delta (Ha): {max_platform_delta:.12e}")
print("model-revision coarse-grid validation: pass")
