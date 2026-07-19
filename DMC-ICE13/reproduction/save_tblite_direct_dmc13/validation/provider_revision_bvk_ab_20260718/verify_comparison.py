#!/usr/bin/env python3
"""Recompute the archived direct-provider BvK comparison from raw JSON."""

from __future__ import annotations

import csv
import hashlib
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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

def verify_full_relative_mesh(mesh: int, table_name: str, summary_name: str) -> None:
    current_errors = []
    author_errors = []
    relative_shifts = []
    phases = []
    with (ROOT / table_name).open(newline="") as handle:
        for row in csv.DictReader(handle):
            phase = row["phase"]
            phases.append(phase)
            if row.get("mesh_n") and int(row["mesh_n"]) != mesh:
                raise SystemExit(f"wrong mesh in {table_name}: {row['mesh_n']}")
            values = {}
            for provider in ("current_save_tblite_cli", "seidler_pbc_cli_linux"):
                values[provider] = (
                    energy(provider, mesh, phase) / N_WATER[phase]
                    - energy(provider, mesh, "Ih") / N_WATER["Ih"]
                ) * HARTREE_TO_KJMOL
            current = values["current_save_tblite_cli"]
            author = values["seidler_pbc_cli_linux"]
            reference = float(row["dmc_reference_kj_mol_per_water"])
            close(current, float(row["current_save_tblite_kj_mol_per_water"]),
                  5.0e-10, f"full current relative {mesh} {phase}")
            close(author, float(row["author_pbc_kj_mol_per_water"]), 5.0e-10,
                  f"full author relative {mesh} {phase}")
            close(author - current,
                  float(row["author_minus_current_kj_mol_per_water"]), 5.0e-10,
                  f"full relative difference {mesh} {phase}")
            close(current - reference,
                  float(row["current_error_kj_mol_per_water"]), 5.0e-10,
                  f"full current error {mesh} {phase}")
            close(author - reference,
                  float(row["author_error_kj_mol_per_water"]), 5.0e-10,
                  f"full author error {mesh} {phase}")
            current_errors.append(abs(current - reference))
            author_errors.append(abs(author - reference))
            relative_shifts.append((abs(author - current), phase))
    expected_phases = set(N_WATER) - {"Ih"}
    if set(phases) != expected_phases or len(phases) != len(expected_phases):
        raise SystemExit(f"incomplete or duplicate phase set in {table_name}: {phases}")

    summary = json.loads((ROOT / summary_name).read_text())
    current_mae = statistics.mean(current_errors)
    author_mae = statistics.mean(author_errors)
    maximum_shift, maximum_phase = max(relative_shifts)
    close(current_mae, summary["current_save_tblite_mae_kj_mol_per_water"],
          5.0e-10, f"full current k{mesh}{mesh}{mesh} MAE")
    close(author_mae, summary["author_pbc_mae_kj_mol_per_water"],
          5.0e-10, f"full author k{mesh}{mesh}{mesh} MAE")
    close(author_mae - current_mae,
          summary["author_minus_current_mae_kj_mol_per_water"], 5.0e-10,
          f"full k{mesh}{mesh}{mesh} MAE shift")
    close(maximum_shift,
          summary["maximum_absolute_relative_energy_shift_kj_mol_per_water"],
          5.0e-10, f"full k{mesh}{mesh}{mesh} maximum relative-energy shift")
    if maximum_phase != summary["maximum_shift_phase"]:
        raise SystemExit(
            f"maximum shift phase at mesh {mesh}: actual={maximum_phase} "
            f"expected={summary['maximum_shift_phase']}"
        )


verify_full_relative_mesh(2, "full_k222_relative_comparison.csv", "full_k222_summary.json")
verify_full_relative_mesh(3, "full_k333_relative_comparison.csv", "full_k333_summary.json")

with (ROOT / "full_k333_absolute_energy_comparison.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
if {row["phase"] for row in rows} != set(N_WATER) or len(rows) != len(N_WATER):
    raise SystemExit("incomplete or duplicate phase set in full k333 absolute comparison")
for row in rows:
    phase = row["phase"]
    if int(row["mesh_n"]) != 3:
        raise SystemExit(f"wrong mesh in full k333 absolute comparison: {row['mesh_n']}")
    current_primitive = energy("current_save_tblite_cli", 3, phase)
    author_primitive = energy("seidler_pbc_cli_linux", 3, phase)
    close(current_primitive * 27, float(row["current_total_Ha"]), 5.0e-10,
          f"full current total 3 {phase}")
    close(author_primitive * 27, float(row["author_total_Ha"]), 5.0e-10,
          f"full author total 3 {phase}")
    close(current_primitive, float(row["current_Ha_per_primitive"]), 5.0e-12,
          f"full current primitive 3 {phase}")
    close(author_primitive, float(row["author_Ha_per_primitive"]), 5.0e-12,
          f"full author primitive 3 {phase}")
    close(author_primitive - current_primitive,
          float(row["author_minus_current_Ha_per_primitive"]), 5.0e-12,
          f"full absolute difference 3 {phase}")

manifest_pairs = set()
with (ROOT / "full_k333_input_manifest.csv").open(newline="") as handle:
    for row in csv.DictReader(handle):
        if int(row["mesh_n"]) != 3:
            raise SystemExit(f"wrong mesh in full k333 manifest: {row['mesh_n']}")
        relative_path = Path(row["path"])
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise SystemExit(f"unsafe path in full k333 manifest: {relative_path}")
        path = ROOT / relative_path
        if sha256(path) != row["sha256"]:
            raise SystemExit(f"hash mismatch in full k333 manifest: {relative_path}")
        manifest_pairs.add((row["phase"], row["implementation"]))
expected_manifest_pairs = {
    (phase, implementation)
    for phase in N_WATER
    for implementation in ("current", "author_pbc")
}
if manifest_pairs != expected_manifest_pairs or len(manifest_pairs) != 2 * len(N_WATER):
    raise SystemExit("incomplete or duplicate full k333 input manifest")

for phase in ("Ih", "VII", "XVII"):
    linux = energy("seidler_pbc_cli_linux", 2, phase) * 8
    macos = energy("seidler_pbc_cli", 2, phase) * 8
    close(linux, macos, 2.0e-11, f"Linux/macOS author pbc 2 {phase}")

print("provider-revision BvK comparison: all archived values verified")
