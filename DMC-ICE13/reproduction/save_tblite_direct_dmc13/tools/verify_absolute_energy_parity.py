#!/usr/bin/env python3
"""Verify current CP2K-native versus direct save_tblite absolute energies."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
REQUIRED = {
    1: PHASES,
    2: PHASES,
    3: PHASES,
    4: ("Ih", "VII", "XVII"),
}
TOLERANCE_HARTREE = 2.0e-7
TABLE_TOTAL_TOLERANCE_HARTREE = 5.0e-10
TABLE_PRIMITIVE_TOLERANCE_HARTREE = 5.0e-12
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(
            f"{label}: actual={actual:.16g} expected={expected:.16g} "
            f"difference={actual - expected:+.3e} tolerance={tolerance:.3e}"
        )


def native_energy(text: str, path: Path) -> float:
    values = [float(match.group(1)) for line in text.splitlines() if (match := ENERGY_RE.match(line))]
    if "PROGRAM ENDED AT" not in text or not values:
        raise AssertionError(f"incomplete CP2K result: {path}")
    return values[-1]


def verify_row(row: dict[str, str]) -> float:
    mesh = int(row["mesh_n"])
    mesh_id = row["mesh_id"]
    phase = row["phase"]
    expected_mesh_id = f"k{mesh}{mesh}{mesh}"
    if mesh_id != expected_mesh_id:
        raise AssertionError(
            f"mesh identifier mismatch: mesh={mesh} actual={mesh_id} expected={expected_mesh_id}"
        )
    native = ROOT / "results" / "current_cp2k_native" / mesh_id / phase / "cp2k.out"
    direct = ROOT / "results" / "current_save_tblite_cli" / mesh_id / phase / "result.json"
    poscar = ROOT / "structures" / mesh_id / phase / "POSCAR"
    if not native.is_file() or not direct.is_file():
        raise AssertionError(f"missing required result: {mesh_id}/{phase}")
    if not poscar.is_file():
        raise AssertionError(f"missing required Cartesian structure: {mesh_id}/{phase}")
    native_text = native.read_text(encoding="utf-8", errors="replace")
    native_value = native_energy(native_text, native)
    direct_payload = json.loads(direct.read_text(encoding="utf-8"))
    direct_total = float(direct_payload["energy"])
    if not math.isfinite(direct_total):
        raise AssertionError(f"invalid direct energy: {direct}")
    replicas = mesh**3
    direct_per_primitive = direct_total / replicas
    computed_delta = native_value - direct_per_primitive
    natom_primitive = int(row["natom_primitive"])
    natom_supercell = int(row["natom_cli_supercell"])
    if natom_supercell != natom_primitive * replicas:
        raise AssertionError(
            f"supercell atom-count mismatch {mesh_id}/{phase}: "
            f"actual={natom_supercell} expected={natom_primitive * replicas}"
        )
    if digest(native) != row["cp2k_output_sha256"]:
        raise AssertionError(f"CP2K hash mismatch: {mesh_id}/{phase}")
    if digest(direct) != row["save_tblite_json_sha256"]:
        raise AssertionError(f"direct hash mismatch: {mesh_id}/{phase}")
    if digest(poscar) != row["poscar_sha256"]:
        raise AssertionError(f"POSCAR hash mismatch: {mesh_id}/{phase}")
    close(
        native_value,
        float(row["cp2k_native_energy_Ha_per_primitive"]),
        TABLE_PRIMITIVE_TOLERANCE_HARTREE,
        f"native table energy {mesh_id}/{phase}",
    )
    close(
        direct_total,
        float(row["save_tblite_cli_energy_Ha_supercell"]),
        TABLE_TOTAL_TOLERANCE_HARTREE,
        f"direct total table energy {mesh_id}/{phase}",
    )
    close(
        direct_per_primitive,
        float(row["save_tblite_cli_energy_Ha_per_primitive"]),
        TABLE_PRIMITIVE_TOLERANCE_HARTREE,
        f"direct primitive table energy {mesh_id}/{phase}",
    )
    close(
        computed_delta,
        float(row["native_minus_cli_per_primitive_Ha"]),
        TABLE_PRIMITIVE_TOLERANCE_HARTREE,
        f"native/direct table difference {mesh_id}/{phase}",
    )
    delta = abs(computed_delta)
    if delta > TOLERANCE_HARTREE:
        raise AssertionError(
            f"parity failure {mesh_id}/{phase}: {delta:.6e} Ha > {TOLERANCE_HARTREE:.6e} Ha"
        )
    return delta


def main() -> None:
    table = ROOT / "tables" / "absolute_energies_vs_mesh.csv"
    with table.open(newline="", encoding="utf-8") as handle:
        rows = {(int(row["mesh_n"]), row["phase"]): row for row in csv.DictReader(handle)}
    for mesh, phases in REQUIRED.items():
        deltas = [verify_row(rows[(mesh, phase)]) for phase in phases]
        rms = math.sqrt(sum(delta * delta for delta in deltas) / len(deltas))
        print(
            f"mesh={mesh} coverage={len(deltas)} max_abs_delta_Ha={max(deltas):.12e} "
            f"rms_delta_Ha={rms:.12e} status=pass"
        )


if __name__ == "__main__":
    main()
