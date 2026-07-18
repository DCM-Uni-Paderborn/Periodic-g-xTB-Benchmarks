#!/usr/bin/env python3
"""Verify current CP2K-native versus direct save_tblite absolute energies."""

from __future__ import annotations

import csv
import hashlib
import json
import math
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


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_row(row: dict[str, str]) -> float:
    mesh = int(row["mesh_n"])
    mesh_id = row["mesh_id"]
    phase = row["phase"]
    native = ROOT / "results" / "current_cp2k_native" / mesh_id / phase / "cp2k.out"
    direct = ROOT / "results" / "current_save_tblite_cli" / mesh_id / phase / "result.json"
    if not native.is_file() or not direct.is_file():
        raise AssertionError(f"missing required result: {mesh_id}/{phase}")
    native_text = native.read_text(encoding="utf-8", errors="replace")
    if "PROGRAM ENDED AT" not in native_text or "ENERGY| Total FORCE_EVAL" not in native_text:
        raise AssertionError(f"incomplete CP2K result: {native}")
    direct_payload = json.loads(direct.read_text(encoding="utf-8"))
    if not math.isfinite(float(direct_payload["energy"])):
        raise AssertionError(f"invalid direct energy: {direct}")
    if digest(native) != row["cp2k_output_sha256"]:
        raise AssertionError(f"CP2K hash mismatch: {mesh_id}/{phase}")
    if digest(direct) != row["save_tblite_json_sha256"]:
        raise AssertionError(f"direct hash mismatch: {mesh_id}/{phase}")
    delta = abs(float(row["native_minus_cli_per_primitive_Ha"]))
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
