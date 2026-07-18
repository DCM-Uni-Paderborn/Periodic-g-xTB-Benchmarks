#!/usr/bin/env python3
"""Collect absolute CP2K-native and direct save_tblite energies."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def cp2k_energy(path: Path) -> float | None:
    if not path.exists():
        return None
    values = re.findall(r"ENERGY\| Total FORCE_EVAL.*?(-?\d+\.\d+)", path.read_text(errors="ignore"))
    return float(values[-1]) if values else None


def cli_energy(path: Path) -> float | None:
    if not path.exists():
        return None
    return float(json.loads(path.read_text())["energy"])


def atom_count(path: Path) -> int:
    lines = path.read_text().splitlines()
    return sum(int(value) for value in lines[6].split())


def main() -> None:
    rows: list[dict[str, object]] = []
    for mesh in (1, 2, 3, 4):
        mesh_id = f"k{mesh}{mesh}{mesh}"
        for phase in PHASES:
            poscar = ROOT / "structures" / mesh_id / phase / "POSCAR"
            native_file = ROOT / "results" / "current_cp2k_native" / mesh_id / phase / "cp2k.out"
            cli_file = ROOT / "results" / "current_save_tblite_cli" / mesh_id / phase / "result.json"
            native = cp2k_energy(native_file)
            cli_supercell = cli_energy(cli_file)
            cli_primitive = cli_supercell / mesh**3 if cli_supercell is not None else None
            delta = native - cli_primitive if native is not None and cli_primitive is not None else None
            rows.append(
                {
                    "mesh_n": mesh,
                    "mesh_id": mesh_id,
                    "phase": phase,
                    "natom_primitive": atom_count(ROOT / "structures" / "k111" / phase / "POSCAR"),
                    "natom_cli_supercell": atom_count(poscar),
                    "cp2k_native_energy_Ha_per_primitive": "" if native is None else f"{native:.15f}",
                    "save_tblite_cli_energy_Ha_supercell": "" if cli_supercell is None else f"{cli_supercell:.15f}",
                    "save_tblite_cli_energy_Ha_per_primitive": "" if cli_primitive is None else f"{cli_primitive:.15f}",
                    "native_minus_cli_per_primitive_Ha": "" if delta is None else f"{delta:+.12e}",
                    "poscar_sha256": digest(poscar),
                    "cp2k_output_sha256": digest(native_file),
                    "save_tblite_json_sha256": digest(cli_file),
                }
            )
    output = ROOT / "tables" / "absolute_energies_vs_mesh.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
