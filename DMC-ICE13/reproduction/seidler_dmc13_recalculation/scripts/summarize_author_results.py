#!/usr/bin/env python3
"""Summarize independent pbc and mstore-inorganic result directories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
HARTREE_TO_KJMOL = 2625.4996394798254


def first_hash(path: Path) -> str:
    return path.read_text(encoding="utf-8").split()[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_branch(root: Path, mesh: int) -> tuple[dict[str, float], str]:
    energies = {}
    hashes = set()
    for phase in PHASES:
        run = root / "runs" / f"k{mesh}{mesh}{mesh}" / phase
        status = run / "exit_status"
        result = run / "result.json"
        if not status.is_file() or status.read_text(encoding="utf-8").strip() != "0":
            raise RuntimeError(f"incomplete run: {run}")
        energy = float(json.loads(result.read_text(encoding="utf-8"))["energy"])
        if not math.isfinite(energy):
            raise RuntimeError(f"non-finite energy: {run}")
        structure = root / "structures" / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
        if not structure.is_file():
            raise RuntimeError(f"missing generated structure: {structure}")
        if first_hash(run / "input.sha256") != sha256(structure):
            raise RuntimeError(f"input hash differs from generated structure: {run}")
        energies[phase] = energy
        hashes.add(first_hash(run / "binary.sha256"))
    if len(hashes) != 1:
        raise RuntimeError(f"mesh {mesh} uses multiple executable hashes in {root}")
    return energies, hashes.pop()


def water_counts() -> dict[str, int]:
    counts = {}
    for phase in PHASES:
        lines = (PACKAGE / "structures" / "primitive" / phase / "POSCAR").read_text(encoding="utf-8").splitlines()
        counts[phase] = sum(int(value) for value in lines[6].split()) // 3
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pbc_root", type=Path)
    parser.add_argument("mstore_root", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--meshes", type=int, nargs="+", default=(2, 3))
    args = parser.parse_args()

    references = {
        row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
        for row in csv.DictReader((PACKAGE / "tables" / "dmc_reference_relative_energies.csv").open())
    }
    waters = water_counts()
    rows = []
    for mesh in args.meshes:
        pbc, pbc_hash = load_branch(args.pbc_root, mesh)
        mstore, mstore_hash = load_branch(args.mstore_root, mesh)
        for phase in PHASES:
            pbc_primitive = pbc[phase] / mesh**3
            mstore_primitive = mstore[phase] / mesh**3
            if phase == "Ih":
                pbc_relative = 0.0
                mstore_relative = 0.0
            else:
                pbc_relative = (
                    pbc_primitive / waters[phase]
                    - pbc["Ih"] / (mesh**3 * waters["Ih"])
                ) * HARTREE_TO_KJMOL
                mstore_relative = (
                    mstore_primitive / waters[phase]
                    - mstore["Ih"] / (mesh**3 * waters["Ih"])
                ) * HARTREE_TO_KJMOL
            reference = 0.0 if phase == "Ih" else references[phase]
            rows.append({
                "mesh_n": mesh,
                "phase": phase,
                "pbc_energy_Ha_per_primitive": f"{pbc_primitive:.15f}",
                "mstore_energy_Ha_per_primitive": f"{mstore_primitive:.15f}",
                "mstore_minus_pbc_Ha_per_primitive": f"{mstore_primitive - pbc_primitive:+.15e}",
                "dmc_reference_kj_mol_per_H2O": f"{reference:.6f}",
                "pbc_relative_kj_mol_per_H2O": f"{pbc_relative:.12f}",
                "mstore_relative_kj_mol_per_H2O": f"{mstore_relative:.12f}",
                "mstore_minus_pbc_kj_mol_per_H2O": f"{mstore_relative - pbc_relative:+.12f}",
                "pbc_binary_sha256": pbc_hash,
                "mstore_binary_sha256": mstore_hash,
            })
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"rows={len(rows)} status=PASS")


if __name__ == "__main__":
    main()
