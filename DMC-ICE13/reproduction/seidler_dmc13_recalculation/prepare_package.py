#!/usr/bin/env python3
"""Regenerate the compact author-facing DMC-ICE13 recalculation package."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "save_tblite_direct_dmc13"
REPOSITORY = HERE.parents[2]
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def compare_poscars(first: Path, second: Path) -> float:
    left = [line.strip() for line in first.read_text(encoding="utf-8").splitlines()]
    right = [line.strip() for line in second.read_text(encoding="utf-8").splitlines()]
    if len(left) != len(right) or left[:2] != right[:2] or left[5:8] != right[5:8]:
        raise AssertionError(f"POSCAR metadata or atom order differs: {first} {second}")
    residual = 0.0
    for index in (*range(2, 5), *range(8, len(left))):
        left_values = [float(value) for value in left[index].split()[:3]]
        right_values = [float(value) for value in right[index].split()[:3]]
        residual = max(
            residual,
            *(abs(a - b) for a, b in zip(left_values, right_values, strict=True)),
        )
    return residual


def main() -> None:
    structure_rows = []
    for phase in PHASES:
        for name in ("POSCAR", "structure.xyz"):
            copy(SOURCE / "structures/k111" / phase / name, HERE / "structures/primitive" / phase / name)
        poscar = HERE / "structures/primitive" / phase / "POSCAR"
        xyz = HERE / "structures/primitive" / phase / "structure.xyz"
        lines = poscar.read_text(encoding="utf-8").splitlines()
        structure_rows.append({
            "phase": phase,
            "atom_count": sum(int(value) for value in lines[6].split()),
            "poscar_sha256": sha256(poscar),
            "extxyz_sha256": sha256(xyz),
        })
    manifest = HERE / "structures/structure_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(structure_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(structure_rows)

    copies = {
        SOURCE / "tables/absolute_energies_vs_mesh.csv": HERE / "tables/current_absolute_energies_by_mesh.csv",
        SOURCE / "tables/dmc_reference_relative_energies.csv": HERE / "tables/dmc_reference_relative_energies.csv",
        REPOSITORY / "validation/three_route_k333_closure_20260719/absolute_energies.csv": HERE / "tables/three_route_absolute_energies_k333.csv",
        REPOSITORY / "validation/three_route_k333_closure_20260719/relative_energies_and_errors.csv": HERE / "tables/three_route_relative_energies_k333.csv",
    }
    for source, target in copies.items():
        copy(source, target)

    author_rows = []
    author_relative_rows = []
    provider = SOURCE / "validation/provider_revision_bvk_ab_20260718"
    generated_k222 = HERE / ".verification/provider-k222"
    subprocess.run(
        [
            sys.executable,
            str(provider / "compare_complete_mesh.py"),
            "--mesh", "2",
            "--current-root", str(SOURCE / "results/current_save_tblite_cli"),
            "--author-root", str(provider / "seidler_pbc_cli_linux"),
            "--reference-csv", str(provider / "full_k222_relative_comparison.csv"),
            "--output-dir", str(generated_k222),
        ],
        check=True,
    )
    for mesh in (2, 3):
        table = (
            generated_k222 / "absolute_energy_comparison.csv"
            if mesh == 2
            else provider / "full_k333_absolute_energy_comparison.csv"
        )
        with table.open(
            newline="", encoding="utf-8"
        ) as handle:
            for row in csv.DictReader(handle):
                author_rows.append({
                    "mesh_n": mesh,
                    "phase": row["phase"],
                    "author_pbc_total_Ha": row["author_total_Ha"],
                    "author_pbc_Ha_per_primitive": row["author_Ha_per_primitive"],
                    "current_cli_total_Ha": row["current_total_Ha"],
                    "current_cli_Ha_per_primitive": row["current_Ha_per_primitive"],
                    "author_minus_current_Ha_per_primitive": row["author_minus_current_Ha_per_primitive"],
                })
        relative_table = (
            generated_k222 / "relative_energy_comparison.csv"
            if mesh == 2
            else provider / "full_k333_relative_comparison.csv"
        )
        with relative_table.open(newline="", encoding="utf-8") as handle:
            author_relative_rows.extend(csv.DictReader(handle))
    author_table = HERE / "tables/author_pbc_absolute_energies.csv"
    author_table.parent.mkdir(parents=True, exist_ok=True)
    with author_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(author_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(author_rows)
    author_relative_table = HERE / "tables/author_pbc_relative_energies.csv"
    with author_relative_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(author_relative_rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(author_relative_rows)

    # Prove that the vendored primitive structures reproduce every archived
    # explicit cell through 4^3 with exact order and sub-picometre residuals.
    builder = HERE / "scripts/build_bvk_from_poscar.py"
    for mesh in (1, 2, 3, 4):
        for phase in PHASES:
            temporary = HERE / ".verification" / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
            subprocess.run(
                [sys.executable, str(builder), str(HERE / "structures/primitive" / phase / "POSCAR"), str(temporary), str(mesh)],
                check=True,
            )
            archived = SOURCE / "structures" / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
            residual = compare_poscars(temporary, archived)
            if residual > 5.0e-12:
                raise AssertionError(
                    f"generated structure differs: mesh={mesh} phase={phase} "
                    f"residual={residual:.6e} Angstrom"
                )

    (HERE / "sources.json").write_text(
        json.dumps({
            "current_save_tblite_source": "15915c9435644eb257178ca8f8bf7220c38b1a84",
            "current_save_tblite_linux_cli_sha256": "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a",
            "current_cp2k_source": "8520b2e592cd04d35081ab4ad46d92c606071e23",
            "current_cp2k_linux_sha256": "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f",
            "author_pbc_source": "c932120d2580811901de6a1fe3f89b943c251766",
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    # Remove transient verifier output and the superseded one-mesh table from
    # early package drafts before constructing the portable manifest.
    shutil.rmtree(HERE / ".verification", ignore_errors=True)
    (HERE / "tables/author_pbc_relative_energies_k333.csv").unlink(missing_ok=True)

    files = sorted(
        path for path in HERE.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and ".verification" not in path.parts
    )
    (HERE / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {path.relative_to(HERE)}\n" for path in files),
        encoding="utf-8",
    )
    print(f"files={len(files)} structures={len(structure_rows)} status=PASS")


if __name__ == "__main__":
    main()
