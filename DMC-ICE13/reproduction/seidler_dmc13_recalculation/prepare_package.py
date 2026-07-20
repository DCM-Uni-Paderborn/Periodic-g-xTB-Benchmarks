#!/usr/bin/env python3
"""Verify the portable DMC-ICE13 comparison package and refresh its hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sidecar_hash(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"missing hash sidecar: {path}")
    fields = path.read_text(encoding="utf-8").split()
    if not fields:
        raise AssertionError(f"empty hash sidecar: {path}")
    return fields[0].lower()


def structure_rows() -> list[dict[str, object]]:
    rows = []
    for phase in PHASES:
        poscar = HERE / "structures" / "primitive" / phase / "POSCAR"
        xyz = HERE / "structures" / "primitive" / phase / "structure.xyz"
        if not poscar.is_file() or not xyz.is_file():
            raise FileNotFoundError(f"missing structure for {phase}")
        lines = poscar.read_text(encoding="utf-8").splitlines()
        atom_count = sum(int(value) for value in lines[6].split())
        xyz_count = int(xyz.read_text(encoding="utf-8").splitlines()[0])
        if xyz_count != atom_count:
            raise AssertionError(f"POSCAR/extxyz atom count differs for {phase}")
        rows.append({
            "phase": phase,
            "atom_count": atom_count,
            "water_molecule_count": atom_count // 3,
            "coordinate_mode": "absolute Cartesian Angstrom",
            "poscar_sha256": sha256(poscar),
            "extxyz_sha256": sha256(xyz),
        })
    return rows


def write_structure_manifest(rows: list[dict[str, object]]) -> None:
    path = HERE / "structures" / "structure_manifest.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_poscar(path: Path) -> tuple[tuple[str, ...], tuple[int, ...], list[list[float]], list[list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    scale = float(lines[1])
    species = tuple(lines[5].split())
    counts = tuple(int(value) for value in lines[6].split())
    cell = [[float(value) * scale for value in lines[index].split()[:3]] for index in (2, 3, 4)]
    start = 8
    coordinates = [
        [float(value) for value in lines[index].split()[:3]]
        for index in range(start, start + sum(counts))
    ]
    return species, counts, cell, coordinates


def maximum_difference(left: list[list[float]], right: list[list[float]]) -> float:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise AssertionError("numeric POSCAR arrays have different shapes")
    return max((abs(a - b) for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b)), default=0.0)


def parse_cp2k_structure(path: Path) -> tuple[tuple[str, ...], list[list[float]], list[list[float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    cell_vectors: dict[str, list[float]] = {}
    coordinate_lines: list[tuple[str, list[float]]] = []
    in_cell = False
    in_coord = False
    scaled = False
    for raw in lines:
        fields = raw.split()
        if not fields:
            continue
        upper = [field.upper() for field in fields]
        if upper[:1] == ["&CELL"]:
            in_cell = True
            continue
        if upper[:2] == ["&END", "CELL"]:
            in_cell = False
            continue
        if upper[:1] == ["&COORD"]:
            in_coord = True
            continue
        if upper[:2] == ["&END", "COORD"]:
            in_coord = False
            continue
        if in_cell and upper[0] in {"A", "B", "C"} and len(fields) >= 4:
            cell_vectors[upper[0]] = [float(value) for value in fields[1:4]]
        if in_coord:
            if upper[0] == "SCALED":
                scaled = True
            elif not fields[0].startswith("#") and len(fields) >= 4:
                coordinate_lines.append(
                    (fields[0], [float(value) for value in fields[1:4]])
                )
    if set(cell_vectors) != {"A", "B", "C"} or not coordinate_lines:
        raise AssertionError(f"incomplete CP2K structure: {path}")
    cell = [cell_vectors[label] for label in ("A", "B", "C")]
    species = tuple(label for label, _ in coordinate_lines)
    coordinates = []
    for _, values in coordinate_lines:
        if scaled:
            coordinates.append([
                sum(values[axis] * cell[axis][component] for axis in range(3))
                for component in range(3)
            ])
        else:
            coordinates.append(values)
    return species, cell, coordinates


def verify_cp2k_input_structures(rows: list[dict[str, object]]) -> int:
    raw = HERE / "raw" / "cp2k_native"
    checked = 0
    for row in rows:
        phase = str(row["phase"])
        primitive = parse_poscar(HERE / "structures" / "primitive" / phase / "POSCAR")
        primitive_species = tuple(
            species
            for species, count in zip(primitive[0], primitive[1])
            for _ in range(count)
        )
        for sidecar in sorted(raw.glob(f"k*-reduced/{phase}/input.sha256")):
            run = sidecar.parent
            cp2k_input = run / "input.inp"
            if not cp2k_input.is_file():
                raise FileNotFoundError(f"missing exact CP2K input: {cp2k_input}")
            if read_sidecar_hash(sidecar) != sha256(cp2k_input):
                raise AssertionError(f"CP2K input hash differs: {cp2k_input}")
            cp2k_species, cp2k_cell, cp2k_coordinates = parse_cp2k_structure(cp2k_input)
            if cp2k_species != primitive_species:
                raise AssertionError(f"CP2K/primitive species order differs: {cp2k_input}")
            cell_difference = maximum_difference(cp2k_cell, primitive[2])
            coordinate_difference = maximum_difference(cp2k_coordinates, primitive[3])
            if max(cell_difference, coordinate_difference) > 5.0e-10:
                raise AssertionError(
                    "CP2K/primitive geometry differs: "
                    f"input={cp2k_input} cell={cell_difference:.3e} "
                    f"coordinates={coordinate_difference:.3e}"
                )
            checked += 1
    return checked


def verify_supercell_builder(rows: list[dict[str, object]]) -> None:
    counts = {str(row["phase"]): int(row["atom_count"]) for row in rows}
    builder = HERE / "scripts" / "build_bvk_from_poscar.py"
    with tempfile.TemporaryDirectory(prefix="gxtb-bvk-check-") as temporary:
        root = Path(temporary)
        for mesh in (1, 2, 3, 4):
            for phase in PHASES:
                target = root / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
                subprocess.run(
                    [
                        sys.executable,
                        str(builder),
                        str(HERE / "structures" / "primitive" / phase / "POSCAR"),
                        str(target),
                        str(mesh),
                    ],
                    check=True,
                )
                lines = target.read_text(encoding="utf-8").splitlines()
                generated_count = sum(int(value) for value in lines[6].split())
                if generated_count != counts[phase] * mesh**3:
                    raise AssertionError(
                        f"generated atom count differs: phase={phase} mesh={mesh}"
                    )
                cli_run = (
                    HERE
                    / "raw"
                    / "current_pbc_cli"
                    / f"cli-k{mesh}{mesh}{mesh}"
                    / phase
                )
                cli_input = cli_run / "POSCAR"
                cli_input_hash = read_sidecar_hash(cli_run / "input.sha256")
                if not cli_input.is_file():
                    raise FileNotFoundError(f"missing archived CLI input: {cli_input}")
                if cli_input_hash != sha256(cli_input):
                    raise AssertionError(
                        "archived CLI input differs from its recorded hash: "
                        f"phase={phase} mesh={mesh}"
                    )
                generated_data = parse_poscar(target)
                cli_data = parse_poscar(cli_input)
                if generated_data[:2] != cli_data[:2]:
                    raise AssertionError(
                        "generated and archived CLI inputs differ in species/counts: "
                        f"phase={phase} mesh={mesh}"
                    )
                cell_difference = maximum_difference(generated_data[2], cli_data[2])
                coordinate_difference = maximum_difference(generated_data[3], cli_data[3])
                if max(cell_difference, coordinate_difference) > 5.0e-12:
                    raise AssertionError(
                        "generated and archived CLI inputs differ numerically: "
                        f"phase={phase} mesh={mesh} "
                        f"cell={cell_difference:.3e} coordinates={coordinate_difference:.3e}"
                    )
                if mesh in (2, 3):
                    archived = HERE / "raw" / "mstore_inorganic_cli" / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
                    if not archived.is_file():
                        raise FileNotFoundError(f"missing archived BvK input: {archived}")
                    archived_data = parse_poscar(archived)
                    if generated_data[:2] != archived_data[:2]:
                        raise AssertionError(f"species/order differs: phase={phase} mesh={mesh}")
                    if maximum_difference(generated_data[2], archived_data[2]) > 5.0e-12:
                        raise AssertionError(f"cell differs: phase={phase} mesh={mesh}")
                    if maximum_difference(generated_data[3], archived_data[3]) > 5.0e-12:
                        raise AssertionError(f"coordinates differ: phase={phase} mesh={mesh}")


def refresh_hash_manifest() -> int:
    files = sorted(
        path
        for path in HERE.rglob("*")
        if path.is_file()
        and path.name not in {"SHA256SUMS", ".DS_Store"}
        and "__pycache__" not in path.parts
    )
    (HERE / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {path.relative_to(HERE)}\n" for path in files),
        encoding="utf-8",
    )
    return len(files)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--refresh-manifest-only",
        action="store_true",
        help="skip table regeneration and refresh only structures and SHA256SUMS",
    )
    args = parser.parse_args()

    rows = structure_rows()
    write_structure_manifest(rows)
    verify_supercell_builder(rows)
    cp2k_input_count = verify_cp2k_input_structures(rows)
    if not args.refresh_manifest_only:
        subprocess.run(
            [sys.executable, str(HERE / "scripts" / "assemble_comparison_tables.py")],
            check=True,
        )
    file_count = refresh_hash_manifest()
    print(
        f"files={file_count} structures={len(rows)} "
        f"cp2k_inputs={cp2k_input_count} status=PASS"
    )


if __name__ == "__main__":
    main()
