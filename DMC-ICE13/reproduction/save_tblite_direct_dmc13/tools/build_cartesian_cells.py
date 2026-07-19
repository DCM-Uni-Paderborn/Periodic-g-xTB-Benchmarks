#!/usr/bin/env python3
"""Build Cartesian POSCAR files from the canonical DMC-ICE13 CP2K inputs."""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import OrderedDict
from itertools import product
from pathlib import Path


PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def parse_input(path: Path) -> tuple[list[list[float]], list[tuple[str, list[float]]]]:
    lines = path.read_text().splitlines()
    cell: dict[str, list[float]] = {}
    coordinates: list[tuple[str, list[float]]] = []
    in_cell = False
    in_coord = False
    scaled = False
    for raw in lines:
        line = raw.strip()
        upper = line.upper()
        if upper.startswith("&CELL"):
            in_cell = True
            continue
        if upper.startswith("&END CELL"):
            in_cell = False
            continue
        if upper.startswith("&COORD"):
            in_coord = True
            continue
        if upper.startswith("&END COORD"):
            in_coord = False
            continue
        if in_cell:
            match = re.match(r"^([ABC])\s+(.+)$", line, flags=re.IGNORECASE)
            if match:
                cell[match.group(1).upper()] = [float(value) for value in match.group(2).split()[:3]]
        elif in_coord:
            if upper == "SCALED":
                scaled = True
                continue
            fields = line.split()
            if len(fields) >= 4 and not line.startswith(("#", "!", "&")):
                coordinates.append((fields[0], [float(value) for value in fields[1:4]]))
    lattice = [cell[key] for key in "ABC"]
    if scaled:
        coordinates = [
            (
                element,
                [sum(frac[j] * lattice[j][i] for j in range(3)) for i in range(3)],
            )
            for element, frac in coordinates
        ]
    return lattice, coordinates


def expanded_atoms(
    mesh: int,
    lattice: list[list[float]],
    atoms: list[tuple[str, list[float]]],
) -> list[tuple[str, list[float]]]:
    by_element: OrderedDict[str, list[list[float]]] = OrderedDict()
    for element, position in atoms:
        by_element.setdefault(element, []).append(position)
    expanded: list[tuple[str, list[float]]] = []
    for element, sites in by_element.items():
        for site in sites:
            for image in product(range(mesh), repeat=3):
                expanded.append(
                    (
                        element,
                        [
                            site[i] + sum(image[j] * lattice[j][i] for j in range(3))
                            for i in range(3)
                        ],
                    )
                )
    return expanded


def write_poscar(path: Path, phase: str, mesh: int, lattice: list[list[float]], atoms: list[tuple[str, list[float]]]) -> None:
    by_element: OrderedDict[str, list[list[float]]] = OrderedDict()
    for element, position in atoms:
        by_element.setdefault(element, []).append(position)
    super_lattice = [[mesh * value for value in vector] for vector in lattice]
    output = [f"DMC-ICE13 ice {phase}; Cartesian {mesh}x{mesh}x{mesh} supercell", "1.0"]
    output.extend(" ".join(f"{value:.15f}" for value in vector) for vector in super_lattice)
    output.append(" ".join(by_element))
    output.append(" ".join(str(len(sites) * mesh**3) for sites in by_element.values()))
    output.append("Cartesian")
    for _element, position in expanded_atoms(mesh, lattice, atoms):
        output.append(" ".join(f"{value:.15f}" for value in position))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(output) + "\n")


def write_extxyz(path: Path, phase: str, mesh: int, lattice: list[list[float]], atoms: list[tuple[str, list[float]]]) -> None:
    super_lattice = [[mesh * value for value in vector] for vector in lattice]
    flat_lattice = " ".join(f"{value:.15f}" for vector in super_lattice for value in vector)
    sites = expanded_atoms(mesh, lattice, atoms)
    output = [
        str(len(sites)),
        (
            f'Lattice="{flat_lattice}" Properties=species:S:1:pos:R:3 '
            f'pbc="T T T" phase="{phase}" mesh="{mesh} {mesh} {mesh}" units="angstrom"'
        ),
    ]
    output.extend(
        f"{element} " + " ".join(f"{value:.15f}" for value in position)
        for element, position in sites
    )
    path.write_text("\n".join(output) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--meshes", type=int, nargs="+", default=(1, 2, 3, 4))
    args = parser.parse_args()
    manifest = ["mesh,phase,natoms,poscar_sha256,extxyz_sha256,source_input_sha256"]
    for mesh in args.meshes:
        for phase in PHASES:
            source = args.inputs / phase / "input.inp"
            lattice, atoms = parse_input(source)
            target = args.output / f"k{mesh}{mesh}{mesh}" / phase / "POSCAR"
            extxyz = target.with_name("structure.xyz")
            write_poscar(target, phase, mesh, lattice, atoms)
            write_extxyz(extxyz, phase, mesh, lattice, atoms)
            manifest.append(
                f"{mesh},{phase},{len(atoms) * mesh**3},"
                f"{hashlib.sha256(target.read_bytes()).hexdigest()},"
                f"{hashlib.sha256(extxyz.read_bytes()).hexdigest()},"
                f"{hashlib.sha256(source.read_bytes()).hexdigest()}"
            )
    (args.output.parent / "provenance" / "structure_hashes.csv").write_text("\n".join(manifest) + "\n")


if __name__ == "__main__":
    main()
