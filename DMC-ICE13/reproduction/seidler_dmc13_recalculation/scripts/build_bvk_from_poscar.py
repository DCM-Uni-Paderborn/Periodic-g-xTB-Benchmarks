#!/usr/bin/env python3
"""Build an explicit Cartesian cubic BvK supercell from a primitive POSCAR."""

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path


def parse(path: Path) -> tuple[str, list[list[float]], list[str], list[int], list[list[list[float]]]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 8:
        raise ValueError(f"incomplete POSCAR: {path}")
    scale = float(lines[1])
    if scale <= 0.0:
        raise ValueError("only positive POSCAR scale factors are supported")
    lattice = [
        [scale * float(value) for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]
    elements = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(elements) != len(counts):
        raise ValueError("element/count mismatch")
    cursor = 7
    if lines[cursor].lower().startswith("s"):
        cursor += 1
    if not lines[cursor].lower().startswith(("c", "k")):
        raise ValueError("the recalculation package requires Cartesian coordinates")
    cursor += 1
    grouped: list[list[list[float]]] = []
    for count in counts:
        sites = []
        for _ in range(count):
            sites.append([scale * float(value) for value in lines[cursor].split()[:3]])
            cursor += 1
        grouped.append(sites)
    return lines[0], lattice, elements, counts, grouped


def build(source: Path, target: Path, mesh: int) -> None:
    if mesh < 1:
        raise ValueError("mesh must be positive")
    title, lattice, elements, counts, grouped = parse(source)
    base_title = title.split("; Cartesian", maxsplit=1)[0]
    output = [f"{base_title}; Cartesian {mesh}x{mesh}x{mesh} supercell", "1.0"]
    output.extend(
        " ".join(f"{mesh * value:.15f}" for value in vector) for vector in lattice
    )
    output.extend((" ".join(elements), " ".join(str(value * mesh**3) for value in counts), "Cartesian"))
    for sites in grouped:
        for site in sites:
            for image in product(range(mesh), repeat=3):
                position = [
                    site[axis]
                    + sum(image[index] * lattice[index][axis] for index in range(3))
                    for axis in range(3)
                ]
                output.append(" ".join(f"{value:.15f}" for value in position))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(output) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("mesh", type=int)
    args = parser.parse_args()
    build(args.source, args.target, args.mesh)


if __name__ == "__main__":
    main()
