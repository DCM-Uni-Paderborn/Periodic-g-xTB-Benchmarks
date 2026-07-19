#!/usr/bin/env python3
"""Verify that a CP2K Gamma-supercell input reproduces an archived POSCAR."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def vector(values: list[str], scale: float = 1.0) -> tuple[float, float, float]:
    if len(values) < 3:
        raise ValueError(f"expected three vector components, got: {values}")
    result = tuple(scale * float(value) for value in values[:3])
    if not all(math.isfinite(value) for value in result):
        raise ValueError(f"non-finite vector: {values}")
    return result


def parse_poscar(path: Path) -> tuple[
    list[tuple[float, float, float]],
    list[str],
    list[tuple[float, float, float]],
]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 8:
        raise ValueError(f"incomplete POSCAR: {path}")
    scale = float(lines[1])
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("only a finite, positive POSCAR scale is supported")
    cell = [vector(lines[index].split(), scale) for index in range(2, 5)]
    species = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(species) != len(counts) or any(count <= 0 for count in counts):
        raise ValueError("invalid POSCAR species/count records")

    cursor = 7
    if lines[cursor].lower().startswith("s"):
        cursor += 1
    coordinate_mode = lines[cursor].lower()
    cursor += 1
    direct = coordinate_mode.startswith("d")
    cartesian = coordinate_mode.startswith(("c", "k"))
    if not (direct or cartesian):
        raise ValueError(f"unsupported POSCAR coordinate mode: {lines[cursor - 1]}")

    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for label, count in zip(species, counts, strict=True):
        for _ in range(count):
            if cursor >= len(lines):
                raise ValueError("POSCAR coordinate block ended early")
            raw = vector(lines[cursor].split(), scale if cartesian else 1.0)
            cursor += 1
            if direct:
                raw = tuple(
                    raw[0] * cell[0][axis]
                    + raw[1] * cell[1][axis]
                    + raw[2] * cell[2][axis]
                    for axis in range(3)
                )
            labels.append(label)
            coordinates.append(raw)
    return cell, labels, coordinates


def section(lines: list[str], name: str) -> list[str]:
    wanted = name.upper()
    start = None
    depth = 0
    result: list[str] = []
    for line in lines:
        words = line.strip().split()
        if not words:
            continue
        token = words[0].upper()
        if start is None:
            if token == f"&{wanted}":
                start = True
                depth = 1
            continue
        if token.startswith("&") and token != "&END":
            depth += 1
        elif token == "&END":
            if len(words) == 1 or words[1].upper() == wanted:
                depth -= 1
                if depth == 0:
                    return result
        if depth == 1:
            result.append(line)
    raise ValueError(f"missing or incomplete CP2K &{wanted} section")


def parse_cp2k(path: Path) -> tuple[
    list[tuple[float, float, float]],
    list[str],
    list[tuple[float, float, float]],
    bool,
    bool,
]:
    lines = path.read_text(encoding="utf-8").splitlines()
    cell_lines = section(lines, "CELL")
    cell_by_name: dict[str, tuple[float, float, float]] = {}
    periodic_xyz = False
    for line in cell_lines:
        words = line.split()
        if not words or words[0].startswith("#"):
            continue
        key = words[0].upper()
        if key in {"A", "B", "C"}:
            cell_by_name[key] = vector(words[1:])
        elif key == "PERIODIC":
            periodic_xyz = len(words) > 1 and words[1].upper() == "XYZ"
    if set(cell_by_name) != {"A", "B", "C"}:
        raise ValueError("CP2K input must define explicit A, B, and C vectors")

    coord_lines = section(lines, "COORD")
    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for line in coord_lines:
        words = line.split()
        if not words or words[0].startswith("#"):
            continue
        if words[0].upper() in {"UNIT", "SCALED"}:
            raise ValueError("the qualification input must use default Cartesian Angstrom coordinates")
        if len(words) < 4:
            raise ValueError(f"invalid CP2K coordinate record: {line}")
        labels.append(words[0])
        coordinates.append(vector(words[1:]))

    has_kpoints = any(
        line.strip().split()[:1]
        and line.strip().split()[0].upper() == "&KPOINTS"
        for line in lines
    )
    cell = [cell_by_name[name] for name in ("A", "B", "C")]
    return cell, labels, coordinates, periodic_xyz, has_kpoints


def maximum_delta(
    left: list[tuple[float, float, float]],
    right: list[tuple[float, float, float]],
) -> float:
    if len(left) != len(right):
        return math.inf
    return max(
        (abs(a - b) for lhs, rhs in zip(left, right, strict=True) for a, b in zip(lhs, rhs, strict=True)),
        default=0.0,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cp2k_input", type=Path)
    parser.add_argument("poscar", type=Path)
    parser.add_argument("--tolerance", type=float, default=5.0e-12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0.0:
        raise ValueError("tolerance must be finite and non-negative")

    poscar_cell, poscar_labels, poscar_coordinates = parse_poscar(args.poscar)
    cp2k_cell, cp2k_labels, cp2k_coordinates, periodic_xyz, has_kpoints = parse_cp2k(
        args.cp2k_input
    )
    cell_delta = maximum_delta(poscar_cell, cp2k_cell)
    coordinate_delta = maximum_delta(poscar_coordinates, cp2k_coordinates)
    labels_equal = cp2k_labels == poscar_labels
    status = (
        "PASS"
        if periodic_xyz
        and not has_kpoints
        and labels_equal
        and cell_delta <= args.tolerance
        and coordinate_delta <= args.tolerance
        else "FAIL"
    )
    result = {
        "atom_count": len(cp2k_coordinates),
        "cell_max_abs_delta_angstrom": cell_delta,
        "coordinate_max_abs_delta_angstrom": coordinate_delta,
        "elements_and_order_equal": labels_equal,
        "explicit_kpoints_section_present": has_kpoints,
        "periodic_xyz": periodic_xyz,
        "status": status,
        "tolerance_angstrom": args.tolerance,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
