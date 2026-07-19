#!/usr/bin/env python3
"""Strict parsing helpers for canonical Gamma-centred BvK inputs."""

from __future__ import annotations

import math
import re
from pathlib import Path


def canonical_bvk_shift(mesh: int) -> float:
    """Return the MacDonald shift for a Gamma-centred BvK mesh."""
    if mesh <= 0:
        raise ValueError("mesh must be positive")
    if mesh % 2:
        return 0.0
    return (mesh - 1) / (2 * mesh)


def format_shift(value: float) -> str:
    return "0.0" if value == 0.0 else repr(value)


def input_mesh_and_water_count(path: Path) -> tuple[int, int]:
    """Parse and validate one cubic canonical MacDonald mesh and O count."""
    mesh_records: list[tuple[tuple[int, int, int], tuple[float, float, float]]] = []
    water_count = 0
    in_coordinates = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        code = re.split(r"[#!]", raw_line, maxsplit=1)[0].split()
        if len(code) >= 2 and tuple(token.upper() for token in code[:2]) == (
            "SCHEME",
            "MACDONALD",
        ):
            if len(code) != 8:
                raise ValueError(f"unsupported MACDONALD syntax in {path}: {raw_line}")
            try:
                dimensions = tuple(int(value) for value in code[2:5])
                shifts = tuple(
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in code[5:8]
                )
            except ValueError as exc:
                raise ValueError(f"invalid MACDONALD mesh in {path}: {raw_line}") from exc
            mesh_records.append((dimensions, shifts))
        line = raw_line.strip()
        upper = line.upper()
        if upper == "&COORD":
            in_coordinates = True
            continue
        if in_coordinates and upper.startswith("&END"):
            in_coordinates = False
            continue
        if in_coordinates and line and not line.startswith(("#", "!")):
            water_count += line.split()[0].upper() == "O"

    if len(mesh_records) != 1:
        raise ValueError(
            f"expected exactly one MACDONALD mesh in {path}, found {len(mesh_records)}"
        )
    dimensions, shifts = mesh_records[0]
    if len(set(dimensions)) != 1 or dimensions[0] <= 0:
        raise ValueError(f"non-positive or anisotropic mesh in {path}: {dimensions}")
    if any(not math.isfinite(value) for value in shifts):
        raise ValueError(f"non-finite MACDONALD shift in {path}: {shifts}")
    mesh = dimensions[0]
    expected_shift = canonical_bvk_shift(mesh)
    if any(not math.isclose(value, expected_shift, rel_tol=0.0, abs_tol=1.0e-14) for value in shifts):
        raise ValueError(
            f"noncanonical Gamma-centred BvK shift in {path}: "
            f"mesh={mesh} actual={shifts} expected={(expected_shift,) * 3}"
        )
    if water_count <= 0:
        raise ValueError(f"no oxygen atoms in coordinate section: {path}")
    return mesh, water_count
