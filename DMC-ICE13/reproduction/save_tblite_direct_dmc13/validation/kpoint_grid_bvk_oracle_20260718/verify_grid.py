#!/usr/bin/env python3
"""Verify that printed CP2K regular meshes are exact BvK reciprocal grids."""

from __future__ import annotations

import itertools
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
POINT = re.compile(
    r"^\s*BRILLOUIN\|\s+(\d+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$"
)
TOLERANCE = 6.0e-6


def canonical(value: float) -> float:
    result = value % 1.0
    if math.isclose(result, 1.0, abs_tol=TOLERANCE):
        result = 0.0
    return result


def parse(path: Path) -> list[tuple[float, tuple[float, float, float]]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "PROGRAM ENDED AT" not in text or "ENERGY| Total FORCE_EVAL" not in text:
        raise AssertionError(f"incomplete CP2K calculation: {path}")
    rows: list[tuple[float, tuple[float, float, float]]] = []
    for line in text.splitlines():
        match = POINT.match(line)
        if match:
            rows.append(
                (
                    float(match.group(2)),
                    tuple(canonical(float(match.group(i))) for i in range(3, 6)),
                )
            )
    return rows


def equivalent(left: tuple[float, ...], right: tuple[float, ...]) -> bool:
    return all(math.isclose(a, b, abs_tol=TOLERANCE) for a, b in zip(left, right))


def verify_mesh(mesh: int) -> None:
    rows = parse(ROOT / f"k{mesh}{mesh}{mesh}" / "cp2k.out")
    expected_points = [
        tuple(component / mesh for component in index)
        for index in itertools.product(range(mesh), repeat=3)
    ]
    if len(rows) != mesh**3:
        raise AssertionError(f"k{mesh}: expected {mesh**3} points, found {len(rows)}")
    expected_weight = 1.0 / mesh**3
    for weight, _ in rows:
        if not math.isclose(weight, expected_weight, abs_tol=TOLERANCE):
            raise AssertionError(
                f"k{mesh}: printed weight {weight} differs from {expected_weight}"
            )
    unmatched = list(expected_points)
    for _, point in rows:
        for index, candidate in enumerate(unmatched):
            if equivalent(point, candidate):
                unmatched.pop(index)
                break
        else:
            raise AssertionError(f"k{mesh}: unexpected reciprocal point {point}")
    if unmatched:
        raise AssertionError(f"k{mesh}: missing reciprocal points {unmatched}")
    print(
        f"mesh={mesh} points={len(rows)} equal_weights=true "
        "bvk_grid_equivalent=true status=pass"
    )


def main() -> None:
    for mesh in (2, 3, 4):
        verify_mesh(mesh)


if __name__ == "__main__":
    main()

