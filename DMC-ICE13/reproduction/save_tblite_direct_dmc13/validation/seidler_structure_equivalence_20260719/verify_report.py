#!/usr/bin/env python3
"""Verify the archived author-structure equivalence report and live inputs."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
STRUCTURE_ROOT = HERE.parents[2] / "seidler_dmc13_recalculation/structures/primitive"
PHASES = {"Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"}
ARCHIVE_SHA256 = "716aeda1d664d6d71d56b3ce1ff9a412d9fac7aab9cea428caa54bed6a9bd600"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    payload = json.loads((HERE / "report.json").read_text(encoding="utf-8"))
    assert payload["status"] == "PASS"
    assert payload["archive_sha256"] == ARCHIVE_SHA256
    assert payload["phase_count"] == len(PHASES)
    rows = payload["rows"]
    assert {row["phase"] for row in rows} == PHASES
    assert len(rows) == len(PHASES)
    metric_tolerance = float(payload["tolerances"]["lattice_metric_angstrom2"])
    atom_tolerance = float(payload["tolerances"]["matched_atom_angstrom"])
    assert metric_tolerance == 1.0e-8
    assert atom_tolerance == 1.0e-8
    for row in rows:
        values = (
            float(row["lattice_metric_max_abs_residual_angstrom2"]),
            float(row["rotation_orthogonality_max_abs_residual"]),
            float(row["volume_abs_difference_angstrom3"]),
            float(row["maximum_matched_atom_distance_angstrom"]),
        )
        assert all(math.isfinite(value) and value >= 0.0 for value in values)
        assert values[0] <= metric_tolerance
        assert values[3] <= atom_tolerance
        poscar = STRUCTURE_ROOT / row["phase"] / "POSCAR"
        assert digest(poscar) == row["production_poscar_sha256"]
        transform = row["lattice_unimodular_transform"]
        determinant = (
            transform[0][0] * (transform[1][1] * transform[2][2] - transform[1][2] * transform[2][1])
            - transform[0][1] * (transform[1][0] * transform[2][2] - transform[1][2] * transform[2][0])
            + transform[0][2] * (transform[1][0] * transform[2][1] - transform[1][1] * transform[2][0])
        )
        assert abs(determinant) == 1
    summary = payload["summary"]
    assert float(summary["maximum_lattice_metric_abs_residual_angstrom2"]) == max(
        float(row["lattice_metric_max_abs_residual_angstrom2"]) for row in rows
    )
    assert float(summary["maximum_matched_atom_distance_angstrom"]) == max(
        float(row["maximum_matched_atom_distance_angstrom"]) for row in rows
    )
    print(
        "author DMC-ICE13 structures: 13/13 equivalent; "
        f"max atom residual={summary['maximum_matched_atom_distance_angstrom']:.12e} Angstrom"
    )


if __name__ == "__main__":
    main()
