#!/usr/bin/env python3
"""Prove that direct-CLI and CP2K-native DMC structures are identical."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "DMC-ICE13" / "reproduction" / "seidler_dmc13_recalculation"
PREPARE = PACKAGE / "prepare_package.py"
TOLERANCE_ANGSTROM = 5.0e-13


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_prepare_module():
    spec = importlib.util.spec_from_file_location("gxtb_prepare_package", PREPARE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {PREPARE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    prepare = load_prepare_module()
    rows: list[dict[str, object]] = []
    maximum_cell = 0.0
    maximum_coordinate = 0.0

    for phase_dir in sorted((PACKAGE / "structures" / "primitive").iterdir()):
        if not phase_dir.is_dir():
            continue
        phase = phase_dir.name
        cp2k_input = PACKAGE / "raw" / "cp2k_native" / "k222-reduced" / phase / "input.inp"
        if not cp2k_input.is_file():
            cp2k_input = PACKAGE / "raw" / "cp2k_native" / "k333-reduced" / phase / "input.inp"
        poscar = phase_dir / "POSCAR"

        cp2k_species, cp2k_cell, cp2k_coordinates = prepare.parse_cp2k_structure(cp2k_input)
        cli_elements, cli_counts, cli_cell, cli_coordinates = prepare.parse_poscar(poscar)
        cli_species = tuple(
            element for element, count in zip(cli_elements, cli_counts) for _ in range(count)
        )
        if tuple(cp2k_species) != cli_species:
            raise RuntimeError(f"{phase}: species order differs")

        cell_difference = prepare.maximum_difference(cp2k_cell, cli_cell)
        coordinate_difference = prepare.maximum_difference(cp2k_coordinates, cli_coordinates)
        maximum_cell = max(maximum_cell, cell_difference)
        maximum_coordinate = max(maximum_coordinate, coordinate_difference)
        rows.append({
            "phase": phase,
            "atom_count": len(cp2k_species),
            "maximum_cell_difference_angstrom": cell_difference,
            "maximum_cartesian_coordinate_difference_angstrom": coordinate_difference,
            "cp2k_input_sha256": sha256(cp2k_input),
            "primitive_poscar_sha256": sha256(poscar),
        })

    maximum = max(maximum_cell, maximum_coordinate)
    status = "PASS" if len(rows) == 13 and maximum <= TOLERANCE_ANGSTROM else "FAIL"
    report = {
        "status": status,
        "phase_count_including_Ih": len(rows),
        "maximum_cell_difference_angstrom": maximum_cell,
        "maximum_cartesian_coordinate_difference_angstrom": maximum_coordinate,
        "tolerance_angstrom": TOLERANCE_ANGSTROM,
        "interpretation": (
            "The primitive structures used to generate every direct-CLI BvK supercell are "
            "identical to the corresponding CP2K-native cells and Cartesian coordinates at "
            "floating-point roundoff. Geometry conversion cannot explain the residual energy "
            "difference."
        ),
        "phases": rows,
    }
    output = Path(__file__).resolve().parent / "verification.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
