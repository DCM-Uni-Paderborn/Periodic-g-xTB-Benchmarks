#!/usr/bin/env python3
"""Verify Cartesian structures and DMC/native/CLI energy tables."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
HARTREE_TO_KJMOL = 2625.499638


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    references = rows(ROOT / "tables" / "dmc_reference_relative_energies.csv")
    assert tuple(row["phase"] for row in references) == PHASES
    reference = {
        row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
        for row in references
    }
    assert reference["Ih"] == 0.0

    manifest = rows(ROOT / "provenance" / "structure_hashes.csv")
    assert len(manifest) == 4 * len(PHASES)
    for row in manifest:
        mesh = int(row["mesh"])
        phase = row["phase"]
        assert phase in PHASES and mesh in (1, 2, 3, 4)
        directory = ROOT / "structures" / f"k{mesh}{mesh}{mesh}" / phase
        poscar = directory / "POSCAR"
        extxyz = directory / "structure.xyz"
        assert poscar.read_text().splitlines()[7].strip().lower().startswith("cart")
        xyz_lines = extxyz.read_text().splitlines()
        assert int(xyz_lines[0]) == int(row["natoms"])
        assert len(xyz_lines) == int(row["natoms"]) + 2
        assert 'pbc="T T T"' in xyz_lines[1] and 'Lattice="' in xyz_lines[1]
        assert digest(poscar) == row["poscar_sha256"]
        assert digest(extxyz) == row["extxyz_sha256"]

    absolute = rows(ROOT / "tables" / "absolute_energies_vs_mesh.csv")
    relative = rows(ROOT / "tables" / "relative_energies_vs_mesh.csv")
    absolute_index = {(row["mesh_id"], row["phase"]): row for row in absolute}
    assert len(absolute_index) == 4 * len(PHASES)
    assert {row["provider"] for row in relative} == {"cp2k_native", "save_tblite_cli"}
    assert all("GFN1" not in str(row) and "GFN2" not in str(row) for row in relative)

    for row in relative:
        source = absolute_index[(row["mesh_id"], row["phase"])]
        column = (
            "cp2k_native_energy_Ha_per_primitive"
            if row["provider"] == "cp2k_native"
            else "save_tblite_cli_energy_Ha_per_primitive"
        )
        hash_column = (
            "cp2k_output_sha256"
            if row["provider"] == "cp2k_native"
            else "save_tblite_json_sha256"
        )
        assert row["absolute_energy_Ha_per_primitive"] == source[column]
        assert row["source_sha256"] == source[hash_column] and row["source_sha256"]
        assert float(row["DMC_reference_kJmol_per_H2O"]) == reference[row["phase"]]
        ih = absolute_index[(row["mesh_id"], "Ih")]
        waters = int(source["natom_primitive"]) / 3
        ih_waters = int(ih["natom_primitive"]) / 3
        expected = HARTREE_TO_KJMOL * (
            float(source[column]) / waters - float(ih[column]) / ih_waters
        )
        observed = float(row["relative_energy_kJmol_per_H2O"])
        assert abs(observed - expected) < 5.0e-10
        assert abs(float(row["signed_error_kJmol_per_H2O"]) - (expected - reference[row["phase"]])) < 5.0e-10

    print(
        f"PASS phases={len(PHASES)} cartesian_structures={len(manifest)} "
        f"absolute_rows={len(absolute)} relative_rows={len(relative)}"
    )


if __name__ == "__main__":
    main()
