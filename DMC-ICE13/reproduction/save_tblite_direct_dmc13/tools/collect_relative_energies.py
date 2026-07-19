#!/usr/bin/env python3
"""Derive Ih-referenced DMC-ICE13 energies from the absolute-energy oracle."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HARTREE_TO_KJMOL = 2625.499638
PROVIDERS = {
    "cp2k_native": "cp2k_native_energy_Ha_per_primitive",
    "save_tblite_cli": "save_tblite_cli_energy_Ha_per_primitive",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    absolute_path = ROOT / "tables" / "absolute_energies_vs_mesh.csv"
    reference_path = ROOT / "tables" / "dmc_reference_relative_energies.csv"
    absolute = load_rows(absolute_path)
    references = {
        row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
        for row in load_rows(reference_path)
    }
    by_mesh: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in absolute:
        by_mesh[row["mesh_id"]].append(row)

    output_rows: list[dict[str, object]] = []
    for mesh_id in sorted(by_mesh, key=lambda value: int(value[1:])):
        mesh_rows = by_mesh[mesh_id]
        by_phase = {row["phase"]: row for row in mesh_rows}
        for provider, energy_column in PROVIDERS.items():
            ih_text = by_phase["Ih"][energy_column]
            if not ih_text:
                continue
            ih_energy_per_water = float(ih_text) / (int(by_phase["Ih"]["natom_primitive"]) / 3)
            for row in mesh_rows:
                energy_text = row[energy_column]
                if not energy_text:
                    continue
                waters = int(row["natom_primitive"]) / 3
                relative = HARTREE_TO_KJMOL * (float(energy_text) / waters - ih_energy_per_water)
                reference = references[row["phase"]]
                source_hash = (
                    row["cp2k_output_sha256"]
                    if provider == "cp2k_native"
                    else row["save_tblite_json_sha256"]
                )
                output_rows.append(
                    {
                        "mesh_n": row["mesh_n"],
                        "mesh_id": mesh_id,
                        "provider": provider,
                        "phase": row["phase"],
                        "n_H2O_primitive": f"{waters:g}",
                        "absolute_energy_Ha_per_primitive": energy_text,
                        "relative_energy_kJmol_per_H2O": f"{relative:.12f}",
                        "DMC_reference_kJmol_per_H2O": f"{reference:.6f}",
                        "signed_error_kJmol_per_H2O": f"{relative - reference:+.12f}",
                        "reference_phase": "Ih",
                        "source_sha256": source_hash,
                        "status": "recorded",
                    }
                )

    output = ROOT / "tables" / "relative_energies_vs_mesh.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
