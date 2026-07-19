#!/usr/bin/env python3
"""Generate the complete model-revision DMC-ICE13 comparison tables."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
    "XIII", "XIV", "XV", "XVII",
)
N_WATER = {
    "Ih": 12, "II": 12, "III": 12, "IV": 16, "VI": 10, "VII": 12,
    "VIII": 8, "IX": 12, "XI": 8, "XIII": 28, "XIV": 12, "XV": 10,
    "XVII": 6,
}
REFERENCE = {
    "II": 0.31, "III": 1.25, "IV": 3.83, "VI": 1.78, "VII": 4.99,
    "VIII": 4.23, "IX": 0.60, "XI": 0.16, "XIII": 2.12, "XIV": 1.70,
    "XV": 1.74, "XVII": 1.75,
}
HARTREE_TO_KJMOL = 2625.499638
MESH_PROVIDERS = {
    1: ("current", "legacy_mstore_inorganic", "gxtb_v201", "dcm_main"),
    2: ("current", "legacy_mstore_inorganic", "gxtb_v201", "dcm_main"),
    3: ("current", "final_pbc", "legacy_mstore_inorganic"),
}


def result_path(provider: str, mesh: int, phase: str) -> Path:
    mesh_id = f"k{mesh}{mesh}{mesh}"
    if provider == "current":
        return PACKAGE / "results/current_save_tblite_cli" / mesh_id / phase / "result.json"
    if provider == "legacy_mstore_inorganic":
        if mesh == 3:
            return ROOT / "raw/authors_exchange_linux_k333" / phase / "result.json"
        return ROOT / "raw/authors_exchange" / mesh_id / phase / "result.json"
    if provider == "final_pbc":
        return (
            ROOT.parent / "provider_revision_bvk_ab_20260718"
            / "seidler_pbc_cli_linux" / mesh_id / phase / "result.json"
        )
    return ROOT / "raw" / provider / mesh_id / phase / "result.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def energy(provider: str, mesh: int, phase: str) -> float:
    return float(json.loads(result_path(provider, mesh, phase).read_text())["energy"]) / mesh**3


def relative(provider: str, mesh: int, phase: str) -> float:
    return (
        energy(provider, mesh, phase) / N_WATER[phase]
        - energy(provider, mesh, "Ih") / N_WATER["Ih"]
    ) * HARTREE_TO_KJMOL


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


mae_rows: list[dict[str, object]] = []
mae: dict[tuple[int, str], float] = {}
for mesh, providers in MESH_PROVIDERS.items():
    for provider in providers:
        value = statistics.mean(
            abs(relative(provider, mesh, phase) - REFERENCE[phase])
            for phase in PHASES[1:]
        )
        mae[(mesh, provider)] = value
        mae_rows.append(
            {
                "mesh": mesh,
                "provider": provider,
                "mae_kj_mol_per_water": f"{value:.12f}",
            }
        )
write_csv(
    ROOT / "coarse_mae_summary.csv",
    ("mesh", "provider", "mae_kj_mol_per_water"),
    mae_rows,
)

for mesh in (2, 3):
    providers = MESH_PROVIDERS[mesh]
    fields = (
        "phase",
        "dmc_reference_kj_mol_per_water",
        *(f"{provider}_kj_mol_per_water" for provider in providers),
    )
    rows: list[dict[str, object]] = []
    for phase in PHASES[1:]:
        row: dict[str, object] = {
            "phase": phase,
            "dmc_reference_kj_mol_per_water": f"{REFERENCE[phase]:.12f}",
        }
        for provider in providers:
            row[f"{provider}_kj_mol_per_water"] = f"{relative(provider, mesh, phase):.12f}"
        rows.append(row)
    write_csv(ROOT / f"relative_energies_k{mesh}{mesh}{mesh}.csv", fields, rows)

manifest_rows: list[dict[str, object]] = []
for provider in MESH_PROVIDERS[3]:
    for phase in PHASES:
        path = result_path(provider, 3, phase)
        manifest_rows.append(
            {
                "mesh": 3,
                "phase": phase,
                "provider": provider,
                "path_from_package_root": path.relative_to(PACKAGE),
                "sha256": sha256(path),
            }
        )
write_csv(
    ROOT / "k333_input_manifest.csv",
    ("mesh", "phase", "provider", "path_from_package_root", "sha256"),
    manifest_rows,
)

shifts = [
    (
        abs(relative("legacy_mstore_inorganic", 3, phase) - relative("final_pbc", 3, phase)),
        phase,
    )
    for phase in PHASES[1:]
]
maximum_shift, maximum_phase = max(shifts)
summary = {
    "current_mae_kj_mol_per_water": mae[(3, "current")],
    "final_pbc_mae_kj_mol_per_water": mae[(3, "final_pbc")],
    "final_pbc_minus_current_mae_kj_mol_per_water": (
        mae[(3, "final_pbc")] - mae[(3, "current")]
    ),
    "legacy_mstore_inorganic_mae_kj_mol_per_water": mae[(3, "legacy_mstore_inorganic")],
    "legacy_minus_final_pbc_mae_kj_mol_per_water": (
        mae[(3, "legacy_mstore_inorganic")] - mae[(3, "final_pbc")]
    ),
    "maximum_legacy_final_relative_energy_shift_kj_mol_per_water": maximum_shift,
    "maximum_shift_phase": maximum_phase,
    "mesh": 3,
    "phase_count_including_Ih": len(PHASES),
    "relative_phase_count": len(PHASES) - 1,
    "status": "PASS",
}
(ROOT / "k333_revision_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(summary, indent=2, sort_keys=True))
