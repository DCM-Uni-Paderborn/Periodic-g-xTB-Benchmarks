#!/usr/bin/env python3
"""Select the first passing adjacent DMC-ICE13 mesh for every phase."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("adaptive_output", type=Path)
    parser.add_argument("statistics_output", type=Path)
    parser.add_argument("--tolerance", type=float, default=0.10)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0.0:
        parser.error("--tolerance must be finite and non-negative")

    by_phase: dict[str, dict[int, dict[str, str]]] = {phase: {} for phase in PHASES}
    for row in read_rows(args.source):
        phase = row["phase"]
        if phase not in by_phase or row["qualification"] != "PASS":
            continue
        mesh = int(row["mesh_n"])
        if mesh in by_phase[phase]:
            raise RuntimeError(f"duplicate qualified endpoint: phase={phase} mesh={mesh}")
        by_phase[phase][mesh] = row

    selected_rows: list[dict[str, object]] = []
    for phase in PHASES:
        endpoints = by_phase[phase]
        if not endpoints:
            raise RuntimeError(f"no qualified endpoint for phase {phase}")
        selected_mesh = max(endpoints)
        adjacent_delta: float | None = None
        converged = False
        for mesh in sorted(endpoints):
            if mesh - 1 not in endpoints:
                continue
            delta = abs(
                float(endpoints[mesh]["relative_energy_kj_mol_per_H2O"])
                - float(endpoints[mesh - 1]["relative_energy_kj_mol_per_H2O"])
            )
            if delta <= args.tolerance:
                selected_mesh = mesh
                adjacent_delta = delta
                converged = True
                break
        if not converged and selected_mesh - 1 in endpoints:
            adjacent_delta = abs(
                float(endpoints[selected_mesh]["relative_energy_kj_mol_per_H2O"])
                - float(endpoints[selected_mesh - 1]["relative_energy_kj_mol_per_H2O"])
            )

        selected = endpoints[selected_mesh]
        relative = float(selected["relative_energy_kj_mol_per_H2O"])
        reference = float(selected["dmc_reference_kj_mol_per_H2O"])
        error = relative - reference
        selected_rows.append({
            "phase": phase,
            "mesh_n": selected_mesh,
            "relative_kj_mol_per_water": f"{relative:.12f}",
            "dmc_reference_kj_mol_per_water": f"{reference:.6f}",
            "error_kj_mol_per_water": f"{error:.12f}",
            "absolute_error_kj_mol_per_water": f"{abs(error):.12f}",
            "absolute_adjacent_delta_kj_mol_per_water": (
                "" if adjacent_delta is None else f"{adjacent_delta:.12f}"
            ),
            "phase_converged": str(converged).lower(),
        })

    errors = [float(row["error_kj_mol_per_water"]) for row in selected_rows]
    converged_count = sum(row["phase_converged"] == "true" for row in selected_rows)
    statistics_rows = [{
        "scope": "current_qualified_phasewise_set",
        "phase_count": len(selected_rows),
        "converged_phase_count": converged_count,
        "largest_mesh_n": max(int(row["mesh_n"]) for row in selected_rows),
        "me_kj_mol_per_water": f"{sum(errors) / len(errors):.12f}",
        "mae_kj_mol_per_water": f"{sum(abs(value) for value in errors) / len(errors):.12f}",
        "rmse_kj_mol_per_water": f"{math.sqrt(sum(value * value for value in errors) / len(errors)):.12f}",
        "maxae_kj_mol_per_water": f"{max(abs(value) for value in errors):.12f}",
        "final_result": str(converged_count == len(selected_rows)).lower(),
    }]

    write_rows(args.adaptive_output, selected_rows, tuple(selected_rows[0]))
    write_rows(args.statistics_output, statistics_rows, tuple(statistics_rows[0]))
    print(
        f"phases={len(selected_rows)} converged={converged_count} "
        f"mae={statistics_rows[0]['mae_kj_mol_per_water']}"
    )


if __name__ == "__main__":
    main()
