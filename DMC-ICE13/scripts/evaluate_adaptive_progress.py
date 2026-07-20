#!/usr/bin/env python3
"""Reconstruct the qualified phase-wise DMC-ICE13 progress snapshots."""

from __future__ import annotations

import csv
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = (
    ROOT
    / "reproduction/seidler_dmc13_recalculation/tables/"
    "cp2k_native_relative_energies_by_mesh.csv"
)
CURRENT_SET = ROOT / "data/dmc_ice13_gxtb_current_adaptive_set.csv"
CURRENT_STATISTICS = ROOT / "data/dmc_ice13_gxtb_current_adaptive_statistics.csv"
OUTPUT = ROOT / "data/dmc_ice13_gxtb_phasewise_progress.csv"
PHASES = (
    "II", "III", "IV", "VI", "VII", "VIII",
    "IX", "XI", "XIII", "XIV", "XV", "XVII",
)
LIMITS = (6, 7, 8)
MAXIMUM_MESH = max(LIMITS)
TOLERANCE = 0.10


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=tuple(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    by_phase: dict[str, dict[int, dict[str, str]]] = {
        phase: {} for phase in PHASES
    }
    for row in read_rows(SOURCE):
        phase = row["phase"]
        if phase not in by_phase or row["qualification"] != "PASS":
            continue
        mesh = int(row["mesh_n"])
        if mesh in by_phase[phase]:
            raise RuntimeError(f"duplicate endpoint: phase={phase} mesh={mesh}")
        by_phase[phase][mesh] = row

    output_rows: list[dict[str, object]] = []
    last_selection: dict[str, tuple[int, float, float | None, bool]] = {}
    for limit in LIMITS:
        selections: dict[str, tuple[int, float, float | None, bool]] = {}
        unresolved: list[str] = []
        pending: list[str] = []
        errors: list[float] = []
        for phase in PHASES:
            endpoints = {
                mesh: row
                for mesh, row in by_phase[phase].items()
                if mesh <= limit
            }
            if not endpoints:
                raise RuntimeError(
                    f"no qualified endpoint through mesh {limit} for phase {phase}"
                )
            selected_mesh = max(endpoints)
            adjacent_delta: float | None = None
            converged = False
            for mesh in sorted(endpoints):
                if mesh - 1 not in endpoints:
                    continue
                delta = abs(
                    float(endpoints[mesh]["relative_energy_kj_mol_per_H2O"])
                    - float(
                        endpoints[mesh - 1]["relative_energy_kj_mol_per_H2O"]
                    )
                )
                if delta <= TOLERANCE:
                    selected_mesh = mesh
                    adjacent_delta = delta
                    converged = True
                    break
            if not converged:
                if selected_mesh - 1 in endpoints:
                    adjacent_delta = abs(
                        float(
                            endpoints[selected_mesh][
                                "relative_energy_kj_mol_per_H2O"
                            ]
                        )
                        - float(
                            endpoints[selected_mesh - 1][
                                "relative_energy_kj_mol_per_H2O"
                            ]
                        )
                    )
                unresolved.append(phase)
                next_mesh = selected_mesh + 1
                if next_mesh <= MAXIMUM_MESH:
                    pending.append(f"{phase}:{next_mesh}")

            error = float(endpoints[selected_mesh]["error_kj_mol_per_H2O"])
            if not math.isfinite(error):
                raise RuntimeError(
                    f"non-finite error: phase={phase} mesh={selected_mesh}"
                )
            selections[phase] = (
                selected_mesh,
                error,
                adjacent_delta,
                converged,
            )
            errors.append(error)

        output_rows.append(
            {
                "mesh_limit_n": limit,
                "mesh_label": f"phase-wise <= {limit}x{limit}x{limit}",
                "phase_count": len(errors),
                "converged_phase_count": sum(
                    selection[3] for selection in selections.values()
                ),
                "me_kj_mol_per_water": f"{sum(errors) / len(errors):.12f}",
                "mae_kj_mol_per_water": (
                    f"{sum(abs(value) for value in errors) / len(errors):.12f}"
                ),
                "rmse_kj_mol_per_water": (
                    f"{math.sqrt(sum(value * value for value in errors) / len(errors)):.12f}"
                ),
                "maxae_kj_mol_per_water": f"{max(abs(value) for value in errors):.12f}",
                "selected_meshes": ";".join(
                    f"{phase}:{selections[phase][0]}" for phase in PHASES
                ),
                "unresolved_phases": ";".join(unresolved),
                "pending_next_endpoints": ";".join(pending),
                "qualification": "PASS",
            }
        )
        last_selection = selections

    current_rows = {row["phase"]: row for row in read_rows(CURRENT_SET)}
    if set(current_rows) != set(PHASES):
        raise RuntimeError("current adaptive set has an invalid phase matrix")
    for phase, selection in last_selection.items():
        mesh, error, delta, converged = selection
        archived = current_rows[phase]
        archived_delta = archived[
            "absolute_adjacent_delta_kj_mol_per_water"
        ].strip()
        if mesh != int(archived["mesh_n"]):
            raise RuntimeError(f"current selected mesh mismatch for phase {phase}")
        if not math.isclose(
            error,
            float(archived["error_kj_mol_per_water"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError(f"current error mismatch for phase {phase}")
        if (delta is None) != (not archived_delta):
            raise RuntimeError(f"current delta availability mismatch for phase {phase}")
        if delta is not None and not math.isclose(
            delta,
            float(archived_delta),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError(f"current delta mismatch for phase {phase}")
        if converged != (archived["phase_converged"].lower() == "true"):
            raise RuntimeError(f"current convergence mismatch for phase {phase}")

    current_statistics = read_rows(CURRENT_STATISTICS)
    if len(current_statistics) != 1:
        raise RuntimeError("current adaptive statistics must contain one row")
    archived_statistics = current_statistics[0]
    latest = output_rows[-1]
    for field in (
        "phase_count",
        "converged_phase_count",
        "me_kj_mol_per_water",
        "mae_kj_mol_per_water",
        "rmse_kj_mol_per_water",
        "maxae_kj_mol_per_water",
    ):
        if not math.isclose(
            float(latest[field]),
            float(archived_statistics[field]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError(f"current aggregate mismatch in {field}")

    write_rows(OUTPUT, output_rows)
    print(
        " ".join(
            f"<={row['mesh_limit_n']}^3:{row['mae_kj_mol_per_water']}"
            for row in output_rows
        )
    )


if __name__ == "__main__":
    main()
