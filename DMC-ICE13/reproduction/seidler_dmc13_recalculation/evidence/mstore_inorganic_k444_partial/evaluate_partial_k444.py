#!/usr/bin/env python3
"""Qualify the completed subset of the historical mstore-inorganic 4^3 matrix.

The XIII calculation was terminated by the operating system before its first
SCC iteration.  This verifier therefore treats the data as a source-state
convergence diagnostic only.  A PASS proves the integrity of all completed
points and of the recorded failed endpoint; it must never be interpreted as a
complete DMC-ICE13 benchmark statistic.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from decimal import Decimal, getcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
if HERE.parent.name == "evidence" and (HERE.parents[1] / "tables").is_dir():
    # Self-contained copy inside the Seidler recalculation package.
    PACKAGE = HERE.parents[1]
    ROOT = PACKAGE.parents[2]
else:
    # Canonical validation copy at repository/validation/<evidence>.
    ROOT = HERE.parents[1]
    PACKAGE = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII",
    "IX", "XI", "XIII", "XIV", "XV", "XVII",
)
COMPLETED_PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII",
    "IX", "XI", "XIV", "XV", "XVII",
)
FAILED_PHASE = "XIII"
EXPECTED_BINARY_SHA256 = (
    "8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa"
)
MESH = 4
REPLICAS = MESH**3
HARTREE_TO_KJ_MOL = Decimal("2625.4996394798254")
EXPECTED_ENERGY_THRESHOLD = Decimal("1e-7")
EXPECTED_DENSITY_THRESHOLD = Decimal("2e-6")
EXPECTED_METRICS = {
    "mstore_partial_mae_k4": Decimal("7.004262971119"),
    "pbc_partial_mae_k4": Decimal("12.217167057179"),
    "mean_abs_branch_gap_k4": Decimal("5.245588734356"),
    "max_abs_branch_gap_k4": Decimal("31.034708379013"),
    "mean_abs_mstore_delta_k3_to_k4": Decimal("12.606362232512"),
}
METRIC_TOLERANCE = Decimal("5e-10")
getcontext().prec = 40


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decimal_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def threshold(text: str, label: str) -> Decimal:
    match = re.search(rf"^{re.escape(label)}\s+([0-9.Ee+-]+)\s+", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {label!r} in process output")
    return Decimal(match.group(1))


def recorded_hash(path: Path) -> str:
    fields = path.read_text(encoding="utf-8").strip().split()
    if not fields or re.fullmatch(r"[0-9a-f]{64}", fields[0]) is None:
        raise RuntimeError(f"invalid SHA-256 record: {path}")
    return fields[0]


def mean(values: list[Decimal]) -> Decimal:
    return sum(values, Decimal(0)) / len(values)


def main() -> None:
    manifest = {
        row["phase"]: row
        for row in read_csv(PACKAGE / "structures/structure_manifest.csv")
    }
    references = {
        row["phase"]: Decimal(row["reference_relative_energy_kJmol_per_H2O"])
        for row in read_csv(PACKAGE / "tables/dmc_reference_relative_energies.csv")
    }
    earlier_rows = read_csv(
        PACKAGE / "tables/mstore_inorganic_relative_energies_by_mesh.csv"
    )
    earlier = {
        (int(row["mesh_n"]), row["phase"]): Decimal(
            row["relative_energy_kj_mol_per_H2O"]
        )
        for row in earlier_rows
    }
    pbc_rows = read_csv(PACKAGE / "tables/cp2k_native_relative_energies_by_mesh.csv")
    pbc = {
        (int(row["mesh_n"]), row["phase"]): Decimal(
            row["relative_energy_kj_mol_per_H2O"]
        )
        for row in pbc_rows
        if row["qualification"] == "PASS"
    }

    checks: dict[str, bool] = {}
    energies: dict[str, Decimal] = {}
    water_counts: dict[str, int] = {}
    raw_hashes: dict[str, dict[str, str]] = {}
    absolute_rows: list[dict[str, str]] = []

    for phase in PHASES:
        run = HERE / "results" / phase
        structure = HERE / "structures" / phase / "POSCAR"
        common = (
            run / "process.out",
            run / "process.err",
            run / "exit_status",
            run / "binary.sha256",
            run / "input.sha256",
            run / "command.json",
            structure,
        )
        if not all(path.is_file() for path in common):
            raise RuntimeError(f"incomplete raw record for {phase}")
        output = (run / "process.out").read_text(encoding="utf-8")
        command = json.loads((run / "command.json").read_text(encoding="utf-8"))
        primitive_atoms = int(manifest[phase]["atom_count"])
        actual_atoms = sum(
            int(value) for value in structure.read_text(encoding="utf-8").splitlines()[6].split()
        )
        phase_checks = {
            "binary_hash": recorded_hash(run / "binary.sha256")
            == EXPECTED_BINARY_SHA256,
            "input_hash": recorded_hash(run / "input.sha256") == sha256(structure),
            "mesh_atom_count": actual_atoms == primitive_atoms * REPLICAS,
            "method": command[command.index("--method") + 1] == "gxtb",
            "accuracy": Decimal(command[command.index("--acc") + 1])
            == Decimal("0.1"),
            "iterations": command[command.index("--iterations") + 1] == "300",
            "restart_disabled": "--no-restart" in command,
            "json_target": command[command.index("--json") + 1] == "result.json",
        }

        result_path = run / "result.json"
        if phase == FAILED_PHASE:
            phase_checks.update(
                {
                    "recorded_exit_minus_nine": (run / "exit_status")
                    .read_text(encoding="utf-8")
                    .strip()
                    == "-9",
                    "no_result_json": not result_path.exists(),
                    "thresholds_reached_before_termination": threshold(
                        output, "energy convergence"
                    )
                    == EXPECTED_ENERGY_THRESHOLD
                    and threshold(output, "density convergence")
                    == EXPECTED_DENSITY_THRESHOLD,
                    "terminated_before_first_scc_result": "cycle        total energy" in output
                    and "total energy                  " not in output,
                }
            )
            raw_files = common
        else:
            if not result_path.is_file():
                raise RuntimeError(f"missing successful result for {phase}")
            result = decimal_json(result_path)
            energy = Decimal(result["energy"])
            primitive_waters = int(manifest[phase]["water_molecule_count"])
            waters = primitive_waters * REPLICAS
            phase_checks.update(
                {
                    "exit_zero": (run / "exit_status")
                    .read_text(encoding="utf-8")
                    .strip()
                    == "0",
                    "energy_threshold": threshold(output, "energy convergence")
                    == EXPECTED_ENERGY_THRESHOLD,
                    "density_threshold": threshold(output, "density convergence")
                    == EXPECTED_DENSITY_THRESHOLD,
                    "normal_scc_completion": "total energy                  " in output
                    and "JSON dump of results written" in output,
                    "json_version": result.get("version") == "0.5.0",
                    "finite_energy": math.isfinite(float(energy)),
                }
            )
            energies[phase] = energy
            water_counts[phase] = waters
            raw_files = common + (result_path,)
            absolute_rows.append(
                {
                    "mesh_n": str(MESH),
                    "phase": phase,
                    "water_molecules_primitive": str(primitive_waters),
                    "water_molecules_supercell": str(waters),
                    "mstore_energy_Ha_supercell": format(energy, ".15f"),
                    "mstore_energy_Ha_per_primitive": format(
                        energy / REPLICAS, ".15f"
                    ),
                    "qualification": "PASS" if all(phase_checks.values()) else "FAIL",
                }
            )

        checks.update({f"{phase}:{name}": value for name, value in phase_checks.items()})
        raw_hashes[phase] = {
            path.relative_to(HERE).as_posix(): sha256(path) for path in raw_files
        }

    completed_set = tuple(phase for phase in PHASES if phase in energies)
    checks["exact_completed_phase_set"] = completed_set == COMPLETED_PHASES
    checks["failed_phase_excluded_from_statistics"] = FAILED_PHASE not in energies
    checks["full_matrix_explicitly_incomplete"] = len(energies) == 12
    checks["ih_reference_completed"] = "Ih" in energies

    ih_per_water = energies["Ih"] / water_counts["Ih"]
    comparison_rows: list[dict[str, str]] = []
    mstore_errors: list[Decimal] = []
    pbc_errors: list[Decimal] = []
    branch_gaps: list[Decimal] = []
    mstore_deltas: list[Decimal] = []
    phase_metrics: dict[str, dict[str, Decimal]] = {}
    for phase in COMPLETED_PHASES[1:]:
        relative = (
            energies[phase] / water_counts[phase] - ih_per_water
        ) * HARTREE_TO_KJ_MOL
        reference = references[phase]
        pbc_relative = pbc[(MESH, phase)]
        mstore_error = relative - reference
        pbc_error = pbc_relative - reference
        branch_gap = relative - pbc_relative
        delta = relative - earlier[(3, phase)]
        mstore_errors.append(mstore_error)
        pbc_errors.append(pbc_error)
        branch_gaps.append(branch_gap)
        mstore_deltas.append(delta)
        phase_metrics[phase] = {
            "mstore_relative": relative,
            "pbc_relative": pbc_relative,
            "reference": reference,
            "mstore_error": mstore_error,
            "pbc_error": pbc_error,
            "branch_gap": branch_gap,
            "mstore_delta_k3_to_k4": delta,
        }
        comparison_rows.append(
            {
                "mesh_n": str(MESH),
                "phase": phase,
                "historical_mstore_relative_kj_mol_per_H2O": format(relative, ".12f"),
                "current_pbc_native_relative_kj_mol_per_H2O": format(
                    pbc_relative, ".12f"
                ),
                "mstore_minus_pbc_kj_mol_per_H2O": format(branch_gap, ".12f"),
                "dmc_reference_kj_mol_per_H2O": format(reference, ".6f"),
                "mstore_absolute_error_kj_mol_per_H2O": format(
                    abs(mstore_error), ".12f"
                ),
                "pbc_absolute_error_kj_mol_per_H2O": format(abs(pbc_error), ".12f"),
                "mstore_delta_k3_to_k4_kj_mol_per_H2O": format(delta, ".12f"),
                "qualification": "PASS",
            }
        )

    metrics = {
        "mstore_partial_mae_k4": mean([abs(value) for value in mstore_errors]),
        "pbc_partial_mae_k4": mean([abs(value) for value in pbc_errors]),
        "mean_abs_branch_gap_k4": mean([abs(value) for value in branch_gaps]),
        "max_abs_branch_gap_k4": max(abs(value) for value in branch_gaps),
        "mean_abs_mstore_delta_k3_to_k4": mean(
            [abs(value) for value in mstore_deltas]
        ),
    }
    max_gap_phase = max(
        phase_metrics,
        key=lambda phase: abs(phase_metrics[phase]["branch_gap"]),
    )
    checks["max_branch_gap_phase_is_VII"] = max_gap_phase == "VII"
    for name, expected in EXPECTED_METRICS.items():
        checks[f"metric_reproduces:{name}"] = abs(metrics[name] - expected) <= METRIC_TOLERANCE

    mesh_trends = []
    for mesh in (2, 3, 4):
        if mesh == 4:
            mstore_values = {
                phase: phase_metrics[phase]["mstore_relative"]
                for phase in COMPLETED_PHASES[1:]
            }
        else:
            mstore_values = {
                phase: earlier[(mesh, phase)] for phase in COMPLETED_PHASES[1:]
            }
        gaps = [mstore_values[phase] - pbc[(mesh, phase)] for phase in mstore_values]
        m_errors = [mstore_values[phase] - references[phase] for phase in mstore_values]
        p_errors = [pbc[(mesh, phase)] - references[phase] for phase in mstore_values]
        gap_phase = max(mstore_values, key=lambda phase: abs(mstore_values[phase] - pbc[(mesh, phase)]))
        mesh_trends.append(
            {
                "mesh_n": mesh,
                "phase_count": len(mstore_values),
                "mstore_partial_mae_kj_mol_per_H2O": float(
                    mean([abs(value) for value in m_errors])
                ),
                "pbc_partial_mae_kj_mol_per_H2O": float(
                    mean([abs(value) for value in p_errors])
                ),
                "mean_abs_branch_gap_kj_mol_per_H2O": float(
                    mean([abs(value) for value in gaps])
                ),
                "max_abs_branch_gap_kj_mol_per_H2O": float(
                    max(abs(value) for value in gaps)
                ),
                "max_gap_phase": gap_phase,
            }
        )
    mean_gaps = [Decimal(str(row["mean_abs_branch_gap_kj_mol_per_H2O"])) for row in mesh_trends]
    checks["branch_gap_decreases_monotonically_k2_to_k4"] = (
        mean_gaps[0] > mean_gaps[1] > mean_gaps[2]
    )
    checks["all_completed_raw_qualifications_pass"] = all(
        row["qualification"] == "PASS" for row in absolute_rows
    )
    checks["all_derived_values_finite"] = all(
        math.isfinite(float(value)) for value in metrics.values()
    )

    absolute_path = HERE / "partial_absolute_energies.csv"
    with absolute_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(absolute_rows[0]))
        writer.writeheader()
        writer.writerows(absolute_rows)
    comparison_path = HERE / "partial_comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0]))
        writer.writeheader()
        writer.writerows(comparison_rows)

    passed = all(checks.values())
    payload = {
        "schema": "periodic-gxtb-mstore-inorganic-k444-partial-v1",
        "status": "PASS" if passed else "FAIL",
        "mesh_n": MESH,
        "completed_phase_count_including_ih": len(energies),
        "completed_benchmark_phase_count": len(mstore_errors),
        "completed_phases": list(COMPLETED_PHASES),
        "failed_phase": FAILED_PHASE,
        "failed_phase_exit_status": -9,
        "full_matrix_complete": False,
        "usable_for_full_benchmark_statistics": False,
        "partial_metrics": {name: float(value) for name, value in metrics.items()},
        "max_branch_gap_phase_k4": max_gap_phase,
        "mesh_trends_same_eleven_phases": mesh_trends,
        "checks": checks,
        "raw_sha256": raw_hashes,
        "derived_sha256": {
            absolute_path.name: sha256(absolute_path),
            comparison_path.name: sha256(comparison_path),
        },
        "interpretation": (
            "Eleven non-reference phases plus Ih form an integrity-qualified historical "
            "mstore-inorganic 4x4x4 subset. XIII was killed before the first SCC result, "
            "so the reported statistics are explicitly same-eleven-phase diagnostics, "
            "not a complete DMC-ICE13 MAE. The mstore/pbc branch gap decreases strongly "
            "from 2x2x2 to 4x4x4, while the historical 3x3x3-to-4x4x4 changes remain "
            "large; the lower sparse-mesh mstore error is therefore finite-size error "
            "cancellation rather than evidence that CP2K-native pbc is incorrect."
        ),
    }
    (HERE / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
