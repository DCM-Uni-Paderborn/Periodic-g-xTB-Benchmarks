#!/usr/bin/env python3
"""Build fail-closed, fixed-mesh DMC-ICE13 MAE data from CP2K runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from bvk_input import input_mesh_and_water_count
from dmc_mixed_mae import (
    HARTREE_TO_KJMOL,
    PHASES,
    SHA256_RE,
    binary_digest,
    energy,
    recorded_digest,
    sha256,
)


def references(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {
            row["phase"]: float(row["DMC_relative_kJmol"])
            for row in csv.DictReader(handle)
        }
    missing = set(PHASES) - set(rows)
    if missing:
        raise RuntimeError(f"reference table lacks phases: {sorted(missing)}")
    return rows


def qualified_energy(
    root: Path,
    mesh: int,
    phase: str,
    expected_digest: str,
) -> tuple[float, int]:
    mesh_name = f"k{mesh}{mesh}{mesh}-reduced"
    input_path = root / "inputs" / mesh_name / phase / "input.inp"
    run_dir = root / "runs" / mesh_name / phase
    value = energy(run_dir / "cp2k.out")
    if value is None:
        raise RuntimeError(f"incomplete CP2K output: mesh={mesh} phase={phase}")
    if binary_digest(run_dir) != expected_digest:
        raise RuntimeError(f"binary mismatch: mesh={mesh} phase={phase}")
    status = run_dir / "exit_status"
    if not status.is_file() or status.read_text(encoding="utf-8").strip() != "0":
        raise RuntimeError(f"nonzero or missing exit status: mesh={mesh} phase={phase}")
    input_digest = run_dir / "input.sha256"
    if (
        not input_path.is_file()
        or not input_digest.is_file()
        or recorded_digest(input_digest) != sha256(input_path)
    ):
        raise RuntimeError(f"input provenance mismatch: mesh={mesh} phase={phase}")
    parsed_mesh, water_count = input_mesh_and_water_count(input_path)
    if parsed_mesh != mesh:
        raise RuntimeError(
            f"directory/input mesh mismatch: directory={mesh} input={parsed_mesh} "
            f"phase={phase}"
        )
    if water_count <= 0:
        raise RuntimeError(f"invalid water count: mesh={mesh} phase={phase}")
    return value, water_count


def mesh_statistics(
    root: Path,
    mesh: int,
    reference: dict[str, float],
    expected_digest: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    ih_energy, ih_water = qualified_energy(root, mesh, "Ih", expected_digest)
    phase_rows: list[dict[str, object]] = []
    errors: list[float] = []
    for phase in PHASES:
        phase_energy, phase_water = qualified_energy(
            root, mesh, phase, expected_digest
        )
        relative = (
            phase_energy / phase_water - ih_energy / ih_water
        ) * HARTREE_TO_KJMOL
        signed_error = relative - reference[phase]
        errors.append(signed_error)
        phase_rows.append(
            {
                "mesh_n": mesh,
                "mesh_label": "Gamma" if mesh == 1 else f"{mesh}x{mesh}x{mesh}",
                "phase": phase,
                "relative_kj_mol_per_water": relative,
                "reference_kj_mol_per_water": reference[phase],
                "error_kj_mol_per_water": signed_error,
                "absolute_error_kj_mol_per_water": abs(signed_error),
            }
        )
    summary = {
        "mesh_n": mesh,
        "mesh_label": "Gamma" if mesh == 1 else f"{mesh}x{mesh}x{mesh}",
        "phase_count": len(errors),
        "me_kj_mol_per_water": sum(errors) / len(errors),
        "mae_kj_mol_per_water": sum(abs(value) for value in errors) / len(errors),
        "rmse_kj_mol_per_water": math.sqrt(
            sum(value * value for value in errors) / len(errors)
        ),
        "maxae_kj_mol_per_water": max(abs(value) for value in errors),
    }
    return summary, phase_rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--meshes", default="1,2,3,4,5")
    parser.add_argument("--require-binary-sha256", required=True)
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--output-phase-csv", type=Path)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    expected_digest = args.require_binary_sha256.lower()
    if not SHA256_RE.fullmatch(expected_digest):
        parser.error(
            "--require-binary-sha256 must be a 64-character hexadecimal digest"
        )
    meshes = tuple(int(value) for value in args.meshes.split(","))
    if not meshes or min(meshes) < 1 or len(set(meshes)) != len(meshes):
        parser.error("--meshes must contain unique positive integers")

    reference = references(args.reference_csv)
    summaries: list[dict[str, object]] = []
    phase_rows: list[dict[str, object]] = []
    for mesh in meshes:
        summary, rows = mesh_statistics(
            root=args.root,
            mesh=mesh,
            reference=reference,
            expected_digest=expected_digest,
        )
        summaries.append(summary)
        phase_rows.extend(rows)

    if args.output_csv:
        write_csv(args.output_csv, summaries)
    if args.output_phase_csv:
        write_csv(args.output_phase_csv, phase_rows)
    payload = {
        "schema_version": 1,
        "benchmark": "DMC-ICE13",
        "method": "g-xTB",
        "status": "complete_fixed_mesh_series",
        "binary_sha256": expected_digest,
        "meshes": summaries,
        "phase_rows": phase_rows,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(text, encoding="utf-8")
    print("mesh\tMAE_kJ_mol_H2O\tRMSE_kJ_mol_H2O\tMaxAE_kJ_mol_H2O")
    for row in summaries:
        print(
            f"{row['mesh_label']}\t{row['mae_kj_mol_per_water']:.12f}"
            f"\t{row['rmse_kj_mol_per_water']:.12f}"
            f"\t{row['maxae_kj_mol_per_water']:.12f}"
        )


if __name__ == "__main__":
    main()
