#!/usr/bin/env python3
"""Select strictly qualified one-step DMC-ICE13 k-point endpoints."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

from bvk_input import input_mesh_and_water_count


HARTREE_TO_KJMOL = 2625.4996394799
PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def mesh_directory(root: Path, mesh: int, area: str) -> Path:
    return root / area / f"k{mesh}{mesh}{mesh}-reduced"


def final_energy(path: Path) -> float:
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values or not math.isfinite(values[-1]):
        raise ValueError(f"incomplete or non-finite CP2K output: {path}")
    return values[-1]


def qualified_energy(
    root: Path, mesh: int, phase: str, required_binary: str
) -> tuple[float, int, str]:
    run = mesh_directory(root, mesh, "runs") / phase
    input_path = mesh_directory(root, mesh, "inputs") / phase / "input.inp"
    exit_status = run / "exit_status"
    if not exit_status.is_file() or exit_status.read_text().strip() != "0":
        raise ValueError(f"missing or nonzero exit status: {exit_status}")
    digest_path = run / "binary.sha256"
    fields = digest_path.read_text(encoding="utf-8", errors="replace").split()
    digest = fields[0].lower() if fields else ""
    if digest != required_binary:
        raise ValueError(
            f"wrong binary at mesh={mesh} phase={phase}: "
            f"actual={digest or 'missing'} required={required_binary}"
        )
    actual_mesh, water_count = input_mesh_and_water_count(input_path)
    if actual_mesh != mesh:
        raise ValueError(
            f"directory/input mesh mismatch for {phase}: "
            f"directory={mesh}, input={actual_mesh}"
        )
    return final_energy(run / "cp2k.out"), water_count, digest


def relative_energy(root: Path, mesh: int, phase: str, required_binary: str) -> float:
    phase_energy, phase_water, _ = qualified_energy(root, mesh, phase, required_binary)
    ih_energy, ih_water, _ = qualified_energy(root, mesh, "Ih", required_binary)
    return (
        phase_energy / phase_water - ih_energy / ih_water
    ) * HARTREE_TO_KJMOL


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = (
        "phase",
        "status",
        "previous_mesh",
        "endpoint_mesh",
        "relative_kj_mol_per_water",
        "reference_kj_mol_per_water",
        "error_kj_mol_per_water",
        "absolute_error_kj_mol_per_water",
        "adjacent_change_kj_mol_per_water",
        "reason",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--meshes", default="1,2,3,4,5,6,7,8")
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--require-binary-sha256", required=True)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    required_binary = args.require_binary_sha256.lower()
    if not SHA256_RE.fullmatch(required_binary):
        parser.error("--require-binary-sha256 must contain 64 hexadecimal characters")
    if not math.isfinite(args.threshold) or args.threshold < 0.0:
        parser.error("--threshold must be finite and non-negative")
    try:
        meshes = tuple(int(value) for value in args.meshes.split(","))
    except ValueError as exc:
        parser.error(f"invalid mesh sequence: {exc}")
    if len(meshes) < 2 or any(mesh <= 0 for mesh in meshes):
        parser.error("at least two positive meshes are required")
    if any(right != left + 1 for left, right in zip(meshes, meshes[1:])):
        parser.error("meshes must be a strictly adjacent ascending sequence")

    with args.reference_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "phase" not in (reader.fieldnames or ()):
            raise ValueError("reference table has no phase column")
        reference_field = next(
            (
                field
                for field in (
                    "DMC_relative_kJmol",
                    "dmc_reference_kj_mol_per_water",
                    "dmc_reference_kJ_mol_per_water",
                )
                if field in (reader.fieldnames or ())
            ),
            None,
        )
        if reference_field is None:
            raise ValueError("reference table has no supported DMC-energy column")
        references = {row["phase"]: float(row[reference_field]) for row in reader}
    missing_references = sorted(set(PHASES) - references.keys())
    if missing_references:
        raise ValueError(f"missing DMC references: {', '.join(missing_references)}")

    rows: list[dict[str, object]] = []
    endpoint_errors: list[float] = []
    incomplete = 0
    unresolved = 0
    for phase in PHASES:
        values: dict[int, float] = {}
        selected: tuple[int, int, float, float] | None = None
        reason = ""
        for previous_mesh, current_mesh in zip(meshes, meshes[1:]):
            try:
                for mesh in (previous_mesh, current_mesh):
                    if mesh not in values:
                        values[mesh] = relative_energy(
                            args.root, mesh, phase, required_binary
                        )
            except (OSError, ValueError) as exc:
                reason = str(exc)
                break
            change = values[current_mesh] - values[previous_mesh]
            if abs(change) <= args.threshold:
                selected = (
                    previous_mesh,
                    current_mesh,
                    values[current_mesh],
                    change,
                )
                break

        if selected is None:
            if reason:
                status = "incomplete"
                incomplete += 1
            else:
                status = "unresolved"
                unresolved += 1
                reason = f"no adjacent change through mesh {meshes[-1]} passed"
            rows.append(
                {
                    "phase": phase,
                    "status": status,
                    "previous_mesh": "",
                    "endpoint_mesh": "",
                    "relative_kj_mol_per_water": "",
                    "reference_kj_mol_per_water": f"{references[phase]:.12f}",
                    "error_kj_mol_per_water": "",
                    "absolute_error_kj_mol_per_water": "",
                    "adjacent_change_kj_mol_per_water": "",
                    "reason": reason,
                }
            )
            continue

        previous_mesh, endpoint_mesh, relative, change = selected
        error = relative - references[phase]
        endpoint_errors.append(error)
        rows.append(
            {
                "phase": phase,
                "status": "converged",
                "previous_mesh": previous_mesh,
                "endpoint_mesh": endpoint_mesh,
                "relative_kj_mol_per_water": f"{relative:.12f}",
                "reference_kj_mol_per_water": f"{references[phase]:.12f}",
                "error_kj_mol_per_water": f"{error:+.12f}",
                "absolute_error_kj_mol_per_water": f"{abs(error):.12f}",
                "adjacent_change_kj_mol_per_water": f"{change:+.12f}",
                "reason": "",
            }
        )

    complete = incomplete == 0 and unresolved == 0 and len(endpoint_errors) == len(PHASES)
    statistics: dict[str, float] | None = None
    if complete:
        statistics = {
            "mae_kj_mol_per_water": sum(map(abs, endpoint_errors)) / len(endpoint_errors),
            "max_abs_error_kj_mol_per_water": max(map(abs, endpoint_errors)),
            "mean_error_kj_mol_per_water": sum(endpoint_errors) / len(endpoint_errors),
            "rmse_kj_mol_per_water": math.sqrt(
                sum(error * error for error in endpoint_errors) / len(endpoint_errors)
            ),
        }
    result = {
        "complete": complete,
        "convergence_rule": "first adjacent phase-local pair; denser endpoint retained",
        "incomplete_phase_count": incomplete,
        "meshes": list(meshes),
        "phase_count": len(PHASES),
        "required_binary_sha256": required_binary,
        "rows": rows,
        "statistics": statistics,
        "threshold_kj_mol_per_water": args.threshold,
        "unresolved_phase_count": unresolved,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.write_text(payload, encoding="utf-8")
    if args.output_csv:
        write_csv(args.output_csv, rows)
    print(payload, end="")
    if incomplete:
        return 2
    if unresolved:
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2) from error
