#!/usr/bin/env python3
"""Independently verify phase-wise adaptive DMC-ICE13 endpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path


HARTREE_TO_KJMOL = 2625.4996394799
PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recorded_digest(path: Path) -> str:
    fields = path.read_text(encoding="utf-8", errors="replace").split()
    digest = fields[0].lower() if fields else ""
    if not SHA256_RE.fullmatch(digest):
        raise ValueError(f"invalid SHA-256 record: {path}")
    return digest


def input_mesh_and_water_count(path: Path) -> tuple[int, int]:
    mesh_values: list[int] = []
    water_count = 0
    in_coordinates = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        code = re.split(r"[#!]", raw_line, maxsplit=1)[0].split()
        if len(code) >= 2 and tuple(token.upper() for token in code[:2]) == (
            "SCHEME",
            "MACDONALD",
        ):
            if len(code) not in (5, 8):
                raise ValueError(f"unsupported MACDONALD syntax in {path}: {raw_line}")
            try:
                dimensions = tuple(int(value) for value in code[2:5])
                shifts = tuple(float(value) for value in code[5:])
            except ValueError as exc:
                raise ValueError(f"invalid MACDONALD mesh in {path}: {raw_line}") from exc
            if len(set(dimensions)) != 1:
                raise ValueError(f"anisotropic mesh in {path}: {dimensions}")
            if any(not math.isfinite(value) for value in shifts):
                raise ValueError(f"non-finite MACDONALD shift in {path}: {shifts}")
            mesh_values.append(dimensions[0])
        line = raw_line.strip()
        upper = line.upper()
        if upper == "&COORD":
            in_coordinates = True
            continue
        if in_coordinates and upper.startswith("&END"):
            in_coordinates = False
            continue
        if in_coordinates and line and not line.startswith(("#", "!")):
            water_count += line.split()[0].upper() == "O"
    if len(mesh_values) != 1:
        raise ValueError(
            f"expected exactly one MACDONALD mesh in {path}, found {len(mesh_values)}"
        )
    if water_count <= 0:
        raise ValueError(f"no oxygen atoms in coordinate section: {path}")
    return mesh_values[0], water_count


def final_energy(path: Path) -> float:
    energies: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            energies.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not energies or not math.isfinite(energies[-1]):
        raise ValueError(f"incomplete or non-finite CP2K output: {path}")
    return energies[-1]


def qualified_energy(
    root: Path, mesh: int, phase: str, required_binary: str
) -> tuple[float, int]:
    mesh_name = f"k{mesh}{mesh}{mesh}-reduced"
    run = root / "runs" / mesh_name / phase
    input_path = root / "inputs" / mesh_name / phase / "input.inp"
    if (run / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise ValueError(f"nonzero exit status: {run}")
    if recorded_digest(run / "binary.sha256") != required_binary:
        raise ValueError(f"wrong execution binary: {run}")
    if recorded_digest(run / "input.sha256") != sha256(input_path):
        raise ValueError(f"input hash mismatch: {run}")
    actual_mesh, water_count = input_mesh_and_water_count(input_path)
    if actual_mesh != mesh:
        raise ValueError(
            f"directory/input mesh mismatch for {phase}: directory={mesh}, input={actual_mesh}"
        )
    return final_energy(run / "cp2k.out"), water_count


def relative_energy(
    root: Path, mesh: int, phase: str, required_binary: str
) -> float:
    phase_energy, phase_waters = qualified_energy(root, mesh, phase, required_binary)
    ih_energy, ih_waters = qualified_energy(root, mesh, "Ih", required_binary)
    return (
        phase_energy / phase_waters - ih_energy / ih_waters
    ) * HARTREE_TO_KJMOL


def load_references(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        reference_field = next(
            (
                candidate
                for candidate in (
                    "DMC_relative_kJmol",
                    "dmc_reference_kj_mol_per_water",
                    "dmc_reference_kJ_mol_per_water",
                )
                if candidate in fields
            ),
            None,
        )
        if "phase" not in fields or reference_field is None:
            raise ValueError("unsupported DMC reference table")
        references = {row["phase"]: float(row[reference_field]) for row in reader}
    if set(PHASES) - references.keys():
        raise ValueError("reference table does not cover all DMC-ICE13 phases")
    return references


def close(actual: float, recorded: object, tolerance: float, label: str) -> None:
    try:
        expected = float(recorded)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid recorded {label}: {recorded!r}") from exc
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(
            f"recorded {label} differs: recomputed={actual:.15g}, recorded={expected:.15g}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("endpoint_json", type=Path)
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--meshes", default="5,6,7,8")
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--require-binary-sha256", required=True)
    parser.add_argument("--numeric-tolerance", type=float, default=5.0e-10)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    required_binary = args.require_binary_sha256.lower()
    if not SHA256_RE.fullmatch(required_binary):
        parser.error("--require-binary-sha256 must contain 64 hexadecimal characters")
    if not math.isfinite(args.threshold) or args.threshold < 0.0:
        parser.error("--threshold must be finite and non-negative")
    if not math.isfinite(args.numeric_tolerance) or args.numeric_tolerance < 0.0:
        parser.error("--numeric-tolerance must be finite and non-negative")
    try:
        meshes = tuple(int(value) for value in args.meshes.split(","))
    except ValueError as exc:
        parser.error(f"invalid mesh sequence: {exc}")
    if len(meshes) < 2 or any(right != left + 1 for left, right in zip(meshes, meshes[1:])):
        parser.error("meshes must be a strictly adjacent ascending sequence")

    endpoint_payload = json.loads(args.endpoint_json.read_text(encoding="utf-8"))
    if endpoint_payload.get("complete") is not True:
        raise ValueError("endpoint file is not marked complete")
    if endpoint_payload.get("required_binary_sha256") != required_binary:
        raise ValueError("endpoint file records a different execution binary")
    close(args.threshold, endpoint_payload.get("threshold_kj_mol_per_water"), 0.0, "threshold")
    endpoint_rows = endpoint_payload.get("rows")
    if not isinstance(endpoint_rows, list):
        raise ValueError("endpoint file has no row list")
    rows_by_phase = {row.get("phase"): row for row in endpoint_rows}
    if (
        len(endpoint_rows) != len(PHASES)
        or len(rows_by_phase) != len(PHASES)
        or set(rows_by_phase) != set(PHASES)
    ):
        raise ValueError("endpoint rows are incomplete or duplicated")

    references = load_references(args.reference_csv)
    verified_rows: list[dict[str, object]] = []
    signed_errors: list[float] = []
    for phase in PHASES:
        row = rows_by_phase[phase]
        if row.get("status") != "converged":
            raise ValueError(f"phase {phase} is not marked converged")
        endpoint_mesh = int(row["endpoint_mesh"])
        previous_mesh = int(row["previous_mesh"])
        if endpoint_mesh != previous_mesh + 1 or endpoint_mesh not in meshes:
            raise ValueError(f"invalid adjacent endpoint for phase {phase}")

        values: dict[int, float] = {}
        first_passing_pair: tuple[int, int] | None = None
        first_change = math.nan
        for left, right in zip(meshes, meshes[1:]):
            if right > endpoint_mesh:
                break
            for mesh in (left, right):
                if mesh not in values:
                    values[mesh] = relative_energy(root=args.root, mesh=mesh, phase=phase, required_binary=required_binary)
            change = values[right] - values[left]
            if abs(change) <= args.threshold:
                first_passing_pair = (left, right)
                first_change = change
                break
        if first_passing_pair != (previous_mesh, endpoint_mesh):
            raise ValueError(
                f"phase {phase} does not record its first passing adjacent pair: "
                f"recomputed={first_passing_pair}, recorded={(previous_mesh, endpoint_mesh)}"
            )

        relative = values[endpoint_mesh]
        error = relative - references[phase]
        signed_errors.append(error)
        close(relative, row.get("relative_kj_mol_per_water"), args.numeric_tolerance, f"{phase} relative energy")
        close(references[phase], row.get("reference_kj_mol_per_water"), args.numeric_tolerance, f"{phase} reference")
        close(error, row.get("error_kj_mol_per_water"), args.numeric_tolerance, f"{phase} signed error")
        close(abs(error), row.get("absolute_error_kj_mol_per_water"), args.numeric_tolerance, f"{phase} absolute error")
        close(first_change, row.get("adjacent_change_kj_mol_per_water"), args.numeric_tolerance, f"{phase} adjacent change")
        verified_rows.append(
            {
                "absolute_error_kj_mol_per_water": abs(error),
                "adjacent_change_kj_mol_per_water": first_change,
                "endpoint_mesh": endpoint_mesh,
                "phase": phase,
                "previous_mesh": previous_mesh,
                "relative_kj_mol_per_water": relative,
            }
        )

    statistics = {
        "mae_kj_mol_per_water": sum(map(abs, signed_errors)) / len(signed_errors),
        "max_abs_error_kj_mol_per_water": max(map(abs, signed_errors)),
        "mean_error_kj_mol_per_water": sum(signed_errors) / len(signed_errors),
        "rmse_kj_mol_per_water": math.sqrt(
            sum(error * error for error in signed_errors) / len(signed_errors)
        ),
    }
    recorded_statistics = endpoint_payload.get("statistics")
    if not isinstance(recorded_statistics, dict):
        raise ValueError("endpoint file has no statistics")
    for key, value in statistics.items():
        close(value, recorded_statistics.get(key), args.numeric_tolerance, key)

    result = {
        "required_binary_sha256": required_binary,
        "rows": verified_rows,
        "statistics": statistics,
        "status": "PASS",
        "threshold_kj_mol_per_water": args.threshold,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2) from error
