#!/usr/bin/env python3
"""Report a phase-wise mixed-grid DMC-ICE13 MAE from completed CP2K runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


HARTREE_TO_KJMOL = 2625.4996394799
PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
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
    if not fields or not SHA256_RE.fullmatch(fields[0].lower()):
        raise RuntimeError(f"invalid SHA-256 provenance: {path}")
    return fields[0].lower()


def binary_digest(run_dir: Path) -> str | None:
    path = run_dir / "binary.sha256"
    if not path.is_file():
        return None
    return recorded_digest(path)


def energy(path: Path) -> float | None:
    if not path.is_file():
        return None
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ENERGY_RE.match(line)
        if match:
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    return values[-1] if ended and values else None


def oxygen_count(path: Path) -> int:
    count = 0
    in_coordinates = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        upper = line.upper()
        if upper == "&COORD":
            in_coordinates = True
            continue
        if in_coordinates and upper.startswith("&END"):
            in_coordinates = False
            continue
        if in_coordinates and line and not line.startswith(("#", "!")):
            count += line.split()[0].upper() == "O"
    if count <= 0:
        raise ValueError(f"no oxygen atoms in {path}")
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--meshes", default="6,5,4,3,2,1")
    parser.add_argument("--paper-json", type=Path)
    parser.add_argument("--require-binary-sha256")
    args = parser.parse_args()
    expected_digest = None
    if args.require_binary_sha256 is not None:
        expected_digest = args.require_binary_sha256.lower()
        if not SHA256_RE.fullmatch(expected_digest):
            parser.error("--require-binary-sha256 must be a 64-character hexadecimal digest")
    meshes = tuple(int(value) for value in args.meshes.split(","))
    with args.reference_csv.open(newline="", encoding="utf-8") as handle:
        references = {
            row["phase"]: float(row["DMC_relative_kJmol"])
            for row in csv.DictReader(handle)
        }

    paper_results = None
    if args.paper_json is not None:
        paper_results = json.loads(args.paper_json.read_text(encoding="utf-8"))["results"]

    rows = []
    for phase in PHASES:
        selected = None
        for mesh in meshes:
            phase_dir = args.root / "runs" / f"k{mesh}{mesh}{mesh}-reduced" / phase
            ih_dir = args.root / "runs" / f"k{mesh}{mesh}{mesh}-reduced" / "Ih"
            phase_output = phase_dir / "cp2k.out"
            ih_output = ih_dir / "cp2k.out"
            phase_input = args.root / "inputs" / f"k{mesh}{mesh}{mesh}-reduced" / phase / "input.inp"
            ih_input = args.root / "inputs" / f"k{mesh}{mesh}{mesh}-reduced" / "Ih" / "input.inp"
            phase_energy = energy(phase_output)
            ih_energy = energy(ih_output)
            if phase_energy is None or ih_energy is None:
                continue
            phase_digest = binary_digest(phase_dir)
            ih_digest = binary_digest(ih_dir)
            if expected_digest is not None:
                if phase_digest != expected_digest or ih_digest != expected_digest:
                    # A denser result from an older build must not mask a
                    # qualified lower-mesh pair.  Search the remaining meshes
                    # and fail below only if no qualified same-mesh pair exists.
                    continue
                qualified = True
                for run_dir, input_path in (
                    (phase_dir, phase_input),
                    (ih_dir, ih_input),
                ):
                    exit_status = run_dir / "exit_status"
                    if (
                        not exit_status.is_file()
                        or exit_status.read_text().strip() != "0"
                    ):
                        qualified = False
                        break
                    input_hash_path = run_dir / "input.sha256"
                    if (
                        not input_hash_path.is_file()
                        or recorded_digest(input_hash_path) != sha256(input_path)
                    ):
                        qualified = False
                        break
                if not qualified:
                    continue
            elif phase_digest is not None or ih_digest is not None:
                if phase_digest is None or ih_digest is None or phase_digest != ih_digest:
                    raise RuntimeError(
                        "phase/reference binary mismatch: "
                        f"phase={phase_digest or 'missing'} Ih={ih_digest or 'missing'} "
                        f"mesh={mesh} phase_name={phase}"
                    )
            relative = (
                phase_energy / oxygen_count(phase_input)
                - ih_energy / oxygen_count(ih_input)
            ) * HARTREE_TO_KJMOL
            selected = (mesh, relative)
            break
        if selected is None:
            qualification = (
                f" for required binary {expected_digest}"
                if expected_digest is not None
                else ""
            )
            raise RuntimeError(
                f"no complete qualified same-mesh result for {phase}{qualification}"
            )
        mesh, relative = selected
        error = relative - references[phase]
        paper_error = None
        paper_mesh = None
        if paper_results is not None:
            for candidate in range(mesh, 0, -1):
                mesh_id = f"k{candidate}{candidate}{candidate}"
                method = paper_results.get(mesh_id, {}).get("GXTB", {})
                relative_values = method.get("relative_kjmol", {})
                if phase in relative_values and "Ih" in relative_values:
                    paper_error = float(relative_values[phase]) - references[phase]
                    paper_mesh = candidate
                    break
            if paper_error is None:
                raise RuntimeError(f"no paper comparator available for {phase} through mesh {mesh}")
        rows.append((phase, mesh, relative, references[phase], error, abs(error), paper_error, paper_mesh))

    print(
        "phase\tmesh\trelative_kj_mol\treference_kj_mol\terror_kj_mol"
        "\tabs_error_kj_mol\tpaper_error_kj_mol\tpaper_mesh"
    )
    for row in rows:
        paper_field = "" if row[6] is None else f"{row[6]:.12f}"
        paper_mesh_field = "" if row[7] is None else str(row[7])
        print(
            f"{row[0]}\t{row[1]}\t{row[2]:.12f}\t{row[3]:.6f}\t{row[4]:.12f}"
            f"\t{row[5]:.12f}\t{paper_field}\t{paper_mesh_field}"
        )
    current_mae = sum(row[5] for row in rows) / len(rows)
    print(f"mixed_mae_kj_mol\t{current_mae:.12f}")
    if paper_results is not None:
        paper_mae = sum(abs(row[6]) for row in rows) / len(rows)
        print(f"paper_comparator_mae_kj_mol\t{paper_mae:.12f}")
        print(f"paper_comparator_all_same_mesh\t{str(all(row[1] == row[7] for row in rows)).lower()}")
        print(f"mae_improvement_kj_mol\t{paper_mae-current_mae:.12f}")
        print(f"mae_improvement_percent\t{100.0*(paper_mae-current_mae)/paper_mae:.12f}")


if __name__ == "__main__":
    main()
