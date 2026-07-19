#!/usr/bin/env python3
"""Test one-step convergence of a DMC relative energy against ice Ih."""

from __future__ import annotations

import argparse
import hashlib
import math
import re
import sys
from pathlib import Path

from bvk_input import input_mesh_and_water_count


HARTREE_TO_KJMOL = 2625.4996394799
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
        raise ValueError(f"invalid SHA-256 provenance: {path}")
    return fields[0].lower()


def binary_digest(run_dir: Path) -> str | None:
    path = run_dir / "binary.sha256"
    if not path.is_file():
        return None
    return recorded_digest(path)


def final_energy(path: Path) -> float:
    values: list[float] = []
    ended = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = ENERGY_RE.match(line)
            if match:
                values.append(float(match.group(1)))
            if "PROGRAM ENDED AT" in line:
                ended = True
    if not ended or not values or not math.isfinite(values[-1]):
        raise ValueError(f"incomplete CP2K output: {path}")
    return values[-1]


def relative_energy(
    root: Path, mesh: int, phase: str, expected_digest: str | None = None
) -> float:
    phase_dir = root / "runs" / f"k{mesh}{mesh}{mesh}-reduced" / phase
    ih_dir = root / "runs" / f"k{mesh}{mesh}{mesh}-reduced" / "Ih"
    phase_input = root / "inputs" / f"k{mesh}{mesh}{mesh}-reduced" / phase / "input.inp"
    ih_input = root / "inputs" / f"k{mesh}{mesh}{mesh}-reduced" / "Ih" / "input.inp"
    phase_digest = binary_digest(phase_dir)
    ih_digest = binary_digest(ih_dir)
    if expected_digest is not None:
        if phase_digest != expected_digest or ih_digest != expected_digest:
            raise ValueError(
                "unqualified phase/reference binary: "
                f"required={expected_digest} phase={phase_digest or 'missing'} "
                f"Ih={ih_digest or 'missing'} mesh={mesh} phase_name={phase}"
            )
    elif phase_digest is not None or ih_digest is not None:
        if phase_digest is None or ih_digest is None or phase_digest != ih_digest:
            raise ValueError(
                "phase/reference binary mismatch: "
                f"phase={phase_digest or 'missing'} Ih={ih_digest or 'missing'} "
                f"mesh={mesh} phase_name={phase}"
            )
    if expected_digest is not None:
        for run_dir, input_path in (
            (phase_dir, phase_input),
            (ih_dir, ih_input),
        ):
            exit_status = run_dir / "exit_status"
            if not exit_status.is_file() or exit_status.read_text().strip() != "0":
                raise ValueError(f"missing or nonzero exit status: {exit_status}")
            input_hash_path = run_dir / "input.sha256"
            recorded_input = recorded_digest(input_hash_path)
            actual_input = sha256(input_path)
            if recorded_input != actual_input:
                raise ValueError(
                    f"input hash mismatch: recorded={recorded_input} "
                    f"actual={actual_input} path={input_path}"
                )
    phase_mesh, phase_water = input_mesh_and_water_count(phase_input)
    ih_mesh, ih_water = input_mesh_and_water_count(ih_input)
    if phase_mesh != mesh or ih_mesh != mesh:
        raise ValueError(
            f"directory/input mesh mismatch: directory={mesh} "
            f"phase={phase_mesh} Ih={ih_mesh} phase_name={phase}"
        )
    phase_per_water = final_energy(phase_dir / "cp2k.out") / phase_water
    ih_per_water = final_energy(ih_dir / "cp2k.out") / ih_water
    return (phase_per_water - ih_per_water) * HARTREE_TO_KJMOL


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("previous_mesh", type=int)
    parser.add_argument("current_mesh", type=int)
    parser.add_argument("phase")
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--require-binary-sha256")
    args = parser.parse_args()
    expected_digest = None
    if args.require_binary_sha256 is not None:
        expected_digest = args.require_binary_sha256.lower()
        if not SHA256_RE.fullmatch(expected_digest):
            parser.error("--require-binary-sha256 must be a 64-character hexadecimal digest")

    try:
        previous = relative_energy(
            args.root, args.previous_mesh, args.phase, expected_digest
        )
        current = relative_energy(
            args.root, args.current_mesh, args.phase, expected_digest
        )
    except (OSError, ValueError) as exc:
        print(f"status=incomplete phase={args.phase} error={exc}", file=sys.stderr)
        return 2

    delta = abs(current - previous)
    converged = delta <= args.threshold
    print(
        f"phase={args.phase}\tprevious_mesh={args.previous_mesh}"
        f"\tcurrent_mesh={args.current_mesh}\trel_previous_kj_mol={previous:.12f}"
        f"\trel_current_kj_mol={current:.12f}\tdelta_kj_mol={delta:.12f}"
        f"\tthreshold_kj_mol={args.threshold:.12f}"
        f"\tstatus={'converged' if converged else 'unresolved'}"
    )
    return 0 if converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
