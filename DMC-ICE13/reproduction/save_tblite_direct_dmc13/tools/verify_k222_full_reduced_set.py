#!/usr/bin/env python3
"""Qualify full-grid versus symmetry-reduced native 2x2x2 DMC-ICE13 energies."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from verify_k222_cli_native_requalification import (
    NATIVE_DISPERSION_RE,
    PHASES,
    SCHEME_RE,
    digest,
    qualify_affinity,
    read_component,
    read_native_energy,
    recorded_digest,
    require_status_zero,
)


SYMMETRY_RE = re.compile(r"^\s*SYMMETRY\s+([TF])\s*$", re.IGNORECASE | re.MULTILINE)
FULL_GRID_RE = re.compile(r"^\s*FULL_GRID\s+([TF])\s*$", re.IGNORECASE | re.MULTILINE)
PROJECT_RE = re.compile(r"^\s*PROJECT\s+.*$", re.IGNORECASE)


def verify_input(path: Path, expected_symmetry: str, expected_full_grid: str) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if SCHEME_RE.search(text) is None:
        raise AssertionError(f"noncanonical 2x2x2 MacDonald mesh: {path}")
    symmetry = SYMMETRY_RE.findall(text)
    full_grid = FULL_GRID_RE.findall(text)
    if symmetry != [expected_symmetry] or full_grid != [expected_full_grid]:
        raise AssertionError(
            f"wrong symmetry/full-grid flags: {path} "
            f"symmetry={symmetry} full_grid={full_grid}"
        )
    invariant_lines = []
    for line in text.splitlines():
        if SYMMETRY_RE.match(line) or FULL_GRID_RE.match(line) or PROJECT_RE.match(line):
            continue
        invariant_lines.append(line.rstrip())
    return digest(path), "\n".join(invariant_lines) + "\n"


def verify_run(
    run_dir: Path,
    input_path: Path,
    expected_binary: str,
    expected_symmetry: str,
    expected_full_grid: str,
) -> dict[str, object]:
    require_status_zero(run_dir / "exit_status", f"run status {run_dir}")
    if recorded_digest(run_dir / "binary.sha256") != expected_binary:
        raise AssertionError(f"binary mismatch: {run_dir}")
    input_hash, invariant = verify_input(
        input_path, expected_symmetry, expected_full_grid
    )
    if recorded_digest(run_dir / "input.sha256") != input_hash:
        raise AssertionError(f"input mismatch: {run_dir}")
    cpu = qualify_affinity(run_dir / "affinity_preexec.txt")
    output = run_dir / "cp2k.out"
    text = output.read_text(encoding="utf-8", errors="replace")
    return {
        "energy_Ha": read_native_energy(output),
        "dispersion_Ha": read_component(
            text,
            NATIVE_DISPERSION_RE,
            f"non-self-consistent dispersion energy: {run_dir}",
        ),
        "output_sha256": digest(output),
        "input_sha256": input_hash,
        "invariant_input": invariant,
        "cpu": cpu,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reduced_run_root", type=Path)
    parser.add_argument("full_run_root", type=Path)
    parser.add_argument("reduced_input_root", type=Path)
    parser.add_argument("full_input_root", type=Path)
    parser.add_argument("--expected-binary", required=True)
    parser.add_argument("--energy-tolerance-ha", type=float, default=5.0e-12)
    parser.add_argument("--dispersion-tolerance-ha", type=float, default=5.0e-12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if re.fullmatch(r"[0-9a-f]{64}", args.expected_binary) is None:
        parser.error("expected binary must be a SHA-256 digest")
    if args.energy_tolerance_ha <= 0.0 or args.dispersion_tolerance_ha <= 0.0:
        parser.error("tolerances must be positive")

    rows = []
    for phase in PHASES:
        reduced_input = args.reduced_input_root / phase / "input.inp"
        full_input = args.full_input_root / phase / "input.inp"
        reduced = verify_run(
            args.reduced_run_root / phase,
            reduced_input,
            args.expected_binary,
            "T",
            "F",
        )
        full = verify_run(
            args.full_run_root / phase,
            full_input,
            args.expected_binary,
            "F",
            "T",
        )
        if reduced["invariant_input"] != full["invariant_input"]:
            raise AssertionError(f"non-symmetry input difference: {phase}")
        energy_delta = float(full["energy_Ha"]) - float(reduced["energy_Ha"])
        dispersion_delta = float(full["dispersion_Ha"]) - float(
            reduced["dispersion_Ha"]
        )
        if abs(energy_delta) > args.energy_tolerance_ha:
            raise AssertionError(
                f"full/reduced energy mismatch {phase}: {energy_delta:+.6e} Ha"
            )
        if abs(dispersion_delta) > args.dispersion_tolerance_ha:
            raise AssertionError(
                f"full/reduced dispersion mismatch {phase}: "
                f"{dispersion_delta:+.6e} Ha"
            )
        rows.append(
            {
                "phase": phase,
                "reduced_energy_Ha": reduced["energy_Ha"],
                "full_energy_Ha": full["energy_Ha"],
                "full_minus_reduced_energy_Ha": energy_delta,
                "reduced_dispersion_Ha": reduced["dispersion_Ha"],
                "full_dispersion_Ha": full["dispersion_Ha"],
                "full_minus_reduced_dispersion_Ha": dispersion_delta,
                "reduced_input_sha256": reduced["input_sha256"],
                "full_input_sha256": full["input_sha256"],
                "reduced_output_sha256": reduced["output_sha256"],
                "full_output_sha256": full["output_sha256"],
                "reduced_cpu": reduced["cpu"],
                "full_cpu": full["cpu"],
            }
        )

    energy_deltas = [float(row["full_minus_reduced_energy_Ha"]) for row in rows]
    dispersion_deltas = [
        float(row["full_minus_reduced_dispersion_Ha"]) for row in rows
    ]
    payload = {
        "status": "PASS",
        "mesh": "2x2x2",
        "phase_count": len(rows),
        "rows": rows,
        "statistics": {
            "max_abs_full_minus_reduced_energy_Ha": max(
                abs(value) for value in energy_deltas
            ),
            "rms_full_minus_reduced_energy_Ha": math.sqrt(
                sum(value * value for value in energy_deltas) / len(energy_deltas)
            ),
            "max_abs_full_minus_reduced_dispersion_Ha": max(
                abs(value) for value in dispersion_deltas
            ),
            "rms_full_minus_reduced_dispersion_Ha": math.sqrt(
                sum(value * value for value in dispersion_deltas)
                / len(dispersion_deltas)
            ),
        },
        "provenance": {"cp2k_binary_sha256": args.expected_binary},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
