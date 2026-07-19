#!/usr/bin/env python3
"""Independently verify the phase-VIII full/reduced component ablation."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
HARTREE_TO_KJMOL = 2625.4996394799
EXPECTED_BINARY = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
CELL_TOLERANCE_HARTREE = 5.0e-10
PER_WATER_TOLERANCE_KJMOL = 1.0e-6
MODES = ("no_exchange", "no_acp", "no_exchange_no_acp")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_run(mode: str, route: str) -> dict[str, object]:
    run = ROOT / "results" / mode / route
    input_path = ROOT / "inputs" / f"{mode}_{route}.inp"
    output_path = run / "cp2k.out"
    output = output_path.read_text(errors="replace")
    input_text = input_path.read_text()

    energy_match = re.search(
        r"^\s*ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$", output, re.MULTILINE
    )
    if energy_match is None:
        raise AssertionError(f"missing total energy in {output_path}")

    oxygen_count = sum(
        1
        for line in input_text.splitlines()
        if re.match(r"^\s*O(?:\s|$)", line)
    )
    binary_hash = (run / "binary.sha256").read_text().split()[0]
    exit_status = int((run / "exit_status").read_text().strip())
    affinity = (run / "affinity_preexec.txt").read_text()
    allowed_match = re.search(r"^Cpus_allowed_list:\s*(\d+)\s*$", affinity, re.MULTILINE)

    return {
        "energy_hartree": float(energy_match.group(1)),
        "water_molecules": oxygen_count,
        "binary_sha256": binary_hash,
        "input_sha256": sha256(input_path),
        "output_sha256": sha256(output_path),
        "exit_status": exit_status,
        "normal_termination": "PROGRAM ENDED AT" in output,
        "singleton_cpu": int(allowed_match.group(1)) if allowed_match else None,
    }


def main() -> int:
    checks: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []

    controller_status = int((ROOT / "controller_exit_status").read_text().strip())
    checks.append({"name": "controller_exit_status", "passed": controller_status == 0})

    for mode in MODES:
        full = parse_run(mode, "full")
        reduced = parse_run(mode, "reduced")
        waters = int(full["water_molecules"])
        if waters <= 0 or waters != reduced["water_molecules"]:
            raise AssertionError(f"invalid molecule count for {mode}")

        delta = float(full["energy_hartree"]) - float(reduced["energy_hartree"])
        delta_per_water = delta * HARTREE_TO_KJMOL / waters
        row_passed = (
            abs(delta) <= CELL_TOLERANCE_HARTREE
            and abs(delta_per_water) <= PER_WATER_TOLERANCE_KJMOL
        )
        rows.append(
            {
                "mode": mode,
                "full": full,
                "reduced": reduced,
                "full_minus_reduced_hartree_per_cell": delta,
                "full_minus_reduced_kjmol_per_water": delta_per_water,
                "within_numerical_equivalence_tolerance": row_passed,
            }
        )
        checks.extend(
            [
                {
                    "name": f"{mode}: successful normal termination",
                    "passed": full["exit_status"] == reduced["exit_status"] == 0
                    and full["normal_termination"]
                    and reduced["normal_termination"],
                },
                {
                    "name": f"{mode}: immutable binary",
                    "passed": full["binary_sha256"]
                    == reduced["binary_sha256"]
                    == EXPECTED_BINARY,
                },
                {
                    "name": f"{mode}: disjoint singleton affinity",
                    "passed": full["singleton_cpu"] is not None
                    and reduced["singleton_cpu"] is not None
                    and full["singleton_cpu"] != reduced["singleton_cpu"],
                },
                {"name": f"{mode}: numerical energy equivalence", "passed": row_passed},
            ]
        )

    magnitudes = [abs(float(row["full_minus_reduced_hartree_per_cell"])) for row in rows]
    all_passed = all(bool(check["passed"]) for check in checks)
    payload = {
        "status": "PASS" if all_passed else "FAIL",
        "phase": "VIII",
        "mesh": "2x2x2",
        "eps_scf": 1.0e-12,
        "tolerances": {
            "hartree_per_cell": CELL_TOLERANCE_HARTREE,
            "kjmol_per_water": PER_WATER_TOLERANCE_KJMOL,
        },
        "checks": checks,
        "rows": rows,
        "maximum_absolute_difference_hartree_per_cell": max(magnitudes),
        "minimum_absolute_difference_hartree_per_cell": min(magnitudes),
        "interpretation": (
            "The sub-nanohartree full/reduced residual remains when exchange, ACP, "
            "or both are disabled. It is therefore not uniquely caused by either "
            "term or their coupling and is numerically negligible."
        ),
    }
    target = ROOT / "verification.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
