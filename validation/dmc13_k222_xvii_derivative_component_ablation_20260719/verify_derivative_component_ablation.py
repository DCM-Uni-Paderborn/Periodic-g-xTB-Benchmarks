#!/usr/bin/env python3
"""Verify full-grid versus SPGLIB-reduced energy, forces, and stress."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_BINARY = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
MODES = ("no_exchange", "no_acp", "no_exchange_no_acp")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse(mode: str, route: str) -> dict[str, object]:
    directory = ROOT / "results" / mode / route
    output_path = directory / "cp2k.out"
    output = output_path.read_text(errors="replace")
    energy = re.findall(
        r"^\s*ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$", output, re.MULTILINE
    )
    if not energy:
        raise AssertionError(f"missing energy in {output_path}")

    force_rows = re.findall(
        r"^\s*FORCES\|\s+\d+\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s*$",
        output,
        re.MULTILINE,
    )
    if not force_rows:
        raise AssertionError(f"missing forces in {output_path}")

    stress_rows = []
    for axis in ("x", "y", "z"):
        match = re.search(
            rf"^\s*STRESS\|\s+{axis}\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
            output,
            re.MULTILINE,
        )
        if match is None:
            raise AssertionError(f"missing stress row {axis} in {output_path}")
        stress_rows.append([float(value) for value in match.groups()])

    affinity = (directory / "affinity_preexec.txt").read_text()
    allowed = re.search(r"^Cpus_allowed_list:\s*(\d+)\s*$", affinity, re.MULTILINE)
    return {
        "energy_hartree": float(energy[-1]),
        "forces_hartree_per_bohr": [[float(value) for value in row] for row in force_rows],
        "stress_bar": stress_rows,
        "binary_sha256": (directory / "binary.sha256").read_text().split()[0],
        "input_sha256": sha256(ROOT / "inputs" / f"{mode}_{route}.inp"),
        "output_sha256": sha256(output_path),
        "exit_status": int((directory / "exit_status").read_text().strip()),
        "normal_termination": "PROGRAM ENDED AT" in output,
        "singleton_cpu": int(allowed.group(1)) if allowed else None,
    }


def maximum_difference(left: list, right: list) -> float:
    if len(left) != len(right):
        raise AssertionError("array shapes differ")
    maximum = 0.0
    for left_row, right_row in zip(left, right):
        if len(left_row) != len(right_row):
            raise AssertionError("array shapes differ")
        maximum = max(
            maximum,
            *(abs(float(a) - float(b)) for a, b in zip(left_row, right_row)),
        )
    return maximum


def main() -> int:
    controller_status = int((ROOT / "controller_exit_status").read_text().strip())
    checks = [{"name": "controller_exit_status", "passed": controller_status == 0}]
    rows = []

    for mode in MODES:
        full = parse(mode, "full")
        reduced = parse(mode, "reduced")
        energy_difference = float(full["energy_hartree"]) - float(reduced["energy_hartree"])
        force_difference = maximum_difference(
            full["forces_hartree_per_bohr"], reduced["forces_hartree_per_bohr"]
        )
        stress_difference = maximum_difference(full["stress_bar"], reduced["stress_bar"])
        row_passed = (
            abs(energy_difference) <= 1.0e-12
            and force_difference <= 1.0e-12
            and stress_difference <= 1.0e-8
        )
        rows.append(
            {
                "mode": mode,
                "full": full,
                "reduced": reduced,
                "full_minus_reduced_energy_hartree": energy_difference,
                "maximum_force_difference_hartree_per_bohr": force_difference,
                "maximum_stress_difference_bar": stress_difference,
                "passed": row_passed,
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
                {"name": f"{mode}: derivative equivalence", "passed": row_passed},
            ]
        )

    passed = all(bool(check["passed"]) for check in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "phase": "XVII",
        "mesh": "2x2x2",
        "checks": checks,
        "rows": rows,
        "interpretation": (
            "For all three component ablations, full-grid and SPGLIB-reduced "
            "energies, forces, and stresses are identical at printed precision."
        ),
    }
    (ROOT / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
