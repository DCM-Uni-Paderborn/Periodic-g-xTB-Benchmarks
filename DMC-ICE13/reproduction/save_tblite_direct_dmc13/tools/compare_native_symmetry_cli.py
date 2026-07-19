#!/usr/bin/env python3
"""Compare full/reduced native-k energies with a direct BvK CLI result."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def cp2k_energy(path: Path, required_binary: str | None) -> float:
    run_dir = path.parent
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values or not math.isfinite(values[-1]):
        raise ValueError(f"incomplete or non-finite CP2K output: {path}")
    exit_status = run_dir / "exit_status"
    if not exit_status.is_file() or exit_status.read_text().strip() != "0":
        raise ValueError(f"missing or nonzero exit status: {exit_status}")
    if required_binary is not None:
        digest_path = run_dir / "binary.sha256"
        fields = digest_path.read_text(encoding="utf-8", errors="replace").split()
        digest = fields[0].lower() if fields else ""
        if digest != required_binary:
            raise ValueError(
                f"wrong CP2K binary for {path}: actual={digest or 'missing'} "
                f"required={required_binary}"
            )
    return values[-1]


def cli_energy(path: Path) -> float:
    value = float(json.loads(path.read_text(encoding="utf-8"))["energy"])
    if not math.isfinite(value):
        raise ValueError(f"non-finite CLI energy: {path}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reduced_cp2k_output", type=Path)
    parser.add_argument("full_cp2k_output", type=Path)
    parser.add_argument("direct_cli_result", type=Path)
    parser.add_argument("--replicas", type=int, required=True)
    parser.add_argument("--symmetry-tolerance", type=float, default=5.0e-12)
    parser.add_argument("--cli-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--require-binary-sha256")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.replicas <= 0:
        raise ValueError("replica count must be positive")
    for name, value in (
        ("symmetry tolerance", args.symmetry_tolerance),
        ("CLI tolerance", args.cli_tolerance),
    ):
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    required_binary = args.require_binary_sha256
    if required_binary is not None:
        required_binary = required_binary.lower()
        if not SHA256_RE.fullmatch(required_binary):
            raise ValueError("required binary digest must contain 64 hexadecimal characters")

    reduced = cp2k_energy(args.reduced_cp2k_output, required_binary)
    full = cp2k_energy(args.full_cp2k_output, required_binary)
    direct_total = cli_energy(args.direct_cli_result)
    direct = direct_total / args.replicas
    deltas = {
        "full_minus_reduced": full - reduced,
        "reduced_minus_cli": reduced - direct,
        "full_minus_cli": full - direct,
    }
    symmetry_pass = abs(deltas["full_minus_reduced"]) <= args.symmetry_tolerance
    cli_pass = max(
        abs(deltas["reduced_minus_cli"]), abs(deltas["full_minus_cli"])
    ) <= args.cli_tolerance
    result = {
        "cli_parity_pass": cli_pass,
        "cli_tolerance_hartree_per_primitive": args.cli_tolerance,
        "deltas_hartree_per_primitive": deltas,
        "direct_cli_hartree_per_primitive": direct,
        "direct_cli_total_hartree": direct_total,
        "full_cp2k_hartree_per_primitive": full,
        "reduced_cp2k_hartree_per_primitive": reduced,
        "replicas": args.replicas,
        "status": "PASS" if symmetry_pass and cli_pass else "FAIL",
        "symmetry_parity_pass": symmetry_pass,
        "symmetry_tolerance_hartree": args.symmetry_tolerance,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
