#!/usr/bin/env python3
"""Compare native k points with CP2K and CLI explicit BvK supercells."""

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


def cp2k_energy(path: Path) -> float:
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values:
        raise ValueError(f"incomplete CP2K output: {path}")
    return values[-1]


def cli_energy(path: Path) -> float:
    value = float(json.loads(path.read_text(encoding="utf-8"))["energy"])
    if not math.isfinite(value):
        raise ValueError(f"non-finite CLI energy: {path}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("native_cp2k_output", type=Path)
    parser.add_argument("gamma_supercell_cp2k_output", type=Path)
    parser.add_argument("direct_cli_result", type=Path)
    parser.add_argument("--replicas", type=int, default=8)
    parser.add_argument("--parity-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--alignment-margin", type=float, default=1.0e-10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.replicas <= 0:
        raise ValueError("replica count must be positive")

    native = cp2k_energy(args.native_cp2k_output)
    gamma_total = cp2k_energy(args.gamma_supercell_cp2k_output)
    cli_total = cli_energy(args.direct_cli_result)
    gamma = gamma_total / args.replicas
    cli = cli_total / args.replicas

    deltas = {
        "native_minus_cli": native - cli,
        "gamma_minus_cli": gamma - cli,
        "native_minus_gamma": native - gamma,
    }
    gamma_to_cli = abs(deltas["gamma_minus_cli"])
    gamma_to_native = abs(deltas["native_minus_gamma"])
    if abs(gamma_to_cli - gamma_to_native) <= args.alignment_margin:
        alignment = "ambiguous_within_margin"
    elif gamma_to_cli < gamma_to_native:
        alignment = "direct_cli"
    else:
        alignment = "native_k_points"

    maximum = max(map(abs, deltas.values()))
    result = {
        "replicas": args.replicas,
        "native_cp2k_hartree": native,
        "gamma_supercell_cp2k_total_hartree": gamma_total,
        "gamma_supercell_cp2k_per_primitive_hartree": gamma,
        "direct_cli_total_hartree": cli_total,
        "direct_cli_per_primitive_hartree": cli,
        "deltas_hartree_per_primitive": deltas,
        "gamma_oracle_alignment": alignment,
        "maximum_absolute_pairwise_delta_hartree": maximum,
        "parity_tolerance_hartree": args.parity_tolerance,
        "status": "PASS" if maximum <= args.parity_tolerance else "FAIL",
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
