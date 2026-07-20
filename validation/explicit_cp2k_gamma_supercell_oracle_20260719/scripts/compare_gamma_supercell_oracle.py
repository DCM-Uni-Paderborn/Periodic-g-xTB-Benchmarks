#!/usr/bin/env python3
"""Compare native k points with CP2K and CLI explicit BvK supercells."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_native_symmetry_cli import SHA256_RE, cp2k_energy, cli_energy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("native_cp2k_output", type=Path)
    parser.add_argument("gamma_supercell_cp2k_output", type=Path)
    parser.add_argument("direct_cli_result", type=Path)
    parser.add_argument("--replicas", type=int, default=8)
    parser.add_argument("--parity-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--alignment-margin", type=float, default=1.0e-10)
    parser.add_argument("--require-binary-sha256")
    parser.add_argument("--require-cli-binary-sha256")
    parser.add_argument("--native-input", type=Path)
    parser.add_argument("--gamma-input", type=Path)
    parser.add_argument("--cli-input", type=Path)
    parser.add_argument("--require-native-input-sha256")
    parser.add_argument("--require-gamma-input-sha256")
    parser.add_argument("--require-cli-input-sha256")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.replicas <= 0:
        raise ValueError("replica count must be positive")

    digests = {
        "CP2K binary": args.require_binary_sha256,
        "CLI binary": args.require_cli_binary_sha256,
        "native input": args.require_native_input_sha256,
        "Gamma input": args.require_gamma_input_sha256,
        "CLI input": args.require_cli_input_sha256,
    }
    for label, digest in tuple(digests.items()):
        if digest is None:
            continue
        normalized = digest.lower()
        if not SHA256_RE.fullmatch(normalized):
            raise ValueError(
                f"required {label} digest must contain 64 hexadecimal characters"
            )
        digests[label] = normalized

    native, native_provenance = cp2k_energy(
        args.native_cp2k_output,
        digests["CP2K binary"],
        args.native_input,
        digests["native input"],
    )
    gamma_total, gamma_provenance = cp2k_energy(
        args.gamma_supercell_cp2k_output,
        digests["CP2K binary"],
        args.gamma_input,
        digests["Gamma input"],
    )
    cli_total, cli_provenance = cli_energy(
        args.direct_cli_result,
        digests["CLI binary"],
        args.cli_input,
        digests["CLI input"],
    )
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
        "provenance": {
            "native_cp2k": native_provenance,
            "gamma_supercell_cp2k": gamma_provenance,
            "direct_cli": cli_provenance,
        },
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
