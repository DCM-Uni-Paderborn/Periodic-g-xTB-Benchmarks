#!/usr/bin/env python3
"""Compare full/reduced native-k energies with a direct BvK CLI result."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def recorded_digest(path: Path, label: str) -> str:
    if not path.is_file():
        raise ValueError(f"missing {label} digest: {path}")
    fields = path.read_text(encoding="utf-8", errors="replace").split()
    digest = fields[0].lower() if fields else ""
    if not SHA256_RE.fullmatch(digest):
        raise ValueError(f"invalid {label} digest: {path}")
    return digest


def file_digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def cp2k_energy(
    path: Path,
    required_binary: str | None,
    input_path: Path | None,
    required_input: str | None,
) -> tuple[float, dict[str, str]]:
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
    provenance: dict[str, str] = {}
    if required_binary is not None:
        digest = recorded_digest(run_dir / "binary.sha256", "CP2K binary")
        if digest != required_binary:
            raise ValueError(
                f"wrong CP2K binary for {path}: actual={digest or 'missing'} "
                f"required={required_binary}"
            )
        provenance["binary_sha256"] = digest
    if required_input is not None:
        if input_path is None or not input_path.is_file():
            raise ValueError(f"missing CP2K input for {path}: {input_path}")
        actual = file_digest(input_path)
        recorded = recorded_digest(run_dir / "input.sha256", "CP2K input")
        if actual != recorded or actual != required_input:
            raise ValueError(
                f"wrong CP2K input for {path}: actual={actual} recorded={recorded} "
                f"required={required_input}"
            )
        provenance["input_sha256"] = actual
    return values[-1], provenance


def cli_energy(
    path: Path,
    required_binary: str | None,
    input_path: Path | None,
    required_input: str | None,
) -> tuple[float, dict[str, str]]:
    run_dir = path.parent
    exit_status = run_dir / "exit_status"
    if not exit_status.is_file() or exit_status.read_text().strip() != "0":
        raise ValueError(f"missing or nonzero CLI exit status: {exit_status}")
    value = float(json.loads(path.read_text(encoding="utf-8"))["energy"])
    if not math.isfinite(value):
        raise ValueError(f"non-finite CLI energy: {path}")
    provenance: dict[str, str] = {}
    if required_binary is not None:
        digest = recorded_digest(run_dir / "binary.sha256", "CLI binary")
        if digest != required_binary:
            raise ValueError(
                f"wrong CLI binary for {path}: actual={digest} required={required_binary}"
            )
        provenance["binary_sha256"] = digest
    if required_input is not None:
        if input_path is None or not input_path.is_file():
            raise ValueError(f"missing CLI input for {path}: {input_path}")
        actual = file_digest(input_path)
        if actual != required_input:
            raise ValueError(
                f"wrong CLI input for {path}: actual={actual} required={required_input}"
            )
        provenance["input_sha256"] = actual
    return value, provenance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reduced_cp2k_output", type=Path)
    parser.add_argument("full_cp2k_output", type=Path)
    parser.add_argument("direct_cli_result", type=Path)
    parser.add_argument("--replicas", type=int, required=True)
    parser.add_argument("--symmetry-tolerance", type=float, default=5.0e-12)
    parser.add_argument("--cli-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--require-binary-sha256")
    parser.add_argument("--require-cli-binary-sha256")
    parser.add_argument("--reduced-input", type=Path)
    parser.add_argument("--full-input", type=Path)
    parser.add_argument("--cli-input", type=Path)
    parser.add_argument("--require-reduced-input-sha256")
    parser.add_argument("--require-full-input-sha256")
    parser.add_argument("--require-cli-input-sha256")
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
    digests = {
        "CP2K binary": args.require_binary_sha256,
        "CLI binary": args.require_cli_binary_sha256,
        "reduced input": args.require_reduced_input_sha256,
        "full input": args.require_full_input_sha256,
        "CLI input": args.require_cli_input_sha256,
    }
    for label, digest in tuple(digests.items()):
        if digest is None:
            continue
        normalized = digest.lower()
        if not SHA256_RE.fullmatch(normalized):
            raise ValueError(f"required {label} digest must contain 64 hexadecimal characters")
        digests[label] = normalized

    reduced, reduced_provenance = cp2k_energy(
        args.reduced_cp2k_output,
        digests["CP2K binary"],
        args.reduced_input,
        digests["reduced input"],
    )
    full, full_provenance = cp2k_energy(
        args.full_cp2k_output,
        digests["CP2K binary"],
        args.full_input,
        digests["full input"],
    )
    direct_total, cli_provenance = cli_energy(
        args.direct_cli_result,
        digests["CLI binary"],
        args.cli_input,
        digests["CLI input"],
    )
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
        "provenance": {
            "reduced_cp2k": reduced_provenance,
            "full_cp2k": full_provenance,
            "direct_cli": cli_provenance,
        },
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
