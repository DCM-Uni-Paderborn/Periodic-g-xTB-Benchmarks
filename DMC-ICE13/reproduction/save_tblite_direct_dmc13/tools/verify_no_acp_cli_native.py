#!/usr/bin/env python3
"""Qualify direct-CLI/native-CP2K energy parity with ACP disabled."""

from __future__ import annotations

import argparse
import json
import math
import re
import tomllib
from pathlib import Path

from verify_absolute_energy_parity import (
    HARTREE_TO_KJ_MOL,
    digest,
    native_energy,
)


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
REVISION_RE = re.compile(r"^[0-9a-f]{40}$")


def recorded_digest(path: Path) -> str:
    fields = path.read_text(encoding="utf-8", errors="replace").split()
    if not fields or not SHA256_RE.fullmatch(fields[0].lower()):
        raise AssertionError(f"invalid SHA-256 record: {path}")
    return fields[0].lower()


def metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def qualify_affinity(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        for field in line.split():
            if "=" in field:
                key, value = field.split("=", 1)
                values[key] = value
    expected = values.get("expected_cpu")
    allowed = values.get("allowed")
    if expected is None or allowed is None or expected != allowed:
        raise AssertionError(
            f"invalid singleton affinity proof: expected={expected} "
            f"allowed={allowed} path={path}"
        )
    if not expected.isdigit():
        raise AssertionError(f"non-singleton CPU affinity: {expected}")
    return {"expected_cpu": expected, "allowed": allowed, "sha256": digest(path)}


def read_direct_energy(path: Path) -> float:
    payload = json.loads(path.read_text(encoding="utf-8"))
    energy = float(payload["energy"])
    if not math.isfinite(energy):
        raise AssertionError(f"non-finite direct energy: {path}")
    return energy


def qualify_no_acp_parameter(path: Path) -> dict[str, int]:
    parameter = tomllib.loads(path.read_text(encoding="utf-8"))
    elements = parameter.get("element")
    if not isinstance(elements, dict):
        raise AssertionError("No-ACP parameter has no element table")
    counts: dict[str, int] = {}
    for symbol in ("H", "O"):
        record = elements.get(symbol)
        acp = record.get("acp") if isinstance(record, dict) else None
        levels = acp.get("acp_levels") if isinstance(acp, dict) else None
        if not isinstance(levels, list) or not levels:
            raise AssertionError(f"No-ACP parameter lacks {symbol} ACP levels")
        values = [float(value) for value in levels]
        if any(value != 0.0 for value in values):
            raise AssertionError(f"No-ACP parameter has nonzero {symbol} ACP levels")
        counts[symbol] = len(values)
    return counts


def water_count(path: Path, replicas: int) -> int:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 7:
        raise AssertionError(f"incomplete POSCAR: {path}")
    atoms = sum(int(value) for value in lines[6].split())
    denominator = 3 * replicas
    if atoms % denominator:
        raise AssertionError(
            f"atom count {atoms} is incompatible with {replicas} replicated waters: {path}"
        )
    return atoms // denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("direct_root", type=Path)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--structure-root", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--controller-exit-status", type=Path, required=True)
    parser.add_argument("--source-identity", type=Path, required=True)
    parser.add_argument("--require-source-revision", required=True)
    parser.add_argument("--require-binary-sha256", required=True)
    parser.add_argument("--require-native-binary-sha256")
    parser.add_argument("--replicas", type=int, default=8)
    parser.add_argument("--absolute-tolerance-ha", type=float, default=2.0e-7)
    parser.add_argument(
        "--relative-tolerance-kj-mol-per-water", type=float, default=5.0e-5
    )
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    expected_binary = args.require_binary_sha256.lower()
    expected_revision = args.require_source_revision.lower()
    expected_native_binary = (
        args.require_native_binary_sha256.lower()
        if args.require_native_binary_sha256 is not None
        else None
    )
    for label, value in (
        ("direct binary", expected_binary),
        ("native binary", expected_native_binary),
    ):
        if value is not None and not SHA256_RE.fullmatch(value):
            parser.error(f"invalid {label} SHA-256 digest")
    if not REVISION_RE.fullmatch(expected_revision):
        parser.error("invalid source revision")
    if args.replicas <= 0:
        parser.error("--replicas must be positive")
    if args.absolute_tolerance_ha <= 0.0:
        parser.error("--absolute-tolerance-ha must be positive")
    if args.relative_tolerance_kj_mol_per_water <= 0.0:
        parser.error("--relative-tolerance-kj-mol-per-water must be positive")

    if args.controller_exit_status.read_text().strip() != "0":
        raise AssertionError("nonzero direct CLI controller exit status")
    source = metadata(args.source_identity)
    if source.get("commit") != expected_revision:
        raise AssertionError("source revision mismatch")
    parameter_hash = digest(args.parameter_file)
    zeroed_acp_levels = qualify_no_acp_parameter(args.parameter_file)

    rows: list[dict[str, object]] = []
    for phase in ("Ih", "XVII"):
        direct_dir = args.direct_root / "k222" / phase
        native_dir = args.native_root / phase
        structure = args.structure_root / phase / "POSCAR"
        direct_json = direct_dir / "result.json"
        native_output = native_dir / "cp2k.out"

        if (direct_dir / "exit_status").read_text().strip() != "0":
            raise AssertionError(f"nonzero direct CLI exit status: {phase}")
        if recorded_digest(direct_dir / "binary.sha256") != expected_binary:
            raise AssertionError(f"direct CLI binary hash mismatch: {phase}")
        if recorded_digest(direct_dir / "input.sha256") != digest(structure):
            raise AssertionError(f"direct CLI input hash mismatch: {phase}")
        if recorded_digest(direct_dir / "parameter.sha256") != parameter_hash:
            raise AssertionError(f"direct CLI parameter hash mismatch: {phase}")
        process_text = (direct_dir / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        if (
            "total energy" not in process_text
            or "JSON dump of results written" not in process_text
        ):
            raise AssertionError(f"incomplete direct CLI process output: {phase}")
        affinity = qualify_affinity(direct_dir / "affinity_preexec.txt")

        if (native_dir / "exit_status").read_text().strip() != "0":
            raise AssertionError(f"nonzero native CP2K exit status: {phase}")
        if expected_native_binary is not None:
            if recorded_digest(native_dir / "binary.sha256") != expected_native_binary:
                raise AssertionError(f"native CP2K binary hash mismatch: {phase}")
        native_text = native_output.read_text(encoding="utf-8", errors="replace")
        if "PROGRAM ENDED AT" not in native_text:
            raise AssertionError(f"incomplete native CP2K output: {phase}")

        direct_supercell = read_direct_energy(direct_json)
        direct_primitive = direct_supercell / args.replicas
        native = native_energy(native_text, native_output)
        delta = native - direct_primitive
        if abs(delta) > args.absolute_tolerance_ha:
            raise AssertionError(
                f"No-ACP native/direct mismatch {phase}: {delta:+.6e} Ha"
            )
        rows.append(
            {
                "phase": phase,
                "water_count_primitive": water_count(structure, args.replicas),
                "direct_supercell_energy_Ha": direct_supercell,
                "direct_primitive_energy_Ha": direct_primitive,
                "cp2k_native_energy_Ha": native,
                "native_minus_direct_Ha": delta,
                "direct_json_sha256": digest(direct_json),
                "native_output_sha256": digest(native_output),
                "structure_sha256": digest(structure),
                "affinity": affinity,
            }
        )

    by_phase = {str(row["phase"]): row for row in rows}
    ih = by_phase["Ih"]
    xvii = by_phase["XVII"]
    native_relative = HARTREE_TO_KJ_MOL * (
        float(xvii["cp2k_native_energy_Ha"]) / int(xvii["water_count_primitive"])
        - float(ih["cp2k_native_energy_Ha"]) / int(ih["water_count_primitive"])
    )
    direct_relative = HARTREE_TO_KJ_MOL * (
        float(xvii["direct_primitive_energy_Ha"]) / int(xvii["water_count_primitive"])
        - float(ih["direct_primitive_energy_Ha"]) / int(ih["water_count_primitive"])
    )
    relative_delta = native_relative - direct_relative
    if abs(relative_delta) > args.relative_tolerance_kj_mol_per_water:
        raise AssertionError(
            "No-ACP native/direct relative mismatch: "
            f"{relative_delta:+.6e} kJ mol-1 per water"
        )

    payload = {
        "status": "PASS",
        "model": "g-xTB without ACP",
        "mesh": "2x2x2",
        "replicas": args.replicas,
        "thresholds": {
            "absolute_Ha": args.absolute_tolerance_ha,
            "relative_kJ_mol_per_water": args.relative_tolerance_kj_mol_per_water,
        },
        "statistics": {
            "max_abs_native_minus_direct_Ha": max(
                abs(float(row["native_minus_direct_Ha"])) for row in rows
            ),
            "native_XVII_minus_Ih_kJ_mol_per_water": native_relative,
            "direct_XVII_minus_Ih_kJ_mol_per_water": direct_relative,
            "native_minus_direct_relative_kJ_mol_per_water": relative_delta,
        },
        "provenance": {
            "direct_binary_sha256": expected_binary,
            "native_binary_sha256": expected_native_binary,
            "source_revision": expected_revision,
            "source_identity_sha256": digest(args.source_identity),
            "parameter_sha256": parameter_hash,
            "zeroed_acp_level_counts": zeroed_acp_levels,
            "controller_exit_status_sha256": digest(args.controller_exit_status),
        },
        "rows": rows,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
