#!/usr/bin/env python3
"""Qualify a fresh all-phase Linux CLI/native 2x2x2 energy comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


PHASES = (
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
)
REPLICAS = 8
HARTREE_TO_KJ_MOL = 2625.4996394799
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
SCHEME_RE = re.compile(
    r"^\s*SCHEME\s+MACDONALD\s+2\s+2\s+2\s+"
    r"0\.25(?:0*)\s+0\.25(?:0*)\s+0\.25(?:0*)\s*$",
    re.IGNORECASE | re.MULTILINE,
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def recorded_digest(path: Path) -> str:
    if not path.is_file():
        raise AssertionError(f"missing digest record: {path}")
    fields = path.read_text(encoding="utf-8").split()
    if not fields or not re.fullmatch(r"[0-9a-f]{64}", fields[0]):
        raise AssertionError(f"invalid digest record: {path}")
    return fields[0]


def require_status_zero(path: Path, label: str) -> None:
    if not path.is_file() or path.read_text(encoding="utf-8").strip() != "0":
        raise AssertionError(f"nonzero or missing {label}: {path}")


def qualify_affinity(path: Path, expected_cpu: int | None = None) -> int:
    if not path.is_file():
        raise AssertionError(f"missing affinity proof: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"^pid=\d+\s+expected_cpu=(\d+)\s+allowed=([^\s]+)$", text, re.MULTILINE)
    if match is None:
        raise AssertionError(f"malformed affinity proof: {path}")
    recorded_cpu = int(match.group(1))
    allowed = match.group(2)
    if allowed != str(recorded_cpu) or (
        expected_cpu is not None and recorded_cpu != expected_cpu
    ):
        raise AssertionError(
            f"non-singleton or wrong affinity: requested={expected_cpu} "
            f"recorded={recorded_cpu} allowed={allowed}"
        )
    return recorded_cpu


def read_cli_energy(path: Path) -> float:
    value = float(json.loads(path.read_text(encoding="utf-8"))["energy"])
    if not math.isfinite(value):
        raise AssertionError(f"non-finite direct energy: {path}")
    return value


def read_native_energy(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    values = [
        float(match.group(1))
        for line in text.splitlines()
        if (match := ENERGY_RE.match(line))
    ]
    if "PROGRAM ENDED AT" not in text or not values:
        raise AssertionError(f"incomplete native output: {path}")
    return values[-1]


def water_count(structure: Path) -> int:
    lines = structure.read_text(encoding="utf-8").splitlines()
    if len(lines) < 7:
        raise AssertionError(f"incomplete POSCAR: {structure}")
    supercell_atoms = sum(int(value) for value in lines[6].split())
    denominator = 3 * REPLICAS
    if supercell_atoms % denominator:
        raise AssertionError(
            f"atom count is incompatible with a 2x2x2 water supercell: {structure}"
        )
    return supercell_atoms // denominator


def verify_native_input(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    if SCHEME_RE.search(text) is None:
        raise AssertionError(f"noncanonical 2x2x2 MacDonald mesh: {path}")
    if re.search(r"^\s*SYMMETRY\s+T\s*$", text, re.IGNORECASE | re.MULTILINE) is None:
        raise AssertionError(f"native symmetry reduction is not enabled: {path}")
    if re.search(r"^\s*FULL_GRID\s+F\s*$", text, re.IGNORECASE | re.MULTILINE) is None:
        raise AssertionError(f"native full-grid materialization is not disabled: {path}")
    return digest(path)


def verify_source_identity(
    path: Path,
    expected_revision: str,
    expected_binary: str,
    expected_native_binary: str,
    expected_native_provider_archive: str,
    expected_native_cp2k_revision: str,
    expected_native_cmake_cache: str,
    expected_native_build_ninja: str,
) -> str:
    if not path.is_file():
        raise AssertionError(f"missing source identity: {path}")
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    if values.get("commit") != expected_revision:
        raise AssertionError("direct provider revision mismatch")
    if values.get("executable_sha256") != expected_binary:
        raise AssertionError("direct provider executable mismatch")
    if values.get("native_provider_commit") != expected_revision:
        raise AssertionError("native provider revision mismatch")
    if values.get("native_provider_archive_sha256") != expected_native_provider_archive:
        raise AssertionError("native provider archive mismatch")
    if values.get("native_cp2k_commit") != expected_native_cp2k_revision:
        raise AssertionError("native CP2K revision mismatch")
    if values.get("native_cp2k_binary_sha256") != expected_native_binary:
        raise AssertionError("native CP2K executable mismatch")
    if values.get("native_cmake_provider") != "SAVE":
        raise AssertionError("native CMake provider is not SAVE")
    if values.get("native_cmake_provider_revision") != expected_revision:
        raise AssertionError("native CMake provider revision mismatch")
    if values.get("native_cmake_cache_sha256") != expected_native_cmake_cache:
        raise AssertionError("native CMake cache mismatch")
    if values.get("native_build_ninja_sha256") != expected_native_build_ninja:
        raise AssertionError("native link plan mismatch")
    return digest(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive_root", type=Path)
    parser.add_argument("direct_root", type=Path)
    parser.add_argument("native_run_root", type=Path)
    parser.add_argument("native_input_root", type=Path)
    parser.add_argument("--direct-controller-status", type=Path, required=True)
    parser.add_argument("--source-identity", type=Path, required=True)
    parser.add_argument("--expected-source-revision", required=True)
    parser.add_argument("--expected-direct-binary", required=True)
    parser.add_argument("--expected-native-binary", required=True)
    parser.add_argument("--expected-native-provider-archive", required=True)
    parser.add_argument("--expected-native-cp2k-revision", required=True)
    parser.add_argument("--expected-native-cmake-cache", required=True)
    parser.add_argument("--expected-native-build-ninja", required=True)
    parser.add_argument("--expected-direct-cpu", type=int, required=True)
    parser.add_argument("--tolerance-ha", type=float, default=2.0e-7)
    parser.add_argument("--relative-tolerance-kj-mol", type=float, default=5.0e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    for label, value, length in (
        ("direct binary", args.expected_direct_binary, 64),
        ("native binary", args.expected_native_binary, 64),
        ("native provider archive", args.expected_native_provider_archive, 64),
        ("native CMake cache", args.expected_native_cmake_cache, 64),
        ("native build.ninja", args.expected_native_build_ninja, 64),
        ("source revision", args.expected_source_revision, 40),
        ("native CP2K revision", args.expected_native_cp2k_revision, 40),
    ):
        if re.fullmatch(rf"[0-9a-f]{{{length}}}", value) is None:
            parser.error(f"invalid {label} digest")
    if args.expected_direct_cpu < 0:
        parser.error("expected direct CPU must be nonnegative")
    if args.tolerance_ha <= 0.0 or args.relative_tolerance_kj_mol <= 0.0:
        parser.error("energy tolerances must be positive")

    require_status_zero(args.direct_controller_status, "direct controller status")
    source_hash = verify_source_identity(
        args.source_identity,
        args.expected_source_revision,
        args.expected_direct_binary,
        args.expected_native_binary,
        args.expected_native_provider_archive,
        args.expected_native_cp2k_revision,
        args.expected_native_cmake_cache,
        args.expected_native_build_ninja,
    )

    rows: list[dict[str, object]] = []
    for phase in PHASES:
        structure = args.archive_root / "structures" / "k222" / phase / "POSCAR"
        direct_dir = args.direct_root / "k222" / phase
        native_dir = args.native_run_root / phase
        native_input = args.native_input_root / phase / "input.inp"
        if not structure.is_file() or not native_input.is_file():
            raise AssertionError(f"missing frozen input for phase {phase}")

        require_status_zero(direct_dir / "exit_status", f"direct phase {phase}")
        require_status_zero(native_dir / "exit_status", f"native phase {phase}")
        if recorded_digest(direct_dir / "binary.sha256") != args.expected_direct_binary:
            raise AssertionError(f"direct binary mismatch: {phase}")
        if recorded_digest(direct_dir / "input.sha256") != digest(structure):
            raise AssertionError(f"direct input mismatch: {phase}")
        if recorded_digest(native_dir / "binary.sha256") != args.expected_native_binary:
            raise AssertionError(f"native binary mismatch: {phase}")
        native_input_hash = verify_native_input(native_input)
        if recorded_digest(native_dir / "input.sha256") != native_input_hash:
            raise AssertionError(f"native input mismatch: {phase}")
        direct_cpu = qualify_affinity(
            direct_dir / "affinity_preexec.txt", args.expected_direct_cpu
        )
        native_cpu = qualify_affinity(native_dir / "affinity_preexec.txt")

        direct_text = (direct_dir / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        if (
            "total energy" not in direct_text
            or "JSON dump of results written" not in direct_text
        ):
            raise AssertionError(f"incomplete direct output: {phase}")
        direct_json = direct_dir / "result.json"
        direct_total = read_cli_energy(direct_json)
        direct_primitive = direct_total / REPLICAS
        native_output = native_dir / "cp2k.out"
        native = read_native_energy(native_output)
        delta = native - direct_primitive
        if abs(delta) > args.tolerance_ha:
            raise AssertionError(
                f"native/direct mismatch {phase}: {delta:+.6e} Ha"
            )
        rows.append(
            {
                "phase": phase,
                "water_count_primitive": water_count(structure),
                "direct_supercell_energy_Ha": direct_total,
                "direct_primitive_energy_Ha": direct_primitive,
                "native_primitive_energy_Ha": native,
                "native_minus_direct_Ha": delta,
                "structure_sha256": digest(structure),
                "direct_json_sha256": digest(direct_json),
                "native_output_sha256": digest(native_output),
                "native_input_sha256": native_input_hash,
                "direct_cpu": direct_cpu,
                "native_cpu": native_cpu,
            }
        )

    by_phase = {str(row["phase"]): row for row in rows}
    ih = by_phase["Ih"]
    ih_delta_per_water = float(ih["native_minus_direct_Ha"]) / int(
        ih["water_count_primitive"]
    )
    relative_deltas = {
        phase: HARTREE_TO_KJ_MOL
        * (
            float(row["native_minus_direct_Ha"])
            / int(row["water_count_primitive"])
            - ih_delta_per_water
        )
        for phase, row in by_phase.items()
        if phase != "Ih"
    }
    maximum_relative = max(abs(value) for value in relative_deltas.values())
    if maximum_relative > args.relative_tolerance_kj_mol:
        raise AssertionError(
            "native/direct relative mismatch: "
            f"{maximum_relative:.6e} kJ mol-1 per water"
        )
    signed = [float(row["native_minus_direct_Ha"]) for row in rows]
    payload = {
        "status": "PASS",
        "mesh": "2x2x2",
        "replicas": REPLICAS,
        "phase_count": len(rows),
        "rows": rows,
        "relative_native_minus_direct_kJ_mol_per_water": relative_deltas,
        "statistics": {
            "max_abs_native_minus_direct_Ha": max(abs(value) for value in signed),
            "rms_native_minus_direct_Ha": math.sqrt(
                sum(value * value for value in signed) / len(signed)
            ),
            "max_abs_relative_native_minus_direct_kJ_mol_per_water": maximum_relative,
        },
        "provenance": {
            "direct_binary_sha256": args.expected_direct_binary,
            "native_binary_sha256": args.expected_native_binary,
            "native_provider_archive_sha256": args.expected_native_provider_archive,
            "native_cp2k_revision": args.expected_native_cp2k_revision,
            "native_cmake_cache_sha256": args.expected_native_cmake_cache,
            "native_build_ninja_sha256": args.expected_native_build_ninja,
            "direct_source_revision": args.expected_source_revision,
            "source_identity_sha256": source_hash,
            "direct_cpu": args.expected_direct_cpu,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
