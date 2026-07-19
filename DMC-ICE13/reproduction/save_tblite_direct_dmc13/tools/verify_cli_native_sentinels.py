#!/usr/bin/env python3
"""Fail-closed comparison of direct save_tblite and CP2K-native sentinels."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


HARTREE_TO_KJ_MOL = 2625.4996394799
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
HEX_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def recorded_digest(path: Path) -> str:
    if not path.is_file():
        raise AssertionError(f"missing digest record: {path}")
    fields = path.read_text(encoding="utf-8").split()
    if not fields or not HEX_RE.fullmatch(fields[0].lower()):
        raise AssertionError(f"invalid digest record: {path}")
    return fields[0].lower()


def require_exit_zero(run: Path) -> None:
    status = run / "exit_status"
    if not status.is_file() or status.read_text(encoding="utf-8").strip() != "0":
        raise AssertionError(f"missing or nonzero exit status: {run}")


def require_singleton_affinity(run: Path) -> int:
    path = run / "affinity_preexec.txt"
    if not path.is_file():
        raise AssertionError(f"missing pre-exec affinity proof: {run}")
    text = path.read_text(encoding="utf-8", errors="replace")
    header = re.search(r"expected_cpu=(\d+)\s+allowed=(\d+)", text)
    allowed = re.search(r"^Cpus_allowed_list:\s*(\d+)\s*$", text, re.MULTILINE)
    if header is None or allowed is None:
        raise AssertionError(f"malformed pre-exec affinity proof: {path}")
    expected_cpu = int(header.group(1))
    if int(header.group(2)) != expected_cpu or int(allowed.group(1)) != expected_cpu:
        raise AssertionError(f"non-singleton or mismatched affinity proof: {path}")
    return expected_cpu


def cp2k_energy(path: Path) -> float:
    if not path.is_file():
        raise AssertionError(f"missing CP2K output: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    values = [
        float(match.group(1))
        for line in text.splitlines()
        if (match := ENERGY_RE.match(line))
    ]
    if "PROGRAM ENDED AT" not in text or not values:
        raise AssertionError(f"incomplete CP2K output: {path}")
    return values[-1]


def direct_energy(path: Path) -> float:
    if not path.is_file():
        raise AssertionError(f"missing direct JSON result: {path}")
    value = float(json.loads(path.read_text(encoding="utf-8"))["energy"])
    if not math.isfinite(value):
        raise AssertionError(f"non-finite direct energy: {path}")
    return value


def poscar_atom_count(path: Path) -> int:
    if not path.is_file():
        raise AssertionError(f"missing POSCAR: {path}")
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 7:
        raise AssertionError(f"truncated POSCAR: {path}")
    for index in (5, 6):
        fields = lines[index]
        if fields and all(field.isdigit() for field in fields):
            return sum(map(int, fields))
    raise AssertionError(f"cannot locate POSCAR atom counts: {path}")


def cp2k_coord_count(path: Path) -> int:
    if not path.is_file():
        raise AssertionError(f"missing CP2K input: {path}")
    inside = False
    count = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        upper = line.upper()
        if upper == "&COORD":
            inside = True
            continue
        if inside and upper.startswith("&END"):
            break
        if (
            inside
            and line
            and not line.startswith(("#", "!"))
            and upper not in {"SCALED", "SCALED T", "SCALED TRUE"}
            and not upper.startswith("UNIT ")
        ):
            count += 1
    if not inside or count == 0:
        raise AssertionError(f"cannot locate CP2K coordinates: {path}")
    return count


def verify_run(
    phase: str,
    mesh: int,
    direct_root: Path,
    native_root: Path,
    direct_input_root: Path,
    native_input_root: Path,
    direct_json_name: str,
    expected_direct_binary: str,
    expected_native_binary: str,
) -> dict[str, object]:
    direct_run = direct_root / phase
    native_run = native_root / phase
    direct_input = direct_input_root / phase / "POSCAR"
    native_input = native_input_root / phase / "input.inp"

    for run in (direct_run, native_run):
        require_exit_zero(run)
    direct_cpu = require_singleton_affinity(direct_run)
    native_cpu = require_singleton_affinity(native_run)

    if recorded_digest(direct_run / "binary.sha256") != expected_direct_binary:
        raise AssertionError(f"direct binary mismatch: {phase}")
    if recorded_digest(native_run / "binary.sha256") != expected_native_binary:
        raise AssertionError(f"native binary mismatch: {phase}")
    if recorded_digest(direct_run / "input.sha256") != sha256(direct_input):
        raise AssertionError(f"direct input digest mismatch: {phase}")
    if recorded_digest(native_run / "input.sha256") != sha256(native_input):
        raise AssertionError(f"native input digest mismatch: {phase}")

    direct_json = direct_run / direct_json_name
    direct_total = direct_energy(direct_json)
    native_total = cp2k_energy(native_run / "cp2k.out")
    replicas = mesh**3
    direct_atoms = poscar_atom_count(direct_input)
    native_atoms = cp2k_coord_count(native_input)
    if direct_atoms != native_atoms * replicas:
        raise AssertionError(
            f"BvK atom-count mismatch {phase}: direct={direct_atoms} "
            f"native={native_atoms} replicas={replicas}"
        )
    if native_atoms % 3:
        raise AssertionError(f"nonintegral primitive water count: {phase}")

    direct_per_primitive = direct_total / replicas
    signed_delta = native_total - direct_per_primitive
    return {
        "phase": phase,
        "mesh": mesh,
        "primitive_atoms": native_atoms,
        "primitive_waters": native_atoms // 3,
        "direct_supercell_atoms": direct_atoms,
        "direct_total_Ha": direct_total,
        "direct_per_primitive_Ha": direct_per_primitive,
        "native_per_primitive_Ha": native_total,
        "native_minus_direct_Ha": signed_delta,
        "direct_cpu": direct_cpu,
        "native_cpu": native_cpu,
        "direct_result_sha256": sha256(direct_json),
        "native_output_sha256": sha256(native_run / "cp2k.out"),
        "direct_input_sha256": sha256(direct_input),
        "native_input_sha256": sha256(native_input),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, required=True)
    parser.add_argument("--phase", action="append", required=True)
    parser.add_argument("--direct-root", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--direct-input-root", type=Path, required=True)
    parser.add_argument("--native-input-root", type=Path, required=True)
    parser.add_argument("--direct-json-name", default="tblite.json")
    parser.add_argument("--direct-binary-sha256", required=True)
    parser.add_argument("--native-binary-sha256", required=True)
    parser.add_argument("--absolute-tolerance-Ha", type=float, default=2.0e-7)
    parser.add_argument(
        "--relative-tolerance-kJ-mol-water", type=float, default=5.0e-5
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mesh < 1:
        raise SystemExit("--mesh must be positive")
    expected_direct = args.direct_binary_sha256.lower()
    expected_native = args.native_binary_sha256.lower()
    for label, value in (("direct", expected_direct), ("native", expected_native)):
        if not HEX_RE.fullmatch(value):
            raise SystemExit(f"invalid {label} binary SHA-256")
    phases = list(dict.fromkeys(args.phase))
    if "Ih" not in phases:
        raise SystemExit("the Ih reference phase is required")

    rows = [
        verify_run(
            phase,
            args.mesh,
            args.direct_root,
            args.native_root,
            args.direct_input_root,
            args.native_input_root,
            args.direct_json_name,
            expected_direct,
            expected_native,
        )
        for phase in phases
    ]
    ih = next(row for row in rows if row["phase"] == "Ih")
    for row in rows:
        relative_delta = HARTREE_TO_KJ_MOL * (
            float(row["native_minus_direct_Ha"]) / int(row["primitive_waters"])
            - float(ih["native_minus_direct_Ha"]) / int(ih["primitive_waters"])
        )
        row["relative_native_minus_direct_kJ_mol_per_water"] = relative_delta

    maximum_absolute = max(abs(float(row["native_minus_direct_Ha"])) for row in rows)
    maximum_relative = max(
        abs(float(row["relative_native_minus_direct_kJ_mol_per_water"]))
        for row in rows
    )
    if maximum_absolute > args.absolute_tolerance_Ha:
        raise AssertionError(
            f"absolute-energy identity failed: {maximum_absolute:.12e} Ha > "
            f"{args.absolute_tolerance_Ha:.12e} Ha"
        )
    if maximum_relative > args.relative_tolerance_kJ_mol_water:
        raise AssertionError(
            f"relative-energy identity failed: {maximum_relative:.12e} kJ mol-1 "
            f"per water > {args.relative_tolerance_kJ_mol_water:.12e}"
        )

    payload = {
        "status": "pass",
        "mesh": args.mesh,
        "phases": phases,
        "thresholds": {
            "absolute_Ha": args.absolute_tolerance_Ha,
            "relative_kJ_mol_per_water": args.relative_tolerance_kJ_mol_water,
        },
        "summary": {
            "max_abs_native_minus_direct_Ha": maximum_absolute,
            "max_abs_relative_native_minus_direct_kJ_mol_per_water": maximum_relative,
        },
        "provenance": {
            "direct_binary_sha256": expected_direct,
            "native_binary_sha256": expected_native,
        },
        "rows": rows,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
