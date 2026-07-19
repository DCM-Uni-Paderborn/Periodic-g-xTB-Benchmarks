#!/usr/bin/env python3
"""Qualify a fresh direct Gamma CLI repetition against both archived oracles."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from verify_absolute_energy_parity import (
    HARTREE_TO_KJ_MOL,
    PHASES,
    RELATIVE_TOLERANCE_KJ_MOL_PER_WATER,
    ROOT,
    TOLERANCE_HARTREE,
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


def read_energy(path: Path) -> float:
    payload = json.loads(path.read_text(encoding="utf-8"))
    energy = float(payload["energy"])
    if not math.isfinite(energy):
        raise AssertionError(f"non-finite direct energy: {path}")
    return energy


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
            f"invalid singleton affinity proof: expected={expected} allowed={allowed} path={path}"
        )
    if not expected.isdigit():
        raise AssertionError(f"non-singleton CPU affinity: {expected}")
    return {"expected_cpu": expected, "allowed": allowed, "sha256": digest(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("new_result_root", type=Path)
    parser.add_argument("--archive-root", type=Path, default=ROOT)
    parser.add_argument("--require-binary-sha256", required=True)
    parser.add_argument("--source-identity", type=Path, required=True)
    parser.add_argument("--require-source-revision", required=True)
    parser.add_argument("--affinity-proof", type=Path, required=True)
    parser.add_argument("--controller-exit-status", type=Path, required=True)
    parser.add_argument("--phases", default=",".join(PHASES))
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    expected_binary = args.require_binary_sha256.lower()
    if not SHA256_RE.fullmatch(expected_binary):
        parser.error("--require-binary-sha256 must be a 64-character hexadecimal digest")
    expected_revision = args.require_source_revision.lower()
    if not REVISION_RE.fullmatch(expected_revision):
        parser.error("--require-source-revision must be a 40-character hexadecimal revision")
    phases = tuple(value.strip() for value in args.phases.split(",") if value.strip())
    if "Ih" not in phases:
        parser.error("--phases must include Ih")
    if len(set(phases)) != len(phases):
        parser.error("--phases contains duplicates")

    if args.controller_exit_status.read_text().strip() != "0":
        raise AssertionError(
            f"nonzero controller exit status: {args.controller_exit_status}"
        )
    affinity = qualify_affinity(args.affinity_proof)
    source_identity = metadata(args.source_identity)
    if source_identity.get("commit") != expected_revision:
        raise AssertionError("source revision mismatch")
    run_metadata = metadata(args.new_result_root / "run_metadata.txt")
    if run_metadata.get("executable_sha256") != expected_binary:
        raise AssertionError("run metadata executable hash mismatch")
    if "1" not in run_metadata.get("meshes", "").split():
        raise AssertionError("Gamma mesh missing from run metadata")
    if run_metadata.get("accuracy") != "0.1" or run_metadata.get("iterations") != "300":
        raise AssertionError("unexpected direct CLI numerical settings")
    metadata_phases = tuple(run_metadata.get("phases", "").split())
    if any(phase not in metadata_phases for phase in phases):
        raise AssertionError("phase coverage missing from run metadata")

    rows: list[dict[str, object]] = []
    for phase in phases:
        current_dir = args.new_result_root / "k111" / phase
        current_json = current_dir / "result.json"
        archived_json = (
            args.archive_root
            / "results"
            / "current_save_tblite_cli"
            / "k111"
            / phase
            / "result.json"
        )
        native_output = (
            args.archive_root
            / "results"
            / "current_cp2k_native"
            / "k111"
            / phase
            / "cp2k.out"
        )
        structure = args.archive_root / "structures" / "k111" / phase / "POSCAR"

        if (current_dir / "exit_status").read_text().strip() != "0":
            raise AssertionError(f"nonzero direct CLI exit status: {phase}")
        if recorded_digest(current_dir / "binary.sha256") != expected_binary:
            raise AssertionError(f"direct CLI binary hash mismatch: {phase}")
        if recorded_digest(current_dir / "input.sha256") != digest(structure):
            raise AssertionError(f"direct CLI input hash mismatch: {phase}")
        process_text = (current_dir / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        if "total energy" not in process_text or "JSON dump of results written" not in process_text:
            raise AssertionError(f"incomplete direct CLI process output: {phase}")

        current = read_energy(current_json)
        archived = read_energy(archived_json)
        native_text = native_output.read_text(encoding="utf-8", errors="replace")
        native = native_energy(native_text, native_output)
        current_minus_archived = current - archived
        native_minus_current = native - current
        if abs(current_minus_archived) > TOLERANCE_HARTREE:
            raise AssertionError(
                f"fresh/archived Gamma mismatch {phase}: {current_minus_archived:+.6e} Ha"
            )
        if abs(native_minus_current) > TOLERANCE_HARTREE:
            raise AssertionError(
                f"native/fresh Gamma mismatch {phase}: {native_minus_current:+.6e} Ha"
            )
        atom_count = sum(
            int(value)
            for value in structure.read_text(encoding="utf-8").splitlines()[6].split()
        )
        if atom_count % 3:
            raise AssertionError(f"nonintegral water count: {phase}")
        rows.append(
            {
                "phase": phase,
                "water_count": atom_count // 3,
                "fresh_cli_energy_Ha": current,
                "archived_cli_energy_Ha": archived,
                "cp2k_native_energy_Ha": native,
                "fresh_minus_archived_Ha": current_minus_archived,
                "native_minus_fresh_Ha": native_minus_current,
                "fresh_json_sha256": digest(current_json),
                "process_output_sha256": digest(current_dir / "process.out"),
                "input_sha256": digest(structure),
            }
        )

    ih = next(row for row in rows if row["phase"] == "Ih")
    maximum_relative = 0.0
    for row in rows:
        relative_delta = HARTREE_TO_KJ_MOL * (
            float(row["native_minus_fresh_Ha"]) / int(row["water_count"])
            - float(ih["native_minus_fresh_Ha"]) / int(ih["water_count"])
        )
        row["native_minus_fresh_relative_kJ_mol_per_water"] = relative_delta
        maximum_relative = max(maximum_relative, abs(relative_delta))
    if maximum_relative > RELATIVE_TOLERANCE_KJ_MOL_PER_WATER:
        raise AssertionError(
            "native/fresh relative Gamma mismatch: "
            f"{maximum_relative:.6e} > {RELATIVE_TOLERANCE_KJ_MOL_PER_WATER:.6e} "
            "kJ mol-1 per water"
        )

    payload = {
        "status": "PASS",
        "coverage": len(rows),
        "thresholds": {
            "absolute_Ha": TOLERANCE_HARTREE,
            "relative_kJ_mol_per_water": RELATIVE_TOLERANCE_KJ_MOL_PER_WATER,
        },
        "statistics": {
            "max_abs_fresh_minus_archived_Ha": max(
                abs(float(row["fresh_minus_archived_Ha"])) for row in rows
            ),
            "max_abs_native_minus_fresh_Ha": max(
                abs(float(row["native_minus_fresh_Ha"])) for row in rows
            ),
            "max_abs_native_minus_fresh_relative_kJ_mol_per_water": maximum_relative,
        },
        "provenance": {
            "binary_sha256": expected_binary,
            "source_revision": expected_revision,
            "source_identity_sha256": digest(args.source_identity),
            "run_metadata_sha256": digest(args.new_result_root / "run_metadata.txt"),
            "controller_exit_status_sha256": digest(args.controller_exit_status),
            "affinity": affinity,
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
