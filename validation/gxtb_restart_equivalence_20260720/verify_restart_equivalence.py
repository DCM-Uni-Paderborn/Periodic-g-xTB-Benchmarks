#!/usr/bin/env python3
"""Verify exact same-mesh and validated cross-mesh g-xTB restarts."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
QUALIFIED_BINARY_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
ENERGY_TOLERANCE_HARTREE = 1.0e-12

CASES = {
    "cold3": {
        "input": "CH4_gxtb_kp_restart_3.inp",
        "mesh": 3,
        "steps": 12,
        "energy": -40.468866070692428,
        "message": None,
    },
    "same3": {
        "input": "CH4_gxtb_kp_restart_3_same.inp",
        "mesh": 3,
        "steps": 1,
        "energy": -40.468866070692428,
        "message": "KPOINT_RESTART| Strict same-mesh restart accepted",
    },
    "transfer4": {
        "input": "CH4_gxtb_kp_restart_4.inp",
        "mesh": 4,
        "steps": 7,
        "energy": -40.468551982577388,
        "message": "KPOINT_RESTART| Validated BvK mesh transfer accepted",
    },
    "cold4": {
        "input": "CH4_gxtb_kp_cold_4.inp",
        "mesh": 4,
        "steps": 11,
        "energy": -40.468551982577395,
        "message": None,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_recorded_hash(path: Path) -> str:
    return path.read_text(encoding="utf-8").split()[0]


def verify_raw_manifest() -> tuple[bool, dict[str, str]]:
    recorded: dict[str, str] = {}
    for line in (HERE / "RAW_SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        relative = relative.lstrip("* ")
        recorded[relative] = digest
    expected_files = sorted(
        path.relative_to(HERE).as_posix()
        for path in RAW.rglob("*")
        if path.is_file()
    )
    return (
        sorted(recorded) == expected_files
        and all(sha256(HERE / relative) == digest for relative, digest in recorded.items()),
        recorded,
    )


def main() -> None:
    raw_manifest_passes, raw_manifest = verify_raw_manifest()
    results: dict[str, dict] = {}

    for name, expected in CASES.items():
        case_dir = RAW / name
        output_path = case_dir / "cp2k.out"
        input_path = RAW / "inputs" / expected["input"]
        output = output_path.read_text(encoding="utf-8")
        input_text = input_path.read_text(encoding="utf-8")

        energy_matches = re.findall(
            r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+([-+0-9.Ee]+)",
            output,
        )
        step_matches = re.findall(r"SCF run converged in\s+([0-9]+) steps", output)
        if len(energy_matches) != 1 or len(step_matches) != 1:
            raise RuntimeError(f"Could not uniquely parse {name}")
        energy = float(energy_matches[0])
        steps = int(step_matches[0])

        mesh = expected["mesh"]
        mesh_declared = bool(
            re.search(
                rf"^\s*SCHEME\s+MONKHORST-PACK\s+{mesh}\s+{mesh}\s+{mesh}\s*$",
                input_text,
                re.MULTILINE,
            )
        )
        affinity = (case_dir / "affinity_preexec.txt").read_text(encoding="utf-8")
        affinity_match = re.search(r"expected_cpu=([0-9]+) allowed=([0-9]+)", affinity)
        singleton_affinity = bool(
            affinity_match
            and affinity_match.group(1) == affinity_match.group(2)
            and re.search(
                rf"^Cpus_allowed_list:\s*{affinity_match.group(1)}\s*$",
                affinity,
                re.MULTILINE,
            )
        )
        message = expected["message"]
        case_checks = {
            "exit_status_zero": (case_dir / "exit_status").read_text().strip() == "0",
            "program_ended": "PROGRAM ENDED AT" in output,
            "scf_converged": "SCF run converged" in output,
            "energy_matches_archived_expectation": math.isclose(
                energy, expected["energy"], rel_tol=0.0, abs_tol=5.0e-15
            ),
            "step_count_matches_archived_expectation": steps == expected["steps"],
            "mesh_declared": mesh_declared,
            "binary_is_qualified": parse_recorded_hash(case_dir / "binary.sha256")
            == QUALIFIED_BINARY_SHA256,
            "input_hash_matches": parse_recorded_hash(case_dir / "input.sha256")
            == sha256(input_path),
            "singleton_cpu_affinity_proven": singleton_affinity,
            "expected_restart_message_present": message is None or message in output,
        }
        results[name] = {
            "status": "PASS" if all(case_checks.values()) else "FAIL",
            "energy_hartree": energy,
            "scf_steps": steps,
            "mesh": [mesh, mesh, mesh],
            "checks": case_checks,
        }

    restart_path = RAW / "work" / "CH4_gxtb_kp_restart_3-RESTART.kp"
    restart_hash_passes = (
        parse_recorded_hash(RAW / "restart_file.sha256") == sha256(restart_path)
    )
    same_mesh_delta = abs(
        results["same3"]["energy_hartree"] - results["cold3"]["energy_hartree"]
    )
    cross_mesh_delta = abs(
        results["transfer4"]["energy_hartree"]
        - results["cold4"]["energy_hartree"]
    )
    equivalence = {
        "same_mesh_abs_energy_difference_hartree": same_mesh_delta,
        "cross_mesh_transfer_vs_cold_abs_energy_difference_hartree": cross_mesh_delta,
        "energy_tolerance_hartree": ENERGY_TOLERANCE_HARTREE,
        "same_mesh_energy_equivalent": same_mesh_delta <= ENERGY_TOLERANCE_HARTREE,
        "cross_mesh_energy_equivalent": cross_mesh_delta <= ENERGY_TOLERANCE_HARTREE,
        "same_mesh_restart_reduces_scf_steps": results["same3"]["scf_steps"]
        < results["cold3"]["scf_steps"],
        "cross_mesh_transfer_reduces_scf_steps": results["transfer4"]["scf_steps"]
        < results["cold4"]["scf_steps"],
        "restart_payload_hash_matches": restart_hash_passes,
        "raw_manifest_complete_and_valid": raw_manifest_passes,
        "raw_file_count": len(raw_manifest),
    }
    passed = all(case["status"] == "PASS" for case in results.values()) and all(
        value
        for key, value in equivalence.items()
        if key
        in {
            "same_mesh_energy_equivalent",
            "cross_mesh_energy_equivalent",
            "same_mesh_restart_reduces_scf_steps",
            "cross_mesh_transfer_reduces_scf_steps",
            "restart_payload_hash_matches",
            "raw_manifest_complete_and_valid",
        }
    )
    payload = {
        "schema": "periodic-gxtb-restart-equivalence-v1",
        "status": "PASS" if passed else "FAIL",
        "qualified_binary_sha256": QUALIFIED_BINARY_SHA256,
        "cases": results,
        "equivalence": equivalence,
        "interpretation": (
            "A strict same-mesh restart is energy-identical to its source and "
            "converges in one SCF step. A validated 3x3x3-to-4x4x4 BvK transfer "
            "reaches the independent cold 4x4x4 energy within the declared "
            "tolerance while reducing the SCF iteration count."
        ),
    }
    (HERE / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
