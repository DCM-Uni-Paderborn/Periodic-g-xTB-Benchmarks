#!/usr/bin/env python3
"""Verify CP2K-native versus direct save_tblite energy/derivative parity."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CP2K_SHA256 = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
CLI_SHA256 = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
ENERGY_LIMIT = 1.0e-9
DERIVATIVE_LIMIT = 1.0e-5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recorded_sha(path: Path) -> str:
    return path.read_text().split()[0]


def verify_raw_manifest() -> tuple[bool, int]:
    entries: dict[str, str] = {}
    for line in (HERE / "RAW_SHA256SUMS").read_text().splitlines():
        digest, relative = line.split(maxsplit=1)
        if relative in entries:
            return False, len(entries)
        entries[relative] = digest
    expected = {
        str(path.relative_to(HERE))
        for parent in (HERE / "inputs", HERE / "raw")
        for path in parent.rglob("*")
        if path.is_file()
    }
    passed = set(entries) == expected and all(
        sha256(HERE / relative) == digest for relative, digest in entries.items()
    )
    return passed, len(entries)


def parse_pair(text: str, label: str) -> tuple[float, float]:
    match = re.search(rf"^\s*{re.escape(label)}:\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)", text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing {label!r}")
    return float(match.group(1)), float(match.group(2))


def parse_energy_triplet(text: str) -> tuple[float, float, float]:
    match = re.search(
        r"^\s*Energy CP2K/CLI/absdiff:\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)",
        text,
        re.MULTILINE,
    )
    if not match:
        raise ValueError("missing in-process energy comparison")
    return tuple(float(match.group(index)) for index in range(1, 4))


def parse_total_energy(text: str) -> float:
    matches = re.findall(
        r"^\s*ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+([+\-0-9.Ee]+)",
        text,
        re.MULTILINE,
    )
    if len(matches) != 1:
        raise ValueError(f"expected exactly one CP2K total energy, found {len(matches)}")
    return float(matches[0])


def verify_case(name: str, input_name: str, require_virial: bool) -> dict:
    case_dir = HERE / "raw" / name
    output_path = case_dir / "cp2k.out"
    input_path = HERE / "inputs" / input_name
    text = output_path.read_text(errors="strict")
    reported_cp_energy, reported_cli_energy, reported_energy_difference = parse_energy_triplet(text)
    cp_energy = parse_total_energy(text)
    cli_payload = json.loads((case_dir / "tblite-reference.json").read_text())
    cli_energy = float(cli_payload["energy"])
    energy_difference = abs(cp_energy - cli_energy)
    gradient_sum, gradient_max = parse_pair(text, "Gradient diff sum/max")
    virial_sum = None
    virial_max = None
    if require_virial:
        virial_sum, virial_max = parse_pair(text, "Virial diff sum/max")

    checks = {
        "exit_zero": (case_dir / "exit_status").read_text().strip() == "0",
        "program_ended": "PROGRAM ENDED AT" in text,
        "cp2k_binary_sha256": recorded_sha(case_dir / "binary.sha256") == CP2K_SHA256,
        "executed_input_sha256": recorded_sha(case_dir / "input.sha256") == sha256(input_path),
        "method_gxtb": bool(re.search(r"^\s*METHOD\s+GXTB\b", input_path.read_text(), re.MULTILINE)),
        "reference_cli_stop_on_error": bool(
            re.search(r"^\s*STOP_ON_ERROR\s+T\b", input_path.read_text(), re.MULTILINE)
        ),
        "energy_difference_within_limit": energy_difference <= ENERGY_LIMIT,
        "gradient_max_within_limit": gradient_max <= DERIVATIVE_LIMIT,
        "virial_max_within_limit": (not require_virial) or (virial_max is not None and virial_max <= DERIVATIVE_LIMIT),
        "reported_energy_difference_consistent": math.isclose(
            energy_difference, reported_energy_difference, rel_tol=0.0, abs_tol=1.0e-22
        ),
        "printed_cp_energy_consistent": math.isclose(
            cp_energy, reported_cp_energy, rel_tol=0.0, abs_tol=5.0e-11
        ),
        "printed_cli_energy_consistent": math.isclose(
            cli_energy, reported_cli_energy, rel_tol=0.0, abs_tol=5.0e-11
        ),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "cp2k_energy_hartree": cp_energy,
        "cli_energy_hartree": cli_energy,
        "absolute_energy_difference_hartree": energy_difference,
        "reported_absolute_energy_difference_hartree": reported_energy_difference,
        "gradient_difference_sum_hartree_per_bohr": gradient_sum,
        "gradient_difference_max_hartree_per_bohr": gradient_max,
        "virial_difference_sum_hartree": virial_sum,
        "virial_difference_max_hartree": virial_max,
        "output_sha256": sha256(output_path),
        "input_sha256": sha256(input_path),
        "cli_json_sha256": sha256(case_dir / "tblite-reference.json"),
        "cli_gradient_file_sha256": sha256(case_dir / "tblite-reference.grad"),
    }


def main() -> None:
    h2o = verify_case(
        "h2o_keep_files", "H2O_gxtb_reference_cli_keep_files.inp", require_virial=True
    )
    o2 = verify_case(
        "o2_keep_files", "O2_gxtb_uks_reference_cli_keep.inp", require_virial=False
    )
    direct_dir = HERE / "raw" / "o2_direct_cli"
    direct_json = json.loads((direct_dir / "result.json").read_text())
    direct_energy = float(direct_json["energy"])
    direct_checks = {
        "exit_zero": (direct_dir / "exit_status").read_text().strip() == "0",
        "cli_binary_sha256": recorded_sha(direct_dir / "binary.sha256") == CLI_SHA256,
        "input_sha256": recorded_sha(direct_dir / "input.sha256") == sha256(direct_dir / "input.gen"),
        "energy_matches_inprocess_cli": math.isclose(
            direct_energy, o2["cli_energy_hartree"], rel_tol=0.0, abs_tol=5.0e-11
        ),
        "result_json_matches_inprocess_cli": (
            sha256(direct_dir / "result.json") == o2["cli_json_sha256"]
        ),
        "gradient_file_matches_inprocess_cli": (
            sha256(direct_dir / "result.grad") == o2["cli_gradient_file_sha256"]
        ),
    }
    program_identity = dict(
        line.split("=", 1)
        for line in (HERE / "cli_program_identity.txt").read_text().splitlines()
        if "=" in line
    )
    program_identity_checks = {
        "program_name_matches_inputs": all(
            program_identity["program_name"] in (HERE / "inputs" / name).read_text()
            for name in (
                "H2O_gxtb_reference_cli_keep_files.inp",
                "O2_gxtb_uks_reference_cli_keep.inp",
            )
        ),
        "resolved_cli_sha256_qualified": program_identity["sha256"] == CLI_SHA256,
    }
    preflight = HERE / "raw" / "preflight_program_path_failure"
    preflight_text = (preflight / "cp2k.out").read_text(errors="strict")
    preflight_checks = {
        "failure_retained": (preflight / "exit_status").read_text().strip() != "0",
        "failure_is_input_length_only": "has more than" in preflight_text and "characters" in preflight_text,
        "failure_precedes_energy_evaluation": "ENERGY| Total FORCE_EVAL" not in preflight_text,
    }
    raw_manifest_passes, raw_manifest_file_count = verify_raw_manifest()
    checks = {
        "h2o_passes": h2o["status"] == "PASS",
        "o2_passes": o2["status"] == "PASS",
        "independent_o2_cli_passes": all(direct_checks.values()),
        "inprocess_cli_program_identity_passes": all(program_identity_checks.values()),
        "preflight_failure_classified": all(preflight_checks.values()),
        "raw_manifest_passes": raw_manifest_passes,
    }
    status = "PASS" if all(checks.values()) else "FAIL"
    result = {
        "schema": "periodic-gxtb-native-cli-inprocess-derivative-parity-v1",
        "status": status,
        "checks": checks,
        "limits": {
            "energy_hartree": ENERGY_LIMIT,
            "gradient_hartree_per_bohr": DERIVATIVE_LIMIT,
            "virial_hartree": DERIVATIVE_LIMIT,
        },
        "qualified_cp2k_binary_sha256": CP2K_SHA256,
        "qualified_direct_cli_binary_sha256": CLI_SHA256,
        "cases": {"h2o_periodic_gamma": h2o, "o2_triplet_molecular": o2},
        "independent_o2_direct_cli": {
            "checks": direct_checks,
            "energy_hartree": direct_energy,
            "result_json_sha256": sha256(direct_dir / "result.json"),
            "gradient_file_sha256": sha256(direct_dir / "result.grad"),
        },
        "inprocess_cli_program_identity": {
            "checks": program_identity_checks,
            "record": program_identity,
            "record_sha256": sha256(HERE / "cli_program_identity.txt"),
        },
        "preflight_program_path_failure": preflight_checks,
        "raw_manifest": {
            "file": "RAW_SHA256SUMS",
            "file_count": raw_manifest_file_count,
            "sha256": sha256(HERE / "RAW_SHA256SUMS"),
        },
        "interpretation": (
            "The qualified CP2K-native path reproduces direct save_tblite energies, Cartesian gradients, "
            "and the tested periodic virial within the predeclared component limits."
        ),
    }
    (HERE / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
