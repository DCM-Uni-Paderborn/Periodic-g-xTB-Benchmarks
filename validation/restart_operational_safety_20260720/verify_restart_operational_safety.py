#!/usr/bin/env python3
"""Verify archived CP2K k-point restart controls and the recovery helpers."""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE / "raw"
QUALIFIED_BINARY_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
REFERENCE_ENERGY_HARTREE = -40.468866070692428
REFERENCE_CHECKPOINT = (
    ROOT
    / "validation/gxtb_restart_equivalence_20260720/raw/work/"
    "CH4_gxtb_kp_restart_3-RESTART.kp"
)
SOURCE_INPUT = (
    ROOT
    / "validation/gxtb_restart_equivalence_20260720/raw/inputs/"
    "CH4_gxtb_kp_restart_3.inp"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_hash(path: Path) -> str:
    return path.read_text(encoding="utf-8").split()[0]


def verify_manifest(path: Path, base: Path, expected: list[str]) -> bool:
    recorded: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        recorded[relative.lstrip("* ")] = digest
    return sorted(recorded) == sorted(expected) and all(
        sha256(base / relative) == digest for relative, digest in recorded.items()
    )


def read_manifest(path: Path) -> dict[str, str]:
    recorded: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, relative = line.split(maxsplit=1)
            recorded[relative.lstrip("* ")] = digest
    return recorded


def parse_output(path: Path) -> tuple[str, int, float]:
    text = path.read_text(encoding="utf-8")
    steps = re.findall(r"SCF run converged in\s+([0-9]+) steps", text)
    energies = re.findall(
        r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+([-+0-9.Ee]+)",
        text,
    )
    if len(steps) != 1 or len(energies) != 1:
        raise RuntimeError(f"could not uniquely parse {path}")
    return text, int(steps[0]), float(energies[0])


def singleton_affinity(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"expected_cpu=([0-9]+) allowed=([0-9]+)", text)
    return bool(
        match
        and match.group(1) == match.group(2)
        and re.search(
            rf"^Cpus_allowed_list:\s*{match.group(1)}\s*$", text, re.MULTILINE
        )
    )


def common_checks(case: Path, input_name: str = "input.inp") -> dict[str, bool]:
    output = (case / "cp2k.out").read_text(encoding="utf-8")
    return {
        "exit_status_zero": (case / "exit_status").read_text().strip() == "0",
        "program_ended": "PROGRAM ENDED AT" in output,
        "qualified_binary": first_hash(case / "binary.sha256")
        == QUALIFIED_BINARY_SHA256,
        "input_hash_matches": first_hash(case / "input.sha256")
        == sha256(case / input_name),
        "singleton_cpu_affinity": singleton_affinity(case / "affinity_preexec.txt"),
    }


def helper_self_tests() -> dict[str, bool]:
    helper = HERE / "prepare_same_mesh_restart.py"
    with tempfile.TemporaryDirectory(prefix="restart-safety-") as temporary:
        work = Path(temporary)
        safe_checkpoint = work / "source-RESTART.kp"
        unsafe_checkpoint = work / "source checkpoint-RESTART.kp"
        shutil.copyfile(REFERENCE_CHECKPOINT, safe_checkpoint)
        shutil.copyfile(REFERENCE_CHECKPOINT, unsafe_checkpoint)

        safe_output = work / "safe.inp"
        safe = subprocess.run(
            [sys.executable, str(helper), str(SOURCE_INPUT), str(safe_checkpoint), str(safe_output)],
            text=True,
            capture_output=True,
            check=False,
        )
        safe_text = safe_output.read_text(encoding="utf-8") if safe_output.exists() else ""
        unsafe = subprocess.run(
            [
                sys.executable,
                str(helper),
                str(SOURCE_INPUT),
                str(unsafe_checkpoint),
                str(work / "unsafe.inp"),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        restart_lines = re.findall(
            r"^\s*WFN_RESTART_FILE_NAME\s+(.+)$", safe_text, re.MULTILINE
        )
        return {
            "safe_input_generated": safe.returncode == 0 and safe_output.is_file(),
            "safe_path_is_unquoted": len(restart_lines) == 1
            and '"' not in restart_lines[0]
            and "'" not in restart_lines[0],
            "scf_guess_changed_to_restart": bool(
                re.search(r"^\s*SCF_GUESS\s+RESTART\s*$", safe_text, re.MULTILINE)
            ),
            "unsafe_path_rejected": unsafe.returncode != 0
            and not (work / "unsafe.inp").exists(),
        }


def pending_input_checks() -> dict[str, bool]:
    pending = RAW / "pending_inputs"
    specifications = {
        "VII-k999.inp": ("VII", 9, "inputs/k999-reduced/VII/input.inp"),
        "XIV-k888.inp": ("XIV", 8, "inputs/k888-reduced/XIV/input.inp"),
        "XI-k777.inp": ("XI", 7, "inputs/k777-reduced/XI/input.inp"),
        "XVII-k777.inp": ("XVII", 7, "inputs/k777-reduced/XVII/input.inp"),
    }
    input_manifest = read_manifest(pending / "input_sha256")
    remote_manifest = read_manifest(pending / "remote_safety_sha256")
    checks: dict[str, bool] = {
        "input_manifest_has_exact_endpoint_set": set(input_manifest)
        == {item[2] for item in specifications.values()},
        "remote_manifest_includes_exact_endpoint_set": set(input_manifest).issubset(
            remote_manifest
        ),
    }
    for local_name, (phase, mesh, remote_name) in specifications.items():
        path = pending / local_name
        text = path.read_text(encoding="utf-8")
        digest = sha256(path)
        checks[f"{phase}_input_hash_matches"] = (
            input_manifest.get(remote_name) == digest
            and remote_manifest.get(remote_name) == digest
        )
        checks[f"{phase}_project_matches_mesh"] = bool(
            re.search(
                rf"^\s*PROJECT\s+ice_{phase}_GXTB_k{mesh}{mesh}{mesh}\s*$",
                text,
                re.MULTILINE,
            )
        )
        checks[f"{phase}_macdonald_mesh_matches"] = bool(
            re.search(
                rf"^\s*SCHEME\s+MACDONALD\s+{mesh}\s+{mesh}\s+{mesh}\b",
                text,
                re.MULTILINE,
            )
        )
        checks[f"{phase}_checkpoint_writing_enabled"] = bool(
            re.search(r"^\s*&RESTART\s+ON\s*$", text, re.MULTILINE)
            and re.search(r"^\s*QS_SCF\s+1\s*$", text, re.MULTILINE)
        )
        checks[f"{phase}_cold_start_remains_explicit"] = bool(
            re.search(r"^\s*SCF_GUESS\s+MOPAC\s*$", text, re.MULTILINE)
            and not re.search(r"^\s*WFN_RESTART_FILE_NAME\b", text, re.MULTILINE)
        )
    tool_mapping = {
        "launch_pinned_cp2k.sh": HERE / "launch_pinned_cp2k.sh",
        "tools/prepare_same_mesh_restart.py": HERE / "prepare_same_mesh_restart.py",
        "tools/resume_pinned_cp2k.sh": HERE / "resume_pinned_cp2k.sh",
    }
    checks["remote_tool_set_matches_archived_tools"] = all(
        remote_manifest.get(remote_name) == sha256(local_path)
        for remote_name, local_path in tool_mapping.items()
    )
    return checks


def main() -> None:
    known = RAW / "known_unquoted"
    quoted = RAW / "quoted_ablation"
    helper = RAW / "resume_helper"

    known_text, known_steps, known_energy = parse_output(known / "cp2k.out")
    quoted_text, quoted_steps, quoted_energy = parse_output(quoted / "cp2k.out")
    helper_text, helper_steps, helper_energy = parse_output(helper / "cp2k.out")
    known_input = (known / "input.inp").read_text(encoding="utf-8")
    quoted_input = (quoted / "input.inp").read_text(encoding="utf-8")
    helper_input = (helper / "input.inp").read_text(encoding="utf-8")

    cases = {
        "known_unquoted_strict_restart": {
            **common_checks(known),
            "acceptance_marker": "KPOINT_RESTART| Strict same-mesh restart accepted"
            in known_text,
            "one_scf_step": known_steps == 1,
            "energy_matches_reference": math.isclose(
                known_energy, REFERENCE_ENERGY_HARTREE, rel_tol=0.0, abs_tol=5.0e-15
            ),
            "unquoted_restart_name": bool(
                re.search(
                    r"^\s*WFN_RESTART_FILE_NAME\s+[^\"']+\.kp\s*$",
                    known_input,
                    re.MULTILINE,
                )
            ),
        },
        "quoted_name_negative_control": {
            **common_checks(quoted),
            "no_acceptance_marker": "KPOINT_RESTART|" not in quoted_text,
            "missing_file_warning": "This file does not exist" in quoted_text,
            "cold_guess_step_count": quoted_steps == 12,
            "energy_still_matches_after_cold_convergence": math.isclose(
                quoted_energy, REFERENCE_ENERGY_HARTREE, rel_tol=0.0, abs_tol=5.0e-15
            ),
            "quoted_restart_name": bool(
                re.search(
                    r'^\s*WFN_RESTART_FILE_NAME\s+"[^\"]+\.kp"\s*$',
                    quoted_input,
                    re.MULTILINE,
                )
            ),
        },
        "recovery_helper_positive_control": {
            **common_checks(helper),
            "wrapper_acceptance_pass": (helper / "restart_acceptance").read_text().strip()
            == "PASS",
            "validated_transfer_marker": (
                "KPOINT_RESTART| Validated BvK mesh transfer accepted" in helper_text
            ),
            "recorded_mode_matches": (helper / "restart_acceptance_mode").read_text().strip()
            == "validated_bvk_transfer",
            "one_scf_step": helper_steps == 1,
            "energy_matches_reference": math.isclose(
                helper_energy, REFERENCE_ENERGY_HARTREE, rel_tol=0.0, abs_tol=5.0e-15
            ),
            "unquoted_restart_name": bool(
                re.search(
                    r"^\s*WFN_RESTART_FILE_NAME\s+[^\"']+\.kp\s*$",
                    helper_input,
                    re.MULTILINE,
                )
            ),
            "restart_writing_remains_enabled": bool(
                re.search(r"^\s*&RESTART\s+ON\s*$", helper_input, re.MULTILINE)
            ),
            "source_checkpoint_identity": first_hash(helper / "source-checkpoint.sha256")
            == sha256(REFERENCE_CHECKPOINT),
        },
    }
    case_results = {
        name: {
            "status": "PASS" if all(checks.values()) else "FAIL",
            "checks": checks,
        }
        for name, checks in cases.items()
    }

    raw_files = sorted(
        path.relative_to(RAW).as_posix() for path in RAW.rglob("*") if path.is_file()
    )
    tool_files = [
        "launch_pinned_cp2k.sh",
        "prepare_same_mesh_restart.py",
        "resume_pinned_cp2k.sh",
    ]
    helper_checks = helper_self_tests()
    pending_checks = pending_input_checks()
    integrity = {
        "raw_manifest_complete_and_valid": verify_manifest(
            HERE / "RAW_SHA256SUMS", RAW, raw_files
        ),
        "tool_manifest_complete_and_valid": verify_manifest(
            HERE / "TOOL_SHA256SUMS", HERE, tool_files
        ),
        "reference_checkpoint_exists": REFERENCE_CHECKPOINT.is_file(),
        "reference_checkpoint_hash_matches_archived_source": (
            first_hash(helper / "source-checkpoint.sha256")
            == sha256(REFERENCE_CHECKPOINT)
        ),
        **helper_checks,
        **pending_checks,
    }
    passed = all(item["status"] == "PASS" for item in case_results.values()) and all(
        integrity.values()
    )
    payload = {
        "schema": "periodic-gxtb-restart-operational-safety-v1",
        "status": "PASS" if passed else "FAIL",
        "qualified_binary_sha256": QUALIFIED_BINARY_SHA256,
        "reference_checkpoint_sha256": sha256(REFERENCE_CHECKPOINT),
        "cases": case_results,
        "integrity": integrity,
        "interpretation": (
            "An unquoted k-point checkpoint is accepted and reproduces the reference "
            "energy in one SCF step. Quoting the same WFN_RESTART_FILE_NAME causes a "
            "silent cold-guess fallback despite exit status zero. The recovery helpers "
            "therefore emit only safe unquoted names, preserve the source checkpoint, "
            "and require an explicit CP2K k-point restart acceptance marker."
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
