#!/usr/bin/env python3
"""Re-run every completed Part-I evidence verifier from a clean checkout."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "requalification.json"
PYTHON = sys.executable

QUALIFIED_CP2K_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
QUALIFIED_CLI_SHA256 = (
    "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
)

STANDARD_VERIFIERS = (
    ("accuracy_equivalence", "validation/accuracy_equivalence_20260720/verify_accuracy_equivalence.py"),
    ("binary_provider_identity", "validation/binary_provider_identity_20260720/verify_binary_provider_identity.py"),
    ("cecl3_tolerance_recheck", "validation/cecl3_tolerance_recheck_20260720/verify_cecl3_tolerance_recheck.py"),
    ("energy_component_ablation", "validation/dmc13_k222_viii_component_ablation_20260719/verify_component_ablation.py"),
    ("derivative_component_ablation", "validation/dmc13_k222_xvii_derivative_component_ablation_20260719/verify_derivative_component_ablation.py"),
    ("geometry_equivalence", "validation/geometry_equivalence_20260720/verify_geometry_equivalence.py"),
    ("restart_equivalence", "validation/gxtb_restart_equivalence_20260720/verify_restart_equivalence.py"),
    ("macdonald_bvk_mesh_equivalence", "validation/macdonald_bvk_mesh_equivalence_20260720/verify_macdonald_bvk_mesh.py"),
    ("mstore_accuracy_equivalence", "validation/mstore_accuracy_equivalence_20260720/verify_mstore_accuracy_equivalence.py"),
    ("native_cli_full_parity", "validation/native_cli_full_parity_20260720/verify_native_cli_full_parity.py"),
    ("native_cli_inprocess_derivatives", "validation/native_cli_inprocess_derivatives_20260720/verify_inprocess_derivative_parity.py"),
    ("pbc_h0_attribution", "validation/pbc_h0_anisotropy_attribution_20260719/verify_h0_attribution.py"),
    ("provider_component_attribution", "validation/provider_component_attribution_20260719/verify_attribution.py"),
    ("qualified_build_head_delta", "validation/qualified_build_head_delta_20260720/verify_qualified_build_head_delta.py"),
    ("relative_energy_postprocessing", "validation/relative_energy_postprocessing_20260720/verify_relative_energy_postprocessing.py"),
    ("restart_operational_safety", "validation/restart_operational_safety_20260720/verify_restart_operational_safety.py"),
    ("save_tblite_periodic_source_tests", "validation/save_tblite_periodic_source_tests_20260719/verify_source_tests.py"),
    ("seidler_package_tracking", "validation/seidler_package_tracking_20260720/verify_package_tracking.py"),
    ("tight_parity_k222", "validation/tight_parity_k222_20260720/verify_tight_parity.py"),
)

MANIFESTS = (
    "validation/explicit_cp2k_gamma_supercell_oracle_20260719/SHA256SUMS",
    "validation/three_route_k333_closure_20260719/SHA256SUMS",
    "validation/gxtb_final_lowk_derivatives_20260719/SHA256SUMS",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tracked_status() -> str:
    completed = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout


def normalized_command(command: list[str]) -> list[str]:
    normalized = []
    for item in command:
        text = str(item)
        if text == PYTHON:
            normalized.append("python3")
            continue
        path = Path(text)
        if path.is_absolute():
            try:
                normalized.append("<repository>/" + path.relative_to(ROOT).as_posix())
                continue
            except ValueError:
                pass
        normalized.append(text)
    return normalized


def run_check(name: str, command: list[str], display: list[str] | None = None) -> dict:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        timeout=180,
    )
    record = {
        "name": name,
        "status": "PASS" if completed.returncode == 0 else "FAIL",
        "returncode": completed.returncode,
        "command": normalized_command(display or command),
        "stdout_sha256": sha256_bytes(completed.stdout),
        "stderr_sha256": sha256_bytes(completed.stderr),
    }
    if completed.returncode != 0:
        sys.stdout.buffer.write(completed.stdout)
        sys.stderr.buffer.write(completed.stderr)
        raise RuntimeError(f"{name} failed with exit status {completed.returncode}")
    return record


def compare_json(actual: Path, expected: Path, label: str) -> dict:
    actual_payload = json.loads(actual.read_text(encoding="utf-8"))
    expected_payload = json.loads(expected.read_text(encoding="utf-8"))
    identical = actual_payload == expected_payload
    result = {
        "name": label,
        "status": "PASS" if identical else "FAIL",
        "actual_sha256": sha256_file(actual),
        "expected_sha256": sha256_file(expected),
        "json_identical": identical,
    }
    if not identical:
        raise RuntimeError(f"{label} did not reproduce its archived JSON")
    return result


def check_manifest(relative: str) -> dict:
    manifest = ROOT / relative
    base = manifest.parent.resolve()
    entries = []
    for number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  ([*]?)(.+)", raw)
        if match is None:
            raise RuntimeError(f"invalid SHA-256 manifest line {relative}:{number}")
        expected, _, written_path = match.groups()
        candidate = (base / written_path).resolve()
        if candidate != base and base not in candidate.parents:
            raise RuntimeError(f"manifest path escapes its evidence directory: {written_path}")
        if not candidate.is_file():
            raise RuntimeError(f"manifest member is missing: {candidate.relative_to(ROOT)}")
        actual = sha256_file(candidate)
        if actual != expected:
            raise RuntimeError(f"manifest mismatch: {candidate.relative_to(ROOT)}")
        entries.append(str(candidate.relative_to(ROOT)))
    return {
        "name": f"manifest:{relative}",
        "status": "PASS",
        "entry_count": len(entries),
        "manifest_sha256": sha256_file(manifest),
    }


def main() -> None:
    initial_status = tracked_status()
    if initial_status:
        raise RuntimeError(
            "refusing requalification because tracked files are already modified:\n"
            + initial_status
        )

    checks = []
    for name, relative in STANDARD_VERIFIERS:
        checks.append(run_check(name, [PYTHON, relative]))

    with tempfile.TemporaryDirectory(prefix="gxtb-part-i-requalification-") as temporary:
        temporary_dir = Path(temporary)

        lowk_output = temporary_dir / "lowk.json"
        checks.append(
            run_check(
                "lowk_derivatives_and_partial_pbc",
                [
                    PYTHON,
                    "validation/gxtb_final_lowk_derivatives_20260719/verify_final_lowk_derivatives.py",
                    "validation/gxtb_final_lowk_derivatives_20260719",
                    "--legacy-partial-manifest",
                    "validation/gxtb_final_lowk_derivatives_20260719/legacy_partial_pbc_manifest.json",
                    "--output",
                    str(lowk_output),
                ],
                [
                    PYTHON,
                    "validation/gxtb_final_lowk_derivatives_20260719/verify_final_lowk_derivatives.py",
                    "validation/gxtb_final_lowk_derivatives_20260719",
                    "--legacy-partial-manifest",
                    "validation/gxtb_final_lowk_derivatives_20260719/legacy_partial_pbc_manifest.json",
                    "--output",
                    "<temporary>/lowk.json",
                ],
            )
        )
        checks.append(
            compare_json(
                lowk_output,
                ROOT / "validation/gxtb_final_lowk_derivatives_20260719/verification.reproduced.json",
                "lowk_reproduced_json",
            )
        )

        checks.append(
            run_check(
                "three_route_k333_closure",
                [PYTHON, "validation/three_route_k333_closure_20260719/generate_and_verify.py"],
            )
        )

        oracle_root = ROOT / "validation/explicit_cp2k_gamma_supercell_oracle_20260719"
        oracle_output = temporary_dir / "oracle.json"
        oracle_command = [
            PYTHON,
            str(oracle_root / "scripts/compare_gamma_supercell_oracle.py"),
            str(oracle_root / "results/native_k222/cp2k.out"),
            str(oracle_root / "results/gamma_supercell_k222/cp2k.out"),
            str(oracle_root / "results/direct_cli_k222/result.json"),
            "--replicas",
            "8",
            "--parity-tolerance",
            "2e-7",
            "--require-binary-sha256",
            QUALIFIED_CP2K_SHA256,
            "--require-cli-binary-sha256",
            QUALIFIED_CLI_SHA256,
            "--native-input",
            str(oracle_root / "inputs/native_k222/input.inp"),
            "--gamma-input",
            str(oracle_root / "inputs/XVII/input.inp"),
            "--cli-input",
            str(oracle_root / "inputs/direct_cli_k222/POSCAR"),
            "--require-native-input-sha256",
            "265692cdca09516e68b10022a5291326c1cc39fc664a0ad60191b15eb0999a1b",
            "--require-gamma-input-sha256",
            "bfd43e54957647f6ad6b8df8c13f2b9dbc34cde41038ea55b4ec5c63c8abbec1",
            "--require-cli-input-sha256",
            "17446bfc858d189bb9e48745c7e14b470ba2b72bc56e7db50afa2864a101b8c7",
            "--output",
            str(oracle_output),
        ]
        oracle_display = [
            "<repository>/" + item.relative_to(ROOT).as_posix()
            if isinstance(item, Path) and item.is_relative_to(ROOT)
            else "<temporary>/oracle.json"
            if str(item) == str(oracle_output)
            else str(item)
            for item in oracle_command
        ]
        checks.append(
            run_check(
                "explicit_gamma_supercell_oracle",
                [str(item) for item in oracle_command],
                oracle_display,
            )
        )
        checks.append(
            compare_json(
                oracle_output,
                oracle_root / "verification.reproduced.json",
                "gamma_supercell_oracle_reproduced_json",
            )
        )

    for relative in MANIFESTS:
        checks.append(check_manifest(relative))

    checks.append(
        run_check(
            "implementation_audit",
            [PYTHON, "validation/implementation_audit_20260720/verify_implementation_audit.py"],
        )
    )
    aggregate = json.loads((HERE / "verification.json").read_text(encoding="utf-8"))
    aggregate_valid = (
        aggregate.get("status") == "PASS"
        and aggregate.get("completed_gate_count") == 23
        and all(
            item.get("passed") is True
            for item in aggregate.get("completed_gates", {}).values()
        )
    )
    checks.append(
        {
            "name": "aggregate_semantics",
            "status": "PASS" if aggregate_valid else "FAIL",
            "completed_gate_count": aggregate.get("completed_gate_count"),
            "all_completed_gates_pass": aggregate_valid,
            "verification_sha256": sha256_file(HERE / "verification.json"),
        }
    )
    if not aggregate_valid:
        raise RuntimeError("the regenerated implementation audit is not fully passing")

    final_status = tracked_status()
    if final_status:
        raise RuntimeError(
            "requalification changed tracked evidence instead of reproducing it:\n"
            + final_status
        )

    payload = {
        "schema": "periodic-gxtb-part-i-completed-evidence-requalification-v1",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "tracked_checkout_clean_after_reproduction": True,
        "interpretation": (
            "Every completed Part-I verifier was executed again from a clean tracked "
            "checkout. Regenerated JSON, derived tables, aggregate gates, and selected "
            "SHA-256 manifests reproduce the archived evidence exactly. The remaining "
            "adaptive DMC-ICE13 endpoints are science calculations, not missing "
            "implementation-verification gates."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
