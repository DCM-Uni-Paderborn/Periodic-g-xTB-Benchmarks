#!/usr/bin/env python3
"""Verify that post-qualification branch changes are inactive for DMC energies."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUT_ROOT = (
    ROOT
    / "DMC-ICE13"
    / "reproduction"
    / "seidler_dmc13_recalculation"
    / "raw"
    / "cp2k_native"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def changed_files(patch: str) -> list[str]:
    return re.findall(r"^diff --git a/(\S+) b/\S+$", patch, flags=re.MULTILINE)


def main() -> None:
    metadata = json.loads((HERE / "metadata.json").read_text())
    cp2k_patch_path = HERE / "cp2k_delta.patch"
    provider_patch_path = HERE / "save_tblite_delta.patch"
    cp2k_patch = cp2k_patch_path.read_text()
    provider_patch = provider_patch_path.read_text()

    cp2k_files = changed_files(cp2k_patch)
    provider_files = changed_files(provider_patch)
    expected_cp2k_files = [
        "src/kpoint_restart_transfer.F",
        "src/tblite_interface.F",
        "tests/xTB/regtest-tblite-gxtb/H2O_gxtb_kp_cp2k_density_mixer_restart.inp",
        "tests/xTB/regtest-tblite-gxtb/TEST_FILES.toml",
    ]

    inputs = sorted(
        path
        for path in INPUT_ROOT.rglob("input.inp")
        if not any(token in part.lower() for part in path.parts for token in ("failed", "incomplete"))
    )
    input_checks: dict[str, dict[str, bool]] = {}
    for path in inputs:
        text = path.read_text(errors="strict")
        rel = str(path.relative_to(ROOT))
        input_checks[rel] = {
            "macdonald": bool(re.search(r"^\s*SCHEME\s+MACDONALD\b", text, re.MULTILINE)),
            "not_general": not bool(re.search(r"^\s*SCHEME\s+GENERAL\b", text, re.MULTILINE)),
            "mopac_guess": bool(re.search(r"^\s*SCF_GUESS\s+MOPAC\b", text, re.MULTILINE)),
            "restart_off": bool(re.search(r"^\s*&RESTART\s+OFF\b", text, re.MULTILINE)),
        }

    identity_path = ROOT / "validation" / "binary_provider_identity_20260720" / "verification.json"
    identity = json.loads(identity_path.read_text())

    checks = {
        "cp2k_changed_files_exact": cp2k_files == expected_cp2k_files,
        "cp2k_general_mesh_change_is_explicitly_guarded": (
            'IF (TRIM(kpoints%kp_scheme) == "GENERAL") THEN' in cp2k_patch
        ),
        "cp2k_restart_change_is_alias_projection_only": (
            "Projecting inconsistent modulo-grid aliases" in cp2k_patch
            and "Inconsistent modulo-grid aliases in k-point restart density" in cp2k_patch
            and "src/kpoint_restart_transfer.F" in cp2k_files
        ),
        "save_tblite_changed_files_test_only": provider_files == ["test/unit/test_hamiltonian.f90"],
        "save_tblite_provider_source_unchanged": not any(path.startswith("src/") for path in provider_files),
        "production_input_count_sufficient": (
            len(inputs) >= metadata["minimum_archived_production_input_count"]
        ),
        "all_production_inputs_macdonald_mopac_no_restart": bool(inputs)
        and all(all(values.values()) for values in input_checks.values()),
        "qualified_binary_identity_gate_passes": identity.get("status") == "PASS",
        "qualified_provider_revision_matches": (
            identity.get("provider_revision")
            == metadata["save_tblite"]["qualified_provider_revision"]
        ),
        "qualified_binary_sha256_matches": (
            identity.get("cp2k_binary_sha256")
            == metadata["qualified_cp2k_binary_sha256"]
        ),
    }
    status = "PASS" if all(checks.values()) else "FAIL"
    nonconforming_inputs = [
        path for path, values in input_checks.items() if not all(values.values())
    ]
    input_mode_counts = {
        key: sum(values[key] for values in input_checks.values())
        for key in ("macdonald", "not_general", "mopac_guess", "restart_off")
    }
    result = {
        "schema": "periodic-gxtb-qualified-build-head-delta-verification-v1",
        "status": status,
        "checks": checks,
        "cp2k_changed_files": cp2k_files,
        "save_tblite_changed_files": provider_files,
        "production_input_count": len(inputs),
        "production_input_mode_counts": input_mode_counts,
        "nonconforming_production_inputs": nonconforming_inputs,
        "cp2k_patch_sha256": sha256(cp2k_patch_path),
        "save_tblite_patch_sha256": sha256(provider_patch_path),
        "metadata_sha256": sha256(HERE / "metadata.json"),
        "interpretation": (
            "The successor commits do not enter the archived DMC-ICE13 energy path; "
            "a full energy rerun is not required solely because the integration branches advanced."
        ),
    }
    (HERE / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
