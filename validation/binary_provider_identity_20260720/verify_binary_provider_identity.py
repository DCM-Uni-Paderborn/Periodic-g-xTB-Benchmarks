#!/usr/bin/env python3
"""Verify the archived Part-I CLI/CP2K provider-identity snapshot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SNAPSHOT = HERE / "remote_snapshot.json"
OUTPUT = HERE / "verification.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    snapshot = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    cp2k = snapshot["cp2k"]
    cli = snapshot["direct_cli"]
    archive = snapshot["provider_archive"]
    runtime = snapshot["runtime"]
    successor = snapshot["successor_scope"]

    checks = {
        "same_provider_revision": (
            cp2k["tblite_revision_from_cmake"] == cli["source_revision"]
        ),
        "save_provider_selected": cp2k["tblite_provider"] == "SAVE",
        "direct_cli_source_clean": cli["source_status_clean"] is True,
        "cli_link_uses_recorded_archive": (
            cli["link_rule_depends_on_build_tree_libtblite_archive"] is True
        ),
        "same_provider_archive_for_build_install_and_cp2k": (
            len(set(archive.values())) == 1
        ),
        "same_compiler_runtime": runtime["cp2k_and_cli_same_gfortran"] is True,
        "same_blas_runtime": runtime["cp2k_and_cli_same_openblas"] is True,
        "later_interface_change_not_entered_by_production_mesh": (
            successor["interface_change_limited_to_general_mesh_inference"] is True
            and successor["production_inputs_use_general_mesh"] is False
            and successor["production_inputs_use_macdonald_mesh"] is True
        ),
        "later_restart_change_not_entered_by_production_runs": (
            successor["restart_change_limited_to_alias_projection_diagnostic"] is True
            and successor["production_inputs_use_restart"] is False
        ),
    }
    passed = all(checks.values())
    output = {
        "schema": "periodic-gxtb-binary-provider-identity-verification-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "cp2k_binary_sha256": cp2k["binary_sha256"],
        "direct_cli_binary_sha256": cli["binary_sha256"],
        "provider_revision": cli["source_revision"],
        "provider_archive_sha256": archive["build_tree_sha256"],
        "snapshot_file": SNAPSHOT.name,
        "snapshot_sha256": sha256(SNAPSHOT),
    }
    OUTPUT.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
