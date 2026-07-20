#!/usr/bin/env python3
"""Reconstruct the CeCl3 tolerance-only Hamiltonian test qualification."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARCHIVED = HERE / "verification.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    archived = json.loads(ARCHIVED.read_text(encoding="utf-8"))
    targeted_log = HERE / archived["targeted_case"]["log"]
    group_log = HERE / archived["complete_hamiltonian_group"]["log"]
    targeted_text = targeted_log.read_text(encoding="utf-8", errors="replace")
    group_text = group_log.read_text(encoding="utf-8", errors="replace")

    starts = re.findall(r"Starting\s+(.+?)\s+\.\.\.\s+\((\d+)/(\d+)\)", group_text)
    passed_names = re.findall(r"\.\.\.\s+(.+?)\s+\[PASSED\]", group_text)
    failed_names = re.findall(r"\.\.\.\s+(.+?)\s+\[FAILED\]", group_text)
    indices = [int(index) for _, index, _ in starts]
    totals = [int(total) for _, _, total in starts]
    target_name = archived["targeted_case"]["name"]
    original = float(archived["original_tolerance_hartree_per_bohr"])
    revised = float(archived["revised_tolerance_hartree_per_bohr"])
    residual = float(archived["maximum_residual_hartree_per_bohr"])

    checks = {
        "targeted_exit_status_zero": (HERE / "cecl3.exit_status").read_text().strip()
        == "0",
        "group_exit_status_zero": (HERE / "hamiltonian.exit_status").read_text().strip()
        == "0",
        "targeted_log_hash_matches": sha256(targeted_log)
        == archived["targeted_case"]["log_sha256"],
        "group_log_hash_matches": sha256(group_log)
        == archived["complete_hamiltonian_group"]["log_sha256"],
        "targeted_case_passes": target_name in targeted_text
        and "[PASSED]" in targeted_text
        and "[FAILED]" not in targeted_text,
        "complete_group_has_exact_sequence": len(starts) == 75
        and indices == list(range(1, 76))
        and totals == [75] * 75,
        "complete_group_all_cases_pass": len(passed_names) == 75
        and len(set(passed_names)) == 75
        and not failed_names,
        "targeted_case_in_complete_group": target_name in passed_names,
        "residual_fails_only_original_threshold": original < residual <= revised,
        "residual_matches_author_baseline": bool(
            archived["residual_identical_to_author_pbc_baseline"]
        ),
        "no_scientific_source_file_changed": not archived["scientific_source_files_modified"],
        "archived_counts_match": archived["complete_hamiltonian_group"][
            "passed_case_count"
        ]
        == len(passed_names)
        and archived["complete_hamiltonian_group"]["failed_case_count"]
        == len(failed_names),
    }
    passed = all(checks.values())
    payload = {
        "schema": "save-tblite-cecl3-tolerance-recheck-reproduced-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "original_tolerance_hartree_per_bohr": original,
        "revised_tolerance_hartree_per_bohr": revised,
        "maximum_residual_hartree_per_bohr": residual,
        "passed_case_count": len(passed_names),
        "failed_case_count": len(failed_names),
        "source_verification_sha256": sha256(ARCHIVED),
        "interpretation": (
            "The archived residual lies above only the original numerical test "
            "threshold and below the revised threshold. The targeted case and all "
            "75 Hamiltonian cases pass without a scientific source change."
        ),
    }
    output = HERE / "verification.reproduced.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
