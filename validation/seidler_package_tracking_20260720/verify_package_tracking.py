#!/usr/bin/env python3
"""Verify that every hashed Seidler-package file is actually tracked by Git."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE_RELATIVE = Path("DMC-ICE13/reproduction/seidler_dmc13_recalculation")
PACKAGE = ROOT / PACKAGE_RELATIVE
MANIFEST = PACKAGE / "SHA256SUMS"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    entries: dict[str, str] = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        entries[relative.lstrip("* ")] = digest

    package_files = {
        path.relative_to(PACKAGE).as_posix()
        for path in PACKAGE.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS"
        and "__pycache__" not in path.parts
    }
    tracked_repo_paths = set(
        subprocess.check_output(
            ["git", "-C", str(ROOT), "ls-files", "--", str(PACKAGE_RELATIVE)],
            text=True,
        ).splitlines()
    )
    tracked_package_paths = {
        Path(path).relative_to(PACKAGE_RELATIVE).as_posix()
        for path in tracked_repo_paths
        if Path(path).name != "SHA256SUMS"
        and "__pycache__" not in Path(path).parts
    }

    missing_on_disk = sorted(set(entries) - package_files)
    unmanifested = sorted(package_files - set(entries))
    untracked = sorted(set(entries) - tracked_package_paths)
    hash_mismatches = sorted(
        relative
        for relative, digest in entries.items()
        if relative in package_files and sha256(PACKAGE / relative) != digest
    )
    raw_text_outputs = sorted(
        relative
        for relative in entries
        if relative.endswith((".out", ".log"))
    )
    checks = {
        "all_manifest_entries_exist": not missing_on_disk,
        "all_package_files_are_manifested": not unmanifested,
        "all_manifest_entries_are_tracked": not untracked,
        "all_hashes_match": not hash_mismatches,
        "raw_text_outputs_present_and_tracked": len(raw_text_outputs) == 176
        and all(relative in tracked_package_paths for relative in raw_text_outputs),
    }
    passed = all(checks.values())
    payload = {
        "schema": "periodic-gxtb-seidler-package-tracking-v1",
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "manifest_entry_count": len(entries),
        "tracked_package_file_count_excluding_manifest": len(tracked_package_paths),
        "raw_text_output_count": len(raw_text_outputs),
        "missing_on_disk": missing_on_disk,
        "unmanifested": unmanifested,
        "untracked": untracked,
        "hash_mismatches": hash_mismatches,
        "interpretation": (
            "Every file named by the Seidler recalculation-package hash manifest, "
            "including all raw CP2K and direct-CLI text outputs, exists, matches its "
            "recorded digest, and is present in the Git index."
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
