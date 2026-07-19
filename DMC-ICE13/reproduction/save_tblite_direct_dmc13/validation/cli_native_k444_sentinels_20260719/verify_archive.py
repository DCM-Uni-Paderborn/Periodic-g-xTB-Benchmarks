#!/usr/bin/env python3
"""Verify the archived 4 x 4 x 4 direct-CLI/native energy sentinels."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
PHASES = ("Ih", "XI", "XIV", "II", "VII")
DIRECT_BINARY = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
NATIVE_BINARY = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_manifest() -> tuple[int, int]:
    manifest = ROOT / "SHA256SUMS"
    entries = 0
    checked_bytes = 0
    for line_number, raw in enumerate(manifest.read_text().splitlines(), 1):
        if not raw.strip():
            continue
        fields = raw.split(None, 1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise AssertionError(f"malformed SHA256SUMS entry at line {line_number}")
        expected, relative_text = fields
        relative = Path(relative_text.lstrip("*").strip())
        if relative.is_absolute() or ".." in relative.parts:
            raise AssertionError(f"nonportable manifest path: {relative}")
        target = ROOT / relative
        if not target.is_file():
            raise AssertionError(f"missing manifest target: {relative}")
        actual = sha256(target)
        if actual != expected.lower():
            raise AssertionError(f"SHA-256 mismatch: {relative}")
        entries += 1
        checked_bytes += target.stat().st_size
    return entries, checked_bytes


def main() -> None:
    verifier = PACKAGE / "tools/verify_cli_native_sentinels.py"
    command = [sys.executable, str(verifier), "--mesh", "4"]
    for phase in PHASES:
        command.extend(("--phase", phase))
    command.extend(
        (
            "--direct-root",
            str(ROOT / "raw/cli"),
            "--native-root",
            str(ROOT / "raw/native"),
            "--direct-input-root",
            str(PACKAGE / "structures/k444"),
            "--native-input-root",
            str(ROOT / "inputs/native"),
            "--direct-binary-sha256",
            DIRECT_BINARY,
            "--native-binary-sha256",
            NATIVE_BINARY,
            "--absolute-tolerance-Ha",
            "2e-7",
            "--relative-tolerance-kJ-mol-water",
            "5e-5",
        )
    )
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise AssertionError(completed.stderr or completed.stdout)
    recomputed = json.loads(completed.stdout)
    archived = json.loads((ROOT / "report.json").read_text())
    if recomputed != archived:
        raise AssertionError("archived comparison differs from recomputed result")

    for phase in PHASES:
        process_text = (ROOT / "raw/cli" / phase / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        if (
            "total energy" not in process_text
            or "JSON dump of results written" not in process_text
        ):
            raise AssertionError(f"incomplete direct output marker: {phase}")

    entries, checked_bytes = verify_manifest()
    print(
        json.dumps(
            {
                "status": "PASS",
                "mesh": 4,
                "phases": list(PHASES),
                "max_abs_native_minus_direct_Ha": archived["summary"][
                    "max_abs_native_minus_direct_Ha"
                ],
                "max_abs_relative_native_minus_direct_kJ_mol_per_water": archived[
                    "summary"
                ]["max_abs_relative_native_minus_direct_kJ_mol_per_water"],
                "manifest_entries": entries,
                "manifest_checked_bytes": checked_bytes,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
