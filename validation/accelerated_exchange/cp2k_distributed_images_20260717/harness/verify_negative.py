#!/usr/bin/env python3
"""Verify frozen distributed-image negative-gate evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = json.loads((ROOT / "negative_cases.json").read_text())["cases"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    binary_hashes = set()
    for case in CASES:
        run_dir = ROOT / "negative_runs" / case["name"]
        metadata_path = run_dir / "run.json"
        output_path = run_dir / "cp2k.out"
        stderr_path = run_dir / "cp2k.err"
        for path in (metadata_path, output_path, stderr_path):
            if not path.is_file():
                raise RuntimeError(f"missing negative evidence: {path}")
        metadata = json.loads(metadata_path.read_text())
        combined = output_path.read_text(errors="replace") + "\n" + stderr_path.read_text(errors="replace")
        if metadata.get("case") != case["name"] or metadata.get("injection") != case["injection"]:
            raise RuntimeError(f"negative metadata mismatch: {case['name']}")
        input_path = ROOT / "inputs" / case["input"]
        if (metadata.get("ranks") != case["ranks"] or
                metadata.get("expected_diagnostic") != case["diagnostic"] or
                metadata.get("input") != case["input"] or
                metadata.get("input_sha256") != sha256(input_path)):
            raise RuntimeError(f"negative input/rank contract mismatch: {case['name']}")
        if metadata.get("returncode") == 0 or metadata.get("timed_out") is not False:
            raise RuntimeError(f"negative termination gate failed: {case['name']}")
        if metadata.get("output_sha256") != sha256(output_path) or metadata.get("stderr_sha256") != sha256(stderr_path):
            raise RuntimeError(f"negative raw-output hash mismatch: {case['name']}")
        if case["diagnostic"] not in combined or "PROGRAM ENDED" in combined:
            raise RuntimeError(f"negative diagnostic gate failed: {case['name']}")
        if (metadata.get("diagnostic_count") != combined.count(case["diagnostic"]) or
                metadata.get("diagnostic_count", 0) < 1):
            raise RuntimeError(f"negative diagnostic count missing: {case['name']}")
        binary_hashes.add(metadata.get("cp2k_sha256"))
    if len(binary_hashes) != 1 or None in binary_hashes:
        raise RuntimeError("negative cases did not use one frozen executable")
    print(f"PASS: {len(CASES)}/{len(CASES)} distributed-image negative gates")


if __name__ == "__main__":
    main()
