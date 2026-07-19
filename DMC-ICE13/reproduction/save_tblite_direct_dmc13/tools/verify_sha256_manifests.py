#!/usr/bin/env python3
"""Verify every portable SHA256SUMS file in the reproduction archive."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()
    root = args.root.resolve()
    manifests = sorted(root.rglob("SHA256SUMS"))
    if not manifests:
        raise ValueError(f"no SHA256SUMS files below {root}")

    failures: list[str] = []
    entry_count = 0
    for manifest in manifests:
        for line_number, raw_line in enumerate(
            manifest.read_text(encoding="utf-8", errors="replace").splitlines(), 1
        ):
            if not raw_line.strip():
                continue
            fields = raw_line.split(None, 1)
            if len(fields) != 2 or len(fields[0]) != 64:
                failures.append(f"{manifest}:{line_number}: malformed entry")
                continue
            expected, recorded_path = fields
            relative = Path(recorded_path.lstrip("*").strip())
            entry_count += 1
            if relative.is_absolute():
                failures.append(
                    f"{manifest}:{line_number}: absolute path {relative}"
                )
                continue
            target = manifest.parent / relative
            if not target.is_file():
                failures.append(f"{manifest}:{line_number}: missing {relative}")
                continue
            actual = sha256(target)
            if actual != expected.lower():
                failures.append(
                    f"{manifest}:{line_number}: hash mismatch {relative}"
                )

    print(
        f"manifests={len(manifests)} entries={entry_count} "
        f"failures={len(failures)}"
    )
    for failure in failures:
        print(failure, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"fatal: {error}", file=sys.stderr)
        raise SystemExit(2) from error
