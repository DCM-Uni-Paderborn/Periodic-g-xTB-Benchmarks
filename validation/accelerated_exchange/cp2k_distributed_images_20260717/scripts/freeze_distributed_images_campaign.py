#!/usr/bin/env python3
"""Freeze a deterministic SHA-256 inventory for the complete campaign."""

from __future__ import annotations

import hashlib
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    paths = sorted(
        path
        for path in CAMPAIGN.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS"
        and path.name != "sha256_check.log"
        and "__pycache__" not in path.parts
    )
    (CAMPAIGN / "SHA256SUMS").write_text(
        "".join(f"{digest(path)}  {path.relative_to(CAMPAIGN)}\n" for path in paths)
    )


if __name__ == "__main__":
    main()
