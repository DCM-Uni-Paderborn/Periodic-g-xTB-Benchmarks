#!/usr/bin/env python3
"""Create deterministic SHA-256 inventory after a completed campaign."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    paths = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and "__pycache__" not in path.parts
    )
    (ROOT / "SHA256SUMS").write_text("".join(
        f"{digest(path)}  {path.relative_to(ROOT)}\n" for path in paths
    ))


if __name__ == "__main__":
    main()
