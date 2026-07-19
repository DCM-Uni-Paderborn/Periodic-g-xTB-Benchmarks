#!/usr/bin/env python3
"""Create a cubic native-k CP2K input by changing only its regular mesh."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


SCHEME_RE = re.compile(
    r"^(?P<indent>\s*)SCHEME\s+MACDONALD\s+"
    r"(?P<n1>\d+)\s+(?P<n2>\d+)\s+(?P<n3>\d+)"
    r"(?P<tail>\s+[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s*)$",
    flags=re.IGNORECASE,
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("mesh", type=int)
    parser.add_argument("--provenance", type=Path)
    args = parser.parse_args()
    if args.mesh <= 0:
        parser.error("mesh must be positive")
    if args.target.exists():
        parser.error(f"target already exists: {args.target}")

    source_text = args.source.read_text(encoding="utf-8")
    lines = source_text.splitlines(keepends=True)
    in_kpoints = False
    replacements = 0
    source_mesh = 0
    for index, raw_line in enumerate(lines):
        stripped = raw_line.strip().upper()
        if stripped == "&KPOINTS":
            if in_kpoints:
                raise ValueError("nested KPOINTS section")
            in_kpoints = True
            continue
        if in_kpoints and stripped.startswith("&END"):
            in_kpoints = False
            continue
        if not in_kpoints:
            continue
        line = raw_line.rstrip("\r\n")
        match = SCHEME_RE.fullmatch(line)
        if not match:
            continue
        mesh_values = tuple(int(match.group(name)) for name in ("n1", "n2", "n3"))
        if len(set(mesh_values)) != 1:
            raise ValueError(f"source mesh is not cubic: {mesh_values}")
        source_mesh = mesh_values[0]
        newline = "\r\n" if raw_line.endswith("\r\n") else "\n" if raw_line.endswith("\n") else ""
        lines[index] = (
            f"{match.group('indent')}SCHEME MACDONALD "
            f"{args.mesh} {args.mesh} {args.mesh}{match.group('tail')}{newline}"
        )
        replacements += 1

    if replacements != 1:
        raise ValueError(f"expected one regular MACDONALD scheme, found {replacements}")
    target_text = "".join(lines)
    if target_text == source_text:
        raise ValueError("target mesh is identical to the source mesh")
    args.target.parent.mkdir(parents=True, exist_ok=True)
    args.target.write_text(target_text, encoding="utf-8")

    changed = [
        index + 1
        for index, (left, right) in enumerate(zip(source_text.splitlines(), target_text.splitlines()))
        if left != right
    ]
    if len(changed) != 1:
        raise AssertionError(f"mesh rewrite changed {len(changed)} lines instead of one")
    result = {
        "changed_line": changed[0],
        "source": str(args.source),
        "source_mesh": source_mesh,
        "source_sha256": digest(args.source),
        "status": "PASS",
        "target": str(args.target),
        "target_mesh": args.mesh,
        "target_sha256": digest(args.target),
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.provenance:
        if args.provenance.exists():
            parser.error(f"provenance target already exists: {args.provenance}")
        args.provenance.parent.mkdir(parents=True, exist_ok=True)
        args.provenance.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
