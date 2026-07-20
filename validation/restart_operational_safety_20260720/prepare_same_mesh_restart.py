#!/usr/bin/env python3
"""Create a fail-closed CP2K input for a strict same-mesh k-point restart."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cp2k_file_name(value: str) -> str:
    # CP2K's WFN_RESTART_FILE_NAME path reaches a legacy file-open path where
    # quoted values are treated as a different (and therefore missing) file.
    # The production paths contain no whitespace; reject unsafe names instead
    # of silently emitting a restart that falls back to an atomic guess.
    if re.search(r"[\s\"']", value):
        raise SystemExit("restart path must not contain whitespace or quotes")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("original_input", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("restart_input", type=Path)
    args = parser.parse_args()

    original = args.original_input.resolve()
    checkpoint = args.checkpoint.resolve()
    output = args.restart_input.resolve()
    if not original.is_file() or not checkpoint.is_file():
        raise SystemExit("original input and checkpoint must both exist")
    if output == original or output.exists():
        raise SystemExit("restart input must be a new path")

    text = original.read_text(encoding="utf-8")
    if re.search(r"^\s*WFN_RESTART_FILE_NAME\b", text, re.MULTILINE):
        raise SystemExit("input already declares WFN_RESTART_FILE_NAME")

    dft_matches = list(re.finditer(r"^(\s*)&DFT\s*$", text, re.MULTILINE | re.IGNORECASE))
    guess_matches = list(
        re.finditer(r"^(\s*)SCF_GUESS\s+\S+\s*$", text, re.MULTILINE | re.IGNORECASE)
    )
    if len(dft_matches) != 1 or len(guess_matches) != 1:
        raise SystemExit("expected exactly one &DFT section and one SCF_GUESS keyword")

    dft = dft_matches[0]
    try:
        # The restart launcher executes CP2K in output.parent.  Prefer a short
        # relative name because CP2K's legacy input strings can truncate long
        # absolute WFN_RESTART_FILE_NAME paths before the file-open call.
        restart_name = str(checkpoint.relative_to(output.parent))
    except ValueError:
        restart_name = str(checkpoint)
    insertion = (
        f"\n{dft.group(1)}  WFN_RESTART_FILE_NAME "
        f"{cp2k_file_name(restart_name)}"
    )
    text = text[: dft.end()] + insertion + text[dft.end() :]
    text, replacements = re.subn(
        r"^(\s*)SCF_GUESS\s+\S+\s*$",
        r"\1SCF_GUESS RESTART",
        text,
        count=1,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    if replacements != 1:
        raise SystemExit("failed to replace SCF_GUESS")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=output.parent, delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, output)

    provenance = {
        "schema": "cp2k-strict-same-mesh-restart-input-v1",
        "original_input": str(original),
        "original_input_sha256": sha256(original),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256(checkpoint),
        "checkpoint_size_bytes": checkpoint.stat().st_size,
        "wfn_restart_file_name": restart_name,
        "restart_input": str(output),
        "restart_input_sha256": sha256(output),
        "scf_guess": "RESTART",
        "restart_mode": "strict_same_mesh",
    }
    provenance_path = output.with_suffix(output.suffix + ".provenance.json")
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
