#!/usr/bin/env python3
"""Fail closed if a CP2K partial-importer path bypasses the partial API."""

from __future__ import annotations

import re
import sys
from pathlib import Path


def subroutine(source: str, name: str) -> str:
    match = re.search(
        rf"(?ims)^\s*SUBROUTINE\s+{re.escape(name)}\b.*?"
        rf"^\s*END\s+SUBROUTINE\s+{re.escape(name)}\s*$",
        source,
    )
    if match is None:
        raise SystemExit(f"FAIL: missing subroutine {name}")
    return match.group(0)


def require(scope: str, needle: str, label: str) -> None:
    if re.search(rf"(?i)\bCALL\s+{re.escape(needle)}\b", scope) is None:
        raise SystemExit(f"FAIL: {label} does not call {needle}")


def forbid(scope: str, needle: str, label: str) -> None:
    if re.search(rf"(?i)\bCALL\s+{re.escape(needle)}\b", scope) is not None:
        raise SystemExit(f"FAIL: {label} bypasses partial API via {needle}")


if len(sys.argv) != 2:
    raise SystemExit("usage: check_partial_api.py /path/to/src/tblite_interface.F")

path = Path(sys.argv[1]).resolve()
text = path.read_text(encoding="utf-8")
forward_name = "tb_build_gxtb_kpoint_exchange_kgroup_partial_root"
reverse_name = "tb_gxtb_kpoint_exchange_gradient_kgroup_partial_root"
forward = subroutine(text, forward_name)
reverse = subroutine(text, reverse_name)

require(forward, "cp2k_exchange_partial_push", "forward partial importer")
forbid(forward, "cp2k_exchange_stream_push", "forward partial importer")
require(reverse, "cp2k_exchange_partial_reverse_push", "reverse partial importer")
forbid(reverse, "cp2k_exchange_stream_reverse_push", "reverse partial importer")

for scope, label in ((forward, "forward"), (reverse, "reverse")):
    require(scope, "cp2k_exchange_partial_set_image_range", label)

print("PASS: partial importers use only the ownership-tracking partial push APIs")
