#!/usr/bin/env python3
"""Compare final CP2K energy, force, and analytical-stress blocks."""

from __future__ import annotations

import re
import sys
from pathlib import Path


FLOAT = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"


def parse(path: Path):
    text = path.read_text(errors="replace")
    energies = [float(value) for value in re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )]
    forces = [
        tuple(float(value) for value in match)
        for match in re.findall(
            rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
            text,
            re.MULTILINE,
        )
    ]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    stress = []
    if stress_blocks:
        stress = [
            tuple(float(value) for value in match)
            for match in re.findall(
                rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
                stress_blocks[-1],
                re.MULTILINE,
            )
        ]
    if not energies or not forces or len(stress) != 3 or "PROGRAM ENDED" not in text:
        raise RuntimeError(f"incomplete CP2K result: {path}")
    natom = len(forces)
    return energies[-1], forces[-natom:], stress


def maximum_delta(left, right):
    return max(abs(x - y) for a, b in zip(left, right) for x, y in zip(a, b))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare.py CANDIDATE ORACLE")
    candidate = parse(Path(sys.argv[1]))
    oracle = parse(Path(sys.argv[2]))
    if len(candidate[1]) != len(oracle[1]):
        raise RuntimeError("force blocks have different atom counts")
    print(
        f"{abs(candidate[0] - oracle[0]):.16e}\t"
        f"{maximum_delta(candidate[1], oracle[1]):.16e}\t"
        f"{maximum_delta(candidate[2], oracle[2]):.16e}"
    )
