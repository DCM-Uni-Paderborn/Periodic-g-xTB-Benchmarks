#!/usr/bin/env python3
"""Move one atom by a complete lattice vector without changing the crystal."""

from __future__ import annotations

import sys
from pathlib import Path


source = Path(sys.argv[1])
target = Path(sys.argv[2])
lines = source.read_text().splitlines()
scale = float(lines[1].split()[0])
lattice_a = [scale * float(value) for value in lines[2].split()[:3]]
coordinate_header = next(
    index for index, line in enumerate(lines) if line.strip().lower().startswith(("cart", "direct"))
)
assert lines[coordinate_header].strip().lower().startswith("cart")
atom_line = coordinate_header + 1
fields = lines[atom_line].split()
xyz = [float(fields[index]) + lattice_a[index] for index in range(3)]
lines[atom_line] = " ".join(f"{value:.15f}" for value in xyz) + (
    " " + " ".join(fields[3:]) if len(fields) > 3 else ""
)
target.write_text("\n".join(lines) + "\n")
