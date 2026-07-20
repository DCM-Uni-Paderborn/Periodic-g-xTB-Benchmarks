#!/usr/bin/env python3
"""Generate and verify an explicit CP2K Gamma-supercell g-xTB input."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


TOLERANCE_ANGSTROM = 5.0e-12


def read_poscar(path: Path) -> tuple[list[list[float]], list[str], list[list[float]]]:
    lines = [line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 8:
        raise RuntimeError(f"POSCAR is too short: {path}")
    scale = float(lines[1].split()[0])
    if scale <= 0.0:
        raise RuntimeError("only a positive scalar POSCAR scale is supported")
    cell = [
        [scale * float(value) for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]
    elements = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(elements) != len(counts) or any(value <= 0 for value in counts):
        raise RuntimeError("invalid POSCAR element/count rows")
    coordinate_header = 7
    if lines[coordinate_header].strip().lower().startswith("s"):
        coordinate_header += 1
    mode = lines[coordinate_header].strip().lower()
    if not mode.startswith(("c", "k")):
        raise RuntimeError("the explicit-BvK POSCAR must use Cartesian coordinates")
    atom_count = sum(counts)
    coordinate_lines = lines[coordinate_header + 1 : coordinate_header + 1 + atom_count]
    if len(coordinate_lines) != atom_count:
        raise RuntimeError("POSCAR atom count and coordinate rows differ")
    coordinates = [
        [scale * float(value) for value in line.split()[:3]]
        for line in coordinate_lines
    ]
    symbols = [element for element, count in zip(elements, counts) for _ in range(count)]
    return cell, symbols, coordinates


def render_input(
    phase: str,
    cell: list[list[float]],
    symbols: list[str],
    coordinates: list[list[float]],
) -> str:
    lines = [
        "# Generated from the archived Cartesian BvK POSCAR.",
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT ice_{phase}_GXTB_k222_gamma_supercell",
        "  RUN_TYPE ENERGY",
        "&END GLOBAL",
        "",
        "&FORCE_EVAL",
        "  METHOD QUICKSTEP",
        "  &DFT",
        "    &QS",
        "      EPS_DEFAULT 1.0E-12",
        "      METHOD xTB",
        "      &XTB",
        "        GFN_TYPE TBLITE",
        "        SCC_MIXER TBLITE",
        "        &TBLITE",
        "          ACCURACY 0.1",
        "          METHOD GXTB",
        "        &END TBLITE",
        "        &TBLITE_MIXER",
        "          ITERATIONS 300",
        "        &END TBLITE_MIXER",
        "      &END XTB",
        "    &END QS",
        "    &SCF",
        "      EPS_SCF 1e-09",
        "      MAX_SCF 300",
        "      SCF_GUESS MOPAC",
        "      &MIXING",
        "        ALPHA 0.2",
        "        METHOD DIRECT_P_MIXING",
        "      &END MIXING",
        "      &PRINT",
        "        &RESTART OFF",
        "        &END RESTART",
        "      &END PRINT",
        "    &END SCF",
        "  &END DFT",
        "  &SUBSYS",
        "    &CELL",
        "      PERIODIC XYZ",
    ]
    for label, vector in zip(("A", "B", "C"), cell):
        lines.append(
            f"      {label} {vector[0]:.15f} {vector[1]:.15f} {vector[2]:.15f}"
        )
    lines.extend(("    &END CELL", "    &COORD"))
    for symbol, coordinate in zip(symbols, coordinates):
        lines.append(
            f"      {symbol}  {coordinate[0]:.15f} "
            f"{coordinate[1]:.15f} {coordinate[2]:.15f}"
        )
    lines.extend(("    &END COORD", "  &END SUBSYS", "&END FORCE_EVAL", ""))
    return "\n".join(lines)


def parse_generated(path: Path) -> tuple[list[list[float]], list[str], list[list[float]]]:
    text = path.read_text(encoding="utf-8")
    cell = []
    for label in ("A", "B", "C"):
        match = re.search(
            rf"^\s*{label}\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)$",
            text,
            flags=re.MULTILINE,
        )
        if match is None:
            raise RuntimeError(f"cannot parse generated cell vector {label}")
        cell.append([float(value) for value in match.groups()])
    section = re.search(
        r"^\s*&COORD\s*$\n(.*?)^\s*&END COORD\s*$",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if section is None:
        raise RuntimeError("cannot parse generated coordinate section")
    symbols = []
    coordinates = []
    for raw in section.group(1).splitlines():
        fields = raw.split()
        if not fields:
            continue
        if len(fields) != 4:
            raise RuntimeError(f"invalid generated coordinate row: {raw}")
        symbols.append(fields[0])
        coordinates.append([float(value) for value in fields[1:]])
    return cell, symbols, coordinates


def maximum_difference(left: list[list[float]], right: list[list[float]]) -> float:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise RuntimeError("matrix dimensions differ")
    return max(
        abs(a - b)
        for left_row, right_row in zip(left, right)
        for a, b in zip(left_row, right_row)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("poscar", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("verification", type=Path)
    parser.add_argument("--phase", required=True)
    args = parser.parse_args()
    cell, symbols, coordinates = read_poscar(args.poscar)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        render_input(args.phase, cell, symbols, coordinates), encoding="utf-8"
    )
    generated_cell, generated_symbols, generated_coordinates = parse_generated(
        args.output
    )
    cell_delta = maximum_difference(cell, generated_cell)
    coordinate_delta = maximum_difference(coordinates, generated_coordinates)
    text = args.output.read_text(encoding="utf-8")
    result = {
        "atom_count": len(symbols),
        "cell_max_abs_delta_angstrom": cell_delta,
        "coordinate_max_abs_delta_angstrom": coordinate_delta,
        "elements_and_order_equal": symbols == generated_symbols,
        "explicit_kpoints_section_present": bool(
            re.search(r"^\s*&KPOINTS\b", text, flags=re.MULTILINE | re.IGNORECASE)
        ),
        "periodic_xyz": bool(
            re.search(r"^\s*PERIODIC\s+XYZ\s*$", text, flags=re.MULTILINE)
        ),
        "tolerance_angstrom": TOLERANCE_ANGSTROM,
    }
    result["status"] = "PASS" if (
        math.isfinite(cell_delta)
        and math.isfinite(coordinate_delta)
        and cell_delta <= TOLERANCE_ANGSTROM
        and coordinate_delta <= TOLERANCE_ANGSTROM
        and result["elements_and_order_equal"]
        and not result["explicit_kpoints_section_present"]
        and result["periodic_xyz"]
    ) else "FAIL"
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.verification.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
