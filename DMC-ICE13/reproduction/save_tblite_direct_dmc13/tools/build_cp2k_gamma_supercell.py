#!/usr/bin/env python3
"""Convert a Cartesian or direct POSCAR into a CP2K g-XTB Gamma input."""

from __future__ import annotations

import argparse
from pathlib import Path


def vector(values: list[str]) -> tuple[float, float, float]:
    if len(values) < 3:
        raise ValueError("expected a three-component vector")
    return tuple(float(value) for value in values[:3])  # type: ignore[return-value]


def cartesian(
    fractional: tuple[float, float, float],
    cell: tuple[tuple[float, float, float], ...],
) -> tuple[float, float, float]:
    return tuple(
        sum(fractional[row] * cell[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def read_poscar(
    path: Path,
) -> tuple[tuple[tuple[float, float, float], ...], list[tuple[str, tuple[float, float, float]]]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 9:
        raise ValueError(f"incomplete POSCAR: {path}")
    scale = float(lines[1])
    if scale <= 0.0:
        raise ValueError("only positive POSCAR scale factors are supported")
    cell = tuple(
        tuple(scale * component for component in vector(lines[index].split()))
        for index in range(2, 5)
    )
    symbols = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(symbols) != len(counts) or any(count <= 0 for count in counts):
        raise ValueError("invalid POSCAR symbol/count header")
    coordinate_line = 7
    if lines[coordinate_line].lower().startswith("s"):
        coordinate_line += 1
    mode = lines[coordinate_line].lower()
    direct = mode.startswith("d")
    if not direct and not mode.startswith(("c", "k")):
        raise ValueError(f"unsupported POSCAR coordinate mode: {lines[coordinate_line]}")
    first = coordinate_line + 1
    natom = sum(counts)
    if len(lines) < first + natom:
        raise ValueError("POSCAR contains fewer coordinates than declared")
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    position = first
    for symbol, count in zip(symbols, counts, strict=True):
        for _ in range(count):
            raw = vector(lines[position].split())
            xyz = cartesian(raw, cell) if direct else tuple(scale * value for value in raw)
            atoms.append((symbol, xyz))
            position += 1
    return cell, atoms


def cp2k_input(
    project: str,
    cell: tuple[tuple[float, float, float], ...],
    atoms: list[tuple[str, tuple[float, float, float]]],
    accuracy: float,
    eps_scf: float,
) -> str:
    cell_lines = "\n".join(
        f"      {label} {values[0]:.15f} {values[1]:.15f} {values[2]:.15f}"
        for label, values in zip(("A", "B", "C"), cell, strict=True)
    )
    coordinate_lines = "\n".join(
        f"      {symbol:<2s} {xyz[0]:.15f} {xyz[1]:.15f} {xyz[2]:.15f}"
        for symbol, xyz in atoms
    )
    return f"""# Generated from the archived Cartesian BvK POSCAR.
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT {project}
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &QS
      EPS_DEFAULT 1.0E-12
      METHOD xTB
      &XTB
        GFN_TYPE TBLITE
        SCC_MIXER TBLITE
        &TBLITE
          ACCURACY {accuracy:.12g}
          METHOD GXTB
        &END TBLITE
        &TBLITE_MIXER
          ITERATIONS 300
        &END TBLITE_MIXER
      &END XTB
    &END QS
    &SCF
      EPS_SCF {eps_scf:.12g}
      MAX_SCF 300
      SCF_GUESS MOPAC
      &MIXING
        ALPHA 0.2
        METHOD DIRECT_P_MIXING
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      PERIODIC XYZ
{cell_lines}
    &END CELL
    &COORD
{coordinate_lines}
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("poscar", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--project", default="gxtb_bvk_gamma_supercell")
    parser.add_argument("--accuracy", type=float, default=0.1)
    parser.add_argument("--eps-scf", type=float, default=1.0e-9)
    args = parser.parse_args()
    cell, atoms = read_poscar(args.poscar)
    args.output.write_text(
        cp2k_input(args.project, cell, atoms, args.accuracy, args.eps_scf),
        encoding="utf-8",
    )
    print(f"wrote {args.output} with {len(atoms)} atoms")


if __name__ == "__main__":
    main()
