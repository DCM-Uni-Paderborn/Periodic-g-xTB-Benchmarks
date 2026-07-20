#!/usr/bin/env python3
"""Create the Part-I DMC-ICE13 uniform/adaptive TikZ figure from CSV data."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected exactly one row in {path}, found {len(rows)}")
    return rows[0]


def read_fixed_mesh_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"fixed-mesh table is empty: {path}")

    previous_mesh = 0
    for row in rows:
        mesh = int(row["mesh_n"])
        phase_count = int(row["phase_count"])
        mae = float(row["mae_kj_mol_per_water"])
        if mesh <= previous_mesh:
            raise ValueError("fixed meshes must be strictly increasing")
        if phase_count != 12:
            raise ValueError(f"mesh {mesh} does not contain twelve phase comparisons")
        if not math.isfinite(mae) or mae <= 0.0:
            raise ValueError(f"mesh {mesh} has an invalid MAE")
        expected_label = "Gamma" if mesh == 1 else f"{mesh}x{mesh}x{mesh}"
        if row["mesh_label"] != expected_label:
            raise ValueError(f"mesh {mesh} has label {row['mesh_label']!r}, expected {expected_label!r}")
        previous_mesh = mesh
    if int(rows[0]["mesh_n"]) != 1:
        raise ValueError("the uniform sequence must begin at Gamma")
    return rows


def tex_float(value: float) -> str:
    return f"{value:.10g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixed_mesh_csv", type=Path)
    parser.add_argument("adaptive_statistics_csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    fixed_rows = read_fixed_mesh_rows(args.fixed_mesh_csv)
    adaptive = read_single_row(args.adaptive_statistics_csv)
    if int(adaptive["phase_count"]) != 12:
        raise ValueError("adaptive statistics do not contain twelve phase comparisons")
    converged = int(adaptive["converged_phase_count"])
    if not 0 <= converged <= 12:
        raise ValueError("invalid adaptive converged-phase count")
    largest_mesh = int(adaptive["largest_mesh_n"])
    adaptive_mae = float(adaptive["mae_kj_mol_per_water"])
    if largest_mesh < 1 or not math.isfinite(adaptive_mae) or adaptive_mae <= 0.0:
        raise ValueError("invalid adaptive mesh or MAE")
    declared_final = adaptive["final_result"].strip().lower()
    if declared_final not in {"true", "false"}:
        raise ValueError("final_result must be true or false")
    if (declared_final == "true") != (converged == 12):
        raise ValueError("final_result and converged_phase_count disagree")

    fixed_maes = [float(row["mae_kj_mol_per_water"]) for row in fixed_rows]
    adaptive_x = len(fixed_rows) + 1
    ticks = ",".join(str(index) for index in range(1, adaptive_x + 1))
    labels = [r"$\Gamma$"] + [rf"${int(row['mesh_n'])}^3$" for row in fixed_rows[1:]]
    labels.append(rf"adaptive $\leq{largest_mesh}^3$")
    ticklabels = ",".join(labels)
    coordinates = "\n".join(
        f"  ({index},{tex_float(mae)})"
        for index, mae in enumerate(fixed_maes, start=1)
    )

    source = rf"""\begin{{tikzpicture}}
\begin{{semilogyaxis}}[
  width=\columnwidth,
  height=0.76\columnwidth,
  xmin=0.65,
  xmax={adaptive_x + 0.35:.2f},
  ymin=1.3,
  ymax=220,
  xtick={{{ticks}}},
  xticklabels={{{ticklabels}}},
  x tick label style={{rotate=45,anchor=north east,font=\scriptsize}},
  y tick label style={{font=\small}},
  label style={{font=\small}},
  xlabel={{k-point mesh}},
  ylabel={{MAE / kJ mol$^{{-1}}$ H$_2$O$^{{-1}}$}},
  grid=both,
  xmajorgrids=false,
  major grid style={{gray!35}},
  minor grid style={{gray!18}},
  axis line style={{black}},
  tick align=outside,
  clip=false,
]
\addplot[
  blue!70!black,
  thick,
  mark=*,
  mark size=2pt,
] coordinates {{
{coordinates}
}};
\node[anchor=west,font=\small] at (axis cs:1.18,{tex_float(fixed_maes[0])}) {{{fixed_maes[0]:.1f}}};
\addplot[
  red!70!black,
  only marks,
  mark=diamond*,
  mark size=3pt,
] coordinates {{({adaptive_x},{tex_float(adaptive_mae)})}};
\node[anchor=south east,font=\small] at (axis cs:{adaptive_x}.0,{tex_float(adaptive_mae * 1.12)}) {{{adaptive_mae:.3f}}};
\end{{semilogyaxis}}
\end{{tikzpicture}}
"""

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(source, encoding="utf-8")
    temporary.replace(args.output)
    print(
        f"fixed_meshes={len(fixed_rows)} adaptive_mesh={largest_mesh} "
        f"converged={converged}/12 mae={adaptive_mae:.12f} final={declared_final}"
    )


if __name__ == "__main__":
    main()
