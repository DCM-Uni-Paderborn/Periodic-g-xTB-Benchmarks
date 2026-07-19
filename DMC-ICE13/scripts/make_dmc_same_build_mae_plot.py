#!/usr/bin/env python3
"""Plot the fail-closed same-build fixed-mesh DMC-ICE13 MAE series."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def read_rows(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle))
    if not raw:
        raise ValueError("fixed-mesh MAE table is empty")
    required = {"mesh_n", "mesh_label", "phase_count", "mae_kj_mol_per_water"}
    if not required.issubset(raw[0]):
        raise ValueError(f"fixed-mesh MAE table lacks columns: {sorted(required - set(raw[0]))}")
    rows: list[dict[str, object]] = []
    previous = 0
    for record in raw:
        mesh = int(record["mesh_n"])
        phase_count = int(record["phase_count"])
        mae = float(record["mae_kj_mol_per_water"])
        if mesh <= previous:
            raise ValueError("mesh sequence must be strictly increasing")
        if phase_count != 12:
            raise ValueError(f"mesh {mesh} does not cover twelve DMC comparisons")
        if mae <= 0.0:
            raise ValueError(f"mesh {mesh} has a nonpositive MAE")
        expected_label = "Gamma" if mesh == 1 else f"{mesh}x{mesh}x{mesh}"
        if record["mesh_label"] != expected_label:
            raise ValueError(f"mesh {mesh} has an inconsistent label")
        rows.append({"mesh_n": mesh, "mesh_label": expected_label, "mae": mae})
        previous = mesh
    if rows[0]["mesh_n"] != 1:
        raise ValueError("the paper curve must start at Gamma")
    return rows


def save_figure(fig: mpl.figure.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        target = output.with_suffix(f".{suffix}")
        options: dict[str, Any] = {"bbox_inches": "tight"}
        if suffix == "pdf":
            options["metadata"] = {
                "CreationDate": None,
                "ModDate": None,
                "Creator": "Periodic-g-xTB-Benchmarks",
            }
        elif suffix == "svg":
            options["metadata"] = {
                "Date": "2026-07-19",
                "Creator": "Periodic-g-xTB-Benchmarks",
            }
        else:
            options["dpi"] = 600
        fig.savefig(target, **options)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixed_mesh_csv", type=Path)
    parser.add_argument("output_stem", type=Path)
    args = parser.parse_args()
    rows = read_rows(args.fixed_mesh_csv)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.8,
            "lines.solid_capstyle": "round",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "dmc-ice13-gxtb-same-build-v1",
        }
    )
    meshes = [int(row["mesh_n"]) for row in rows]
    values = [float(row["mae"]) for row in rows]
    labels = [r"$\Gamma$"] + [rf"${mesh}^3$" for mesh in meshes[1:]]

    fig, ax = plt.subplots(figsize=(3.37, 2.55), constrained_layout=True)
    ax.plot(meshes, values, color="#0072B2", marker="o", markersize=4.0, linewidth=1.6)
    ax.set_yscale("log")
    ax.set_xticks(meshes, labels)
    ax.set_xlabel("k-point mesh")
    ax.set_ylabel(r"DMC-ICE13 MAE (kJ mol$^{-1}$ H$_2$O$^{-1}$)")
    ax.grid(axis="y", which="major", color="0.86", linewidth=0.6)
    ax.grid(axis="y", which="minor", color="0.93", linewidth=0.4)
    ax.spines[["top", "right"]].set_visible(False)
    for mesh, value in ((meshes[0], values[0]), (meshes[-1], values[-1])):
        ax.annotate(
            f"{value:.3g}",
            xy=(mesh, value),
            xytext=(4 if mesh == meshes[0] else -4, 7),
            textcoords="offset points",
            ha="left" if mesh == meshes[0] else "right",
            va="bottom",
        )
    save_figure(fig, args.output_stem)
    plt.close(fig)


if __name__ == "__main__":
    main()
