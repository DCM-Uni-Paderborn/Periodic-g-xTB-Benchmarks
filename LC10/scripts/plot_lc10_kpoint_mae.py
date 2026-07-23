#!/usr/bin/env python3
"""Plot LC10 lattice-constant and cohesive-energy k-point convergence."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    data_path = root / "data/lc10_gxtb_uniform_mesh_mae.csv"
    with data_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    meshes = [int(row["mesh_n"]) for row in rows]
    lattice_mae = [float(row["lattice_mae_A"]) for row in rows]
    cohesive_mae = [float(row["cohesive_mae_eV_per_atom"]) for row in rows]
    tick_labels = [r"$\Gamma$"] + [rf"${mesh}^3$" for mesh in meshes[1:]]

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
            "svg.hashsalt": "lc10-gxtb-kpoint-mae-v1",
        }
    )

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(3.37, 4.25),
        sharex=True,
        constrained_layout=True,
    )
    color = "#0072B2"
    for ax, values in zip(axes, (lattice_mae, cohesive_mae), strict=True):
        ax.plot(
            meshes,
            values,
            color=color,
            marker="o",
            markersize=4.0,
            linewidth=1.6,
        )
        ax.plot(
            meshes[0],
            values[0],
            color=color,
            marker="o",
            markerfacecolor="white",
            markersize=4.0,
            linestyle="none",
        )
        ax.grid(axis="y", color="0.88", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(0.7, 7.3)

    axes[0].set_ylabel(r"Lattice-constant MAE / \AA")
    axes[1].set_ylabel(r"Cohesive-energy MAE / eV atom$^{-1}$")
    axes[1].set_xlabel("k-point mesh")
    axes[1].set_xticks(meshes, tick_labels)
    axes[0].text(
        0.02,
        0.95,
        "(a)",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontweight="bold",
    )
    axes[1].text(
        0.02,
        0.95,
        "(b)",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontweight="bold",
    )

    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    stem = figures / "lc10_gxtb_uniform_mesh_mae"
    for suffix in ("pdf", "svg", "png"):
        options: dict[str, object] = {"bbox_inches": "tight"}
        if suffix == "pdf":
            options["metadata"] = {"CreationDate": None, "ModDate": None}
        elif suffix == "svg":
            options["metadata"] = {"Date": "2026-07-23"}
        else:
            options["dpi"] = 600
        fig.savefig(stem.with_suffix(f".{suffix}"), **options)
    plt.close(fig)


if __name__ == "__main__":
    main()
