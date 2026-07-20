#!/usr/bin/env python3
"""Render the Part-I DMC-ICE13 uniform and phase-wise progress figure."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty table: {path}")
    return rows


def save_figure(figure: mpl.figure.Figure, output: Path) -> None:
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
                "Date": "2026-07-20",
                "Creator": "Periodic-g-xTB-Benchmarks",
            }
        else:
            options["dpi"] = 600
        figure.savefig(target, **options)
        if suffix == "svg":
            normalized = "\n".join(
                line.rstrip()
                for line in target.read_text(encoding="utf-8").splitlines()
            )
            target.write_text(normalized + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixed_mesh_csv", type=Path)
    parser.add_argument("progress_csv", type=Path)
    parser.add_argument("output_stem", type=Path)
    args = parser.parse_args()

    fixed_rows = read_rows(args.fixed_mesh_csv)
    progress_rows = read_rows(args.progress_csv)
    fixed_meshes = [int(row["mesh_n"]) for row in fixed_rows]
    fixed_values = [float(row["mae_kj_mol_per_water"]) for row in fixed_rows]
    progress_meshes = [int(row["mesh_limit_n"]) for row in progress_rows]
    progress_values = [
        float(row["mae_kj_mol_per_water"]) for row in progress_rows
    ]
    all_meshes = fixed_meshes + progress_meshes
    if all_meshes != list(range(1, all_meshes[-1] + 1)):
        raise ValueError("the combined figure mesh sequence must be contiguous")
    if any(int(row["phase_count"]) != 12 for row in fixed_rows + progress_rows):
        raise ValueError("every plotted statistic must contain twelve phases")
    if any(row.get("qualification") != "PASS" for row in progress_rows):
        raise ValueError("every phase-wise progress point must be qualified")

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
            "svg.hashsalt": "dmc-ice13-gxtb-phasewise-progress-v1",
        }
    )

    figure, axis = plt.subplots(figsize=(3.37, 2.55), constrained_layout=True)
    axis.plot(
        fixed_meshes,
        fixed_values,
        color="#0072B2",
        marker="o",
        markersize=4.0,
        linewidth=1.6,
    )
    latest_progress = progress_rows[-1]
    provisional = int(latest_progress["converged_phase_count"]) < int(
        latest_progress["phase_count"]
    )
    blue_progress_count = (
        len(progress_meshes) - 1 if provisional else len(progress_meshes)
    )
    axis.plot(
        [fixed_meshes[-1], *progress_meshes[:blue_progress_count]],
        [fixed_values[-1], *progress_values[:blue_progress_count]],
        color="#0072B2",
        marker="D",
        markersize=4.0,
        linewidth=1.6,
    )
    if provisional:
        previous_mesh = (
            progress_meshes[-2] if len(progress_meshes) > 1 else fixed_meshes[-1]
        )
        previous_value = (
            progress_values[-2] if len(progress_values) > 1 else fixed_values[-1]
        )
        axis.plot(
            [previous_mesh, progress_meshes[-1]],
            [previous_value, progress_values[-1]],
            color="#D55E00",
            marker="D",
            markevery=[1],
            markersize=4.0,
            linewidth=1.4,
            linestyle="--",
        )
    axis.set_yscale("log")
    axis.set_xlim(0.65, all_meshes[-1] + 0.35)
    axis.set_ylim(1.3, 220.0)
    tick_labels = [r"$\Gamma$"] + [rf"${mesh}^3$" for mesh in all_meshes[1:]]
    if provisional:
        tick_labels[-1] = rf"$\leq {all_meshes[-1]}^3$"
    axis.set_xticks(
        all_meshes,
        tick_labels,
        rotation=45,
        ha="right",
    )
    axis.set_xlabel("k-point mesh")
    axis.set_ylabel(r"MAE / kJ mol$^{-1}$ H$_2$O$^{-1}$")
    axis.grid(axis="y", which="major", color="0.86", linewidth=0.6)
    axis.grid(axis="y", which="minor", color="0.93", linewidth=0.4)
    axis.spines[["top", "right"]].set_visible(False)
    axis.annotate(
        f"{fixed_values[0]:.1f}",
        xy=(fixed_meshes[0], fixed_values[0]),
        xytext=(4, 5),
        textcoords="offset points",
        ha="left",
        va="bottom",
    )
    for mesh, value in zip(progress_meshes, progress_values, strict=True):
        axis.annotate(
            f"{value:.3f}",
            xy=(mesh, value),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
        )
    save_figure(figure, args.output_stem)
    plt.close(figure)


if __name__ == "__main__":
    main()
