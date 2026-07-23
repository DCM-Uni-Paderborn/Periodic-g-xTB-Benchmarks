#!/usr/bin/env python3
"""Render the DMC-ICE13 MAE ranking used in the Supporting Information."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
CURRENT_MODELS = {
    r"g-xTB (adaptive $\leq 8^3$)": 1.74098,
    r"GFN2-xTB (adaptive $\leq 4^3$)": 3.46143,
    r"GFN1-xTB (adaptive $\leq 4^3$)": 8.00642,
}


def published_maes(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {row["method"]: row for row in csv.DictReader(handle)}
    reference = rows.pop("DMC")
    dmc = {phase: float(reference[phase]) - float(reference["Ih"]) for phase in PHASES}
    result: dict[str, float] = {}
    for method, row in rows.items():
        errors = [
            (float(row[phase]) - float(row["Ih"])) - dmc[phase]
            for phase in PHASES
        ]
        result[method] = sum(abs(value) for value in errors) / len(errors)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("published_energies", type=Path)
    parser.add_argument("output_stem", type=Path)
    args = parser.parse_args()

    values = published_maes(args.published_energies)
    values.update(CURRENT_MODELS)
    ranked = sorted(values.items(), key=lambda item: item[1])
    methods = [item[0] for item in ranked]
    maes = [item[1] for item in ranked]
    y = list(range(len(ranked)))

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 6.7,
            "axes.labelsize": 7.6,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.0,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "dmc-ice13-method-mae-ranking-v1",
        }
    )
    fig = plt.figure(figsize=(7.3, 10.0), facecolor="white")
    ax = fig.add_axes((0.225, 0.072, 0.72, 0.835), facecolor="#fafbfc")

    neutral = "#c8cdd3"
    current_colors = {
        r"g-xTB (adaptive $\leq 8^3$)": "#0072B2",
        r"GFN2-xTB (adaptive $\leq 4^3$)": "#6C8EBF",
        r"GFN1-xTB (adaptive $\leq 4^3$)": "#C94F45",
    }
    comparison_colors = {
        "B3LYP": "#D6B24C",
        "B3LYP-D4": "#D6B24C",
        "PBE0-D4": "#8F80C2",
        "PBE0": "#B9A7D8",
        "PBE-D4": "#55A868",
        "PBE": "#77B77A",
    }
    for yi, (method, value) in enumerate(ranked):
        color = current_colors.get(method, comparison_colors.get(method, neutral))
        linewidth = 2.4 if method in current_colors else 1.2
        if method in comparison_colors:
            linewidth = 1.8
        marker_size = 4.8 if method in current_colors else 3.0
        ax.hlines(yi, 0.3, value, color=color, linewidth=linewidth, zorder=2)
        ax.plot(value, yi, "o", ms=marker_size, mfc="white", mec=color, mew=1.0, zorder=3)
        ax.annotate(f"{value:.2f}", (value, yi), xytext=(4, 0), textcoords="offset points", va="center", fontsize=5.8, color="#344154")

    ax.set_xscale("log")
    ax.set_xlim(0.3, 15.0)
    ax.set_xticks([0.3, 0.5, 1, 2, 5, 10, 15], ["0.3", "0.5", "1", "2", "5", "10", "15"])
    ax.set_yticks(y, methods)
    ax.invert_yaxis()
    ax.set_xlabel(r"Relative-energy MAE vs DMC / kJ mol$^{-1}$ H$_2$O$^{-1}$")
    ax.grid(axis="x", which="major", color="#e5e8ec", linewidth=0.55)
    ax.axvline(1.0, color="#91a3ba", linestyle=(0, (3, 3)), linewidth=0.65)
    ax.axvline(10.0, color="#91a3ba", linestyle=(0, (3, 3)), linewidth=0.65)
    for spine in ax.spines.values():
        spine.set_color("#d9dee5")
        spine.set_linewidth(0.6)
    ax.tick_params(axis="y", length=0)
    for tick in ax.get_yticklabels():
        if tick.get_text() in current_colors:
            tick.set_fontweight("bold")

    ax.text(1.0, 0.996, r"1 kJ mol$^{-1}$", transform=ax.get_xaxis_transform(),
            fontsize=5.8, color="#536176", va="top", ha="left")
    ax.text(10.0, 0.996, r"10 kJ mol$^{-1}$", transform=ax.get_xaxis_transform(),
            fontsize=5.8, color="#536176", va="top", ha="left")

    fig.text(0.055, 0.962, "DMC-ICE13 relative-energy MAE ranking", fontsize=13.0, fontweight="bold", color="#202936")
    fig.text(
        0.055,
        0.938,
        "All values are relative to ice Ih; MAE is evaluated over the 12 non-reference polymorphs.",
        fontsize=6.8,
        color="#657184",
    )
    fig.text(0.055, 0.914, "Method", fontsize=7.3, fontweight="bold", color="#334055")
    fig.text(0.945, 0.914, "MAE", fontsize=7.3, fontweight="bold", color="#334055", ha="right")

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        target = args.output_stem.with_suffix(f".{suffix}")
        options: dict[str, object] = {"bbox_inches": "tight"}
        if suffix == "png":
            options["dpi"] = 400
        elif suffix == "pdf":
            options["metadata"] = {"CreationDate": None, "ModDate": None}
        else:
            options["metadata"] = {"Date": "2026-07-22"}
        fig.savefig(target, **options)
    plt.close(fig)


if __name__ == "__main__":
    main()
