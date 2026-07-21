#!/usr/bin/env python3
"""Generate the LC10 method-context figure used by the Part-I manuscript."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "lc10_method_comparison.csv"
OUT = ROOT / "figures" / "lc10_method_comparison"
ORDER = ("HF", "MP2", "SCS-MP2", "SOS-MP2", "GFN1-xTB", "GFN2-xTB", "g-xTB")
LABELS = {
    "HF": "HF",
    "MP2": "MP2",
    "SCS-MP2": "SCS-MP2",
    "SOS-MP2": "SOS-MP2",
    "GFN1-xTB": "GFN1-xTB ($4^3$)",
    "GFN2-xTB": "GFN2-xTB ($4^3$)",
    "g-xTB": "g-xTB current ($7^3$)",
}


def add_panel(ax, methods, values, xlabel, xmax, decimals):
    y = np.arange(len(methods))
    colors = ["#d9d9d9"] * 4 + ["#a6a6a6"] * 2 + ["#2166ac"]
    bars = ax.barh(y, values, color=colors, edgecolor="#404040", linewidth=0.7)
    ax.set_yticks(y, methods)
    ax.set_xlim(0.0, xmax)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#d0d0d0", linewidth=0.6, alpha=0.75)
    ax.set_axisbelow(True)
    pad = 0.012 * xmax
    for bar, value in zip(bars, values, strict=True):
        ax.text(
            value + pad,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.{decimals}f}",
            va="center",
            ha="left",
            fontsize=8.5,
        )


def main():
    with DATA.open(newline="") as handle:
        by_method = {row["method"]: row for row in csv.DictReader(handle)}
    methods = [LABELS[name] for name in ORDER]
    lattice_mae = np.array([float(by_method[name]["a0_mae_A"]) for name in ORDER])
    cohesive_mae = np.array(
        [float(by_method[name]["ecoh_mae_eV_per_atom"]) for name in ORDER]
    )
    plt.rcParams.update(
        {
            "font.size": 9.5,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "font.family": "DejaVu Sans",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.25), constrained_layout=True)
    add_panel(axes[0], methods, lattice_mae, r"MAE($a_0$) / $\AA$", 0.165, 3)
    add_panel(
        axes[1],
        methods,
        cohesive_mae,
        r"MAE($E_{\mathrm{coh}}$) / eV atom$^{-1}$",
        1.82,
        3,
    )
    axes[0].text(0.015, 0.975, "(a)", transform=axes[0].transAxes, va="top", fontweight="bold")
    axes[1].text(0.015, 0.975, "(b)", transform=axes[1].transAxes, va="top", fontweight="bold")
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=240, bbox_inches="tight")


if __name__ == "__main__":
    main()
