#!/usr/bin/env python3
"""Plot the final LC10 uniform and adaptive-endpoint MAE sequence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LABELS = (
    r"$3^3$",
    r"$4^3$",
    r"$5^3$",
    r"$6^3$",
    r"$\leq 7^3$",
    r"$\leq 8^3$",
    r"$\leq 9^3$",
)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(LABELS):
        raise ValueError("LC10 paper curve requires exactly seven stages")
    if [int(row["stage"]) for row in rows] != list(range(len(LABELS))):
        raise ValueError("LC10 stages must be ordered from zero through six")
    return rows


def plot(rows: list[dict[str, str]], output_pdf: Path, output_png: Path) -> None:
    x = list(range(len(rows)))
    a0 = [float(row["a0_mae_A"]) for row in rows]
    ecoh = [float(row["ecoh_mae_eV_per_atom"]) for row in rows]

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 5.6), sharex=True)
    series = (
        (axes[0], a0, "MAE($a_0$) / Å", "o", "#1f1f1f", 3),
        (axes[1], ecoh, r"MAE($E_{\rm coh}$) / eV atom$^{-1}$", "s", "#245b8a", 3),
    )
    for ax, values, ylabel, marker, color, decimals in series:
        ax.plot(x, values, color=color, marker=marker, linewidth=1.5, markersize=4.5)
        ax.axvline(4.5, color="0.65", linewidth=0.9, linestyle="--")
        ax.grid(axis="y", color="0.88", linewidth=0.6)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.25, len(rows) - 0.75)
        spread = max(values) - min(values)
        pad = max(0.08 * spread, 0.004 if ax is axes[0] else 0.015)
        ax.set_ylim(min(values) - pad, max(values) + 2.1 * pad)
        for idx, value in enumerate(values):
            ax.annotate(
                f"{value:.{decimals}f}",
                (idx, value),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.2,
            )
    axes[0].text(2.0, axes[0].get_ylim()[1], "uniform meshes", ha="center", va="top", fontsize=8)
    axes[0].text(5.0, axes[0].get_ylim()[1], "adaptive endpoints", ha="center", va="top", fontsize=8)
    axes[1].set_xticks(x, LABELS)
    axes[1].set_xlabel("LC10 k-point stage")
    fig.tight_layout(pad=0.8)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight")
    fig.savefig(output_png, dpi=240, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=root / "data" / "lc10_adaptive_progress_mae.csv",
    )
    parser.add_argument(
        "--pdf",
        type=Path,
        default=root / "figures" / "lc10_gxtb_adaptive_mae.pdf",
    )
    parser.add_argument(
        "--png",
        type=Path,
        default=root / "figures" / "lc10_gxtb_adaptive_mae.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot(load_rows(args.input), args.pdf, args.png)


if __name__ == "__main__":
    main()
