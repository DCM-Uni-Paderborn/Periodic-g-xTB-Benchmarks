#!/usr/bin/env python3
"""Plot the LC10 uniform-to-adaptive g-xTB k-point convergence history."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STAGES = (
    "k333",
    "k444",
    "k555",
    "k666",
    "k777_cleaned",
    "mixed_k888",
    "mixed_k999",
    "adaptive_final",
)

LABELS = {
    "k333": r"$3^3$",
    "k444": r"$4^3$",
    "k555": r"$5^3$",
    "k666": r"$6^3$",
    "k777_cleaned": r"$7^3$",
    "mixed_k888": r"$\leq 8^3$",
    "mixed_k999": r"$\leq 9^3$",
    "adaptive_final": "adaptive\nfinal",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        by_stage = {row["stage"]: row for row in csv.DictReader(handle)}
    missing = [stage for stage in STAGES if stage not in by_stage]
    if missing:
        raise ValueError(f"missing LC10 convergence stages: {', '.join(missing)}")
    return [by_stage[stage] for stage in STAGES]


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
        (axes[0], a0, r"MAE($a_0$) / $\mathrm{\AA}$", "o", "#1f1f1f", 3),
        (
            axes[1],
            ecoh,
            r"MAE($E_{\rm coh}$) / eV atom$^{-1}$",
            "s",
            "#245b8a",
            3,
        ),
    )
    for axis, values, ylabel, marker, color, decimals in series:
        axis.plot(
            x,
            values,
            color=color,
            marker=marker,
            linewidth=1.5,
            markersize=4.5,
        )
        axis.axvline(4.5, color="0.65", linewidth=0.9, linestyle="--")
        axis.grid(axis="y", color="0.88", linewidth=0.6)
        axis.set_ylabel(ylabel)
        axis.set_xlim(-0.25, len(rows) - 0.75)
        spread = max(values) - min(values)
        pad = max(0.08 * spread, 0.004 if axis is axes[0] else 0.015)
        axis.set_ylim(min(values) - pad, max(values) + 2.1 * pad)
        for index, value in enumerate(values):
            axis.annotate(
                f"{value:.{decimals}f}",
                (index, value),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.2,
            )
    axes[0].text(
        2.0,
        axes[0].get_ylim()[1],
        "uniform meshes",
        ha="center",
        va="top",
        fontsize=8,
    )
    axes[0].text(
        6.0,
        axes[0].get_ylim()[1],
        "adaptive selections",
        ha="center",
        va="top",
        fontsize=8,
    )
    axes[1].set_xticks(x, [LABELS[row["stage"]] for row in rows])
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
        default=root / "data" / "lc10_gxtb_adaptive_kmesh_mae.csv",
    )
    parser.add_argument(
        "--pdf",
        type=Path,
        default=root / "figures" / "lc10_gxtb_adaptive_kmesh_mae.pdf",
    )
    parser.add_argument(
        "--png",
        type=Path,
        default=root / "figures" / "lc10_gxtb_adaptive_kmesh_mae.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot(load_rows(args.input), args.pdf, args.png)


if __name__ == "__main__":
    main()
