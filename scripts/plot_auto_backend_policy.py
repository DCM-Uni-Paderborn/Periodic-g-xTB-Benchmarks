#!/usr/bin/env python3
"""Render the deterministic periodic g-xTB AUTO backend policy."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation" / "accelerated_exchange" / "figures"

BLUE = "#1769aa"
BLUE_LIGHT = "#e9f2f9"
BLUE_MID = "#b9d7ea"
INK = "#17202a"
GRAY = "#5f6b73"
FIXED_DATE = datetime(2026, 7, 22, tzinfo=timezone.utc)


def box(ax, xy, width, height, text, *, fill=BLUE_LIGHT, linewidth=1.0):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        facecolor=fill,
        edgecolor=BLUE,
        linewidth=linewidth,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=8.2,
        color=INK,
    )
    return patch


def decision(ax, center, width, height, text):
    x, y = center
    vertices = [
        (x, y + height / 2),
        (x + width / 2, y),
        (x, y - height / 2),
        (x - width / 2, y),
    ]
    patch = Polygon(vertices, closed=True, facecolor=BLUE_MID, edgecolor=BLUE)
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=8.3, color=INK)
    return patch


def arrow(ax, start, end, label=None, label_offset=(0.0, 0.0)):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.0,
        color=GRAY,
        shrinkA=2,
        shrinkB=2,
        connectionstyle="arc3,rad=0",
    )
    ax.add_patch(patch)
    if label:
        ax.text(
            (start[0] + end[0]) / 2 + label_offset[0],
            (start[1] + end[1]) / 2 + label_offset[1],
            label,
            ha="center",
            va="center",
            fontsize=7.5,
            color=GRAY,
        )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams["svg.hashsalt"] = "periodic-gxtb-part-II-auto-policy"
    fig, ax = plt.subplots(figsize=(10.0, 4.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.6)
    ax.axis("off")

    box(
        ax,
        (0.25, 1.75),
        1.55,
        0.85,
        "Validated complete\nregular k mesh",
        fill="#f5f8fa",
    )
    decision(ax, (2.75, 2.18), 1.45, 1.2, "$N_k=1$?")
    box(
        ax,
        (4.0, 3.35),
        2.05,
        0.9,
        "Complete array / dense reverse\nDense transform / dense star",
    )
    decision(ax, (4.7, 2.18), 1.55, 1.25, "$N_{red}<N_k$?")
    box(
        ax,
        (6.1, 3.05),
        2.35,
        1.2,
        "Symmetry fused / streamed reverse\nMixed-radix FFT\nStreamed star with native mixer",
    )
    decision(ax, (6.45, 1.5), 1.65, 1.25, "MPI ranks\n$>1$?")
    box(
        ax,
        (8.25, 1.65),
        1.5,
        1.05,
        "Bounded batch\nStreamed reverse\nDense transform",
    )
    box(
        ax,
        (8.25, 0.25),
        1.5,
        1.05,
        "Distributed partial images\nStreamed reverse\nDense transform",
    )

    arrow(ax, (1.8, 2.18), (2.02, 2.18))
    arrow(ax, (3.12, 2.72), (4.0, 3.55), "yes", (0.0, 0.12))
    arrow(ax, (3.48, 2.18), (3.92, 2.18), "no", (0.0, 0.14))
    arrow(ax, (5.17, 2.69), (6.1, 3.2), "yes", (0.0, 0.12))
    arrow(ax, (5.48, 1.85), (5.75, 1.68), "no", (-0.02, 0.14))
    arrow(ax, (7.28, 1.5), (8.25, 2.05), "no", (0.0, 0.13))
    arrow(ax, (6.82, 0.98), (8.25, 0.78), "yes", (0.0, 0.14))

    ax.text(
        0.25,
        0.28,
        "Every branch: bounded ACP cache enabled; implicit qualification disabled; dense oracle retained in MANUAL mode.",
        ha="left",
        va="bottom",
        fontsize=8.1,
        color=INK,
    )

    fig.tight_layout(pad=0.15)
    metadata = {
        "pdf": {
            "Title": "Periodic g-xTB automatic exact backend policy",
            "Author": "Thomas D. Kühne",
            "Subject": "Part-II reproducibility figure",
            "Creator": "scripts/plot_auto_backend_policy.py",
            "CreationDate": FIXED_DATE,
            "ModDate": FIXED_DATE,
        },
        "svg": {
            "Title": "Periodic g-xTB automatic exact backend policy",
            "Creator": "scripts/plot_auto_backend_policy.py",
            "Description": "Part-II reproducibility figure",
            "Date": "2026-07-22",
        },
    }
    for suffix in ("pdf", "svg"):
        fig.savefig(
            OUTPUT / f"auto_backend_policy.{suffix}",
            bbox_inches="tight",
            pad_inches=0.03,
            metadata=metadata[suffix],
        )
    plt.close(fig)

    svg_path = OUTPUT / "auto_backend_policy.svg"
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
