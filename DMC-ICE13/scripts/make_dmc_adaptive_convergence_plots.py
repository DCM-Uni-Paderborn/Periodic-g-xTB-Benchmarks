#!/usr/bin/env python3
"""Build the paper and SI plots for adaptive DMC-ICE13 k convergence."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PHASES = (
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
)
MESH_EDGES = tuple(range(1, 12))
THRESHOLD_KJMOL_PER_H2O = 0.05


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return value


def mesh_id(edge: int) -> str:
    if edge < 1:
        raise ValueError("mesh edge must be positive")
    return "k" + str(edge) * 3


def error_statistics(errors: list[float]) -> dict[str, float]:
    if len(errors) != len(PHASES):
        raise ValueError("DMC-ICE13 statistics require twelve non-reference phases")
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(value) for value in errors) / len(errors),
        "RMSE": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "MaxAE": max(abs(value) for value in errors),
    }


def build_data(root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    raw_path = root / "data/dmc_ice13_gxtb_spglib_kpoint_results.json"
    phasewise_path = root / "data/dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json"
    raw = load_object(raw_path)
    phasewise = load_object(phasewise_path)

    selected = phasewise["methods"]["GXTB"]["phase_convergence"]
    if set(selected) != set(PHASES):
        raise ValueError("phase-wise source does not contain the canonical twelve phases")

    frontier: list[dict[str, Any]] = []
    for edge in MESH_EDGES:
        mesh = mesh_id(edge)
        errors: list[float] = []
        newly_frozen: list[str] = []
        for phase in PHASES:
            record = selected[phase]
            selected_edge = int(record["mesh_n"])
            if selected_edge <= edge:
                relative = float(record["relative_energy_kjmol_per_h2o"])
            else:
                try:
                    relative = float(raw["results"][mesh]["GXTB"]["relative_kjmol"][phase])
                except KeyError as error:
                    raise ValueError(f"missing active value for {phase} at {mesh}") from error
            reference = float(record["dmc_relative_kjmol_per_h2o"])
            errors.append(relative - reference)
            if selected_edge == edge:
                newly_frozen.append(phase)

        statistics = error_statistics(errors)
        frontier.append(
            {
                "mesh_n": edge,
                "mesh_id": mesh,
                "mesh_label": "Gamma" if edge == 1 else f"{edge}x{edge}x{edge}",
                "frontier_me_kjmol_h2o": statistics["ME"],
                "frontier_mae_kjmol_h2o": statistics["MAE"],
                "frontier_rmse_kjmol_h2o": statistics["RMSE"],
                "frontier_maxae_kjmol_h2o": statistics["MaxAE"],
                "newly_frozen_phases": ";".join(newly_frozen),
                "active_phases_after_mesh": sum(
                    int(selected[phase]["mesh_n"]) > edge for phase in PHASES
                ),
            }
        )

    vii_reference = float(selected["VII"]["dmc_relative_kjmol_per_h2o"])
    phase_vii: list[dict[str, Any]] = []
    previous: float | None = None
    for edge in MESH_EDGES:
        mesh = mesh_id(edge)
        relative = float(raw["results"][mesh]["GXTB"]["relative_kjmol"]["VII"])
        phase_vii.append(
            {
                "mesh_n": edge,
                "mesh_id": mesh,
                "mesh_label": "Gamma" if edge == 1 else f"{edge}x{edge}x{edge}",
                "relative_energy_kjmol_h2o": relative,
                "dmc_reference_kjmol_h2o": vii_reference,
                "error_kjmol_h2o": relative - vii_reference,
                "absolute_adjacent_delta_kjmol_h2o": (
                    None if previous is None else abs(relative - previous)
                ),
            }
        )
        previous = relative

    if not math.isclose(frontier[-1]["frontier_mae_kjmol_h2o"], 3.906294510467664, abs_tol=1.0e-12):
        raise ValueError("unexpected final adaptive-frontier MAE")
    if not math.isclose(
        phase_vii[-1]["absolute_adjacent_delta_kjmol_h2o"],
        0.031136099878160906,
        abs_tol=1.0e-14,
    ):
        raise ValueError("unexpected final ice-VII adjacent-mesh change")
    if phase_vii[-2]["absolute_adjacent_delta_kjmol_h2o"] <= THRESHOLD_KJMOL_PER_H2O:
        raise ValueError("ice VII unexpectedly passes before 11x11x11")
    if phase_vii[-1]["absolute_adjacent_delta_kjmol_h2o"] > THRESHOLD_KJMOL_PER_H2O:
        raise ValueError("ice VII does not satisfy the declared final threshold")

    manifest = {
        "schema_version": 1,
        "benchmark": "DMC-ICE13",
        "method": "g-xTB",
        "mesh_sequence": [mesh_id(edge) for edge in MESH_EDGES],
        "first_mesh_presentation": (
            "explicit unshifted 1x1x1 endpoint labelled Gamma after a separate "
            "implicit-Gamma/explicit-1x1x1 equivalence check"
        ),
        "phase_convergence_threshold_kjmol_per_h2o": THRESHOLD_KJMOL_PER_H2O,
        "frontier_definition": (
            "At mesh N, retain the accepted denser endpoint for phases already "
            "converged and use the current N value for unresolved phases."
        ),
        "frontier_is_stopping_gate": False,
        "sources": {
            raw_path.relative_to(root).as_posix(): sha256(raw_path),
            phasewise_path.relative_to(root).as_posix(): sha256(phasewise_path),
        },
        "adaptive_frontier": frontier,
        "phase_vii": phase_vii,
    }
    return frontier, phase_vii, manifest


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig: mpl.figure.Figure, figures: Path, stem: str) -> None:
    for suffix in ("pdf", "svg", "png"):
        options: dict[str, Any] = {"bbox_inches": "tight"}
        if suffix == "pdf":
            options["metadata"] = {
                "CreationDate": None,
                "ModDate": None,
                "Creator": "Periodic-g-xTB-Benchmarks",
            }
        elif suffix == "svg":
            options["metadata"] = {
                "Date": "2026-07-15",
                "Creator": "Periodic-g-xTB-Benchmarks",
            }
        elif suffix == "png":
            options["dpi"] = 600
        fig.savefig(figures / f"{stem}.{suffix}", **options)
    svg = figures / f"{stem}.svg"
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
    )


def tick_labels(meshes: np.ndarray) -> list[str]:
    return [r"$\Gamma$"] + [rf"${int(edge)}^3$" for edge in meshes[1:]]


def plot_frontier(frontier: list[dict[str, Any]], figures: Path) -> None:
    mesh = np.asarray([row["mesh_n"] for row in frontier], dtype=float)
    mae = np.asarray([row["frontier_mae_kjmol_h2o"] for row in frontier])
    fig, ax = plt.subplots(figsize=(3.37, 2.55), constrained_layout=True)
    ax.plot(mesh, mae, color="#0072B2", marker="o", markersize=4.0, linewidth=1.6)
    ax.set_yscale("log")
    ax.set_xlim(0.65, 11.35)
    ax.set_ylim(3.2, 220.0)
    ax.set_xticks(mesh)
    labels = ax.set_xticklabels(tick_labels(mesh), rotation=45, ha="right")
    labels[0].set_rotation(0)
    labels[0].set_ha("center")
    ax.set_xlabel("k-point mesh")
    ax.set_ylabel("Adaptive-frontier MAE\n" + r"(kJ mol$^{-1}$ H$_2$O$^{-1}$)")
    ax.grid(axis="y", which="major", color="0.86", linewidth=0.6)
    ax.grid(axis="y", which="minor", color="0.93", linewidth=0.4)
    ax.annotate(
        "163.8",
        xy=(mesh[0], mae[0]),
        xytext=(4, -12),
        textcoords="offset points",
        ha="left",
        va="top",
    )
    ax.annotate(
        "3.906",
        xy=(mesh[-1], mae[-1]),
        xytext=(-4, 8),
        textcoords="offset points",
        ha="right",
        va="bottom",
    )
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, figures, "dmc_ice13_gxtb_adaptive_frontier_mae")
    plt.close(fig)


def plot_phase_vii(rows: list[dict[str, Any]], figures: Path) -> None:
    mesh = np.asarray([row["mesh_n"] for row in rows], dtype=float)
    relative = np.asarray([row["relative_energy_kjmol_h2o"] for row in rows])
    reference = float(rows[0]["dmc_reference_kjmol_h2o"])
    delta = np.asarray([row["absolute_adjacent_delta_kjmol_h2o"] for row in rows[1:]])

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(6.6, 4.8),
        sharex=True,
        gridspec_kw={"height_ratios": (1.25, 1.0)},
        constrained_layout=True,
    )
    top, bottom = axes
    top.plot(mesh, relative, color="#0072B2", marker="o", markersize=4.5, linewidth=1.6, label="g-xTB")
    top.axhline(reference, color="#D55E00", linestyle="--", linewidth=1.4, label="DMC reference")
    top.set_ylabel(r"$R_{\mathrm{VII}}(N)$ (kJ mol$^{-1}$ H$_2$O$^{-1}$)")
    top.legend(frameon=False, loc="lower right", ncol=2)
    top.grid(axis="y", color="0.90", linewidth=0.5)
    top.annotate(
        r"$R_{\mathrm{VII}}(11)=5.032$",
        xy=(mesh[-1], relative[-1]),
        xytext=(-8, 18),
        textcoords="offset points",
        ha="right",
        va="bottom",
        arrowprops={"arrowstyle": "-", "color": "0.35", "linewidth": 0.7},
    )

    bottom.plot(mesh[1:], delta, color="#0072B2", marker="o", markersize=4.5, linewidth=1.6)
    bottom.axhline(
        THRESHOLD_KJMOL_PER_H2O,
        color="#D55E00",
        linestyle="--",
        linewidth=1.4,
        label=r"criterion: $0.05$",
    )
    bottom.set_yscale("log")
    bottom.set_ylabel(r"$|R(N)-R(N-1)|$" + "\n" + r"(kJ mol$^{-1}$ H$_2$O$^{-1}$)")
    bottom.set_xlabel("k-point mesh")
    bottom.set_xticks(mesh)
    bottom.set_xticklabels(tick_labels(mesh))
    bottom.grid(axis="y", which="major", color="0.86", linewidth=0.6)
    bottom.grid(axis="y", which="minor", color="0.93", linewidth=0.4)
    bottom.legend(frameon=False, loc="upper right")
    bottom.annotate(
        "0.0311",
        xy=(mesh[-1], delta[-1]),
        xytext=(-7, 9),
        textcoords="offset points",
        ha="right",
        va="bottom",
    )
    for ax in axes:
        ax.set_xlim(0.65, 11.35)
        ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, figures, "dmc_ice13_gxtb_phase_vii_kpoint_convergence")
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data"
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.8,
            "lines.solid_capstyle": "round",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "dmc-ice13-gxtb-adaptive-v1",
        }
    )

    frontier, phase_vii, manifest = build_data(root)
    write_csv(data_dir / "dmc_ice13_gxtb_adaptive_frontier.csv", frontier)
    write_csv(data_dir / "dmc_ice13_gxtb_phase_vii_kpoint_convergence.csv", phase_vii)
    (data_dir / "dmc_ice13_gxtb_adaptive_frontier.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    plot_frontier(frontier, figures)
    plot_phase_vii(phase_vii, figures)


if __name__ == "__main__":
    main()
