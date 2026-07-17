#!/usr/bin/env python3
"""Summarize the controlled DMC-ICE13 0D/3D dilation diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np


HARTREE_TO_KJMOL = 2625.4996394799
BUILD = {
    "cp2k_revision": "28df9380abb327d56bbf216d2469a1fd8c953fc0",
    "save_tblite_revision": "257ba442684c39454175e5192c8a2342b4c6380f",
    "cp2k_executable_sha256": "86949402a326e3118551fea079dd2b79df331087b47cb355b8803964f15deb89",
    "save_tblite_library_sha256": "8ac8c98f462c6b29a2350ed341bf310addbc4692b0f6339a28ecb26c996c13a4",
}
MODELS = {
    "full": None,
    "no_exchange": "c2da7bf8fd7157227433cce3f2f222eb22d4ccfe109fbdce9baff810f1671357",
    "frozen_qvszp": "321e9f68609f1ad48ff961da6e77e82eec128ed4e44568105a225dd3b13f9e56",
    "no_anisotropic_multipole": "bdf99c6bda2cbd6b815ef86e31dd4a213b1e526081cf11a8d337408bfe477333",
    "no_acp": "d191752b3811ee5fcf2132bf088bf67358b399ff2c4f0540dcbde2a451235cef",
}
COMPONENT_LABELS = {
    "core": "Core Hamiltonian energy",
    "repulsive": "Repulsive potential energy",
    "electrostatic": "Electrostatic energy",
    "dispersion_sc": "Self-consistent dispersion energy",
    "dispersion_non_sc": "Non-self consistent dispersion energy",
    "halogen": "Correction for halogen bonding",
    "electronic_entropy": "Electronic entropic energy",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_case_name(name: str) -> tuple[str, str, float]:
    match = re.fullmatch(r"dmc_(Ih|VII|XIV)_GXTB_(mol|pbc)_s(1p5|2|3|4)_k222", name)
    if not match:
        raise RuntimeError(f"unexpected case name: {name}")
    return match.group(1), match.group(2), float(match.group(3).replace("p", "."))


def parse_input(path: Path) -> dict:
    text = path.read_text(errors="replace")
    cell_match = re.search(r"(?ms)&CELL\s*(.*?)&END CELL", text)
    coord_match = re.search(r"(?ms)&COORD\s*(.*?)&END COORD", text)
    if not cell_match or not coord_match:
        raise RuntimeError(f"CELL/COORD missing: {path}")
    cell_rows = []
    for key in ("A", "B", "C"):
        match = re.search(
            rf"(?m)^\s*{key}\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
            cell_match.group(1),
        )
        if not match:
            raise RuntimeError(f"cell vector {key} missing: {path}")
        cell_rows.append([float(match.group(i)) for i in range(1, 4)])
    cell = np.asarray(cell_rows, dtype=float)
    lengths = np.linalg.norm(cell, axis=1)

    coords = []
    symbols = []
    for line in coord_match.group(1).splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0] in {"H", "O"}:
            symbols.append(fields[0])
            coords.append([float(value) for value in fields[1:]])
    if len(coords) != 36 or symbols != [item for _ in range(12) for item in ("O", "H", "H")]:
        raise RuntimeError(f"expected 12 ordered waters in {path}")

    periodic_tokens = re.findall(r"(?m)^\s*PERIODIC\s+(\S+)\s*$", text)
    poisson_match = re.search(r"(?m)^\s*POISSON_SOLVER\s+(\S+)\s*$", text)
    kmesh_match = re.search(r"(?m)^\s*SCHEME\s+MONKHORST-PACK\s+(\d+)\s+(\d+)\s+(\d+)\s*$", text)
    return {
        "cell": cell,
        "cell_lengths_angstrom": lengths.tolist(),
        "max_cell_vector_angstrom": float(max(lengths)),
        "symbols": symbols,
        "coords": np.asarray(coords, dtype=float),
        "periodic_tokens": periodic_tokens,
        "poisson_solver": poisson_match.group(1) if poisson_match else None,
        "kmesh": [int(kmesh_match.group(i)) for i in range(1, 4)] if kmesh_match else None,
        "has_kpoints_section": "&KPOINTS" in text,
        "full_grid": bool(re.search(r"(?m)^\s*FULL_GRID\s+(?:ON|T|TRUE)\s*$", text, re.I)),
        "symmetry_off": bool(re.search(r"(?m)^\s*SYMMETRY\s+(?:OFF|F|FALSE)\s*$", text, re.I)),
        "input_sha256": sha256(path),
    }


def parse_output(path: Path) -> dict:
    text = path.read_text(errors="replace")
    energy_matches = re.findall(r"(?m)^\s*ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$", text)
    energy_is_force_eval = bool(energy_matches)
    if not energy_matches:
        energy_matches = re.findall(r"(?m)^\s*Total energy:\s+([-+0-9.Ee]+)\s*$", text)
    if not energy_matches:
        raise RuntimeError(f"final or last-iteration energy missing: {path}")
    total = float(energy_matches[-1])
    steps_match = re.search(r"SCF run converged in\s+(\d+) steps", text)
    iteration_lines = [
        line.split()
        for line in text.splitlines()
        if re.match(r"^\s*\d+\s+GXTB-", line)
    ]
    final_residual = float(iteration_lines[-1][4]) if iteration_lines else math.nan
    components = {}
    for key, label in COMPONENT_LABELS.items():
        matches = re.findall(rf"(?m)^\s*{re.escape(label)}:\s+([-+0-9.Ee]+)\s*$", text)
        components[key] = float(matches[-1]) if matches else 0.0
    components["other_interactions_residual"] = total - sum(components.values())
    return {
        "total_energy_eh": total,
        "scf_converged": steps_match is not None,
        "scf_steps": int(steps_match.group(1)) if steps_match else None,
        "final_scc_residual": final_residual,
        "program_ended": "PROGRAM ENDED AT" in text,
        "energy_is_completed_force_eval": energy_is_force_eval,
        "components_eh": components,
        "output_sha256": sha256(path),
    }


def geometry_audit(pbc: dict, mol: dict) -> dict:
    cell = pbc["cell"]
    pcoords = pbc["coords"]
    mcoords = mol["coords"]
    global_shift = mcoords[0] - pcoords[0]
    max_alignment = 0.0
    max_intra_delta = 0.0
    lattice_shifts = []
    for water in range(12):
        start = 3 * water
        diff = mcoords[start] - pcoords[start] - global_shift
        fractional = np.linalg.solve(cell.T, diff)
        integer = np.rint(fractional).astype(int)
        residual = diff - cell.T @ integer
        max_alignment = max(max_alignment, float(np.linalg.norm(residual)))
        lattice_shifts.append(integer.tolist())
        for h_offset in (1, 2):
            p_dist = np.linalg.norm(pcoords[start + h_offset] - pcoords[start])
            m_dist = np.linalg.norm(mcoords[start + h_offset] - mcoords[start])
            max_intra_delta = max(max_intra_delta, abs(float(p_dist - m_dist)))
            aligned = pcoords[start + h_offset] + global_shift + cell.T @ integer
            max_alignment = max(max_alignment, float(np.linalg.norm(aligned - mcoords[start + h_offset])))
    return {
        "max_coordinate_alignment_error_angstrom": max_alignment,
        "max_intramolecular_oh_distance_delta_angstrom": max_intra_delta,
        "global_translation_angstrom": global_shift.tolist(),
        "per_water_lattice_shifts": lattice_shifts,
    }


def load_old(path: Path) -> dict[tuple[str, float], dict[str, float]]:
    result = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["method"] == "GXTB" and row["phase"] in {"Ih", "VII", "XIV"}:
                key = (row["phase"], float(row["scale"]))
                result[key] = {
                    "delta_kjmol_per_h2o": float(row["pbc_minus_mol_kjmol_per_h2o"]),
                    "pbc_energy_per_h2o_eh": float(row["pbc_energy_per_h2o_hartree"]),
                    "mol_energy_per_h2o_eh": float(row["mol_energy_per_h2o_hartree"]),
                }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-root", type=Path, required=True)
    parser.add_argument("--component-root", type=Path, required=True)
    parser.add_argument("--old-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    roots = {"full": args.full_root}
    roots.update({model: args.component_root / model for model in MODELS if model != "full"})
    old = load_old(args.old_csv)
    rows = []
    geometry_records = []

    for model, root in roots.items():
        parsed = {}
        for run in sorted((root / "runs").iterdir()):
            if not run.is_dir():
                continue
            phase, mode, scale = parse_case_name(run.name)
            inp = parse_input(run / "input.inp")
            out = parse_output(run / "cp2k.out")
            rc = int((run / "returncode.txt").read_text().strip())
            parsed[(phase, scale, mode)] = {"input": inp, "output": out, "returncode": rc}

        for phase, scale in sorted({(key[0], key[1]) for key in parsed}):
            pbc = parsed[(phase, scale, "pbc")]
            mol = parsed[(phase, scale, "mol")]
            pbc_in, mol_in = pbc["input"], mol["input"]
            pbc_out, mol_out = pbc["output"], mol["output"]
            old_point = old.get((phase, scale)) if model == "full" else None
            audit = geometry_audit(pbc_in, mol_in)
            geometry_records.append({"model": model, "phase": phase, "scale": scale, **audit})

            if max(pbc_in["cell_lengths_angstrom"] + mol_in["cell_lengths_angstrom"]) > 40.0 + 1.0e-8:
                raise RuntimeError(f">40 A point entered analysis: {model} {phase} {scale}")
            if pbc_in["periodic_tokens"] != ["XYZ", "XYZ"]:
                raise RuntimeError(f"unexpected PBC mask: {model} {phase} {scale}: {pbc_in['periodic_tokens']}")
            if mol_in["periodic_tokens"] != ["NONE", "NONE"]:
                raise RuntimeError(f"unexpected 0D mask: {model} {phase} {scale}: {mol_in['periodic_tokens']}")
            if pbc_in["poisson_solver"] != "PERIODIC" or mol_in["poisson_solver"] != "ANALYTIC":
                raise RuntimeError(f"unexpected Poisson route: {model} {phase} {scale}")
            if pbc_in["kmesh"] != [2, 2, 2] or mol_in["kmesh"] is not None:
                raise RuntimeError(f"unexpected k mesh: {model} {phase} {scale}")

            factor = HARTREE_TO_KJMOL / 12.0
            valid_pair = all(
                (
                    pbc["returncode"] == 0,
                    mol["returncode"] == 0,
                    pbc_out["program_ended"],
                    mol_out["program_ended"],
                    pbc_out["scf_converged"],
                    mol_out["scf_converged"],
                )
            )
            component_delta = {
                key: ((pbc_out["components_eh"][key] - mol_out["components_eh"][key]) * factor if valid_pair else None)
                for key in pbc_out["components_eh"]
            }
            delta = ((pbc_out["total_energy_eh"] - mol_out["total_energy_eh"]) * factor if valid_pair else None)
            row = {
                "model": model,
                "phase": phase,
                "scale": scale,
                "n_h2o": 12,
                "pbc_max_cell_vector_A": pbc_in["max_cell_vector_angstrom"],
                "zero_d_cell_vector_A": mol_in["max_cell_vector_angstrom"],
                "pbc_minus_0d_kJ_mol_H2O": delta,
                "valid_converged_pair": valid_pair,
                "old_full_pbc_minus_0d_kJ_mol_H2O": old_point["delta_kjmol_per_h2o"] if old_point else None,
                "current_minus_old_pbc_kJ_mol_H2O": ((pbc_out["total_energy_eh"] / 12.0 - old_point["pbc_energy_per_h2o_eh"]) * HARTREE_TO_KJMOL if old_point else None),
                "current_minus_old_zero_d_kJ_mol_H2O": ((mol_out["total_energy_eh"] / 12.0 - old_point["mol_energy_per_h2o_eh"]) * HARTREE_TO_KJMOL if old_point else None),
                "pbc_scc_steps": pbc_out["scf_steps"],
                "zero_d_scc_steps": mol_out["scf_steps"],
                "pbc_final_scc_residual": pbc_out["final_scc_residual"],
                "zero_d_final_scc_residual": mol_out["final_scc_residual"],
                "pbc_returncode": pbc["returncode"],
                "zero_d_returncode": mol["returncode"],
                "pbc_program_ended": pbc_out["program_ended"],
                "zero_d_program_ended": mol_out["program_ended"],
                "geometry_alignment_max_error_A": audit["max_coordinate_alignment_error_angstrom"],
                "intramolecular_OH_max_delta_A": audit["max_intramolecular_oh_distance_delta_angstrom"],
                "parameter_sha256": MODELS[model],
                "build_id": f"CP2K {BUILD['cp2k_revision'][:12]} / save_tblite {BUILD['save_tblite_revision'][:12]}",
                "pbc_input_sha256": pbc_in["input_sha256"],
                "zero_d_input_sha256": mol_in["input_sha256"],
                "pbc_output_sha256": pbc_out["output_sha256"],
                "zero_d_output_sha256": mol_out["output_sha256"],
            }
            for key, value in component_delta.items():
                row[f"delta_{key}_kJ_mol_H2O"] = value
            rows.append(row)

    csv_path = args.output_dir / "dmc_dilation_component_diagnostics.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_point = {(row["model"], row["phase"], row["scale"]): row for row in rows}
    points = sorted({(row["phase"], row["scale"]) for row in rows})
    compact_rows = []
    for phase, scale in points:
        full = by_point[("full", phase, scale)]
        compact = {
            "phase": phase,
            "scale": scale,
            "max_periodic_cell_vector_A": full["pbc_max_cell_vector_A"],
            "full_deltaE_kJ_mol_H2O": full["pbc_minus_0d_kJ_mol_H2O"],
            "no_exchange_deltaE_kJ_mol_H2O": by_point[("no_exchange", phase, scale)]["pbc_minus_0d_kJ_mol_H2O"],
            "frozen_qvszp_deltaE_kJ_mol_H2O": by_point[("frozen_qvszp", phase, scale)]["pbc_minus_0d_kJ_mol_H2O"],
            "no_anisotropic_multipole_deltaE_kJ_mol_H2O": by_point[("no_anisotropic_multipole", phase, scale)]["pbc_minus_0d_kJ_mol_H2O"],
            "no_acp_deltaE_kJ_mol_H2O": by_point[("no_acp", phase, scale)]["pbc_minus_0d_kJ_mol_H2O"],
            "full_pbc_zero_d_scc_steps": f"{full['pbc_scc_steps']}/{full['zero_d_scc_steps']}",
            "full_max_final_scc_residual": max(full["pbc_final_scc_residual"], full["zero_d_final_scc_residual"]),
        }
        compact_rows.append(compact)
    compact_csv = args.output_dir / "dmc_dilation_energy_ablation_table.csv"
    with compact_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(compact_rows[0]))
        writer.writeheader()
        writer.writerows(compact_rows)

    component_rows = []
    for phase, scale in points:
        full = by_point[("full", phase, scale)]
        component_rows.append(
            {
                "phase": phase,
                "scale": scale,
                "delta_total_kJ_mol_H2O": full["pbc_minus_0d_kJ_mol_H2O"],
                "delta_core_kJ_mol_H2O": full["delta_core_kJ_mol_H2O"],
                "delta_repulsive_kJ_mol_H2O": full["delta_repulsive_kJ_mol_H2O"],
                "delta_electrostatic_kJ_mol_H2O": full["delta_electrostatic_kJ_mol_H2O"],
                "delta_dispersion_total_kJ_mol_H2O": full["delta_dispersion_sc_kJ_mol_H2O"] + full["delta_dispersion_non_sc_kJ_mol_H2O"],
                "delta_other_interactions_residual_kJ_mol_H2O": full["delta_other_interactions_residual_kJ_mol_H2O"],
                "pbc_zero_d_scc_steps": f"{full['pbc_scc_steps']}/{full['zero_d_scc_steps']}",
                "max_final_scc_residual": max(full["pbc_final_scc_residual"], full["zero_d_final_scc_residual"]),
            }
        )
    component_csv = args.output_dir / "dmc_dilation_full_model_component_table.csv"
    with component_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(component_rows[0]))
        writer.writeheader()
        writer.writerows(component_rows)

    legacy_rows = []
    for phase, scale in points:
        full = by_point[("full", phase, scale)]
        legacy_rows.append(
            {
                "phase": phase,
                "scale": scale,
                "legacy_folded_deltaE_kJ_mol_H2O": full["old_full_pbc_minus_0d_kJ_mol_H2O"],
                "current_coupled_deltaE_kJ_mol_H2O": full["pbc_minus_0d_kJ_mol_H2O"],
                "current_minus_legacy_pbc_kJ_mol_H2O": full["current_minus_old_pbc_kJ_mol_H2O"],
                "current_minus_legacy_zero_d_kJ_mol_H2O": full["current_minus_old_zero_d_kJ_mol_H2O"],
                "current_pbc_zero_d_scc_steps": f"{full['pbc_scc_steps']}/{full['zero_d_scc_steps']}",
            }
        )
    legacy_csv = args.output_dir / "dmc_dilation_full_vs_legacy_table.csv"
    with legacy_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(legacy_rows[0]))
        writer.writeheader()
        writer.writerows(legacy_rows)

    def fmt(value: float | None) -> str:
        return "--" if value is None else f"{value:.3f}"

    energy_tex = args.output_dir / "dmc_dilation_energy_ablation_table.tex"
    energy_lines = [
        r"\begin{table*}",
        r"\caption{Controlled DMC--ICE13 periodic-minus-0D energies per water molecule. All cell-vector lengths are at most 40~\AA. Values in kJ\,mol$^{-1}$\,H$_2$O$^{-1}$. The parameter deletions are diagnostics, not reparameterized physical models. Build: CP2K 28df9380abb3; save\_tblite 257ba442684c.}",
        r"\begin{tabular}{lrrrrrrl}",
        r"\hline\hline",
        r"Phase & $s$ & Full & no exchange & frozen q-vSZP & no anisotropic multipole & no ACP & SCC steps (3D/0D) \\",
        r"\hline",
    ]
    for row in compact_rows:
        energy_lines.append(
            f"{row['phase']} & {row['scale']:g} & {fmt(row['full_deltaE_kJ_mol_H2O'])} & "
            f"{fmt(row['no_exchange_deltaE_kJ_mol_H2O'])} & {fmt(row['frozen_qvszp_deltaE_kJ_mol_H2O'])} & "
            f"{fmt(row['no_anisotropic_multipole_deltaE_kJ_mol_H2O'])} & {fmt(row['no_acp_deltaE_kJ_mol_H2O'])} & "
            f"{row['full_pbc_zero_d_scc_steps']} \\\\"
        )
    energy_lines.extend([r"\hline\hline", r"\end{tabular}", r"\end{table*}", ""])
    energy_tex.write_text("\n".join(energy_lines))

    component_tex = args.output_dir / "dmc_dilation_full_model_component_table.tex"
    component_lines = [
        r"\begin{table*}",
        r"\caption{Full-model component differences for the controlled DMC--ICE13 dilation test in kJ\,mol$^{-1}$\,H$_2$O$^{-1}$. The residual is total energy minus separately printed CP2K components; it contains exchange and other tblite interactions and is not pure exchange. Core and residual terms are representation-dependent and strongly cancel.}",
        r"\begin{tabular}{lrrrrrrl}",
        r"\hline\hline",
        r"Phase & $s$ & Total & Core & Repulsion & Electrostatic & Dispersion & Other residual / SCC steps \\",
        r"\hline",
    ]
    for row in component_rows:
        component_lines.append(
            f"{row['phase']} & {row['scale']:g} & {row['delta_total_kJ_mol_H2O']:.3f} & "
            f"{row['delta_core_kJ_mol_H2O']:.1f} & {row['delta_repulsive_kJ_mol_H2O']:.3f} & "
            f"{row['delta_electrostatic_kJ_mol_H2O']:.3f} & {row['delta_dispersion_total_kJ_mol_H2O']:.3f} & "
            f"{row['delta_other_interactions_residual_kJ_mol_H2O']:.1f} / {row['pbc_zero_d_scc_steps']} \\\\"
        )
    component_lines.extend([r"\hline\hline", r"\end{tabular}", r"\end{table*}", ""])
    component_tex.write_text("\n".join(component_lines))

    legacy_tex = args.output_dir / "dmc_dilation_full_vs_legacy_table.tex"
    legacy_lines = [
        r"\begin{table}",
        r"\caption{Change from the legacy folded 2x2x2 route to the current coupled Brillouin-zone implementation. Energies are in kJ\,mol$^{-1}$\,H$_2$O$^{-1}$. The near-invariance of the 0D energy localizes the large change to the periodic route.}",
        r"\begin{tabular}{lrrrr}",
        r"\hline\hline",
        r"Phase & $s$ & Legacy $\Delta E$ & Current $\Delta E$ & 3D shift / 0D shift \\",
        r"\hline",
    ]
    for row in legacy_rows:
        legacy_lines.append(
            f"{row['phase']} & {row['scale']:g} & {row['legacy_folded_deltaE_kJ_mol_H2O']:.3f} & "
            f"{row['current_coupled_deltaE_kJ_mol_H2O']:.3f} & {row['current_minus_legacy_pbc_kJ_mol_H2O']:.3f} / "
            f"{row['current_minus_legacy_zero_d_kJ_mol_H2O']:.3f} \\\\"
        )
    legacy_lines.extend([r"\hline\hline", r"\end{tabular}", r"\end{table}", ""])
    legacy_tex.write_text("\n".join(legacy_lines))

    failure_audit = {
        "excluded_from_energy_tables": [
            {
                "campaign": "initial full-model legacy SCC_MIXER location",
                "count": 22,
                "reason": "post-#5582 input parser rejects SCC_MIXER inside TBLITE; corrected production is separate and complete",
            },
            {
                "campaign": "initial component smokes with long PARAM path",
                "count": 3,
                "reason": "CP2K input string exceeded 80 characters; identical parameters under short paths subsequently passed",
            },
            {
                "campaign": "no_acp",
                "case": "XIV scale 4 periodic",
                "reason": "TBLITE potential mixer did not converge in 250 iterations; no energy difference is reported",
            },
            {
                "campaign": "no_exchange_no_acp",
                "count": 11,
                "scope": "all periodic cases; all 11 corresponding 0D cases ended normally",
                "reason": "CP2K aborts before SCC step 1 with 'Missing complete g-xTB k-mesh map before potential mixing'; this is an integration-path assumption, not a model result",
            },
        ]
    }
    failure_json = args.output_dir / "campaign_failure_audit.json"
    failure_json.write_text(json.dumps(failure_audit, indent=2, sort_keys=True) + "\n")

    metadata = {
        "build": BUILD,
        "models": {
            "full": "unmodified exported g-xTB model",
            "no_exchange": "top-level exchange block absent; diagnostic only",
            "frozen_qvszp": "all coeffs_env entries zero; q-vSZP environment response frozen; diagnostic only",
            "no_anisotropic_multipole": "top-level multipole block absent; anisotropic multipole diagnostic only; periodic electrostatics/images remain",
            "no_acp": "top-level ACP model marker absent; ACP-free diagnostic only",
        },
        "selection": {
            "phases": ["Ih", "VII", "XIV"],
            "scales": {"Ih": [1.5, 2, 3, 4], "VII": [1.5, 2, 3], "XIV": [1.5, 2, 3, 4]},
            "upper_cell_vector_limit_angstrom": 40.0,
            "excluded": "VII scale 4: largest periodic vector 42.2194773541 A",
        },
        "route_audit": {
            "periodic": "CELL/POISSON PERIODIC XYZ/PERIODIC, Gamma-centered MP 2x2x2, FULL_GRID ON, SYMMETRY OFF",
            "zero_d": "CELL/POISSON PERIODIC NONE/NONE with ANALYTIC solver and no KPOINTS section; CP2K passes periodic=[False,False,False] to save_tblite, so save_tblite has no periodic lattice images",
        },
        "geometry_audit": geometry_records,
        "row_count": len(rows),
    }
    json_path = args.output_dir / "dmc_dilation_component_diagnostics.json"
    json_path.write_text(json.dumps({"metadata": metadata, "rows": rows}, indent=2, sort_keys=True) + "\n")

    readme = args.output_dir / "README.md"
    readme.write_text(
        "# DMC-ICE13 controlled 0D/3D cell-dilation diagnostic\n\n"
        "This archive compares the same 12 rigid water molecules as a true nonperiodic 0D cluster and as a 3D-periodic primitive cell. "
        "The 0D inputs use `PERIODIC NONE`, the analytic Poisson solver and no k-point section; CP2K therefore passes three false periodic flags to save_tblite. "
        "The periodic inputs use `PERIODIC XYZ`, periodic Poisson, and an unreduced Gamma-centered 2x2x2 mesh.\n\n"
        "All reported cell-vector lengths are <=40 A. VII at scale 4 is deliberately excluded because its largest periodic vector is 42.2194773541 A. "
        "The geometry audit aligns each periodic water by an integer lattice translation plus one global cluster translation and checks every O-H distance.\n\n"
        "The three modified parameter files are component-deletion diagnostics, not physical replacements or reparameterized models. In particular, removing the anisotropic multipole block does not disable all periodic electrostatics or image interactions. "
        "The residual component in the CSV is the total energy minus all CP2K components printed separately; in the full model it contains exchange and other tblite interactions and must not be labeled as pure exchange.\n\n"
        f"Build: CP2K `{BUILD['cp2k_revision']}`, save_tblite `{BUILD['save_tblite_revision']}`.\n"
    )

    manifest = args.output_dir / "SHA256SUMS"
    files = [
        csv_path,
        json_path,
        readme,
        compact_csv,
        component_csv,
        energy_tex,
        component_tex,
        legacy_csv,
        legacy_tex,
        failure_json,
    ]
    manifest.write_text("".join(f"{sha256(path)}  {path.name}\n" for path in sorted(files)))


if __name__ == "__main__":
    main()
