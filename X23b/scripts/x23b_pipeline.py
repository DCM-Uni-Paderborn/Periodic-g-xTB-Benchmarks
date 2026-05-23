#!/usr/bin/env python3
"""Prepare and analyse CP2K/tblite calculations for the X23b benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
HARTREE_TO_KJMOL = 2625.499638
BOHR_TO_ANG = 0.529177210903

SYSTEMS = [
    {"id": "14-cyclohexanedione", "label": "1,4-Cyclohexanedione", "ref_energy": 90.0, "ref_volume": 262.5},
    {"id": "acetic_acid", "label": "Acetic acid", "ref_energy": 73.6, "ref_volume": 288.8},
    {"id": "adamantane", "label": "Adamantane", "ref_energy": 71.8, "ref_volume": 357.6},
    {"id": "ammonia", "label": "Ammonia", "ref_energy": 38.7, "ref_volume": 121.5},
    {"id": "anthracene", "label": "Anthracene", "ref_energy": 110.4, "ref_volume": 441.2},
    {"id": "benzene", "label": "Benzene", "ref_energy": 54.8, "ref_volume": 444.3},
    {"id": "co2", "label": "Carbon dioxide", "ref_energy": 29.4, "ref_volume": 164.8},
    {"id": "cyanamide", "label": "Cyanamide", "ref_energy": 81.5, "ref_volume": 407.9},
    {"id": "cytosine", "label": "Cytosine", "ref_energy": 163.5, "ref_volume": 440.3},
    {"id": "ethylcarbamate", "label": "Ethyl carbamate", "ref_energy": 88.2, "ref_volume": 231.2},
    {"id": "formamide", "label": "Formamide", "ref_energy": 81.1, "ref_volume": 211.9},
    {"id": "hexamine", "label": "Hexamine", "ref_energy": 84.1, "ref_volume": 321.6, "ref_volume_molecules": 2},
    {"id": "imidazole", "label": "Imidazole", "ref_energy": 90.4, "ref_volume": 336.4},
    {"id": "naphthalene", "label": "Naphthalene", "ref_energy": 81.3, "ref_volume": 329.7},
    {"id": "oxalic_acid_alpha", "label": "Oxalic acid alpha", "ref_energy": 98.8, "ref_volume": 293.2},
    {"id": "oxalic_acid_beta", "label": "Oxalic acid beta", "ref_energy": 96.8, "ref_volume": 150.5},
    {"id": "pyrazine", "label": "Pyrazine", "ref_energy": 64.3, "ref_volume": 189.6},
    {"id": "pyrazole", "label": "Pyrazole", "ref_energy": 78.8, "ref_volume": 662.5},
    {"id": "succinic_acid", "label": "Succinic acid", "ref_energy": 130.1, "ref_volume": 233.3},
    {"id": "triazine", "label": "s-Triazine", "ref_energy": 62.6, "ref_volume": 528.0},
    {"id": "trioxane", "label": "s-Trioxane", "ref_energy": 64.6, "ref_volume": 580.7},
    {"id": "uracil", "label": "Uracil", "ref_energy": 136.2, "ref_volume": 442.0},
    {"id": "urea", "label": "Urea", "ref_energy": 102.1, "ref_volume": 140.8},
]

METHODS = ["GFN1", "GFN2"]
MESHES = [
    {"id": "gamma", "label": "Gamma", "scheme": None},
    {"id": "k111", "label": "1x1x1", "scheme": "MACDONALD 1 1 1 0.0 0.0 0.0"},
    {"id": "k222", "label": "2x2x2", "scheme": "MACDONALD 2 2 2 0.25 0.25 0.25"},
    {"id": "k333", "label": "3x3x3", "scheme": "MACDONALD 3 3 3 0.0 0.0 0.0"},
]

PLOT_LABELS = {
    "14-cyclohexanedione": "1,4-CHD",
    "acetic_acid": "AcOH",
    "adamantane": "Adam",
    "ammonia": "NH3",
    "anthracene": "Anth",
    "benzene": "Benz",
    "co2": "CO2",
    "cyanamide": "Cyan",
    "cytosine": "Cyt",
    "ethylcarbamate": "EtCarb",
    "formamide": "Form",
    "hexamine": "Hex",
    "imidazole": "Imid",
    "naphthalene": "Naph",
    "oxalic_acid_alpha": "Ox-a",
    "oxalic_acid_beta": "Ox-b",
    "pyrazine": "Pyz",
    "pyrazole": "Pyr",
    "succinic_acid": "Succ",
    "triazine": "Triaz",
    "trioxane": "Triox",
    "uracil": "Ura",
    "urea": "Urea",
}

DMC_X23 = {
    "14-cyclohexanedione": (88.3, 1.0),
    "acetic_acid": (71.7, 0.6),
    "adamantane": (61.0, 2.3),
    "ammonia": (38.2, 0.1),
    "anthracene": (100.2, 0.5),
    "benzene": (49.8, 0.2),
    "co2": (29.4, 0.2),
    "cyanamide": (83.6, 0.4),
    "cytosine": (156.2, 1.0),
    "ethylcarbamate": (84.2, 1.3),
    "formamide": (81.0, 1.0),
    "hexamine": (86.2, 0.6),
    "imidazole": (88.2, 0.8),
    "naphthalene": (75.5, 0.5),
    "oxalic_acid_alpha": (102.6, 1.4),
    "oxalic_acid_beta": (102.3, 0.6),
    "pyrazine": (61.1, 1.1),
    "pyrazole": (77.3, 0.5),
    "succinic_acid": (125.2, 0.5),
    "triazine": (60.5, 0.6),
    "trioxane": (62.1, 1.9),
    "uracil": (134.3, 0.7),
    "urea": (108.5, 0.3),
}


def clean_element(element: str) -> str:
    element = re.sub(r"[^A-Za-z]", "", element)
    return element[:1].upper() + element[1:].lower()


def cell_vectors(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> list[list[float]]:
    ar, br, gr = [math.radians(x) for x in (alpha, beta, gamma)]
    avec = [a, 0.0, 0.0]
    bvec = [b * math.cos(gr), b * math.sin(gr), 0.0]
    cx = c * math.cos(br)
    cy = c * (math.cos(ar) - math.cos(br) * math.cos(gr)) / math.sin(gr)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return [avec, bvec, [cx, cy, cz]]


def parse_cif(path: Path) -> dict[str, object]:
    values: dict[str, float] = {}
    atoms: list[dict[str, object]] = []
    lines = path.read_text().splitlines()
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0].startswith("_cell_"):
            try:
                values[parts[0]] = float(parts[1].strip("'\""))
            except ValueError:
                pass
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            headers: list[str] = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip())
                i += 1
            if "_atom_site_type_symbol" in headers and "_atom_site_fract_x" in headers:
                while i < len(lines):
                    line = lines[i].strip()
                    if not line or line.startswith("_") or line == "loop_" or line.startswith("data_"):
                        break
                    parts = line.split()
                    if len(parts) >= len(headers):
                        row = dict(zip(headers, parts))
                        atoms.append(
                            {
                                "element": clean_element(row["_atom_site_type_symbol"]),
                                "frac": [
                                    float(row["_atom_site_fract_x"]),
                                    float(row["_atom_site_fract_y"]),
                                    float(row["_atom_site_fract_z"]),
                                ],
                            }
                        )
                    i += 1
                continue
        i += 1
    cell = cell_vectors(
        values["_cell_length_a"],
        values["_cell_length_b"],
        values["_cell_length_c"],
        values["_cell_angle_alpha"],
        values["_cell_angle_beta"],
        values["_cell_angle_gamma"],
    )
    return {"cell": cell, "atoms": atoms, "volume": values.get("_cell_volume")}


def parse_qe_molecule(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    lines = path.read_text().splitlines()
    i = next(i for i, line in enumerate(lines) if line.strip().lower().startswith("atomic_positions"))
    i += 1
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.upper().startswith("K_POINTS") or line.startswith("&"):
            break
        parts = line.split()
        if len(parts) >= 4:
            atoms.append(
                {
                    "element": clean_element(parts[0]),
                    "coord": [float(parts[1]), float(parts[2]), float(parts[3])],
                }
            )
        i += 1
    return atoms


def parse_qe_crystal(path: Path) -> dict[str, object]:
    lines = path.read_text().splitlines()
    celldm = 1.0
    for line in lines:
        match = re.search(r"celldm\(1\)\s*=\s*([0-9.EDed+-]+)", line)
        if match:
            celldm = float(match.group(1).replace("D", "E").replace("d", "e"))
            break

    i = next(i for i, line in enumerate(lines) if line.strip().lower().startswith("atomic_positions"))
    atoms: list[dict[str, object]] = []
    i += 1
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.upper().startswith("K_POINTS"):
            break
        parts = line.split()
        if len(parts) >= 4:
            atoms.append(
                {
                    "element": clean_element(parts[0]),
                    "frac": [float(parts[1]) % 1.0, float(parts[2]) % 1.0, float(parts[3]) % 1.0],
                }
            )
        i += 1

    i = next(i for i, line in enumerate(lines) if line.strip().lower().startswith("cell_parameters"))
    header = lines[i].strip().lower()
    factor = 1.0
    if "bohr" in header or "cubic" in header or "alat" in header:
        factor = celldm * BOHR_TO_ANG
    cell = []
    for line in lines[i + 1 : i + 4]:
        cell.append([float(value) * factor for value in line.split()[:3]])
    return {"cell": cell, "atoms": atoms, "volume": abs(determinant_3x3(cell))}


def determinant_3x3(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def write_p1_cif(path: Path, title: str, geom: dict[str, object]) -> None:
    lines = [
        f"data_{title.replace('-', '_')}",
        f"_cell_volume           {float(geom['volume']):.6f}",
        "_symmetry_space_group_name_H-M 'P 1'",
        "_symmetry_Int_Tables_number 1",
        "loop_",
        "_symmetry_equiv_pos_site_id",
        "_symmetry_equiv_pos_as_xyz",
        "1 x,y,z",
    ]
    cell = geom["cell"]
    a = math.sqrt(sum(x * x for x in cell[0]))
    b = math.sqrt(sum(x * x for x in cell[1]))
    c = math.sqrt(sum(x * x for x in cell[2]))

    def angle(u: list[float], v: list[float]) -> float:
        dot = sum(x * y for x, y in zip(u, v))
        nu = math.sqrt(sum(x * x for x in u))
        nv = math.sqrt(sum(x * x for x in v))
        return math.degrees(math.acos(max(min(dot / (nu * nv), 1.0), -1.0)))

    lines += [
        f"_cell_length_a         {a:.10f}",
        f"_cell_length_b         {b:.10f}",
        f"_cell_length_c         {c:.10f}",
        f"_cell_angle_alpha      {angle(cell[1], cell[2]):.6f}",
        f"_cell_angle_beta       {angle(cell[0], cell[2]):.6f}",
        f"_cell_angle_gamma      {angle(cell[0], cell[1]):.6f}",
        "_cell_formula_units_Z 1",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    element_counts: dict[str, int] = {}
    for atom in geom["atoms"]:
        element = str(atom["element"])
        element_counts[element] = element_counts.get(element, 0) + 1
        x, y, z = atom["frac"]
        lines.append(f"{element}{element_counts[element]} {element} {x:.14f} {y:.14f} {z:.14f}")
    path.write_text("\n".join(lines) + "\n")


def count_by_element(atoms: list[dict[str, object]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for atom in atoms:
        element = str(atom["element"])
        counts[element] = counts.get(element, 0) + 1
    return counts


def molecules_per_cell(crystal_atoms: list[dict[str, object]], molecule_atoms: list[dict[str, object]]) -> int:
    crystal_counts = count_by_element(crystal_atoms, "frac")
    molecule_counts = count_by_element(molecule_atoms, "coord")
    ratios = []
    for element, count in molecule_counts.items():
        if element == "H" and element not in crystal_counts:
            continue
        if element not in crystal_counts or crystal_counts[element] % count != 0:
            raise ValueError(f"Cannot infer molecule count for {element}.")
        ratios.append(crystal_counts[element] // count)
    if len(set(ratios)) != 1:
        raise ValueError(f"Inconsistent molecule count ratios: {ratios}")
    return ratios[0]


def cp2k_header(project: str, run_type: str) -> list[str]:
    return [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        f"  RUN_TYPE {run_type}",
        "&END GLOBAL",
        "",
    ]


def cp2k_dft(method: str, mesh: dict[str, object] | None = None, periodic: bool = True) -> list[str]:
    lines = [
        "  &DFT",
        "    &QS",
        "      EPS_DEFAULT 1.0E-12",
        "      METHOD xTB",
        "      &XTB",
        "        GFN_TYPE TBLITE",
        "        &TBLITE",
        f"          METHOD {method}",
        "          ACCURACY 0.1",
        "        &END TBLITE",
        "      &END XTB",
        "    &END QS",
    ]
    if mesh and mesh.get("scheme"):
        lines += [
            "    &KPOINTS",
            f"      SCHEME {mesh['scheme']}",
            "      FULL_GRID T",
            "    &END KPOINTS",
        ]
    if not periodic:
        lines += [
            "    &POISSON",
            "      PERIODIC NONE",
            "    &END POISSON",
        ]
    lines += [
        "    &SCF",
        "      EPS_SCF 1.0E-9",
        "      MAX_SCF 300",
        "      SCF_GUESS MOPAC",
        "      &MIXING",
        "        METHOD DIRECT_P_MIXING",
        "        ALPHA 0.2",
        "      &END MIXING",
        "    &END SCF",
        "  &END DFT",
    ]
    return lines


def crystal_input(system: dict[str, object], geom: dict[str, object], method: str, mesh: dict[str, object], run_type: str) -> str:
    project = f"{system['id']}_{method}_{mesh['id']}_{run_type.lower()}".replace("-", "_")
    lines = cp2k_header(project, run_type)
    lines += [
        "&FORCE_EVAL",
        "  METHOD Quickstep",
        "  STRESS_TENSOR ANALYTICAL",
    ]
    lines += cp2k_dft(method, mesh if run_type == "ENERGY" else None, periodic=True)
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      PERIODIC XYZ",
    ]
    for name, vec in zip(("A", "B", "C"), geom["cell"]):
        lines.append(f"      {name} {vec[0]:.12f} {vec[1]:.12f} {vec[2]:.12f}")
    lines += [
        "    &END CELL",
        "    &COORD",
        "      SCALED",
    ]
    for atom in geom["atoms"]:
        x, y, z = atom["frac"]
        lines.append(f"      {atom['element']:<2} {x: .12f} {y: .12f} {z: .12f}")
    lines += [
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    if run_type == "CELL_OPT":
        lines += [
            "",
            "&MOTION",
        "  &CELL_OPT",
        "    OPTIMIZER BFGS",
        "    MAX_ITER 800",
        "    EXTERNAL_PRESSURE 0.0",
        "  &END CELL_OPT",
            "&END MOTION",
        ]
    return "\n".join(lines) + "\n"


def molecule_input(system: dict[str, object], atoms: list[dict[str, object]], method: str) -> str:
    project = f"{system['id']}_{method}_mol_geoopt".replace("-", "_")
    lines = cp2k_header(project, "GEO_OPT")
    lines += [
        "&FORCE_EVAL",
        "  METHOD Quickstep",
    ]
    lines += cp2k_dft(method, periodic=False)
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      ABC 30.0 30.0 30.0",
        "      PERIODIC NONE",
        "    &END CELL",
        "    &COORD",
    ]
    for atom in atoms:
        x, y, z = atom["coord"]
        lines.append(f"      {atom['element']:<2} {x: .12f} {y: .12f} {z: .12f}")
    lines += [
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
        "",
        "&MOTION",
        "  &GEO_OPT",
        "    OPTIMIZER BFGS",
        "    MAX_ITER 200",
        "  &END GEO_OPT",
        "&END MOTION",
    ]
    return "\n".join(lines) + "\n"


def prepare(refdata: Path) -> None:
    expt = refdata / "25_x23" / "expt"
    qe = refdata / "25_x23" / "b86bpbe-xdm"
    for directory in [DATA, FIGURES, ROOT / "structures" / "cif", ROOT / "structures" / "molecules_xyz"]:
        directory.mkdir(parents=True, exist_ok=True)
    metadata = {"systems": [], "methods": METHODS, "meshes": MESHES}
    for system in SYSTEMS:
        sid = str(system["id"])
        cif_path = expt / f"{sid}.cif"
        mol_path = qe / f"mol_{sid}.scf.in"
        if not mol_path.exists() and sid in {"oxalic_acid_alpha", "oxalic_acid_beta"}:
            mol_path = qe / "mol_oxalic_acid.scf.in"
        if sid == "hexamine":
            geom = parse_qe_crystal(qe / "hexamine.scf.in")
            source_note = "complete refdata X23 Quantum ESPRESSO crystal input"
        else:
            geom = parse_cif(cif_path)
            source_note = "refdata X23 experimental CIF"
        mol_atoms = parse_qe_molecule(mol_path)
        n_mol = molecules_per_cell(geom["atoms"], mol_atoms)
        ref_volume_molecules = int(system.get("ref_volume_molecules", n_mol))
        ref_volume_cell = float(system["ref_volume"]) * n_mol / ref_volume_molecules
        write_p1_cif(ROOT / "structures" / "cif" / cif_path.name, sid, geom)
        xyz = [str(len(mol_atoms)), f"{system['label']} gas-phase starting geometry from refdata"]
        xyz += [
            f"{atom['element']:<2} {atom['coord'][0]: .12f} {atom['coord'][1]: .12f} {atom['coord'][2]: .12f}"
            for atom in mol_atoms
        ]
        (ROOT / "structures" / "molecules_xyz" / f"{sid}.xyz").write_text("\n".join(xyz) + "\n")
        metadata["systems"].append(
            {
                **system,
                "n_atoms_crystal": len(geom["atoms"]),
                "n_atoms_molecule": len(mol_atoms),
                "molecules_per_cell": n_mol,
                "input_volume": geom["volume"],
                "x23b_reported_ref_volume": system["ref_volume"],
                "x23b_reported_ref_volume_molecules": ref_volume_molecules,
                "x23b_same_cell_ref_volume": ref_volume_cell,
                "structure_source": source_note,
            }
        )
        for method in METHODS:
            mol_dir = ROOT / "inputs" / "molecule_geoopt" / method
            mol_dir.mkdir(parents=True, exist_ok=True)
            (mol_dir / f"{sid}_{method}_mol_geoopt.inp").write_text(molecule_input(system, mol_atoms, method))
            for mesh in MESHES:
                sp_dir = ROOT / "inputs" / "crystal_sp" / str(mesh["id"]) / method
                sp_dir.mkdir(parents=True, exist_ok=True)
                (sp_dir / f"{sid}_{method}_{mesh['id']}_sp.inp").write_text(
                    crystal_input(system, geom, method, mesh, "ENERGY")
                )
            cell_dir = ROOT / "inputs" / "cellopt_gamma" / method
            cell_dir.mkdir(parents=True, exist_ok=True)
            (cell_dir / f"{sid}_{method}_gamma_cellopt.inp").write_text(
                crystal_input(system, geom, method, MESHES[0], "CELL_OPT")
            )
    (DATA / "metadata.json").write_text(json.dumps(metadata, indent=2))
    write_reference_csv(metadata)


def write_reference_csv(metadata: dict[str, object]) -> None:
    with (DATA / "x23b_reference.csv").open("w", newline="") as handle:
        fields = [
            "system",
            "label",
            "molecules_per_cell",
            "input_volume_A3",
            "x23b_reported_ref_volume_A3",
            "x23b_reported_ref_volume_molecules",
            "x23b_same_cell_ref_volume_A3",
            "x23b_ref_lattice_energy_kJmol",
            "structure_source",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for system in metadata["systems"]:
            writer.writerow(
                {
                    "system": system["id"],
                    "label": system["label"],
                    "molecules_per_cell": system["molecules_per_cell"],
                    "input_volume_A3": f"{float(system['input_volume']):.6f}",
                    "x23b_reported_ref_volume_A3": f"{float(system['x23b_reported_ref_volume']):.6f}",
                    "x23b_reported_ref_volume_molecules": system["x23b_reported_ref_volume_molecules"],
                    "x23b_same_cell_ref_volume_A3": f"{float(system['x23b_same_cell_ref_volume']):.6f}",
                    "x23b_ref_lattice_energy_kJmol": f"{float(system['ref_energy']):.6f}",
                    "structure_source": system["structure_source"],
                }
            )


def parse_energy(output: Path) -> float | None:
    if not output.exists():
        return None
    energy = None
    for line in output.read_text(errors="ignore").splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def completed_optimization(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text:
        return False
    return "GEOMETRY OPTIMIZATION COMPLETED" in text or "CELL OPTIMIZATION COMPLETED" in text


def completed_cp2k_run(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "ABORT" in text or re.search(r"SCF.*NOT|NOT.*SCF|DID NOT CONVERGE|convergence failure", text, re.I):
        return False
    return "PROGRAM ENDED" in text


def parse_last_volume(output: Path) -> float | None:
    if not output.exists():
        return None
    volume = None
    patterns = ["CELL| Volume [angstrom^3]", "CELL| Volume"]
    for line in output.read_text(errors="ignore").splitlines():
        if any(pattern in line for pattern in patterns):
            try:
                volume = float(line.split()[-1])
            except ValueError:
                pass
    return volume


def stats(errors: list[float]) -> dict[str, float]:
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(e) for e in errors) / len(errors),
        "RMSE": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "MaxAE": max(abs(e) for e in errors),
    }


def analyse() -> dict[str, object]:
    metadata = json.loads((DATA / "metadata.json").read_text())
    systems = metadata["systems"]
    rows_energy: list[dict[str, object]] = []
    rows_volume: list[dict[str, object]] = []
    results: dict[str, object] = {"methods": METHODS, "meshes": MESHES, "systems": systems, "energy_rows": [], "volume_rows": []}

    gas: dict[tuple[str, str], float | None] = {}
    for system in systems:
        for method in METHODS:
            stem = f"{system['id']}_{method}_mol_geoopt"
            output = ROOT / "runs" / "molecule_geoopt" / method / stem / f"{stem}.out"
            gas[(system["id"], method)] = parse_energy(output) if completed_optimization(output) else None

    for system in systems:
        n_mol = int(system["molecules_per_cell"])
        ref_energy = float(system["ref_energy"])
        for method in METHODS:
            gas_energy = gas[(system["id"], method)]
            for mesh in MESHES:
                stem = f"{system['id']}_{method}_{mesh['id']}_sp"
                output = ROOT / "runs" / "crystal_sp" / str(mesh["id"]) / method / stem / f"{stem}.out"
                crystal_energy = parse_energy(output) if completed_cp2k_run(output) else None
                complete = crystal_energy is not None and gas_energy is not None
                lattice = None
                error = None
                if complete:
                    lattice = (float(gas_energy) - float(crystal_energy) / n_mol) * HARTREE_TO_KJMOL
                    error = lattice - ref_energy
                rows_energy.append(
                    {
                        "calculation": "single_point",
                        "mesh": mesh["id"],
                        "method": f"{method}-xTB",
                        "system": system["id"],
                        "label": system["label"],
                        "complete": complete,
                        "lattice_energy_kJmol": "" if lattice is None else f"{lattice:.6f}",
                        "x23b_ref_lattice_energy_kJmol": f"{ref_energy:.6f}",
                        "error_kJmol": "" if error is None else f"{error:.6f}",
                    }
                )
            stem = f"{system['id']}_{method}_gamma_cellopt"
            output = ROOT / "runs" / "cellopt_gamma" / method / stem / f"{stem}.out"
            continuation = output.parent / "continue_800.out"
            if completed_optimization(continuation):
                output = continuation
            crystal_energy = parse_energy(output)
            volume = parse_last_volume(output)
            complete = completed_optimization(output) and crystal_energy is not None and gas_energy is not None and volume is not None
            lattice = None
            energy_error = None
            volume_error = None
            if complete:
                lattice = (float(gas_energy) - float(crystal_energy) / n_mol) * HARTREE_TO_KJMOL
                energy_error = lattice - ref_energy
                ref_volume = float(system["x23b_same_cell_ref_volume"])
                volume_error = 100.0 * (float(volume) - ref_volume) / ref_volume
            rows_energy.append(
                {
                    "calculation": "cell_opt",
                    "mesh": "gamma",
                    "method": f"{method}-xTB",
                    "system": system["id"],
                    "label": system["label"],
                    "complete": complete,
                    "lattice_energy_kJmol": "" if lattice is None else f"{lattice:.6f}",
                    "x23b_ref_lattice_energy_kJmol": f"{ref_energy:.6f}",
                    "error_kJmol": "" if energy_error is None else f"{energy_error:.6f}",
                }
            )
            rows_volume.append(
                {
                    "calculation": "cell_opt",
                    "mesh": "gamma",
                    "method": f"{method}-xTB",
                    "system": system["id"],
                    "label": system["label"],
                    "complete": complete,
                    "volume_A3": "" if volume is None else f"{float(volume):.6f}",
                    "x23b_same_cell_ref_volume_A3": f"{float(system['x23b_same_cell_ref_volume']):.6f}",
                    "x23b_reported_ref_volume_A3": f"{float(system['x23b_reported_ref_volume']):.6f}",
                    "volume_error_percent": "" if volume_error is None else f"{volume_error:.6f}",
                }
            )

    write_csv(
        DATA / "x23b_lattice_energies.csv",
        rows_energy,
        [
            "calculation",
            "mesh",
            "method",
            "system",
            "label",
            "complete",
            "lattice_energy_kJmol",
            "x23b_ref_lattice_energy_kJmol",
            "error_kJmol",
        ],
    )
    write_csv(
        DATA / "x23b_cell_volumes.csv",
        rows_volume,
        [
            "calculation",
            "mesh",
            "method",
            "system",
            "label",
            "complete",
            "volume_A3",
            "x23b_same_cell_ref_volume_A3",
            "x23b_reported_ref_volume_A3",
            "volume_error_percent",
        ],
    )
    summaries = summarize(rows_energy, rows_volume)
    write_csv(DATA / "x23b_summary.csv", summaries, ["quantity", "calculation", "mesh", "method", "ME", "MAE", "RMSE", "MaxAE"])
    results["energy_rows"] = rows_energy
    results["volume_rows"] = rows_volume
    results["summary"] = summaries
    (DATA / "x23b_results.json").write_text(json.dumps(results, indent=2))
    make_plots(summaries, rows_energy)
    return results


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows_energy: list[dict[str, object]], rows_volume: list[dict[str, object]]) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for calculation in ["single_point", "cell_opt"]:
        meshes = sorted({str(row["mesh"]) for row in rows_energy if row["calculation"] == calculation})
        for mesh in meshes:
            for method in [f"{name}-xTB" for name in METHODS]:
                errors = [
                    float(row["error_kJmol"])
                    for row in rows_energy
                    if row["calculation"] == calculation
                    and row["mesh"] == mesh
                    and row["method"] == method
                    and row["error_kJmol"] != ""
                ]
                if errors:
                    summaries.append(
                        {
                            "quantity": "lattice_energy_kJmol",
                            "calculation": calculation,
                            "mesh": mesh,
                            "method": method,
                            **{key: f"{value:.6f}" for key, value in stats(errors).items()},
                        }
                    )
    for method in [f"{name}-xTB" for name in METHODS]:
        errors = [
            float(row["volume_error_percent"])
            for row in rows_volume
            if row["method"] == method and row["volume_error_percent"] != ""
        ]
        if errors:
            summaries.append(
                {
                    "quantity": "volume_error_percent",
                    "calculation": "cell_opt",
                    "mesh": "gamma",
                    "method": method,
                    **{key: f"{value:.6f}" for key, value in stats(errors).items()},
                }
            )
    return summaries


def make_plots(summaries: list[dict[str, object]], rows_energy: list[dict[str, object]]) -> None:
    if not summaries or shutil.which("gnuplot") is None:
        return
    FIGURES.mkdir(exist_ok=True)
    dat = DATA / "x23b_summary_for_plot.dat"
    order = [
        ("single_point", "gamma", "GFN1-xTB"),
        ("single_point", "gamma", "GFN2-xTB"),
        ("single_point", "k111", "GFN1-xTB"),
        ("single_point", "k111", "GFN2-xTB"),
        ("single_point", "k222", "GFN1-xTB"),
        ("single_point", "k222", "GFN2-xTB"),
        ("single_point", "k333", "GFN1-xTB"),
        ("single_point", "k333", "GFN2-xTB"),
        ("cell_opt", "gamma", "GFN1-xTB"),
        ("cell_opt", "gamma", "GFN2-xTB"),
    ]
    lookup = {
        (row["quantity"], row["calculation"], row["mesh"], row["method"]): row
        for row in summaries
    }
    with dat.open("w") as handle:
        handle.write("# index label lattice_MAE volume_MAE\n")
        for index, (calc, mesh, method) in enumerate(order, start=1):
            energy = lookup.get(("lattice_energy_kJmol", calc, mesh, method))
            volume = lookup.get(("volume_error_percent", calc, mesh, method))
            label = f"{method.replace('-xTB','')} {mesh}" if calc == "single_point" else f"{method.replace('-xTB','')} opt"
            handle.write(f'{index} "{label}" {energy["MAE"] if energy else "NaN"} {volume["MAE"] if volume else "NaN"}\n')
    script = f"""
set terminal svg enhanced font 'Helvetica,12' size 920,520
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set output '{FIGURES / 'x23b_mae_summary.svg'}'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set key top right spacing 1.2 samplen 2
set ylabel 'MAE'
set yrange [0:*]
set xtics rotate by -35
set style data histogram
set style histogram clustered gap 1
set style fill solid 0.85 border -1
set boxwidth 0.75
plot '{dat}' using 3:xtic(2) lc rgb '#4c72b0' title 'Lattice energy / kJ mol^{{-1}}', \\
     '' using 4 lc rgb '#dd8452' title 'Cell volume / %'
"""
    subprocess.run(["gnuplot"], input=script.encode(), check=True)
    svg = FIGURES / "x23b_mae_summary.svg"
    if shutil.which("rsvg-convert") is not None:
        subprocess.run(["rsvg-convert", str(svg), "-o", str(svg.with_suffix(".png"))], check=True)
        subprocess.run(["rsvg-convert", "-f", "pdf", str(svg), "-o", str(svg.with_suffix(".pdf"))], check=True)
    make_prl_style_plot(rows_energy)


def make_prl_style_plot(rows_energy: list[dict[str, object]]) -> None:
    lookup = {
        (row["system"], row["method"]): float(row["error_kJmol"])
        for row in rows_energy
        if row["calculation"] == "cell_opt"
        and row["mesh"] == "gamma"
        and str(row["complete"]) == "True"
        and row["error_kJmol"] != ""
    }
    systems = sorted(SYSTEMS, key=lambda item: float(item["ref_energy"]))
    dat = DATA / "x23b_prl_style.dat"
    with dat.open("w") as handle:
        handle.write("# index label x23b_ref dmc dmc_err dmc_minus_x23b gfn1_minus_x23b gfn2_minus_x23b\n")
        for index, system in enumerate(systems, start=1):
            system_id = str(system["id"])
            ref = float(system["ref_energy"])
            dmc, dmc_err = DMC_X23[system_id]
            gfn1 = lookup.get((system_id, "GFN1-xTB"), float("nan"))
            gfn2 = lookup.get((system_id, "GFN2-xTB"), float("nan"))
            label = PLOT_LABELS[system_id]
            handle.write(
                f'{index} "{label}" {ref:.6f} {dmc:.6f} {dmc_err:.6f} '
                f"{dmc - ref:.6f} {gfn1:.6f} {gfn2:.6f}\n"
            )

    svg = FIGURES / "x23b_lattice_energy_prl_style.svg"
    script = f"""
set terminal svg enhanced font 'Helvetica,12' size 1120,760
set output '{svg}'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d8d8d8' lw 0.6
set xrange [0.5:23.5]
set lmargin 10
set rmargin 3
set multiplot layout 2,1

set tmargin 2
set bmargin 1
set object 100 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
unset xtics
set ylabel 'Lattice-energy magnitude / kJ mol^{-1}'
set yrange [20:170]
set key top left spacing 1.15 samplen 2
plot '{dat}' using 1:3 with linespoints lt 1 lw 1.4 pt 7 ps 0.7 lc rgb '#222222' title 'X23b reference', \\
     '' using 1:4:5 with yerrorbars pt 7 ps 0.75 lw 1.2 lc rgb '#4C78A8' title 'DMC X23'

unset object 100
set tmargin 1
set bmargin 8
set xtics rotate by -55 font 'Helvetica,10'
set ylabel 'Deviation from X23b / kJ mol^{-1}'
set yrange [-70:235]
set ytics 50
set yzeroaxis lt -1 lw 1.0 lc rgb '#555555'
set object 1 rectangle from graph 0, first -4.184 to graph 1, first 4.184 fillcolor rgb '#e6e6e6' behind
set key top left spacing 1.1 samplen 1.8
plot '{dat}' using 1:6:5 with yerrorbars pt 7 ps 0.75 lw 1.1 lc rgb '#4C78A8' title 'DMC X23 - X23b', \\
     '' using 1:7:xtic(2) with linespoints lt 1 lw 1.4 pt 5 ps 0.8 lc rgb '#E45756' title 'GFN1-xTB opt - X23b', \\
     '' using 1:8 with linespoints lt 1 lw 1.4 pt 9 ps 0.8 lc rgb '#54A24B' title 'GFN2-xTB opt - X23b'
unset multiplot
"""
    subprocess.run(["gnuplot"], input=script.encode(), check=True)
    if shutil.which("rsvg-convert") is not None:
        subprocess.run(["rsvg-convert", str(svg), "-o", str(svg.with_suffix(".png"))], check=True)
        subprocess.run(["rsvg-convert", "-f", "pdf", str(svg), "-o", str(svg.with_suffix(".pdf"))], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prepare", "analyse", "all"], nargs="?", default="all")
    parser.add_argument("--refdata", type=Path, default=Path(os.environ.get("REFDATA_X23", "../refdata")))
    args = parser.parse_args()
    if args.command in {"prepare", "all"}:
        prepare(args.refdata)
    if args.command in {"analyse", "all"}:
        analyse()


if __name__ == "__main__":
    main()
