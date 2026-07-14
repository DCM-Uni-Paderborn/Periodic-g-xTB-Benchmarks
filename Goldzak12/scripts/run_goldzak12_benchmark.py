#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GXTB_CAMPAIGN_MANIFEST = (
    ROOT.parent / "campaigns" / "gxtb-pbc-v1-20260714" / "build_manifest.json"
)
DEFAULT_CP2K = Path(os.environ.get("CP2K", "cp2k.ssmp"))
DEFAULT_TBLITE = Path(os.environ.get("TBLITE", "tblite"))
DEFAULT_SAVE_TBLITE = Path(os.environ.get("SAVE_TBLITE", os.environ.get("TBLITE", "tblite")))
DEFAULT_CP2K_SOURCE = Path(os.environ.get("CP2K_SOURCE", "../cp2k"))
DEFAULT_TBLITE_SOURCE = Path(os.environ.get("TBLITE_SOURCE", "../tblite"))
DEFAULT_SAVE_TBLITE_SOURCE = Path(os.environ.get("SAVE_TBLITE_SOURCE", "../save_tblite"))
HARTREE_TO_EV = 27.211386245988
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
CAMPAIGN_SCHEMA = 1
CAMPAIGN_IDENTITY_FIELDS = (
    "schema",
    "campaign_id",
    "cp2k_executable_sha256",
    "cp2k_loaded_library_sha256",
    "cp2k_cmake_cache_sha256",
    "cp2k_embedded_source_revision",
    "cp2k_source_revision",
    "save_tblite_executable_sha256",
    "save_tblite_source_revision",
    "save_tblite_library_sha256",
    "save_tblite_cmake_cache_sha256",
    "dependency_lock_sha256",
)


@dataclass(frozen=True)
class Reference:
    solid: str
    structure: str
    formula: tuple[tuple[str, int], ...]
    a_exp: float
    a_hf: float
    a_mp2: float
    a_scs_mp2: float
    a_sos_mp2: float
    ecoh_exp: float
    ecoh_hf: float
    ecoh_mp2: float
    ecoh_scs_mp2: float
    ecoh_sos_mp2: float


REFERENCES: tuple[Reference, ...] = (
    Reference("C", "diamond", (("C", 1),), 3.553, 3.547, 3.540, 3.550, 3.554, 7.55, 5.38, 7.98, 7.65, 7.50),
    Reference("Si", "diamond", (("Si", 1),), 5.421, 5.508, 5.399, 5.425, 5.437, 4.70, 3.03, 4.97, 4.69, 4.56),
    Reference("SiC", "zincblende", (("Si", 1), ("C", 1)), 4.347, 4.371, 4.350, 4.358, 4.362, 6.47, 4.53, 6.79, 6.49, 6.35),
    Reference("BN", "zincblende", (("B", 1), ("N", 1)), 3.593, 3.596, 3.596, 3.603, 3.606, 6.76, 4.78, 7.13, 6.92, 6.82),
    Reference("BP", "zincblende", (("B", 1), ("P", 1)), 4.525, 4.584, 4.495, 4.517, 4.528, 5.14, 3.42, 5.58, 5.30, 5.16),
    Reference("AlN", "zincblende", (("Al", 1), ("N", 1)), 4.368, 4.365, 4.388, 4.389, 4.389, 5.85, 3.86, 6.00, 5.85, 5.78),
    Reference("AlP", "zincblende", (("Al", 1), ("P", 1)), 5.448, 5.542, 5.444, 5.465, 5.475, 4.31, 2.71, 4.42, 4.23, 4.14),
    Reference("MgO", "rocksalt", (("Mg", 1), ("O", 1)), 4.189, 4.176, 4.227, 4.224, 4.223, 5.19, 3.62, 5.37, 5.18, 5.09),
    Reference("MgS", "rocksalt", (("Mg", 1), ("S", 1)), 5.188, 5.281, 5.171, 5.191, 5.201, 4.04, 2.78, 4.20, 3.97, 3.86),
    Reference("LiH", "rocksalt", (("Li", 1), ("H", 1)), 3.979, 4.094, 3.996, 4.009, 4.015, 2.49, 1.85, 2.41, 2.45, 2.48),
    Reference("LiF", "rocksalt", (("Li", 1), ("F", 1)), 3.973, 3.964, 3.990, 3.992, 3.993, 4.46, 3.41, 4.58, 4.49, 4.44),
    Reference("LiCl", "rocksalt", (("Li", 1), ("Cl", 1)), 5.072, 5.253, 5.021, 5.059, 5.078, 3.58, 2.73, 3.69, 3.58, 3.52),
)

METHODS = ("GFN1", "GFN2", "GXTB")
LEGACY_METHODS = ("GFN1", "GFN2")
METHOD_COLORS = {"GFN1": "#4C78A8", "GFN2": "#F58518", "GXTB": "#54A24B"}
KPOINT_MESH_CONTRACT = (
    "CP2K native Bloch MACDONALD meshes with SPGLIB reduction for GFN1, GFN2, and GXTB "
    "(SYMMETRY T, FULL_GRID F, SYMMETRY_BACKEND SPGLIB, "
    "SYMMETRY_REDUCTION_METHOD SPGLIB); the GXTB interface expands the irreducible data "
    "to the complete mesh before save_tblite and folds the response back afterwards"
)
LEGACY_GXTB_FULL_GRID_POLICY = (
    "pre-SPGLIB GXTB full-grid inputs and outputs are diagnostics only and are never "
    "accepted as LC12 production results"
)
GXTB_ENERGY_STRESS_POLICY = (
    "LC12 GXTB EOS, final, and isolated-atom ENERGY inputs do not request a stress tensor; "
    "GFN1/GFN2 frozen inputs retain STRESS_TENSOR ANALYTICAL"
)
GXTB_ATOM_SCF_POLICY = (
    "LC12 cohesive energies use save_tblite CLI atom energies. The independent CP2K/CLI "
    "interface gate uses CP2K's supported nonperiodic Gamma/no-smear OT path, where "
    "SCC_MIXER is ignored. This avoids CP2K's otherwise forced 300 K tblite smearing, which "
    "changes several isolated-atom states and cannot be represented at all when a minimal "
    "atomic basis has no virtual MO. Li is the sole documented exception: OT converges to "
    "an atomic state 0.077782 Eh above the save_tblite CLI state, while native g-XTB FDIIS "
    "with diagonalization reproduces the CLI energy."
)
ELEMENT_MULTIPLICITY = {
    "H": 2,
    "Li": 2,
    "B": 2,
    "C": 3,
    "N": 4,
    "O": 3,
    "F": 2,
    "Mg": 1,
    "Al": 2,
    "Si": 3,
    "P": 4,
    "S": 3,
    "Cl": 2,
}
GXTB_NO_SMEAR_OT_ATOMS = frozenset(ELEMENT_MULTIPLICITY) - {"Li"}


def selected_methods(requested: list[str] | tuple[str, ...] | None) -> tuple[str, ...]:
    """Validate and de-duplicate a CLI method selection in paper order."""
    if not requested:
        return LEGACY_METHODS
    normalized = {method.upper() for method in requested}
    unknown = sorted(normalized - set(METHODS))
    if unknown:
        raise ValueError(f"Unknown method(s): {', '.join(unknown)}")
    return tuple(method for method in METHODS if method in normalized)


def method_cli_name(method: str) -> str:
    return method.lower()


def fcc_sites() -> list[tuple[float, float, float]]:
    return [
        (0.0, 0.0, 0.0),
        (0.0, 0.5, 0.5),
        (0.5, 0.0, 0.5),
        (0.5, 0.5, 0.0),
    ]


def frac_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple((a[i] + b[i]) % 1.0 for i in range(3))  # type: ignore[return-value]


def conventional_cell_atoms(ref: Reference) -> list[tuple[str, float, float, float]]:
    sites = fcc_sites()
    if ref.structure == "diamond":
        element = ref.formula[0][0]
        return [(element, *p) for p in sites] + [(element, *frac_add(p, (0.25, 0.25, 0.25))) for p in sites]
    if ref.structure == "zincblende":
        a, b = ref.formula[0][0], ref.formula[1][0]
        return [(a, *p) for p in sites] + [(b, *frac_add(p, (0.25, 0.25, 0.25))) for p in sites]
    if ref.structure == "rocksalt":
        a, b = ref.formula[0][0], ref.formula[1][0]
        return [(a, *p) for p in sites] + [(b, *frac_add(p, (0.5, 0.0, 0.0))) for p in sites]
    raise ValueError(ref.structure)


def atom_counts(ref: Reference) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, *_ in conventional_cell_atoms(ref):
        counts[element] = counts.get(element, 0) + 1
    return counts


def kpoint_block(mesh: str, method: str) -> list[str]:
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}")
    if not mesh.startswith("k") or not mesh[1:].isdigit():
        raise ValueError(f"Bad mesh {mesh!r}; expected k333, k444, ...")
    digits = mesh[1:]
    if len(digits) == 3 and len(set(digits)) == 1:
        n = int(digits[0])
    else:
        n = int(digits)
    shift = 0.0 if n % 2 else 1.0 / (2.0 * n)
    lines = [
        "    &KPOINTS",
        f"      SCHEME MACDONALD {n} {n} {n} {shift:.10g} {shift:.10g} {shift:.10g}",
        "      SYMMETRY T",
        "      FULL_GRID F",
        "      SYMMETRY_BACKEND SPGLIB",
        "      SYMMETRY_REDUCTION_METHOD SPGLIB",
    ]
    return lines + ["    &END KPOINTS"]


def validate_method_input(
    text: str, method: str, *, gxtb_atom_reference: bool = False
) -> None:
    """Reject an LC12 input that can silently use the wrong model or k mesh."""
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}")
    if gxtb_atom_reference and method != "GXTB":
        raise ValueError("the GXTB isolated-atom SCF exception is only valid for METHOD GXTB")
    tblite_methods = re.findall(r"^\s*METHOD\s+(GFN1|GFN2|GXTB)\s*$", text, flags=re.I | re.M)
    if not tblite_methods or tblite_methods[-1].upper() != method:
        raise ValueError(f"input does not select METHOD {method}")
    if method != "GXTB":
        return
    global_section = re.search(
        r"^\s*&GLOBAL\b.*?^\s*&END\s+GLOBAL\s*$", text, flags=re.I | re.M | re.S
    )
    if global_section is not None and re.search(
        r"^\s*BACKUP_COPIES\b", global_section.group(0), flags=re.I | re.M
    ):
        raise ValueError(
            "BACKUP_COPIES is not valid in GLOBAL; place it in the relevant PRINT/RESTART section"
        )
    mixers = [
        value.upper()
        for value in re.findall(r"^\s*SCC_MIXER\s+(\S+)\s*$", text, flags=re.I | re.M)
    ]
    expected_mixer = "CP2K" if gxtb_atom_reference else "TBLITE"
    if mixers != [expected_mixer]:
        raise ValueError(
            f"GXTB {'isolated-atom' if gxtb_atom_reference else 'production'} inputs require "
            f"exactly one 'SCC_MIXER {expected_mixer}'"
        )
    if re.search(r"^\s*&TBLITE_MIXER\b", text, flags=re.I | re.M):
        raise ValueError("GXTB production inputs must not override the native save_tblite FDIIS mixer")
    run_types = [
        value.upper()
        for value in re.findall(r"^\s*RUN_TYPE\s+(\S+)\s*$", text, flags=re.I | re.M)
    ]
    stress_requests = [
        value.upper()
        for value in re.findall(r"^\s*STRESS_TENSOR\s+(\S+)\s*$", text, flags=re.I | re.M)
    ]
    if run_types == ["ENERGY"] and any(value != "NONE" for value in stress_requests):
        raise ValueError("GXTB LC12 ENERGY inputs must not request analytical stress")
    if gxtb_atom_reference:
        atom_contract = (
            run_types == ["ENERGY"],
            re.search(r"^\s*&OT\b", text, flags=re.I | re.M) is not None,
            re.search(r"^\s*&SMEAR\s+OFF\b", text, flags=re.I | re.M) is not None,
            re.search(r"^\s*PERIODIC\s+NONE\s*$", text, flags=re.I | re.M) is not None,
            re.search(r"^\s*&KPOINTS\b", text, flags=re.I | re.M) is None,
        )
        if not all(atom_contract):
            raise ValueError(
                "GXTB isolated-atom references require ENERGY, nonperiodic Gamma-point OT "
                "without smearing"
            )
    kpoints = re.search(
        r"^\s*&KPOINTS\b.*?^\s*&END\s+KPOINTS\s*$", text, flags=re.I | re.M | re.S
    )
    if kpoints is not None:
        block = kpoints.group(0)
        symmetries = re.findall(r"^\s*SYMMETRY\s+(\S+)\s*$", block, flags=re.I | re.M)
        full_grids = re.findall(r"^\s*FULL_GRID\s+(\S+)\s*$", block, flags=re.I | re.M)
        backends = re.findall(r"^\s*SYMMETRY_BACKEND\s+(\S+)\s*$", block, flags=re.I | re.M)
        reductions = re.findall(
            r"^\s*SYMMETRY_REDUCTION_METHOD\s+(\S+)\s*$", block, flags=re.I | re.M
        )
        required = (
            [value.upper() for value in symmetries] == ["T"],
            [value.upper() for value in full_grids] == ["F"],
            [value.upper() for value in backends] == ["SPGLIB"],
            [value.upper() for value in reductions] == ["SPGLIB"],
        )
        if not all(required):
            raise ValueError(
                "GXTB k-point production inputs require the SPGLIB-reduced mesh contract "
                "(SYMMETRY T, FULL_GRID F, SPGLIB backend and reduction method)"
            )


def quickstep_block(
    method: str,
    periodic: bool,
    mesh: str | None = None,
    multiplicity: int | None = None,
    added_mos: int | None = None,
    smear_off: bool = False,
    request_stress: bool = True,
    gxtb_atom_reference: bool = False,
) -> list[str]:
    if gxtb_atom_reference and (
        method != "GXTB" or periodic or mesh is not None or not smear_off or added_mos is not None
    ):
        raise ValueError(
            "the GXTB isolated-atom SCF exception requires nonperiodic Gamma-point GXTB, "
            "SMEAR OFF, and no ADDED_MOS"
        )
    lines = ["  METHOD Quickstep"]
    if request_stress:
        lines.append("  STRESS_TENSOR ANALYTICAL")
    lines.append("  &DFT")
    if multiplicity is not None:
        lines += [
            "    UKS T",
            f"    MULTIPLICITY {multiplicity}",
        ]
    lines += [
        "    &QS",
        "      EPS_DEFAULT 1.0E-12",
        "      METHOD xTB",
        "      &XTB",
        "        GFN_TYPE TBLITE",
    ]
    if method == "GXTB":
        # AUTO currently resolves to the same provider-native path, but the
        # explicit selection is part of the paper protocol and provenance.
        lines.append(
            "        SCC_MIXER CP2K"
            if gxtb_atom_reference
            else "        SCC_MIXER TBLITE"
        )
    lines += [
        "        &TBLITE",
        f"          METHOD {method}",
        "          ACCURACY 0.05",
        "        &END TBLITE",
        "      &END XTB",
        "    &END QS",
    ]
    if not periodic:
        lines += [
            "    &POISSON",
            "      PERIODIC NONE",
            "    &END POISSON",
        ]
    if mesh:
        lines += kpoint_block(mesh, method)
    lines += [
        "    &SCF",
        "      EPS_SCF 1.0E-9",
        "      MAX_SCF 300",
        "      SCF_GUESS MOPAC",
    ]
    if added_mos is not None:
        lines.append(f"      ADDED_MOS {added_mos}")
    if smear_off:
        lines += [
            "      &SMEAR OFF",
            "      &END SMEAR",
        ]
    if gxtb_atom_reference:
        # CP2K otherwise re-enables its 300 K tblite smearing after parsing
        # SMEAR OFF.  This changes several isolated-atom states, while fully
        # occupied minimal bases (H/N/O/F) have no extra MO and abort outright.
        # Gamma-point OT is CP2K's supported no-smear tblite path; SCC_MIXER
        # CP2K prevents the diagonalization-only native g-XTB potential mixer
        # from being selected for this atomic exception.
        lines += [
            "      &OT",
            "        MINIMIZER DIIS",
            "        PRECONDITIONER FULL_SINGLE_INVERSE",
            "      &END OT",
        ]
    else:
        lines += [
            "      &MIXING",
            "        METHOD DIRECT_P_MIXING",
            "        ALPHA 0.2",
            "      &END MIXING",
        ]
    lines += [
        "      &PRINT",
        "        &RESTART OFF",
    ]
    if method == "GXTB":
        lines.append("          BACKUP_COPIES 0")
    lines += [
        "        &END RESTART",
        "      &END PRINT",
        "    &END SCF",
        "  &END DFT",
    ]
    return lines


def solid_input(ref: Reference, method: str, run_type: str, mesh: str, lattice_a: float, project: str) -> str:
    atoms = conventional_cell_atoms(ref)
    lines = [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        f"  RUN_TYPE {run_type}",
    ]
    lines += ["&END GLOBAL", "", "&FORCE_EVAL"]
    lines += quickstep_block(
        method,
        periodic=True,
        mesh=mesh,
        request_stress=not (method == "GXTB" and run_type == "ENERGY"),
    )
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      CANONICALIZE TRUE",
        f"      ABC {lattice_a:.12f} {lattice_a:.12f} {lattice_a:.12f}",
        "      PERIODIC XYZ",
        "      SYMMETRY CUBIC",
        "    &END CELL",
        "    &COORD",
        "      SCALED",
    ]
    for element, x, y, z in atoms:
        lines.append(f"      {element:<2} {x: .12f} {y: .12f} {z: .12f}")
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
            "    MAX_ITER 160",
            "    EXTERNAL_PRESSURE [bar] 0.0",
            "    KEEP_ANGLES T",
            "    KEEP_SYMMETRY T",
            "    MAX_DR 2.0E-3",
            "    RMS_DR 1.0E-3",
            "    MAX_FORCE 6.0E-4",
            "    RMS_FORCE 3.0E-4",
            "    PRESSURE_TOLERANCE [bar] 150.0",
            "    &BFGS",
            "      TRUST_RADIUS [angstrom] 0.05",
            "    &END BFGS",
            "  &END CELL_OPT",
            "&END MOTION",
        ]
    text = "\n".join(lines) + "\n"
    validate_method_input(text, method)
    return text


def atom_input(element: str, method: str) -> str:
    multiplicity = ELEMENT_MULTIPLICITY[element]
    project = f"atom_{element}_{method}"
    gxtb_atom_reference = method == "GXTB" and element in GXTB_NO_SMEAR_OT_ATOMS
    lines = [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        "  RUN_TYPE ENERGY",
    ]
    lines += ["&END GLOBAL", "", "&FORCE_EVAL"]
    lines += quickstep_block(
        method,
        periodic=False,
        multiplicity=multiplicity,
        smear_off=True,
        request_stress=method != "GXTB",
        gxtb_atom_reference=gxtb_atom_reference,
    )
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      ABC 30.0 30.0 30.0",
        "      PERIODIC NONE",
        "    &END CELL",
        "    &COORD",
        f"      {element:<2} 0.0 0.0 0.0",
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    text = "\n".join(lines) + "\n"
    validate_method_input(text, method, gxtb_atom_reference=gxtb_atom_reference)
    return text


def write_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def project_name(solid: str, method: str, kind: str, mesh: str) -> str:
    clean_solid = solid.replace("/", "_")
    return f"{clean_solid}_{method}_{kind}_{mesh}"


def run_cp2k(cp2k: Path, inp: Path, out: Path, threads: int) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "false"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    # CP2K/OpenMP owns the outer parallelism.  Nested BLAS threads otherwise
    # oversubscribe every concurrent LC12 job and make timings/non-convergence
    # unnecessarily machine-dependent.
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["VECLIB_MAXIMUM_THREADS"] = "1"
    main_log = inp.parent / "mainLog.out"
    if main_log.exists():
        main_log.unlink()
    proc = subprocess.run(
        [str(cp2k), "-i", str(inp.name), "-o", str(out.name)],
        cwd=inp.parent,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env=env,
    )
    if main_log.exists() and (not out.exists() or "PROGRAM ENDED" not in out.read_text(errors="ignore")):
        shutil.copyfile(main_log, out)
    return proc.returncode


def output_ok(output: Path, require_opt: bool = False) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "PROGRAM ENDED" not in text:
        return False
    bad = ("ABORT", "DID NOT CONVERGE", "SCF run NOT converged")
    if any(token in text for token in bad):
        return False
    if require_opt and "CELL OPTIMIZATION COMPLETED" not in text and "GEOMETRY OPTIMIZATION COMPLETED" not in text:
        return False
    return True


GXTB_TRANSIENT_PATTERNS = (
    "*-RESTART.kp",
    "*-RESTART.kp.bak*",
    "*-RESTART.wfn",
    "*-RESTART.wfn.bak*",
    "*.wfn.bak*",
)


def prune_gxtb_transients(search_roots: tuple[Path, ...] | None = None) -> tuple[int, int]:
    """Remove only large GXTB restart transients beside a validated output."""
    runs_root = (ROOT / "runs").resolve()
    roots = search_roots or (
        ROOT / "runs" / "eos" / "GXTB",
        ROOT / "runs" / "eos_final_sp" / "GXTB",
        ROOT / "runs" / "single_point" / "GXTB",
    )
    removed_files = 0
    removed_bytes = 0
    for root in roots:
        resolved_root = root.resolve()
        if not resolved_root.is_relative_to(runs_root) or "GXTB" not in resolved_root.parts:
            raise ValueError(f"Refusing to prune outside a GXTB run tree: {root}")
        if not root.exists():
            continue
        for directory in sorted({path.parent for path in root.rglob("*.out")}):
            if not any(output_ok(output) for output in directory.glob("*.out")):
                continue
            candidates: set[Path] = set()
            for pattern in GXTB_TRANSIENT_PATTERNS:
                candidates.update(path for path in directory.glob(pattern) if path.is_file())
            for candidate in sorted(candidates):
                relative = candidate.resolve().relative_to(runs_root)
                if "GXTB" not in relative.parts:
                    raise ValueError(f"Refusing to prune non-GXTB file: {candidate}")
                removed_bytes += candidate.stat().st_size
                candidate.unlink()
                removed_files += 1
    return removed_files, removed_bytes


def parse_energy(output: Path) -> float | None:
    extrapolated_energy = None
    force_eval_energy = None
    if not output.exists():
        return None
    for line in output.read_text(errors="ignore").splitlines():
        if "Total energy (extrapolated to T->0)" in line:
            extrapolated_energy = float(line.split()[-1])
        if "ENERGY| Total FORCE_EVAL" in line:
            force_eval_energy = float(line.split()[-1])
    if extrapolated_energy is not None:
        return extrapolated_energy
    return force_eval_energy


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_executable(executable: Path) -> Path | None:
    """Resolve an executable without pretending an unhashable command is safe."""
    if executable.is_file():
        return executable.resolve()
    located = shutil.which(str(executable))
    if located:
        candidate = Path(located).resolve()
        if candidate.is_file():
            return candidate
    return None


def executable_fingerprint(executable: Path) -> dict[str, str]:
    resolved = resolve_executable(executable)
    return {
        "path": str(resolved) if resolved is not None else str(executable),
        "sha256": sha256(resolved) if resolved is not None else "",
    }


def controlled_subprocess_env(omp_threads: int = 1) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["OMP_PROC_BIND"] = "false"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["VECLIB_MAXIMUM_THREADS"] = "1"
    return env


def job_stamp_path(result: Path) -> Path:
    return result.with_suffix(result.suffix + ".job.json")


def job_signature(
    executable: Path,
    input_path: Path,
    *,
    command_contract: dict[str, object] | None = None,
    executable_identity: dict[str, str] | None = None,
    campaign_fingerprint: dict[str, object] | None = None,
) -> dict[str, object]:
    if campaign_fingerprint is not None:
        validate_campaign_identity(campaign_fingerprint)
    fingerprint = executable_identity or executable_fingerprint(executable)
    signature: dict[str, object] = {
        "schema_version": 1,
        "executable": fingerprint["path"],
        "executable_sha256": fingerprint["sha256"],
        "input": str(input_path.resolve()),
        "input_sha256": sha256(input_path),
        "command_contract": command_contract or {},
    }
    if campaign_fingerprint is None:
        # Preserve the frozen GFN1/GFN2 stamp format.
        signature["campaign_fingerprint"] = {}
    else:
        signature["campaign_identity"] = campaign_fingerprint
    return signature


def write_job_stamp(
    result: Path,
    signature: dict[str, object],
    *,
    completed: bool,
    return_code: int,
) -> None:
    payload = dict(signature)
    payload.update({"completed": completed, "return_code": return_code})
    path = job_stamp_path(result)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def job_stamp_matches(result: Path, signature: dict[str, object]) -> bool:
    """A resume is valid only for the exact input and a hashable executable."""
    if not signature.get("executable_sha256"):
        return False
    path = job_stamp_path(result)
    if not path.exists():
        return False
    try:
        recorded = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    campaign_key = (
        "campaign_identity" if "campaign_identity" in signature else "campaign_fingerprint"
    )
    return bool(recorded.get("completed")) and all(
        recorded.get(key) == signature.get(key)
        for key in (
            "schema_version",
            "executable",
            "executable_sha256",
            "input",
            "input_sha256",
            "command_contract",
            campaign_key,
        )
    )


def campaign_fingerprint_sha256(fingerprint: dict[str, object]) -> str:
    canonical = dict(fingerprint)
    canonical.pop("fingerprint_sha256", None)
    payload = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def make_campaign_identity(
    *,
    campaign_id: str,
    cp2k_executable_sha256: str,
    cp2k_loaded_library_sha256: str,
    cp2k_cmake_cache_sha256: str,
    cp2k_embedded_source_revision: str,
    cp2k_source_revision: str,
    save_tblite_executable_sha256: str,
    save_tblite_source_revision: str,
    save_tblite_library_sha256: str,
    save_tblite_cmake_cache_sha256: str,
    dependency_lock_sha256: str,
) -> dict[str, object]:
    """Return the path-independent build identity shared by all GXTB campaigns."""
    identity: dict[str, object] = {
        "schema": CAMPAIGN_SCHEMA,
        "campaign_id": campaign_id,
        "cp2k_executable_sha256": cp2k_executable_sha256.lower(),
        "cp2k_loaded_library_sha256": cp2k_loaded_library_sha256.lower(),
        "cp2k_cmake_cache_sha256": cp2k_cmake_cache_sha256.lower(),
        "cp2k_embedded_source_revision": cp2k_embedded_source_revision.lower(),
        "cp2k_source_revision": cp2k_source_revision.lower(),
        "save_tblite_executable_sha256": save_tblite_executable_sha256.lower(),
        "save_tblite_source_revision": save_tblite_source_revision.lower(),
        "save_tblite_library_sha256": save_tblite_library_sha256.lower(),
        "save_tblite_cmake_cache_sha256": save_tblite_cmake_cache_sha256.lower(),
        "dependency_lock_sha256": dependency_lock_sha256.lower(),
    }
    identity["fingerprint_sha256"] = campaign_fingerprint_sha256(identity)
    validate_campaign_identity(identity)
    return identity


def validate_campaign_identity(identity: Mapping[str, object]) -> None:
    missing = [field for field in CAMPAIGN_IDENTITY_FIELDS if not identity.get(field)]
    if missing:
        raise ValueError("incomplete GXTB campaign identity: " + ", ".join(missing))
    if int(identity["schema"]) != CAMPAIGN_SCHEMA:
        raise ValueError(f"unsupported GXTB campaign schema: {identity['schema']}")
    core = {field: identity[field] for field in CAMPAIGN_IDENTITY_FIELDS}
    if identity.get("fingerprint_sha256") != campaign_fingerprint_sha256(core):
        raise ValueError("GXTB campaign fingerprint is internally inconsistent")


def campaign_identity_from_manifest(
    manifest: Mapping[str, object], manifest_path: Path
) -> dict[str, object]:
    """Read the frozen build declarations which are the campaign source of truth."""
    if manifest.get("campaign_state") != "production_ready":
        raise ValueError(
            f"GXTB campaign {manifest.get('campaign_id', '<unknown>')} is not production_ready"
        )
    cp2k = manifest.get("cp2k")
    save = manifest.get("save_tblite")
    dependencies = manifest.get("fetched_dependencies")
    if not isinstance(cp2k, Mapping) or not isinstance(save, Mapping):
        raise ValueError(f"invalid CP2K/save_tblite records in {manifest_path}")
    if not isinstance(dependencies, Mapping) or not dependencies:
        raise ValueError(f"missing fetched-dependency lock in {manifest_path}")
    return make_campaign_identity(
        campaign_id=str(manifest.get("campaign_id", "")),
        cp2k_executable_sha256=str(cp2k.get("binary_sha256", "")),
        cp2k_loaded_library_sha256=str(cp2k.get("loaded_library_sha256", "")),
        cp2k_cmake_cache_sha256=str(cp2k.get("cmake_cache_sha256", "")),
        cp2k_embedded_source_revision=str(cp2k.get("reported_revision", "")),
        cp2k_source_revision=str(cp2k.get("revision", "")),
        save_tblite_executable_sha256=str(save.get("cli_sha256", "")),
        save_tblite_source_revision=str(save.get("revision", "")),
        save_tblite_library_sha256=str(save.get("static_library_sha256", "")),
        save_tblite_cmake_cache_sha256=str(save.get("cmake_cache_sha256", "")),
        dependency_lock_sha256=campaign_fingerprint_sha256(dict(dependencies)),
    )


def completed_stamp_campaign_issue(
    result: Path,
    campaign_fingerprint: dict[str, object],
    *,
    executable_role: str,
    require_completed: bool = True,
) -> str | None:
    """Validate a GXTB result against the complete, immutable campaign identity."""
    try:
        validate_campaign_identity(campaign_fingerprint)
    except (TypeError, ValueError) as exc:
        return f"invalid campaign identity for {result}: {exc}"
    path = job_stamp_path(result)
    if not path.is_file():
        return f"missing campaign stamp for {result}"
    try:
        recorded = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid campaign stamp for {result}: {exc}"
    if require_completed and not recorded.get("completed"):
        return f"incomplete campaign stamp for {result}"
    if recorded.get("campaign_identity") != campaign_fingerprint:
        return f"campaign identity mismatch for {result}"
    input_path = Path(str(recorded.get("input", "")))
    if not input_path.is_file() or recorded.get("input_sha256") != sha256(input_path):
        return f"input fingerprint mismatch for {result}"
    expected_fields = {
        "cp2k": "cp2k_executable_sha256",
        "save_tblite": "save_tblite_executable_sha256",
    }
    field = expected_fields.get(executable_role)
    if field is None:
        return f"unknown executable role {executable_role} for {result}"
    expected_executable_sha = campaign_fingerprint.get(field)
    if not expected_executable_sha or recorded.get("executable_sha256") != expected_executable_sha:
        return f"{executable_role} executable fingerprint mismatch for {result}"
    return None


def command_output(command: list[str], allow_empty: bool = False) -> str:
    try:
        proc = subprocess.run(command, text=True, capture_output=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"unavailable: {exc}"
    text = (proc.stdout + proc.stderr).strip()
    if allow_empty and not text and proc.returncode == 0:
        return ""
    return text if text else f"exit status {proc.returncode}"


def shared_library_hashes(executable: Path) -> dict[str, str]:
    roots = (executable.resolve().parent, executable.resolve().parent.parent / "lib")
    patterns = (
        "libcp2k*.dylib",
        "libcp2k*.so*",
        "libtblite*.dylib",
        "libtblite*.so*",
        "libsave_tblite*.dylib",
        "libsave_tblite*.so*",
    )
    libraries: dict[str, str] = {}
    seen: set[Path] = set()
    for root in roots:
        if not root.is_dir():
            continue
        for pattern in patterns:
            for candidate in sorted(root.glob(pattern)):
                resolved = candidate.resolve()
                if not resolved.is_file() or resolved in seen:
                    continue
                seen.add(resolved)
                libraries[resolved.name] = sha256(resolved)
    return libraries


def version_summary(executable: Path) -> str:
    lines = command_output([str(executable), "--version"]).splitlines()
    summary = []
    for line in lines:
        if line.strip().lower().startswith("compiler options:"):
            break
        summary.append(line.rstrip())
    return "\n".join(summary)


def git_metadata(source: Path) -> dict[str, object]:
    if not source.exists():
        return {"available": False}
    revision = command_output(["git", "-C", str(source), "rev-parse", "HEAD"])
    branch = command_output(["git", "-C", str(source), "branch", "--show-current"])
    status = command_output(["git", "-C", str(source), "status", "--short"], allow_empty=True)
    diff = command_output(["git", "-C", str(source), "diff", "--binary"], allow_empty=True)
    return {
        "path": str(source.resolve()),
        "available": True,
        "revision": revision,
        "branch": branch,
        "dirty": bool(status),
        "working_tree_diff_sha256": hashlib.sha256(diff.encode()).hexdigest(),
    }


def _expand_macho_path(value: str, executable: Path) -> Path:
    loader = executable.resolve().parent
    expanded = value.replace("@loader_path", str(loader)).replace(
        "@executable_path", str(loader)
    )
    return Path(expanded)


def loaded_cp2k_library(executable: Path) -> Path:
    """Resolve the libcp2k image that the supplied launcher actually loads."""
    executable = executable.resolve(strict=True)
    system = platform.system()
    if system == "Darwin":
        dependencies = command_output(["otool", "-L", str(executable)])
        dependency = next(
            (
                line.strip().split(" (", 1)[0]
                for line in dependencies.splitlines()[1:]
                if Path(line.strip().split(" (", 1)[0]).name.startswith("libcp2k")
            ),
            "",
        )
        if not dependency:
            raise ValueError(f"No dynamically loaded libcp2k dependency found for {executable}")
        if dependency.startswith("@rpath/"):
            load_commands = command_output(["otool", "-l", str(executable)])
            rpaths: list[str] = []
            waiting_for_path = False
            for line in load_commands.splitlines():
                stripped = line.strip()
                if stripped == "cmd LC_RPATH":
                    waiting_for_path = True
                elif waiting_for_path and stripped.startswith("path "):
                    rpaths.append(stripped.split()[1])
                    waiting_for_path = False
            suffix = dependency.removeprefix("@rpath/")
            candidates = [
                _expand_macho_path(rpath, executable) / suffix for rpath in rpaths
            ]
        else:
            candidates = [_expand_macho_path(dependency, executable)]
    else:
        dependencies = command_output(["ldd", str(executable)])
        candidates = []
        for line in dependencies.splitlines():
            if "libcp2k" not in line:
                continue
            target = line.split("=>", 1)[1].strip().split()[0] if "=>" in line else line.split()[0]
            candidates.append(Path(target))
    resolved = [candidate.resolve() for candidate in candidates if candidate.is_file()]
    if len(set(resolved)) != 1:
        raise ValueError(
            f"Could not resolve one unambiguous loaded libcp2k for {executable}: "
            + ", ".join(str(path) for path in resolved)
        )
    return resolved[0]


def _validated_clean_source(source: Path, label: str) -> dict[str, object]:
    metadata = git_metadata(source.resolve(strict=True))
    revision = str(metadata.get("revision", ""))
    if not re.fullmatch(r"[0-9a-fA-F]{40}", revision):
        raise ValueError(f"Cannot determine exact {label} source revision from {source}")
    if metadata.get("dirty"):
        raise ValueError(f"{label} source checkout is dirty: {source}")
    return metadata


def validated_gxtb_campaign_fingerprint(
    cp2k: Path,
    cp2k_library: Path,
    cp2k_source: Path,
    save_tblite: Path,
    save_tblite_library: Path,
    save_tblite_source: Path,
) -> dict[str, object]:
    """Observe and validate the campaign artifacts available on this host."""
    cp2k = cp2k.resolve(strict=True)
    cp2k_library = cp2k_library.resolve(strict=True)
    save_tblite = save_tblite.resolve(strict=True)
    save_tblite_library = save_tblite_library.resolve(strict=True)
    if not os.access(cp2k, os.X_OK) or not os.access(save_tblite, os.X_OK):
        raise ValueError("CP2K and save_tblite campaign executables must be executable")
    actual_library = loaded_cp2k_library(cp2k).resolve(strict=True)
    if actual_library != cp2k_library:
        raise ValueError(
            f"CP2K loaded-library mismatch: launcher resolves {actual_library}, "
            f"but --cp2k-library is {cp2k_library}"
        )
    cp2k_source_record = _validated_clean_source(cp2k_source, "CP2K")
    save_source_record = _validated_clean_source(save_tblite_source, "save_tblite")
    cp2k_version = command_output([str(cp2k), "--version"])
    embedded = re.search(
        r"^\s*Source code revision\s+([0-9a-fA-F]+)\s*$", cp2k_version, flags=re.M
    )
    if embedded is None:
        raise ValueError("CP2K executable does not report an embedded source revision")
    embedded_revision = embedded.group(1).lower()
    source_revision = str(cp2k_source_record["revision"]).lower()
    if not source_revision.startswith(embedded_revision):
        raise ValueError(
            "CP2K executable/source revision mismatch: binary embeds "
            f"{embedded_revision}, source checkout is {source_revision}"
        )
    save_version = command_output([str(save_tblite), "--version"])
    observations: dict[str, object] = {
        "schema_version": 1,
        "cp2k": {
            "launcher": str(cp2k),
            "launcher_sha256": sha256(cp2k),
            "loaded_libcp2k": str(cp2k_library),
            "loaded_libcp2k_sha256": sha256(cp2k_library),
            "embedded_source_revision": embedded_revision,
            "source": cp2k_source_record,
            "version_output_sha256": hashlib.sha256(cp2k_version.encode()).hexdigest(),
        },
        "save_tblite": {
            "cli": str(save_tblite),
            "cli_sha256": sha256(save_tblite),
            "libtblite_a": str(save_tblite_library),
            "libtblite_a_sha256": sha256(save_tblite_library),
            "source": save_source_record,
            "version_output_sha256": hashlib.sha256(
                save_version.encode()
            ).hexdigest(),
            "version_output": save_version,
        },
    }
    return observations


def validated_gxtb_campaign_from_manifest(
    manifest_path: Path,
    cp2k_source: Path,
    save_tblite_source: Path,
    *,
    cp2k_override: Path | None = None,
    cp2k_library_override: Path | None = None,
    save_tblite_override: Path | None = None,
    save_tblite_library_override: Path | None = None,
) -> tuple[dict[str, object], dict[str, Path]]:
    """Resolve and verify all GXTB build artifacts from the central frozen manifest."""
    manifest_path = manifest_path.resolve(strict=True)
    try:
        manifest = json.loads(manifest_path.read_text())
        cp2k_record = manifest["cp2k"]
        save_record = manifest["save_tblite"]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise ValueError(f"Invalid GXTB campaign manifest {manifest_path}: {exc}") from exc
    declared_identity = campaign_identity_from_manifest(manifest, manifest_path)
    if cp2k_record.get("source_clean") is not True:
        raise ValueError("campaign manifest does not certify a clean CP2K source")
    if save_record.get("source_clean") is not True:
        raise ValueError("campaign manifest does not certify a clean save_tblite source")
    paths = {
        "cp2k": Path(str(cp2k_record["binary"])).resolve(strict=True),
        "cp2k_library": Path(str(cp2k_record["loaded_library"])).resolve(strict=True),
        "save_tblite": Path(str(save_record["cli"])).resolve(strict=True),
        "save_tblite_library": Path(str(save_record["static_library"])).resolve(strict=True),
    }
    overrides = {
        "cp2k": cp2k_override,
        "cp2k_library": cp2k_library_override,
        "save_tblite": save_tblite_override,
        "save_tblite_library": save_tblite_library_override,
    }
    for key, override in overrides.items():
        if override is not None and override.resolve(strict=True) != paths[key]:
            raise ValueError(
                f"{key} override {override.resolve()} differs from campaign manifest {paths[key]}"
            )
    expected_hashes = {
        "cp2k": str(cp2k_record["binary_sha256"]),
        "cp2k_library": str(cp2k_record["loaded_library_sha256"]),
        "save_tblite": str(save_record["cli_sha256"]),
        "save_tblite_library": str(save_record["static_library_sha256"]),
    }
    for key, expected in expected_hashes.items():
        observed = sha256(paths[key])
        if observed != expected:
            raise ValueError(
                f"{key} hash differs from campaign manifest: expected {expected}, observed {observed}"
            )
    if paths["save_tblite_library"].name != "libtblite.a":
        raise ValueError(
            "GXTB campaign must explicitly freeze save_tblite's libtblite.a archive; "
            f"got {paths['save_tblite_library']}"
        )
    observations = validated_gxtb_campaign_fingerprint(
        paths["cp2k"],
        paths["cp2k_library"],
        cp2k_source,
        paths["save_tblite"],
        paths["save_tblite_library"],
        save_tblite_source,
    )
    cp_observation = observations["cp2k"]  # type: ignore[assignment]
    save_observation = observations["save_tblite"]  # type: ignore[assignment]
    cp_source = cp_observation["source"]  # type: ignore[index]
    save_source = save_observation["source"]  # type: ignore[index]
    if cp_source["revision"] != cp2k_record["revision"]:  # type: ignore[index]
        raise ValueError("CP2K source HEAD differs from campaign manifest revision")
    if save_source["revision"] != save_record["revision"]:  # type: ignore[index]
        raise ValueError("save_tblite source HEAD differs from campaign manifest revision")
    if cp_observation["embedded_source_revision"] != cp2k_record["reported_revision"]:  # type: ignore[index]
        raise ValueError("CP2K embedded revision differs from campaign manifest")
    reported_save_version = str(save_record.get("reported_version", ""))
    if reported_save_version and reported_save_version not in str(
        save_observation.get("version_output", "")  # type: ignore[union-attr]
    ):
        raise ValueError("save_tblite CLI version differs from campaign manifest")
    observed_identity = make_campaign_identity(
        campaign_id=str(manifest.get("campaign_id", "")),
        cp2k_executable_sha256=expected_hashes["cp2k"],
        cp2k_loaded_library_sha256=expected_hashes["cp2k_library"],
        cp2k_cmake_cache_sha256=str(cp2k_record.get("cmake_cache_sha256", "")),
        cp2k_embedded_source_revision=str(cp_observation["embedded_source_revision"]),  # type: ignore[index]
        cp2k_source_revision=str(cp_source["revision"]),  # type: ignore[index]
        save_tblite_executable_sha256=expected_hashes["save_tblite"],
        save_tblite_source_revision=str(save_source["revision"]),  # type: ignore[index]
        save_tblite_library_sha256=expected_hashes["save_tblite_library"],
        save_tblite_cmake_cache_sha256=str(save_record.get("cmake_cache_sha256", "")),
        dependency_lock_sha256=campaign_fingerprint_sha256(
            dict(manifest["fetched_dependencies"])
        ),
    )
    if observed_identity != declared_identity:
        raise ValueError("observed GXTB build artifacts differ from the declared campaign identity")
    return observed_identity, paths


def repository_patch_metadata() -> dict[str, dict[str, str]]:
    patches = {
        "cp2k": ROOT.parent / "patches" / "cp2k_trunk_tblite_full_symmetry_scc.patch",
        "tblite": ROOT.parent / "patches" / "tblite_main_pr350_wsc_derivatives.patch",
    }
    return {
        name: {"path": f"../../patches/{path.name}", "sha256": sha256(path)}
        for name, path in patches.items()
        if path.is_file()
    }


def write_build_provenance(
    cp2k: Path,
    tblite: Path,
    cp2k_source: Path,
    tblite_source: Path,
    protocol: dict[str, object],
) -> None:
    payload = {
        "cp2k": {
            "executable": cp2k.name,
            "sha256": sha256(cp2k),
            "shared_library_sha256": shared_library_hashes(cp2k),
            "version": version_summary(cp2k),
            "source": git_metadata(cp2k_source),
        },
        "tblite": {
            "executable": tblite.name,
            "sha256": sha256(tblite),
            "shared_library_sha256": shared_library_hashes(tblite),
            "version": version_summary(tblite),
            "source": git_metadata(tblite_source),
        },
        "repository_patches": repository_patch_metadata(),
        "protocol": protocol,
    }
    path = ROOT / "data" / "build_provenance.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def executable_provenance(executable: Path, source: Path) -> dict[str, object]:
    full_version = command_output([str(executable), "--version"])
    return {
        "executable": str(executable.resolve()),
        "sha256": sha256(executable),
        "shared_library_sha256": shared_library_hashes(executable),
        "version": version_summary(executable),
        "version_and_build_flags": full_version,
        "source": git_metadata(source),
    }


def write_gxtb_build_provenance(
    cp2k: Path,
    save_tblite: Path,
    cp2k_source: Path,
    save_tblite_source: Path,
    protocol: dict[str, object],
    campaign_fingerprint: dict[str, object],
    campaign_manifest: Path,
) -> None:
    """Write GXTB provenance separately from the frozen GFN1/GFN2 record."""
    validate_campaign_identity(campaign_fingerprint)
    manifest_path = campaign_manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text())
    if campaign_identity_from_manifest(manifest, manifest_path) != campaign_fingerprint:
        raise ValueError("current campaign manifest build identity differs from LC12 provenance")
    payload = {
        "cp2k": executable_provenance(cp2k, cp2k_source),
        "save_tblite": executable_provenance(save_tblite, save_tblite_source),
        "campaign_identity": campaign_fingerprint,
        "campaign_manifest": {
            "path": str(manifest_path),
            "file_sha256": sha256(manifest_path),
            "campaign_id": str(manifest.get("campaign_id", "")),
            "campaign_state": str(manifest.get("campaign_state", "")),
            "authority": "single source of truth for the frozen build artifacts",
        },
        "protocol": protocol,
    }
    path = ROOT / "data" / "build_provenance_gxtb.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parse_cell_from_restart(path: Path) -> tuple[list[float], list[float], list[float]] | None:
    if not path.exists():
        return None
    lines = path.read_text(errors="ignore").splitlines()
    in_cell = False
    vectors: dict[str, list[float]] = {}
    abc: list[float] | None = None
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.split()[0] == "&CELL":
            in_cell = True
            continue
        if in_cell and upper in {"&END CELL", "&END"}:
            break
        if not in_cell:
            continue
        parts = stripped.split()
        if len(parts) >= 4 and parts[0].upper() in {"A", "B", "C"}:
            vectors[parts[0].upper()] = [float(parts[1]), float(parts[2]), float(parts[3])]
        elif len(parts) >= 4 and parts[0].upper() == "ABC":
            abc = [float(parts[1]), float(parts[2]), float(parts[3])]
    if {"A", "B", "C"} <= vectors.keys():
        return vectors["A"], vectors["B"], vectors["C"]
    if abc is not None:
        return [abc[0], 0.0, 0.0], [0.0, abc[1], 0.0], [0.0, 0.0, abc[2]]
    return None


def latest_restart(run_dir: Path, project: str) -> Path | None:
    matches = list(run_dir.glob(f"{project}*.restart"))
    if not matches:
        return None
    return max(matches, key=lambda p: p.stat().st_mtime)


def norm(vec: list[float]) -> float:
    return math.sqrt(sum(x * x for x in vec))


def volume(a: list[float], b: list[float], c: list[float]) -> float:
    return abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def optimized_lattice(run_dir: Path, project: str) -> float | None:
    restart = latest_restart(run_dir, project)
    if restart is None:
        return None
    cell = parse_cell_from_restart(restart)
    if cell is None:
        return None
    return sum(norm(v) for v in cell) / 3.0


def reference_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for ref in REFERENCES:
        rows.append(
            {
                "solid": ref.solid,
                "structure": ref.structure,
                "formula": "".join(f"{el}{n if n != 1 else ''}" for el, n in ref.formula),
                "a_exp_A": ref.a_exp,
                "a_HF_A": ref.a_hf,
                "a_MP2_A": ref.a_mp2,
                "a_SCS_MP2_A": ref.a_scs_mp2,
                "a_SOS_MP2_A": ref.a_sos_mp2,
                "ecoh_exp_eV_per_atom": ref.ecoh_exp,
                "ecoh_HF_eV_per_atom": ref.ecoh_hf,
                "ecoh_MP2_eV_per_atom": ref.ecoh_mp2,
                "ecoh_SCS_MP2_eV_per_atom": ref.ecoh_scs_mp2,
                "ecoh_SOS_MP2_eV_per_atom": ref.ecoh_sos_mp2,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def merge_method_rows(
    path: Path,
    rows: list[dict[str, object]],
    methods: tuple[str, ...],
    sort_key: Callable[[dict[str, object]], object] | None = None,
) -> list[dict[str, object]]:
    """Replace only selected-method records and preserve frozen legacy rows."""
    preserved: list[dict[str, object]] = [
        dict(row) for row in read_csv(path) if row.get("method", "") not in methods
    ]
    merged = preserved + rows
    if sort_key is not None:
        merged.sort(key=sort_key)
    write_csv(path, merged)
    return merged


def setup_inputs(cell_mesh: str, energy_meshes: list[str], methods: tuple[str, ...] = METHODS) -> None:
    reference_path = ROOT / "data" / "reference_goldzak2022.csv"
    if not reference_path.exists():
        write_csv(reference_path, reference_rows())
    for ref in REFERENCES:
        for method in methods:
            project = project_name(ref.solid, method, "cellopt", cell_mesh)
            text = solid_input(ref, method, "CELL_OPT", cell_mesh, ref.a_exp, project)
            write_file(ROOT / "inputs" / "cellopt" / method / ref.solid / f"{project}.inp", text)
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                text = solid_input(ref, method, "ENERGY", mesh, ref.a_exp, sp_project)
                write_file(ROOT / "inputs" / "single_point_initial" / method / ref.solid / f"{sp_project}.inp", text)
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    for method in methods:
        for element in elements:
            write_file(ROOT / "inputs" / "atoms" / method / f"atom_{element}_{method}.inp", atom_input(element, method))


def run_jobs(
    job_specs: list[tuple[str, Path, Path, bool]],
    cp2k: Path,
    jobs: int,
    threads: int,
    force: bool,
    campaign_fingerprint: dict[str, object] | None = None,
) -> None:
    cp2k_identity = executable_fingerprint(cp2k)

    def method_of(label: str) -> str:
        fields = label.split()
        return fields[1] if len(fields) > 1 and fields[1] in METHODS else ""

    def is_gxtb_no_smear_atom(label: str) -> bool:
        fields = label.split()
        return (
            len(fields) == 3
            and fields[0] == "atom"
            and fields[1] == "GXTB"
            and fields[2] in GXTB_NO_SMEAR_OT_ATOMS
        )

    def signature(inp: Path) -> dict[str, object]:
        return job_signature(
            cp2k,
            inp,
            command_contract={"driver": "cp2k", "omp_threads": threads},
            executable_identity=cp2k_identity,
            campaign_fingerprint=campaign_fingerprint,
        )

    pending: list[tuple[str, Path, Path, bool]] = []
    for label, inp, out, require_opt in job_specs:
        method = method_of(label)
        if method == "GXTB" and campaign_fingerprint is None:
            raise ValueError("GXTB CP2K jobs require a validated campaign fingerprint")
        if method:
            validate_method_input(
                inp.read_text(),
                method,
                gxtb_atom_reference=is_gxtb_no_smear_atom(label),
            )
        if not force and output_ok(out, require_opt=require_opt):
            if method != "GXTB" or job_stamp_matches(out, signature(inp)):
                continue
        pending.append((label, inp, out, require_opt))
    if not pending:
        print("No jobs pending.")
        return

    def worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool]:
        label, inp, out, require_opt = spec
        code = run_cp2k(cp2k, inp, out, threads)
        ok = output_ok(out, require_opt=require_opt)
        if method_of(label) == "GXTB":
            write_job_stamp(out, signature(inp), completed=ok, return_code=code)
        return label, code, ok

    print(f"Running {len(pending)} jobs with {jobs} worker(s), OMP_NUM_THREADS={threads}.")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in pending}
        done = 0
        failed: list[str] = []
        for future in concurrent.futures.as_completed(futures):
            label, code, ok = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(pending):3d}] {status:14s} {label}", flush=True)
            if not ok:
                failed.append(label)
    if failed:
        raise RuntimeError(
            f"{len(failed)} CP2K job(s) failed after preserving completed jobs: "
            + ", ".join(sorted(failed))
        )


def atom_job_specs(methods: tuple[str, ...] = METHODS) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    for method in methods:
        for element in elements:
            inp = ROOT / "inputs" / "atoms" / method / f"atom_{element}_{method}.inp"
            out = ROOT / "runs" / "atoms" / method / element / f"atom_{element}_{method}.out"
            run_inp = out.parent / inp.name
            if not run_inp.exists() or run_inp.read_text() != inp.read_text():
                write_file(run_inp, inp.read_text())
            specs.append((f"atom {method} {element}", run_inp, out, False))
    return specs


def run_tblite_atom_jobs(
    tblite: Path,
    jobs: int,
    force: bool,
    methods: tuple[str, ...] = METHODS,
    save_tblite: Path | None = None,
    campaign_fingerprint: dict[str, object] | None = None,
) -> None:
    if "GXTB" in methods and campaign_fingerprint is None:
        raise ValueError("GXTB save_tblite jobs require a validated campaign fingerprint")
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    identities = {
        "legacy": executable_fingerprint(tblite),
        "gxtb": executable_fingerprint(save_tblite or tblite),
    }
    specs: list[tuple[str, str, Path, Path, dict[str, object]]] = []
    for method in methods:
        for element in elements:
            run_dir = ROOT / "runs" / "atoms_cli" / method / element
            json_path = run_dir / f"atom_{element}_{method}.json"
            xyz_path = run_dir / f"atom_{element}.xyz"
            write_file(xyz_path, f"1\n{element} atom\n{element} 0.0 0.0 0.0\n")
            spin = ELEMENT_MULTIPLICITY[element] - 1
            executable = save_tblite if method == "GXTB" and save_tblite is not None else tblite
            signature = job_signature(
                executable,
                xyz_path,
                command_contract={
                    "driver": "tblite_run",
                    "method": method_cli_name(method),
                    "spin_2S": spin,
                    "accuracy": 0.05,
                    "restart": False,
                },
                executable_identity=identities["gxtb" if method == "GXTB" else "legacy"],
                campaign_fingerprint=campaign_fingerprint if method == "GXTB" else None,
            )
            if not force and json_path.exists() and parse_tblite_json_energy(json_path) is not None:
                if method != "GXTB" or job_stamp_matches(json_path, signature):
                    continue
            specs.append((method, element, run_dir, executable, signature))
    if not specs:
        print("No tblite atom jobs pending.")
        return

    def worker(spec: tuple[str, str, Path, Path, dict[str, object]]) -> tuple[str, int, bool]:
        method, element, run_dir, executable, signature = spec
        spin = ELEMENT_MULTIPLICITY[element] - 1
        xyz_path = run_dir / f"atom_{element}.xyz"
        json_path = run_dir / f"atom_{element}_{method}.json"
        out_path = run_dir / f"atom_{element}_{method}.out"
        cmd = [
            str(executable),
            "run",
            "--method",
            method_cli_name(method),
            "--spin",
            str(spin),
            "--acc",
            "0.05",
            "--json",
            json_path.name,
            "--no-restart",
            xyz_path.name,
        ]
        if json_path.exists():
            json_path.unlink()
        with out_path.open("w") as handle:
            proc = subprocess.run(
                cmd,
                cwd=run_dir,
                stdout=handle,
                stderr=subprocess.STDOUT,
                env=controlled_subprocess_env(1),
            )
        ok = proc.returncode == 0 and parse_tblite_json_energy(json_path) is not None
        if method == "GXTB":
            write_job_stamp(
                json_path,
                signature,
                completed=ok,
                return_code=proc.returncode,
            )
        return f"atom-cli {method} {element}", proc.returncode, ok

    print(f"Running {len(specs)} tblite atom jobs with {jobs} worker(s).")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in specs}
        done = 0
        failed: list[str] = []
        for future in concurrent.futures.as_completed(futures):
            label, code, ok = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(specs):3d}] {status:14s} {label}", flush=True)
            if not ok:
                failed.append(label)
    if failed:
        raise RuntimeError(
            f"{len(failed)} tblite/save_tblite atom job(s) failed after preserving completed jobs: "
            + ", ".join(sorted(failed))
        )


def parse_tblite_json_energy(path: Path) -> float | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    energy = data.get("energy")
    return float(energy) if energy is not None else None


def cellopt_job_specs(
    cell_mesh: str, methods: tuple[str, ...] = METHODS
) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in REFERENCES:
        for method in methods:
            project = project_name(ref.solid, method, "cellopt", cell_mesh)
            inp = ROOT / "inputs" / "cellopt" / method / ref.solid / f"{project}.inp"
            out = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh / f"{project}.out"
            # CP2K writes restart files into the input directory, so run from an isolated copy.
            run_inp = out.parent / inp.name
            if not run_inp.exists() or run_inp.read_text() != inp.read_text():
                write_file(run_inp, inp.read_text())
            specs.append((f"cellopt {method} {ref.solid} {cell_mesh}", run_inp, out, True))
    return specs


def generate_final_sp_inputs(
    cell_mesh: str, energy_meshes: list[str], methods: tuple[str, ...] = METHODS
) -> None:
    missing: list[str] = []
    for ref in REFERENCES:
        for method in methods:
            cell_project = project_name(ref.solid, method, "cellopt", cell_mesh)
            run_dir = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh
            opt_output = run_dir / f"{cell_project}.out"
            if not output_ok(opt_output, require_opt=True):
                missing.append(f"{method} {ref.solid}")
                continue
            a_opt = optimized_lattice(run_dir, cell_project)
            if a_opt is None:
                missing.append(f"{method} {ref.solid}")
                continue
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                text = solid_input(ref, method, "ENERGY", mesh, a_opt, sp_project)
                write_file(ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{sp_project}.inp", text)
    if missing:
        print("Missing optimized cells for:", ", ".join(missing), file=sys.stderr)


def sp_job_specs(
    energy_meshes: list[str], methods: tuple[str, ...] = METHODS
) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in REFERENCES:
        for method in methods:
            for mesh in energy_meshes:
                project = project_name(ref.solid, method, "sp", mesh)
                inp = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{project}.inp"
                if not inp.exists():
                    continue
                out = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{project}.out"
                specs.append((f"sp {method} {ref.solid} {mesh}", inp, out, False))
    return specs


def atom_energies(
    methods: tuple[str, ...] = METHODS,
    campaign_fingerprint: dict[str, object] | None = None,
) -> dict[tuple[str, str], float]:
    energies: dict[tuple[str, str], float] = {}
    for method in methods:
        atom_root = ROOT / "runs" / "atoms_cli" / method
        for element_dir in atom_root.glob("*"):
            if not element_dir.is_dir():
                continue
            out = element_dir / f"atom_{element_dir.name}_{method}.json"
            if method == "GXTB":
                if campaign_fingerprint is None:
                    raise ValueError("GXTB atom collection requires a campaign fingerprint")
                issue = completed_stamp_campaign_issue(
                    out, campaign_fingerprint, executable_role="save_tblite"
                )
                if issue:
                    raise RuntimeError(issue)
            energy = parse_tblite_json_energy(out)
            if energy is not None:
                energies[(method, element_dir.name)] = energy
    rows = [
            {
                "method": method,
                "element": element,
                "energy_hartree": f"{energy:.12f}",
                "source": "save_tblite_cli" if method == "GXTB" else "tblite_cli",
                "multiplicity": ELEMENT_MULTIPLICITY[element],
                "spin_2S": ELEMENT_MULTIPLICITY[element] - 1,
            }
            for (method, element), energy in sorted(energies.items())
        ]
    legacy = tuple(method for method in methods if method in LEGACY_METHODS)
    if legacy:
        merge_method_rows(
            ROOT / "data" / "atom_energies_tblite_cli.csv",
            [row for row in rows if row["method"] in legacy],
            legacy,
            sort_key=lambda row: (METHODS.index(str(row["method"])), str(row["element"])),
        )
    if "GXTB" in methods:
        write_csv(
            ROOT / "data" / "atom_energies_save_tblite_cli_gxtb.csv",
            [row for row in rows if row["method"] == "GXTB"],
        )
    return energies


def validate_atom_reference_agreement(
    methods: tuple[str, ...],
    tolerance_hartree: float = 1.0e-6,
    campaign_fingerprint: dict[str, object] | None = None,
) -> list[dict[str, object]]:
    """Compare CP2K and matching CLI isolated-atom energies method by method."""
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for method in methods:
        for element in sorted(ELEMENT_MULTIPLICITY):
            cp2k_path = ROOT / "runs" / "atoms" / method / element / f"atom_{element}_{method}.out"
            cli_path = ROOT / "runs" / "atoms_cli" / method / element / f"atom_{element}_{method}.json"
            cp2k_energy = parse_energy(cp2k_path)
            cli_energy = parse_tblite_json_energy(cli_path)
            stamp_issue = None
            if method == "GXTB":
                if campaign_fingerprint is None:
                    raise ValueError("GXTB atom validation requires a campaign fingerprint")
                stamp_issue = completed_stamp_campaign_issue(
                    cp2k_path, campaign_fingerprint, executable_role="cp2k"
                ) or completed_stamp_campaign_issue(
                    cli_path, campaign_fingerprint, executable_role="save_tblite"
                )
            delta = cp2k_energy - cli_energy if cp2k_energy is not None and cli_energy is not None else None
            passed = (
                output_ok(cp2k_path)
                and cli_energy is not None
                and delta is not None
                and abs(delta) <= tolerance_hartree
                and stamp_issue is None
            )
            rows.append(
                {
                    "method": method,
                    "element": element,
                    "multiplicity": ELEMENT_MULTIPLICITY[element],
                    "spin_2S": ELEMENT_MULTIPLICITY[element] - 1,
                    "cp2k_energy_hartree": f"{cp2k_energy:.12f}" if cp2k_energy is not None else "",
                    "cli_energy_hartree": f"{cli_energy:.12f}" if cli_energy is not None else "",
                    "delta_cp2k_minus_cli_hartree": f"{delta:.12e}" if delta is not None else "",
                    "tolerance_hartree": f"{tolerance_hartree:.3e}",
                    "passed": passed,
                    "cli_provider": "save_tblite" if method == "GXTB" else "tblite",
                    "cp2k_scf_contract": (
                        "gamma_no_smear_ot_interface_gate"
                        if method == "GXTB" and element in GXTB_NO_SMEAR_OT_ATOMS
                        else (
                            "li_native_gxtb_fdiis_diagonalization_300K_interface_gate"
                            if method == "GXTB"
                            else "legacy_atom_protocol"
                        )
                    ),
                    "cohesive_energy_atom_reference": (
                        "save_tblite_cli_only" if method == "GXTB" else "tblite_cli"
                    ),
                    "campaign_stamp_issue": stamp_issue or "",
                }
            )
            if not passed:
                problems.append(f"{method}/{element}")
    legacy = tuple(method for method in methods if method in LEGACY_METHODS)
    if legacy:
        merge_method_rows(
            ROOT / "data" / "atom_reference_cp2k_vs_tblite_cli.csv",
            [row for row in rows if row["method"] in legacy],
            legacy,
            sort_key=lambda row: (METHODS.index(str(row["method"])), str(row["element"])),
        )
    if "GXTB" in methods:
        write_csv(
            ROOT / "data" / "atom_reference_cp2k_vs_save_tblite_gxtb.csv",
            [row for row in rows if row["method"] == "GXTB"],
        )
    if problems:
        raise RuntimeError(
            "CP2K/CLI isolated-atom validation failed for "
            + ", ".join(problems)
            + f" (tolerance {tolerance_hartree:.3e} Eh)"
        )
    return rows


def analyse(
    cell_mesh: str,
    energy_meshes: list[str],
    result_mesh: str,
    methods: tuple[str, ...] = METHODS,
) -> None:
    atom_e = atom_energies(methods)
    rows: list[dict[str, object]] = []
    for ref in REFERENCES:
        n_atoms = len(conventional_cell_atoms(ref))
        counts = atom_counts(ref)
        for method in methods:
            cell_project = project_name(ref.solid, method, "cellopt", cell_mesh)
            cell_run = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh
            opt_out = cell_run / f"{cell_project}.out"
            cellopt_completed = output_ok(opt_out, require_opt=True)
            a_opt = optimized_lattice(cell_run, cell_project) if cellopt_completed else None
            opt_energy = parse_energy(opt_out)
            atom_sum = None
            if all((method, el) in atom_e for el in counts):
                atom_sum = sum(atom_e[(method, el)] * count for el, count in counts.items())
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                sp_out = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{sp_project}.out"
                sp_energy = parse_energy(sp_out)
                ecoh = None
                if atom_sum is not None and sp_energy is not None:
                    ecoh = (atom_sum - sp_energy) * HARTREE_TO_EV / n_atoms
                rows.append(
                    {
                        "solid": ref.solid,
                        "structure": ref.structure,
                        "method": method,
                        "cell_mesh": cell_mesh,
                        "energy_mesh": mesh,
                        "cellopt_completed": cellopt_completed,
                        "sp_completed": output_ok(sp_out, require_opt=False),
                        "a_calc_A": f"{a_opt:.8f}" if a_opt is not None else "",
                        "a_ref_exp_A": ref.a_exp,
                        "a_error_A": f"{(a_opt - ref.a_exp):.8f}" if a_opt is not None else "",
                        "a_abs_error_A": f"{abs(a_opt - ref.a_exp):.8f}" if a_opt is not None else "",
                        "ecoh_calc_eV_per_atom": f"{ecoh:.8f}" if ecoh is not None else "",
                        "ecoh_ref_exp_eV_per_atom": ref.ecoh_exp,
                        "ecoh_error_eV_per_atom": f"{(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                        "ecoh_abs_error_eV_per_atom": f"{abs(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                        "solid_energy_hartree": f"{sp_energy:.12f}" if sp_energy is not None else "",
                        "cellopt_last_energy_hartree": f"{opt_energy:.12f}" if opt_energy is not None else "",
                        "n_atoms_conventional_cell": n_atoms,
                        "atom_reference_source": "save_tblite_cli" if method == "GXTB" else "tblite_cli",
                    }
                )
    rows = merge_method_rows(
        ROOT / "data" / "results.csv",
        rows,
        methods,
        sort_key=lambda row: (
            [ref.solid for ref in REFERENCES].index(str(row["solid"])),
            METHODS.index(str(row["method"])),
            str(row["energy_mesh"]),
        ),
    )

    final_rows = [r for r in rows if r["energy_mesh"] == result_mesh]
    summary: list[dict[str, object]] = []
    available_methods = tuple(method for method in METHODS if any(r["method"] == method for r in final_rows))
    for method in available_methods:
        method_rows = [
            r
            for r in final_rows
            if r["method"] == method and truth(r["cellopt_completed"]) and truth(r["sp_completed"])
        ]
        a_err = [float(r["a_error_A"]) for r in method_rows if r["a_error_A"] != ""]
        e_err = [float(r["ecoh_error_eV_per_atom"]) for r in method_rows if r["ecoh_error_eV_per_atom"] != ""]
        summary.append(
            {
                "method": method,
                "n_complete": len(method_rows),
                "result_mesh": result_mesh,
                "cell_mesh": cell_mesh,
                "a_ME_A": f"{sum(a_err) / len(a_err):.8f}" if a_err else "",
                "a_MAE_A": f"{sum(abs(x) for x in a_err) / len(a_err):.8f}" if a_err else "",
                "a_RMSE_A": f"{math.sqrt(sum(x * x for x in a_err) / len(a_err)):.8f}" if a_err else "",
                "ecoh_ME_eV_per_atom": f"{sum(e_err) / len(e_err):.8f}" if e_err else "",
                "ecoh_MAE_eV_per_atom": f"{sum(abs(x) for x in e_err) / len(e_err):.8f}" if e_err else "",
                "ecoh_RMSE_eV_per_atom": f"{math.sqrt(sum(x * x for x in e_err) / len(e_err)):.8f}" if e_err else "",
            }
        )
    write_csv(ROOT / "data" / "summary.csv", summary)
    write_markdown(final_rows, summary, result_mesh)
    plot(final_rows, result_mesh)


def write_markdown(rows: list[dict[str, object]], summary: list[dict[str, object]], result_mesh: str) -> None:
    order = [ref.solid for ref in REFERENCES]
    by_key = {(r["solid"], r["method"]): r for r in rows}
    lines = [
        f"# Goldzak12 CP2K/tblite results ({result_mesh} final energies)",
        "",
        "All GFN values use native Bloch k-points in CP2K. Cohesive energies are in eV per atom.",
        "",
        "## Summary",
        "",
        "| method | n | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['method']} | {row['n_complete']} | {row['a_ME_A']} | {row['a_MAE_A']} | "
            f"{row['ecoh_ME_eV_per_atom']} | {row['ecoh_MAE_eV_per_atom']} |"
        )
    methods = tuple(method for method in METHODS if any(r["method"] == method for r in rows))
    header = ["solid", "a exp"]
    for method in methods:
        header += [f"a {method}", f"da {method}"]
    header.append("Ecoh exp")
    for method in methods:
        header += [f"Ecoh {method}", f"dE {method}"]
    lines += ["", "## Per-system comparison to experiment", "", "| " + " | ".join(header) + " |"]
    lines.append("|---|" + "---:|" * (len(header) - 1))
    refs = {ref.solid: ref for ref in REFERENCES}
    for solid in order:
        ref = refs[solid]
        cells = [solid, f"{ref.a_exp:.3f}"]
        for method in methods:
            row = by_key.get((solid, method), {})
            cells += [fmt(row.get("a_calc_A"), 4), fmt(row.get("a_error_A"), 4)]
        cells.append(f"{ref.ecoh_exp:.2f}")
        for method in methods:
            row = by_key.get((solid, method), {})
            cells += [fmt(row.get("ecoh_calc_eV_per_atom"), 3), fmt(row.get("ecoh_error_eV_per_atom"), 3)]
        lines.append("| " + " | ".join(cells) + " |")
    (ROOT / "data").mkdir(parents=True, exist_ok=True)
    (ROOT / "data" / "results.md").write_text("\n".join(lines) + "\n")


def fmt(value: object, digits: int) -> str:
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def plot(rows: list[dict[str, object]], result_mesh: str) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    complete = [r for r in rows if truth(r["cellopt_completed"]) and truth(r["sp_completed"])]
    if not complete:
        return
    solids = [ref.solid for ref in REFERENCES]
    x = np.arange(len(solids))
    methods = tuple(method for method in METHODS if any(r["method"] == method for r in complete))
    width = min(0.8 / max(len(methods), 1), 0.36)
    colors = METHOD_COLORS
    for prop, ylabel, filename, key in [
        ("lattice", "lattice-constant error (A)", "goldzak12_lattice_errors", "a_error_A"),
        ("cohesive", "cohesive-energy error (eV/atom)", "goldzak12_cohesive_errors", "ecoh_error_eV_per_atom"),
    ]:
        fig, ax = plt.subplots(figsize=(10.5, 4.6))
        for idx, method in enumerate(methods):
            values = []
            for solid in solids:
                row = next((r for r in complete if r["solid"] == solid and r["method"] == method), None)
                values.append(float(row[key]) if row and row[key] != "" else np.nan)
            offset = (idx - (len(methods) - 1) / 2.0) * width
            ax.bar(x + offset, values, width, label=method, color=colors[method])
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(solids, rotation=45, ha="right")
        ax.set_title(f"Goldzak12 CP2K/tblite native-Bloch {result_mesh}: {prop} errors")
        ax.legend(frameon=False)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
        fig.tight_layout()
        out_base = ROOT / "figures" / filename
        out_base.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_base.with_suffix(".png"), dpi=220)
        fig.savefig(out_base.with_suffix(".pdf"))
        plt.close(fig)

    ref_by_solid = {ref.solid: ref for ref in REFERENCES}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    for method in methods:
        a_ref = []
        a_calc = []
        e_ref = []
        e_calc = []
        for solid in solids:
            row = next((r for r in complete if r["solid"] == solid and r["method"] == method), None)
            if not row:
                continue
            a_ref.append(ref_by_solid[solid].a_exp)
            a_calc.append(float(row["a_calc_A"]))
            e_ref.append(ref_by_solid[solid].ecoh_exp)
            e_calc.append(float(row["ecoh_calc_eV_per_atom"]))
        axes[0].scatter(a_ref, a_calc, label=method, color=colors[method], s=44)
        axes[1].scatter(e_ref, e_calc, label=method, color=colors[method], s=44)
    for ax, label in zip(axes, ["lattice constant (A)", "cohesive energy (eV/atom)"]):
        lo, hi = ax.get_xlim()
        ylo, yhi = ax.get_ylim()
        mn, mx = min(lo, ylo), max(hi, yhi)
        ax.plot([mn, mx], [mn, mx], color="black", linewidth=0.8)
        ax.set_xlim(mn, mx)
        ax.set_ylim(mn, mx)
        ax.set_xlabel(f"experiment {label}")
        ax.set_ylabel(f"CP2K/tblite {label}")
        ax.grid(color="#d0d0d0", linewidth=0.6, alpha=0.7)
    axes[0].legend(frameon=False)
    fig.suptitle(f"Goldzak12 CP2K/tblite native-Bloch {result_mesh} vs experiment")
    fig.tight_layout()
    out_base = ROOT / "figures" / "goldzak12_scatter"
    fig.savefig(out_base.with_suffix(".png"), dpi=220)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--cell-mesh", default="k444")
        p.add_argument("--energy-mesh", action="append", default=[])
        p.add_argument(
            "--method",
            action="append",
            choices=METHODS,
            help="method to prepare/run/analyse; repeat as needed (default: GFN1 and GFN2)",
        )

    p_setup = sub.add_parser("setup")
    common(p_setup)

    p_run = sub.add_parser("run")
    common(p_run)
    p_run.add_argument("--cp2k", type=Path)
    p_run.add_argument("--tblite", type=Path, default=DEFAULT_TBLITE)
    p_run.add_argument("--save-tblite", type=Path)
    p_run.add_argument("--campaign-manifest", type=Path, default=DEFAULT_GXTB_CAMPAIGN_MANIFEST)
    p_run.add_argument("--cp2k-library", type=Path)
    p_run.add_argument("--save-tblite-library", type=Path)
    p_run.add_argument("--cp2k-source", type=Path, default=DEFAULT_CP2K_SOURCE)
    p_run.add_argument("--save-tblite-source", type=Path, default=DEFAULT_SAVE_TBLITE_SOURCE)
    p_run.add_argument("--jobs", type=int, default=4)
    p_run.add_argument("--threads", type=int, default=1)
    p_run.add_argument("--force", action="store_true")
    p_run.add_argument("--prune-transients", action="store_true")

    p_sp = sub.add_parser("single-points")
    common(p_sp)
    p_sp.add_argument("--cp2k", type=Path)
    p_sp.add_argument("--save-tblite", type=Path)
    p_sp.add_argument("--campaign-manifest", type=Path, default=DEFAULT_GXTB_CAMPAIGN_MANIFEST)
    p_sp.add_argument("--cp2k-library", type=Path)
    p_sp.add_argument("--save-tblite-library", type=Path)
    p_sp.add_argument("--cp2k-source", type=Path, default=DEFAULT_CP2K_SOURCE)
    p_sp.add_argument("--save-tblite-source", type=Path, default=DEFAULT_SAVE_TBLITE_SOURCE)
    p_sp.add_argument("--jobs", type=int, default=4)
    p_sp.add_argument("--threads", type=int, default=1)
    p_sp.add_argument("--force", action="store_true")

    p_atom_check = sub.add_parser(
        "atom-check",
        help="run and compare CP2K and method-matched CLI isolated atoms",
    )
    common(p_atom_check)
    p_atom_check.add_argument("--cp2k", type=Path)
    p_atom_check.add_argument("--tblite", type=Path, default=DEFAULT_TBLITE)
    p_atom_check.add_argument("--save-tblite", type=Path)
    p_atom_check.add_argument("--campaign-manifest", type=Path, default=DEFAULT_GXTB_CAMPAIGN_MANIFEST)
    p_atom_check.add_argument("--cp2k-library", type=Path)
    p_atom_check.add_argument("--save-tblite-library", type=Path)
    p_atom_check.add_argument("--cp2k-source", type=Path, default=DEFAULT_CP2K_SOURCE)
    p_atom_check.add_argument("--save-tblite-source", type=Path, default=DEFAULT_SAVE_TBLITE_SOURCE)
    p_atom_check.add_argument("--jobs", type=int, default=4)
    p_atom_check.add_argument("--threads", type=int, default=1)
    p_atom_check.add_argument("--force", action="store_true")
    p_atom_check.add_argument("--tolerance-hartree", type=float, default=1.0e-6)

    p_analyse = sub.add_parser("analyse")
    common(p_analyse)
    p_analyse.add_argument("--result-mesh", default="")

    args = parser.parse_args()
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]
    methods = selected_methods(args.method)
    campaign_fingerprint: dict[str, object] | None = None
    if "GXTB" in methods and args.command in {"run", "single-points", "atom-check"}:
        try:
            campaign_fingerprint, campaign_paths = validated_gxtb_campaign_from_manifest(
                args.campaign_manifest,
                args.cp2k_source,
                args.save_tblite_source,
                cp2k_override=args.cp2k,
                cp2k_library_override=args.cp2k_library,
                save_tblite_override=args.save_tblite,
                save_tblite_library_override=args.save_tblite_library,
            )
            args.cp2k = campaign_paths["cp2k"]
            args.cp2k_library = campaign_paths["cp2k_library"]
            args.save_tblite = campaign_paths["save_tblite"]
            args.save_tblite_library = campaign_paths["save_tblite_library"]
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
    elif args.command in {"run", "single-points", "atom-check"}:
        args.cp2k = args.cp2k or DEFAULT_CP2K
        args.save_tblite = args.save_tblite or DEFAULT_SAVE_TBLITE

    if args.command == "setup":
        setup_inputs(args.cell_mesh, energy_meshes, methods)
        return 0

    if args.command == "run":
        setup_inputs(args.cell_mesh, energy_meshes, methods)
        run_tblite_atom_jobs(
            args.tblite,
            args.jobs,
            args.force,
            methods,
            args.save_tblite,
            campaign_fingerprint,
        )
        run_jobs(
            cellopt_job_specs(args.cell_mesh, methods),
            args.cp2k,
            args.jobs,
            args.threads,
            args.force,
            campaign_fingerprint,
        )
        generate_final_sp_inputs(args.cell_mesh, energy_meshes, methods)
        run_jobs(
            sp_job_specs(energy_meshes, methods),
            args.cp2k,
            args.jobs,
            args.threads,
            args.force,
            campaign_fingerprint,
        )
        analyse(args.cell_mesh, energy_meshes, energy_meshes[-1], methods)
        if args.prune_transients and "GXTB" in methods:
            count, size = prune_gxtb_transients()
            print(f"Pruned {count} validated GXTB transient file(s), {size} byte(s).")
        return 0

    if args.command == "single-points":
        generate_final_sp_inputs(args.cell_mesh, energy_meshes, methods)
        run_jobs(
            sp_job_specs(energy_meshes, methods),
            args.cp2k,
            args.jobs,
            args.threads,
            args.force,
            campaign_fingerprint,
        )
        analyse(args.cell_mesh, energy_meshes, energy_meshes[-1], methods)
        return 0

    if args.command == "atom-check":
        if args.tolerance_hartree <= 0.0:
            parser.error("--tolerance-hartree must be positive")
        setup_inputs(args.cell_mesh, energy_meshes, methods)
        run_tblite_atom_jobs(
            args.tblite,
            args.jobs,
            args.force,
            methods,
            args.save_tblite,
            campaign_fingerprint,
        )
        run_jobs(
            atom_job_specs(methods),
            args.cp2k,
            args.jobs,
            args.threads,
            args.force,
            campaign_fingerprint,
        )
        validate_atom_reference_agreement(
            methods, args.tolerance_hartree, campaign_fingerprint
        )
        return 0

    if args.command == "analyse":
        result_mesh = args.result_mesh or energy_meshes[-1]
        analyse(args.cell_mesh, energy_meshes, result_mesh, methods)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
