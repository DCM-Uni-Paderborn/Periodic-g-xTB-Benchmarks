#!/usr/bin/env python3
"""Prepare and analyse k-point dependent DMC-ICE13 GFN benchmarks."""

from __future__ import annotations

import argparse
import csv
import hashlib
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

PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]
BASELINE_METHODS = ["GFN1", "GFN2"]
METHODS = [*BASELINE_METHODS, "GXTB"]
METHOD_LABELS = {"GFN1": "GFN1-xTB", "GFN2": "GFN2-xTB", "GXTB": "g-xTB"}
METHOD_COLORS = {"GFN1": "#c44e52", "GFN2": "#4c72b0", "GXTB": "#55a868"}
GXTB_PRODUCTION_PREFIX = "gxtb_spglib"
GXTB_PROTOCOL_ID = "dmc13-gxtb-spglib-reduced-v1"

DMC_ABS_KJMOL = {
    "Ih": -59.45,
    "II": -59.14,
    "III": -58.20,
    "IV": -55.62,
    "VI": -57.67,
    "VII": -54.46,
    "VIII": -55.22,
    "IX": -58.85,
    "XI": -59.29,
    "XIII": -57.33,
    "XIV": -57.75,
    "XV": -57.71,
    "XVII": -57.70,
}

# Table I of Della Pia et al. reports the relative XI value explicitly as
# 0.15 kJ/mol/H2O.  The historic benchmark derives 0.16 from the rounded
# absolute values above.  Keep the latter for byte-for-byte comparability with
# the published GFN1/GFN2 numbers and write the primary-table sensitivity as a
# separate result.
DMC_REL_PRIMARY_KJMOL = {
    phase: DMC_ABS_KJMOL[phase] - DMC_ABS_KJMOL["Ih"] for phase in PHASES
}
DMC_REL_PRIMARY_KJMOL["XI"] = 0.15

# The six-member production core is frozen for direct comparability with the
# GFN1/GFN2 campaign.  Denser meshes are opt-in convergence extensions: adding
# them must not silently change the default 78-job g-xTB production matrix.
MESHES = [
    {
        "id": "gamma",
        "label": "Gamma",
        "nk_total": 1,
        "scheme": "GAMMA",
        "shift": "",
        "della_pia_role": "Gamma-only reference",
    },
    {
        "id": "k111",
        "label": "1x1x1",
        "nk_total": 1,
        "scheme": "MACDONALD 1 1 1 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "explicit one-point MacDonald mesh",
    },
    {
        "id": "k222",
        "label": "2x2x2",
        "nk_total": 8,
        "scheme": "MACDONALD 2 2 2 0.25 0.25 0.25",
        "shift": "0.25 0.25 0.25",
        "della_pia_role": "coarse Gamma-centered convergence check",
    },
    {
        "id": "k333",
        "label": "3x3x3",
        "nk_total": 27,
        "scheme": "MACDONALD 3 3 3 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "GGA, meta-GGA, and vdW single points",
    },
    {
        "id": "k444",
        "label": "4x4x4",
        "nk_total": 64,
        "scheme": "MACDONALD 4 4 4 0.375 0.375 0.375",
        "shift": "0.375 0.375 0.375",
        "della_pia_role": "hybrid-XC single points",
    },
    {
        "id": "k555",
        "label": "5x5x5",
        "nk_total": 125,
        "scheme": "MACDONALD 5 5 5 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "dense convergence check",
    },
]

DENSE_EXTENSION_MESHES = [
    {
        "id": "k666",
        "label": "6x6x6",
        "nk_total": 216,
        "scheme": (
            "MACDONALD 6 6 6 0.4166666666666667 "
            "0.4166666666666667 0.4166666666666667"
        ),
        "shift": "0.4166666666666667 0.4166666666666667 0.4166666666666667",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k777",
        "label": "7x7x7",
        "nk_total": 343,
        "scheme": "MACDONALD 7 7 7 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k888",
        "label": "8x8x8",
        "nk_total": 512,
        "scheme": "MACDONALD 8 8 8 0.4375 0.4375 0.4375",
        "shift": "0.4375 0.4375 0.4375",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k999",
        "label": "9x9x9",
        "nk_total": 729,
        "scheme": "MACDONALD 9 9 9 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k101010",
        "label": "10x10x10",
        "nk_total": 1000,
        "scheme": "MACDONALD 10 10 10 0.45 0.45 0.45",
        "shift": "0.45 0.45 0.45",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k111111",
        "label": "11x11x11",
        "nk_total": 1331,
        "scheme": "MACDONALD 11 11 11 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k121212",
        "label": "12x12x12",
        "nk_total": 1728,
        "scheme": (
            "MACDONALD 12 12 12 0.4583333333333333 "
            "0.4583333333333333 0.4583333333333333"
        ),
        "shift": (
            "0.4583333333333333 0.4583333333333333 0.4583333333333333"
        ),
        "della_pia_role": "optional dense convergence extension",
    },
    {
        "id": "k131313",
        "label": "13x13x13",
        "nk_total": 2197,
        "scheme": "MACDONALD 13 13 13 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "optional dense convergence extension",
    },
]
SUPPORTED_MESHES = [*MESHES, *DENSE_EXTENSION_MESHES]
MESH_BY_ID = {str(mesh["id"]): mesh for mesh in SUPPORTED_MESHES}
CONVERGENCE_MESH_IDS = [
    "k111",
    "k222",
    "k333",
    "k444",
    "k555",
    "k666",
    "k777",
    "k888",
    "k999",
    "k101010",
    "k111111",
    "k121212",
    "k131313",
]
CONVERGENCE_MAX_ABS_KJMOL = 0.10
CONVERGENCE_RMS_KJMOL = 0.05
CONVERGENCE_REQUIRED_CONSECUTIVE_PAIRS = 2
PHASEWISE_KPOINT_MAX_ABS_KJMOL = 0.05
BUILD_IDENTITY_FIELDS = (
    "cp2k_sha256",
    "cp2k_library_sha256",
    "tblite_static_library_sha256",
    "cp2k_source_revision",
    "tblite_source_revision",
)
QUALIFICATION_EVIDENCE_SCHEMA_VERSION = 3
MAX_TOTAL_ENERGY_TOLERANCE_HARTREE = 1.0e-10
MAX_RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O = 1.0e-3


def insert_kpoints(
    input_text: str,
    mesh: dict[str, object],
    *,
    implicit_gamma: bool = False,
) -> str:
    # Keep the implicit Gamma input free of a KPOINTS section.  The explicit
    # k111 member intentionally takes the independent MacDonald 1x1x1 path.
    if mesh["id"] == "gamma" and implicit_gamma:
        return input_text

    block = [
        "    &KPOINTS",
        f"      SCHEME {mesh['scheme']}",
    ]
    if mesh["id"] != "gamma":
        block += [
            "      SYMMETRY T",
            "      FULL_GRID F",
            "      SYMMETRY_BACKEND SPGLIB",
            "      SYMMETRY_REDUCTION_METHOD SPGLIB",
        ]
    block.append("    &END KPOINTS")
    return input_text.replace("    &END QS\n", "    &END QS\n" + "\n".join(block) + "\n", 1)


def remove_kpoints(input_text: str) -> str:
    """Remove one CP2K KPOINTS section from a reusable g-xTB template."""
    return re.sub(
        r"(?ms)^    &KPOINTS\n.*?^    &END KPOINTS\n",
        "",
        input_text,
        count=1,
    )


def enable_cell_canonicalization(input_text: str) -> str:
    if "      CANONICALIZE TRUE\n" in input_text:
        return input_text
    return input_text.replace("    &CELL\n", "    &CELL\n      CANONICALIZE TRUE\n", 1)


def method_label(method: str) -> str:
    return METHOD_LABELS[method]


def gxtb_from_gfn2_template(text: str, phase: str) -> str:
    text = text.replace(f"PROJECT ice_{phase}_GFN2", f"PROJECT ice_{phase}_GXTB", 1)
    text = text.replace("          METHOD GFN2\n", "          METHOD GXTB\n", 1)
    marker = "        &END TBLITE\n      &END XTB"
    replacement = """        &END TBLITE
        SCC_MIXER TBLITE
        &TBLITE_MIXER
          ITERATIONS 300
        &END TBLITE_MIXER
      &END XTB"""
    if marker not in text:
        raise ValueError(f"Cannot select the native g-xTB mixer for ice {phase}")
    text = text.replace(marker, replacement, 1)
    return f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}\n" + text


def prepare_inputs(
    methods: list[str] | None = None,
    gxtb_input_root: Path | None = None,
    mesh_ids: list[str] | None = None,
    phase_ids: list[str] | None = None,
) -> None:
    # This repository owns g-xTB inputs.  Canonical GFN1/GFN2 inputs live in
    # DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks and are never regenerated
    # implicitly here.
    selected = methods or ["GXTB"]
    selected_meshes = (
        [MESH_BY_ID[mesh_id] for mesh_id in mesh_ids]
        if mesh_ids
        else MESHES
    )
    selected_phases = phase_ids or PHASES
    gxtb_input_root = gxtb_input_root or ROOT / f"{GXTB_PRODUCTION_PREFIX}_inputs"
    for mesh in selected_meshes:
        mesh_id = str(mesh["id"])
        for method in selected:
            out_root = gxtb_input_root if method == "GXTB" else ROOT / "kpoint_inputs"
            out_dir = out_root / mesh_id
            out_dir.mkdir(parents=True, exist_ok=True)
            for phase in selected_phases:
                if method == "GXTB":
                    base_path = (
                        ROOT
                        / "kpoint_inputs"
                        / "gamma"
                        / f"ice_{phase}_GXTB_gamma.inp"
                    )
                    text = remove_kpoints(base_path.read_text())
                    text = re.sub(
                        rf"(?m)^  PROJECT ice_{re.escape(phase)}_GXTB_gamma$",
                        f"  PROJECT ice_{phase}_GXTB",
                        text,
                        count=1,
                    )
                    if not text.startswith(f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}\n"):
                        text = f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}\n" + text
                else:
                    base_path = ROOT / "inputs" / f"ice_{phase}_{method}.inp"
                    if not base_path.exists():
                        raise FileNotFoundError(
                            f"{base_path} is not stored in the g-xTB repository; "
                            "use DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks "
                            "for GFN1/GFN2 input generation"
                        )
                    text = base_path.read_text()
                text = enable_cell_canonicalization(text)
                project = f"ice_{phase}_{method}_{mesh_id}"
                text = text.replace(f"PROJECT ice_{phase}_{method}", f"PROJECT {project}")
                text = insert_kpoints(
                    text,
                    mesh,
                    implicit_gamma=(method == "GXTB"),
                )
                (out_dir / f"{project}.inp").write_text(text)


def parse_energy(output: Path) -> float | None:
    if not output.exists():
        return None
    text = output.read_text(errors="ignore")
    if (
        "PROGRAM ENDED" not in text
        or "SCF run converged" not in text
        or "SCF run NOT converged" in text
        or "ABORT" in text
    ):
        return None
    energy = None
    for line in text.splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def output_path(
    mesh_id: str,
    method: str,
    phase: str,
    gxtb_run_root: Path | None = None,
) -> Path:
    if method == "GXTB" and gxtb_run_root is not None:
        return gxtb_run_root / mesh_id / phase / f"ice_{phase}_{method}_{mesh_id}.out"
    if mesh_id == "gamma":
        return ROOT / "runs" / method / phase / f"ice_{phase}_{method}.out"
    return ROOT / "runs_kpoints" / mesh_id / method / phase / f"ice_{phase}_{method}_{mesh_id}.out"


def data_output_path(stem: str, output_prefix: str | None = None) -> Path:
    if output_prefix:
        return DATA / f"dmc_ice13_{output_prefix}_{stem}"
    if stem == "kpoint_results.json":
        return DATA / stem
    return DATA / f"dmc_ice13_{stem}"


def stats(errors: list[float]) -> dict[str, float]:
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(e) for e in errors) / len(errors),
        "RMSE": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "MaxAE": max(abs(e) for e in errors),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def frozen_build_identity_from_manifest(
    manifest: dict[str, object],
) -> dict[str, str]:
    cp2k = manifest.get("cp2k")
    save_tblite = manifest.get("save_tblite")
    if not isinstance(cp2k, dict) or not isinstance(save_tblite, dict):
        raise ValueError("campaign manifest lacks frozen build identity")
    identity = {
        "cp2k_sha256": str(cp2k.get("binary_sha256")),
        "cp2k_library_sha256": str(cp2k.get("loaded_library_sha256")),
        "tblite_static_library_sha256": str(
            save_tblite.get("static_library_sha256")
        ),
        "cp2k_source_revision": str(cp2k.get("revision")),
        "tblite_source_revision": str(save_tblite.get("revision")),
    }
    for field in BUILD_IDENTITY_FIELDS[:3]:
        if not re.fullmatch(r"[0-9a-f]{64}", identity[field]):
            raise ValueError(f"campaign manifest has invalid frozen SHA256: {field}")
    for field in BUILD_IDENTITY_FIELDS[3:]:
        if not re.fullmatch(r"[0-9a-f]{40}", identity[field]):
            raise ValueError(f"campaign manifest has invalid frozen revision: {field}")
    return identity


def atomic_write_bytes(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _validation_path(value: object, schema_version: int) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("validation record contains an empty artifact path")
    path = Path(value)
    if schema_version == 1:
        return (path if path.is_absolute() else ROOT / path).resolve()
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"schema-v2 validation path is not relative and safe: {value}")
    root = ROOT.resolve()
    resolved = (root / path).resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError(f"schema-v2 validation path escapes artifact root: {value}")
    return resolved


def _read_validation_artifact_bytes(
    value: object,
    expected_hash: object,
    schema_version: int,
    label: str,
) -> tuple[Path, bytes]:
    path = _validation_path(value, schema_version)
    if not isinstance(expected_hash, str) or not re.fullmatch(
        r"[0-9a-f]{64}", expected_hash
    ):
        raise ValueError(f"validation record {label} has invalid SHA256")
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(f"validation record {label} cannot be read: {error}") from error
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError(f"validation record {label} has invalid hash")
    return path, content


def validation_build_id(identity: dict[str, object]) -> str:
    canonical = {
        field: str(identity.get(field)) for field in BUILD_IDENTITY_FIELDS
    }
    content = json.dumps(
        canonical, sort_keys=True, separators=(",", ":")
    ).encode()
    digest = hashlib.sha256(b"dmc13-execution-build-v1\0" + content).hexdigest()
    return f"sha256:{digest}"


def _finite_evidence_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"qualification evidence {label} is not numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"qualification evidence {label} is not finite")
    return number


def _read_hashed_evidence_artifact(
    path_value: object,
    expected_hash: object,
    label: str,
) -> tuple[Path, bytes]:
    if not isinstance(expected_hash, str) or not re.fullmatch(
        r"[0-9a-f]{64}", expected_hash
    ):
        raise ValueError(f"qualification evidence {label} has invalid SHA256")
    path = _validation_path(path_value, 2)
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(
            f"qualification evidence {label} cannot be read: {error}"
        ) from error
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError(f"qualification evidence {label} hash mismatch")
    return path, content


def _gxtb_evidence_input_contract_errors(text: str, mesh: str) -> list[str]:
    normalised_lines = [
        line.strip().upper()
        for line in text.splitlines()
        if line.strip()
    ]
    lines = set(normalised_lines)
    errors: list[str] = []
    for required in (
        f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}".upper(),
        "METHOD XTB",
        "METHOD GXTB",
        "ACCURACY 0.1",
        "SCC_MIXER TBLITE",
        "ITERATIONS 300",
        "EPS_SCF 1.0E-9",
        "METHOD DIRECT_P_MIXING",
        "ALPHA 0.2",
        "CANONICALIZE TRUE",
    ):
        count = normalised_lines.count(required)
        if count == 0:
            errors.append(f"missing {required}")
        elif count != 1:
            errors.append(f"duplicate critical setting {required}")
    for forbidden_method in ("METHOD GFN1", "METHOD GFN2"):
        if forbidden_method in lines:
            errors.append(f"conflicting tblite method {forbidden_method}")
    if mesh == "gamma":
        if "&KPOINTS" in lines:
            errors.append("Gamma production input must use implicit Gamma without &KPOINTS")
    else:
        for required in (
            "&KPOINTS",
            f"SCHEME {MESH_BY_ID[mesh]['scheme']}".upper(),
            "SYMMETRY T",
            "FULL_GRID F",
            "SYMMETRY_BACKEND SPGLIB",
            "SYMMETRY_REDUCTION_METHOD SPGLIB",
        ):
            count = normalised_lines.count(required)
            if count == 0:
                errors.append(f"missing {required}")
            elif count != 1:
                errors.append(f"duplicate critical setting {required}")
        for forbidden in ("SYMMETRY F", "FULL_GRID T"):
            if forbidden in lines:
                errors.append(f"forbidden legacy setting {forbidden}")
    return errors


def _cp2k_input_water_count(text: str, label: str) -> int:
    in_coordinates = False
    oxygen_count = 0
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        upper = line.upper()
        if upper.startswith("&COORD"):
            in_coordinates = True
            continue
        if in_coordinates and upper.startswith("&END"):
            in_coordinates = False
            continue
        if in_coordinates and line and line.split()[0].upper() == "O":
            oxygen_count += 1
    if oxygen_count <= 0:
        raise ValueError(
            f"qualification evidence {label} input lacks explicit water oxygens"
        )
    return oxygen_count


def _validate_evidence_input(
    path: Path,
    content: bytes,
    *,
    mesh: str,
    phase: str,
    label: str,
) -> int:
    expected_project = f"ice_{phase}_GXTB_{mesh}"
    if path.name != f"{expected_project}.inp":
        raise ValueError(
            f"qualification evidence {label} input filename mismatch"
        )
    text = content.decode(errors="ignore")
    lines = {line.strip() for line in text.splitlines() if line.strip()}
    if f"PROJECT {expected_project}" not in lines:
        raise ValueError(f"qualification evidence {label} project mismatch")
    errors = _gxtb_evidence_input_contract_errors(text, mesh)
    if errors:
        raise ValueError(
            f"qualification evidence {label} violates input contract: "
            + "; ".join(errors)
        )
    return _cp2k_input_water_count(text, label)


def _evidence_total_energy(
    content: bytes,
    label: str,
    *,
    expected_project: str,
    expected_source_revision: str,
    expected_tblite_source_revision: str,
    allow_unknown_tblite_revision: bool,
) -> float:
    text = content.decode(errors="ignore")
    if (
        "PROGRAM ENDED" not in text
        or "ENERGY| Total FORCE_EVAL" not in text
        or "SCF run converged" not in text
        or "SCF run NOT converged" in text
        or "ABORT" in text
    ):
        raise ValueError(f"qualification evidence {label} is not a completed output")
    input_matches = re.findall(
        r"^\s*CP2K\|\s+Input file name\s+(\S+)\s*$", text, re.MULTILINE
    )
    project_matches = re.findall(
        r"^\s*GLOBAL\|\s+Project name\s+(\S+)\s*$", text, re.MULTILINE
    )
    revision_matches = re.findall(
        r"^\s*CP2K\|\s+source code revision number:\s*"
        r"([0-9a-fA-F]{7,40})\s*$",
        text,
        re.MULTILINE,
    )
    if len(input_matches) != 1 or input_matches[0] != f"{expected_project}.inp":
        raise ValueError(f"qualification evidence {label} input header mismatch")
    if len(project_matches) != 1 or project_matches[0] != expected_project:
        raise ValueError(f"qualification evidence {label} project header mismatch")
    if "tblite_gxtb" not in text.lower():
        raise ValueError(f"qualification evidence {label} lacks tblite_gxtb header")
    if (
        len(revision_matches) != 1
        or not expected_source_revision.startswith(revision_matches[0].lower())
    ):
        raise ValueError(f"qualification evidence {label} source revision mismatch")
    tblite_revisions = re.findall(
        r"^\s*tblite source revision:\s*(\S+)\s*$",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if len(tblite_revisions) != 1:
        raise ValueError(
            f"qualification evidence {label} must contain exactly one tblite "
            "source revision"
        )
    tblite_revision = tblite_revisions[0].lower()
    if tblite_revision != expected_tblite_source_revision.lower() and not (
        allow_unknown_tblite_revision and tblite_revision == "unknown"
    ):
        raise ValueError(
            f"qualification evidence {label} tblite source revision mismatch"
        )
    energy_lines = [
        line for line in text.splitlines() if "ENERGY| Total FORCE_EVAL" in line
    ]
    if len(energy_lines) != 1:
        raise ValueError(
            f"qualification evidence {label} must contain exactly one total energy"
        )
    try:
        energy = float(energy_lines[0].split()[-1])
    except (ValueError, IndexError) as error:
        raise ValueError(
            f"qualification evidence {label} has an invalid energy"
        ) from error
    if not math.isfinite(energy):
        raise ValueError(f"qualification evidence {label} lacks a finite energy")
    return energy


def _evidence_execution_environment(
    content: bytes, label: str
) -> tuple[str, str, str, str]:
    text = content.decode(errors="ignore")
    fields = {
        "start time": re.findall(
            r"PROGRAM STARTED AT\s+(.+?)\s*$", text, re.MULTILINE
        ),
        "host": re.findall(
            r"PROGRAM STARTED ON\s+(\S+)\s*$", text, re.MULTILINE
        ),
        "compiled host": re.findall(
            r"^\s*CP2K\|\s+Program compiled on\s*(.*?)\s*$",
            text,
            re.MULTILINE,
        ),
        "platform": re.findall(
            r"^\s*CP2K\|\s+Program compiled for\s+(\S+)\s*$",
            text,
            re.MULTILINE,
        ),
    }
    for field, matches in fields.items():
        if len(matches) != 1 or (
            field != "compiled host" and not matches[0].strip()
        ):
            raise ValueError(
                f"qualification evidence {label} lacks unique execution {field}"
            )
    return (
        fields["host"][0].strip(),
        fields["compiled host"][0].strip(),
        fields["platform"][0].strip(),
        fields["start time"][0].strip(),
    )


def _validate_qualification_run_stamp(
    content: bytes,
    label: str,
    *,
    campaign_id: str,
    identity: dict[str, object],
    identity_id: str,
    mesh: str,
    phase: str,
    input_name: str,
    input_sha256: str,
    output_name: str,
    output_sha256: str,
    require_schema_v2: bool,
) -> None:
    try:
        stamp = json.loads(content)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"qualification evidence {label} stamp is invalid: {error}"
        ) from error
    if not isinstance(stamp, dict):
        raise ValueError(f"qualification evidence {label} stamp is not an object")
    expected = {
        "campaign_id": campaign_id,
        "method": "GXTB",
        "mesh": mesh,
        "phase": phase,
        "input": input_name,
        "input_sha256": input_sha256,
        "output": output_name,
        "output_sha256": output_sha256,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "input_contract_valid": True,
        "adopted_existing_output": False,
        **{field: identity[field] for field in BUILD_IDENTITY_FIELDS},
    }
    if require_schema_v2:
        expected.update(
            {
                "schema_version": 2,
                "build_id": identity_id,
                "frozen_input": input_name,
                "frozen_input_sha256": input_sha256,
            }
        )
    else:
        stamp_identity = {
            field: str(stamp.get(field)) for field in BUILD_IDENTITY_FIELDS
        }
        if validation_build_id(stamp_identity) != identity_id:
            raise ValueError(
                f"qualification evidence {label} derived build mismatch"
            )
        if stamp.get("build_id", identity_id) != identity_id:
            raise ValueError(f"qualification evidence {label} build mismatch")
    mismatches = [
        field for field, expected_value in expected.items()
        if stamp.get(field) != expected_value
    ]
    if mismatches:
        raise ValueError(
            f"qualification evidence {label} stamp mismatch: "
            + ", ".join(mismatches)
        )


def _validate_qualification_evidence(
    qualification: dict[str, object],
    *,
    campaign_id: str,
    remote_identity: dict[str, object],
    reference_identity: dict[str, object],
    reference_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    remote_build_id = validation_build_id(remote_identity)
    reference_build_id = validation_build_id(reference_identity)
    remote_cp2k_source_revision = str(remote_identity["cp2k_source_revision"])
    reference_cp2k_source_revision = str(
        reference_identity["cp2k_source_revision"]
    )
    if qualification.get("status") != "passed":
        raise ValueError("qualification.status must be 'passed'")
    if qualification.get("evidence_schema_version") != QUALIFICATION_EVIDENCE_SCHEMA_VERSION:
        raise ValueError("qualification evidence_schema_version mismatch")
    total_tolerance = _finite_evidence_number(
        qualification.get("total_energy_tolerance_hartree"),
        "total_energy_tolerance_hartree",
    )
    relative_tolerance = _finite_evidence_number(
        qualification.get("relative_energy_tolerance_kjmol_per_h2o"),
        "relative_energy_tolerance_kjmol_per_h2o",
    )
    if not 0.0 < total_tolerance <= MAX_TOTAL_ENERGY_TOLERANCE_HARTREE:
        raise ValueError("qualification total-energy tolerance is looser than 1e-10 Eh")
    if not 0.0 < relative_tolerance <= MAX_RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O:
        raise ValueError(
            "qualification relative-energy tolerance is looser than 0.001 kJ/mol/H2O"
        )
    observed_total = _finite_evidence_number(
        qualification.get("observed_max_abs_total_energy_delta_hartree"),
        "observed_max_abs_total_energy_delta_hartree",
    )
    observed_relative = _finite_evidence_number(
        qualification.get(
            "observed_max_abs_relative_energy_delta_kjmol_per_h2o"
        ),
        "observed_max_abs_relative_energy_delta_kjmol_per_h2o",
    )
    sentinels = qualification.get("same_mesh_dense_relative_sentinels")
    count = qualification.get("same_mesh_dense_relative_sentinel_count")
    if (
        not isinstance(sentinels, list)
        or not sentinels
        or isinstance(count, bool)
        or not isinstance(count, int)
        or count != len(sentinels)
    ):
        raise ValueError(
            "qualification requires at least one counted same-mesh dense sentinel"
        )
    total_deltas: list[float] = []
    relative_deltas: list[float] = []
    for index, value in enumerate(sentinels):
        label = f"sentinel[{index}]"
        if not isinstance(value, dict):
            raise ValueError(f"qualification evidence {label} is not an object")
        if value.get("kind") != "same_mesh_dense_relative_energy":
            raise ValueError(f"qualification evidence {label} kind mismatch")
        if value.get("mesh") not in {
            str(mesh["id"]) for mesh in DENSE_EXTENSION_MESHES
        }:
            raise ValueError(f"qualification evidence {label} is not a dense mesh")
        if value.get("phase") not in PHASES[1:]:
            raise ValueError(f"qualification evidence {label} lacks a non-Ih phase")
        if value.get("remote_build_id") != remote_build_id:
            raise ValueError(f"qualification evidence {label} remote build mismatch")
        if value.get("reference_build_id") != reference_build_id:
            raise ValueError(
                f"qualification evidence {label} reference build mismatch"
            )
        if remote_build_id == reference_build_id:
            raise ValueError(
                f"qualification evidence {label} does not compare distinct builds"
            )
        resolved: dict[str, Path] = {}
        artifact_content: dict[str, bytes] = {}
        artifact_hashes: dict[str, str] = {}
        for artifact in (
            "phase_input",
            "ih_input",
            "remote_phase_output",
            "remote_ih_output",
            "reference_phase_output",
            "reference_ih_output",
            "remote_phase_stamp",
            "remote_ih_stamp",
            "reference_phase_stamp",
            "reference_ih_stamp",
        ):
            artifact_hash = value.get(f"{artifact}_sha256")
            artifact_path, content = _read_hashed_evidence_artifact(
                value.get(artifact), artifact_hash, f"{label}/{artifact}"
            )
            resolved[artifact] = artifact_path
            artifact_content[artifact] = content
            artifact_hashes[artifact] = str(artifact_hash)
        if (
            resolved["phase_input"] == resolved["ih_input"]
            or artifact_hashes["phase_input"] == artifact_hashes["ih_input"]
        ):
            raise ValueError(
                f"qualification evidence {label} phase/Ih inputs are not distinct"
            )
        output_artifacts = (
            "remote_phase_output",
            "remote_ih_output",
            "reference_phase_output",
            "reference_ih_output",
        )
        if len({resolved[name] for name in output_artifacts}) != len(
            output_artifacts
        ):
            raise ValueError(
                f"qualification evidence {label} output paths are not distinct"
            )
        for side in ("remote", "reference"):
            if (
                artifact_hashes[f"{side}_phase_output"]
                == artifact_hashes[f"{side}_ih_output"]
            ):
                raise ValueError(
                    f"qualification evidence {label} {side} phase/Ih outputs "
                    "are not distinct"
                )
        for system in ("phase", "ih"):
            if (
                artifact_hashes[f"remote_{system}_output"]
                == artifact_hashes[f"reference_{system}_output"]
            ):
                raise ValueError(
                    f"qualification evidence {label} copied reference {system} "
                    "output"
                )
        stamp_artifacts = (
            "remote_phase_stamp",
            "remote_ih_stamp",
            "reference_phase_stamp",
            "reference_ih_stamp",
        )
        if len({resolved[name] for name in stamp_artifacts}) != len(
            stamp_artifacts
        ):
            raise ValueError(
                f"qualification evidence {label} stamp paths are not distinct"
            )
        execution_environments = {
            (side, system): _evidence_execution_environment(
                artifact_content[f"{side}_{system}_output"],
                f"{label}/{side}_{system}_output",
            )
            for side in ("remote", "reference")
            for system in ("phase", "ih")
        }
        for side in ("remote", "reference"):
            if execution_environments[(side, "phase")][:3] != (
                execution_environments[(side, "ih")][:3]
            ):
                raise ValueError(
                    f"qualification evidence {label} {side} execution "
                    "environment mismatch"
                )
        if (
            execution_environments[("remote", "phase")][0]
            == execution_environments[("reference", "phase")][0]
            or execution_environments[("remote", "phase")][2]
            == execution_environments[("reference", "phase")][2]
        ):
            raise ValueError(
                f"qualification evidence {label} remote/reference execution "
                "environments are not distinct"
            )
        expected_remote_environment = qualification.get(
            "remote_execution_environment"
        )
        if not isinstance(expected_remote_environment, dict):
            raise ValueError(
                "qualification lacks remote_execution_environment"
            )
        expected_environment = (
            expected_remote_environment.get("program_started_on"),
            expected_remote_environment.get("program_compiled_on"),
            expected_remote_environment.get("program_compiled_for"),
        )
        if any(not isinstance(item, str) or not item for item in expected_environment):
            raise ValueError(
                "qualification remote_execution_environment is invalid"
            )
        if execution_environments[("remote", "phase")][:3] != expected_environment:
            raise ValueError(
                f"qualification evidence {label} remote execution environment "
                "does not match manifest"
            )
        counts: dict[str, int] = {}
        phase_name = str(value["phase"])
        mesh = str(value["mesh"])
        for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
            count_value = value.get(f"{system}_water_count")
            if (
                isinstance(count_value, bool)
                or not isinstance(count_value, int)
                or count_value <= 0
            ):
                raise ValueError(
                    f"qualification evidence {label}/{system}_water_count is invalid"
                )
            derived_count = _validate_evidence_input(
                resolved[f"{system}_input"],
                artifact_content[f"{system}_input"],
                mesh=mesh,
                phase=system_phase,
                label=f"{label}/{system}_input",
            )
            if count_value != derived_count:
                raise ValueError(
                    f"qualification evidence {label}/{system}_water_count "
                    "does not match input"
                )
            counts[system] = count_value
        for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
            reference_record = reference_records.get((mesh, system_phase))
            if reference_record is None:
                raise ValueError(
                    f"qualification evidence {label} lacks trusted reference "
                    f"record {mesh}/{system_phase}"
                )
            if reference_record.get("build_id") != reference_build_id:
                raise ValueError(
                    f"qualification evidence {label} reference record build mismatch"
                )
            input_hash = artifact_hashes[f"{system}_input"]
            if (
                reference_record.get("input_sha256") != input_hash
                or reference_record.get("frozen_input_sha256") != input_hash
            ):
                raise ValueError(
                    f"qualification evidence {label} {system} input is not "
                    "canonical reference bytes"
                )
            expected_project = f"ice_{system_phase}_GXTB_{mesh}"
            expected_input_name = f"{expected_project}.inp"
            expected_output_name = f"{expected_project}.out"
            remote_output = resolved[f"remote_{system}_output"]
            remote_stamp = resolved[f"remote_{system}_stamp"]
            if remote_output.name != expected_output_name:
                raise ValueError(
                    f"qualification evidence {label} remote {system} output "
                    "filename mismatch"
                )
            if resolved[f"{system}_input"] != remote_output.parent / expected_input_name:
                raise ValueError(
                    f"qualification evidence {label} remote {system} input path mismatch"
                )
            if remote_stamp != remote_output.with_suffix(".run.json"):
                raise ValueError(
                    f"qualification evidence {label} remote {system} stamp path mismatch"
                )
            reference_output_path = _validation_path(
                reference_record.get("output"), 2
            )
            reference_stamp_path = _validation_path(
                reference_record.get("stamp"), 2
            )
            if (
                resolved[f"reference_{system}_output"] != reference_output_path
                or artifact_hashes[f"reference_{system}_output"]
                != reference_record.get("output_sha256")
            ):
                raise ValueError(
                    f"qualification evidence {label} reference {system} output "
                    "does not match trusted record"
                )
            if (
                resolved[f"reference_{system}_stamp"] != reference_stamp_path
                or artifact_hashes[f"reference_{system}_stamp"]
                != reference_record.get("stamp_sha256")
            ):
                raise ValueError(
                    f"qualification evidence {label} reference {system} stamp "
                    "does not match trusted record"
                )
            _validate_qualification_run_stamp(
                artifact_content[f"remote_{system}_stamp"],
                f"{label}/remote_{system}_stamp",
                campaign_id=campaign_id,
                identity=remote_identity,
                identity_id=remote_build_id,
                mesh=mesh,
                phase=system_phase,
                input_name=expected_input_name,
                input_sha256=input_hash,
                output_name=expected_output_name,
                output_sha256=artifact_hashes[f"remote_{system}_output"],
                require_schema_v2=True,
            )
            _validate_qualification_run_stamp(
                artifact_content[f"reference_{system}_stamp"],
                f"{label}/reference_{system}_stamp",
                campaign_id=campaign_id,
                identity=reference_identity,
                identity_id=reference_build_id,
                mesh=mesh,
                phase=system_phase,
                input_name=expected_input_name,
                input_sha256=input_hash,
                output_name=expected_output_name,
                output_sha256=artifact_hashes[f"reference_{system}_output"],
                require_schema_v2=False,
            )
        conversion = _finite_evidence_number(
            value.get("hartree_to_kjmol"), f"{label}/hartree_to_kjmol"
        )
        if conversion != HARTREE_TO_KJMOL:
            raise ValueError(
                f"qualification evidence {label} Hartree conversion mismatch"
            )
        numeric_energies: dict[str, float] = {}
        for side, source_revision in (
            ("remote", remote_cp2k_source_revision),
            ("reference", reference_cp2k_source_revision),
        ):
            for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
                output = f"{side}_{system}"
                numeric_energies[output] = _evidence_total_energy(
                    artifact_content[f"{output}_output"],
                    f"{label}/{output}_output",
                    expected_project=f"ice_{system_phase}_GXTB_{mesh}",
                    expected_source_revision=source_revision,
                    expected_tblite_source_revision=str(
                        (remote_identity if side == "remote" else reference_identity)[
                            "tblite_source_revision"
                        ]
                    ),
                    allow_unknown_tblite_revision=(side == "reference"),
                )
        phase_total_delta = abs(
            numeric_energies["remote_phase"]
            - numeric_energies["reference_phase"]
        )
        ih_total_delta = abs(
            numeric_energies["remote_ih"] - numeric_energies["reference_ih"]
        )
        for system, computed_total in (
            ("phase", phase_total_delta),
            ("ih", ih_total_delta),
        ):
            declared_total = _finite_evidence_number(
                value.get(f"{system}_total_energy_delta_hartree"),
                f"{label}/{system}_total_energy_delta_hartree",
            )
            if declared_total < 0.0 or not math.isclose(
                declared_total,
                computed_total,
                rel_tol=1.0e-12,
                abs_tol=1.0e-14,
            ):
                raise ValueError(
                    f"qualification evidence {label} {system} total-energy mismatch"
                )
            if declared_total > total_tolerance:
                raise ValueError(f"qualification evidence {label} exceeds tolerance")
            total_deltas.append(declared_total)
        computed_remote_relative = conversion * (
            numeric_energies["remote_phase"] / counts["phase"]
            - numeric_energies["remote_ih"] / counts["ih"]
        )
        computed_reference_relative = conversion * (
            numeric_energies["reference_phase"] / counts["phase"]
            - numeric_energies["reference_ih"] / counts["ih"]
        )
        declared_remote_relative = _finite_evidence_number(
            value.get("remote_relative_energy_kjmol_per_h2o"),
            f"{label}/remote_relative_energy_kjmol_per_h2o",
        )
        declared_reference_relative = _finite_evidence_number(
            value.get("reference_relative_energy_kjmol_per_h2o"),
            f"{label}/reference_relative_energy_kjmol_per_h2o",
        )
        for side, declared, computed in (
            ("remote", declared_remote_relative, computed_remote_relative),
            ("reference", declared_reference_relative, computed_reference_relative),
        ):
            if not math.isclose(declared, computed, rel_tol=1.0e-12, abs_tol=1.0e-9):
                raise ValueError(
                    f"qualification evidence {label} {side} relative-energy mismatch"
                )
        declared_relative = _finite_evidence_number(
            value.get("relative_energy_delta_kjmol_per_h2o"),
            f"{label}/relative_energy_delta_kjmol_per_h2o",
        )
        if declared_relative < 0.0 or not math.isclose(
            declared_relative,
            abs(computed_remote_relative - computed_reference_relative),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        ):
            raise ValueError(f"qualification evidence {label} relative-energy mismatch")
        if declared_relative > relative_tolerance:
            raise ValueError(f"qualification evidence {label} exceeds tolerance")
        relative_deltas.append(declared_relative)
    if observed_total < 0.0 or not math.isclose(
        observed_total, max(total_deltas), rel_tol=1.0e-12, abs_tol=1.0e-14
    ):
        raise ValueError("qualification observed total-energy maximum mismatch")
    if observed_relative < 0.0 or not math.isclose(
        observed_relative,
        max(relative_deltas),
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ValueError("qualification observed relative-energy maximum mismatch")
    if observed_total > total_tolerance or observed_relative > relative_tolerance:
        raise ValueError("qualification observed maximum exceeds tolerance")


def _validate_execution_manifest_payload(
    manifest: dict[str, object],
    campaign_id: str,
    campaign_manifest_sha256: str,
    identity_id: str,
    identity: dict[str, object],
    reference_identity: dict[str, object],
    reference_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    manifest_identity = manifest.get("build_identity")
    qualification = manifest.get("qualification")
    if not isinstance(manifest_identity, dict) or not isinstance(qualification, dict):
        raise ValueError("execution-build manifest lacks identity/qualification")
    checks = {
        "schema_version": (manifest.get("schema_version"), 1),
        "campaign_id": (manifest.get("campaign_id"), campaign_id),
        "campaign_manifest_sha256": (
            manifest.get("campaign_manifest_sha256"),
            campaign_manifest_sha256,
        ),
        "gxtb_protocol_id": (manifest.get("gxtb_protocol_id"), GXTB_PROTOCOL_ID),
        "build_id": (manifest.get("build_id"), identity_id),
    }
    for field in BUILD_IDENTITY_FIELDS:
        checks[field] = (manifest_identity.get(field), identity.get(field))
    mismatches = [key for key, pair in checks.items() if pair[0] != pair[1]]
    cli_hash = manifest_identity.get("tblite_cli_sha256")
    if not isinstance(cli_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", cli_hash):
        mismatches.append("tblite_cli_sha256")
    if mismatches:
        raise ValueError("execution-build manifest mismatch: " + ", ".join(mismatches))
    _validate_qualification_evidence(
        qualification,
        campaign_id=campaign_id,
        remote_identity=identity,
        reference_identity=reference_identity,
        reference_records=reference_records,
    )


def read_validation_index(
    path: Path,
    *,
    campaign_manifest_path: Path | None = None,
) -> dict[str, object]:
    """Read and independently verify the runner's hash-valid phase index."""
    try:
        index_bytes = path.read_bytes()
        payload = json.loads(index_bytes)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read validation index {path}: {error}") from error
    if not isinstance(payload, dict) or payload.get("schema_version") not in (1, 2):
        raise ValueError("validation index must be a schema-version-1/2 JSON object")
    if payload.get("benchmark") != "DMC-ICE13" or payload.get("method") != "GXTB":
        raise ValueError("validation index benchmark/method mismatch")
    source_schema_version = int(payload["schema_version"])
    if payload.get("gxtb_protocol_id") != GXTB_PROTOCOL_ID:
        raise ValueError("validation index g-xTB protocol mismatch")
    records = payload.get("records")
    if not isinstance(records, list):
        raise ValueError("validation index lacks records")
    pending_execution_manifests: list[
        tuple[dict[str, object], str, dict[str, object]]
    ] = []
    if source_schema_version == 1:
        legacy_identity = payload.get("build_identity")
        if not isinstance(legacy_identity, dict):
            raise ValueError("schema-v1 validation index lacks build_identity")
        identity = {
            field: str(legacy_identity.get(field))
            for field in BUILD_IDENTITY_FIELDS
        }
        for field in BUILD_IDENTITY_FIELDS[:3]:
            if not re.fullmatch(r"[0-9a-f]{64}", identity[field]):
                raise ValueError(f"schema-v1 build identity has invalid SHA256: {field}")
        for field in BUILD_IDENTITY_FIELDS[3:]:
            if not re.fullmatch(r"[0-9a-f]{40}", identity[field]):
                raise ValueError(f"schema-v1 build identity has invalid revision: {field}")
        identity_id = validation_build_id(identity)
        build_identities: dict[str, dict[str, object]] = {identity_id: identity}
        source_identity = {
            "cp2k_source_revision": identity["cp2k_source_revision"],
            "tblite_source_revision": identity["tblite_source_revision"],
        }
    else:
        raw_identities = payload.get("build_identities")
        source_identity = payload.get("source_identity")
        campaign_manifest_hash = payload.get("campaign_manifest_sha256")
        if not isinstance(raw_identities, dict) or not isinstance(source_identity, dict):
            raise ValueError("schema-v2 validation index lacks build/source identity")
        if not isinstance(campaign_manifest_hash, str) or not re.fullmatch(
            r"[0-9a-f]{64}", campaign_manifest_hash
        ):
            raise ValueError("schema-v2 campaign_manifest_sha256 is invalid")
        local_campaign_manifest = (
            campaign_manifest_path
            if campaign_manifest_path is not None
            else ROOT.parent
            / "campaigns"
            / str(payload.get("campaign_id"))
            / "build_manifest.json"
        )
        if not local_campaign_manifest.is_file():
            raise ValueError(
                "schema-v2 trusted campaign manifest does not exist: "
                f"{local_campaign_manifest}"
            )
        trusted_manifest_bytes = local_campaign_manifest.read_bytes()
        if campaign_manifest_hash != hashlib.sha256(trusted_manifest_bytes).hexdigest():
            raise ValueError("schema-v2 campaign manifest hash mismatch")
        try:
            trusted_manifest = json.loads(trusted_manifest_bytes)
        except json.JSONDecodeError as error:
            raise ValueError(
                f"cannot read campaign manifest {local_campaign_manifest}: {error}"
            ) from error
        if not isinstance(trusted_manifest, dict):
            raise ValueError(
                f"campaign manifest {local_campaign_manifest} is not a JSON object"
            )
        if trusted_manifest.get("campaign_id") != payload.get("campaign_id"):
            raise ValueError("schema-v2 trusted campaign manifest campaign mismatch")
        frozen_identity = frozen_build_identity_from_manifest(trusted_manifest)
        if {
            "cp2k_source_revision": str(
                source_identity.get("cp2k_source_revision")
            ),
            "tblite_source_revision": str(
                source_identity.get("tblite_source_revision")
            ),
        } != {
            "cp2k_source_revision": frozen_identity["cp2k_source_revision"],
            "tblite_source_revision": frozen_identity["tblite_source_revision"],
        }:
            raise ValueError("schema-v2 source revisions differ from campaign manifest")
        build_identities = {}
        frozen_identity_ids: set[str] = set()
        for identity_id, value in raw_identities.items():
            if not isinstance(identity_id, str) or not isinstance(value, dict):
                raise ValueError("validation build identity is malformed")
            if validation_build_id(value) != identity_id:
                raise ValueError(f"validation build identity digest mismatch: {identity_id}")
            for field in BUILD_IDENTITY_FIELDS[:3]:
                if not re.fullmatch(r"[0-9a-f]{64}", str(value.get(field))):
                    raise ValueError(
                        f"validation build identity has invalid SHA256: {identity_id}/{field}"
                    )
            for field in BUILD_IDENTITY_FIELDS[3:]:
                if not re.fullmatch(r"[0-9a-f]{40}", str(value.get(field))):
                    raise ValueError(
                        f"validation build identity has invalid revision: {identity_id}/{field}"
                    )
            for field in ("cp2k_source_revision", "tblite_source_revision"):
                if value.get(field) != source_identity.get(field):
                    raise ValueError(
                        f"validation build identity source mismatch: {identity_id}/{field}"
                    )
            manifest_path_value = value.get("execution_build_manifest")
            manifest_hash = value.get("execution_build_manifest_sha256")
            if (manifest_path_value is None) != (manifest_hash is None):
                raise ValueError(
                    f"validation build identity has incomplete execution manifest: {identity_id}"
                )
            is_frozen_identity = all(
                value.get(field) == frozen_identity[field]
                for field in BUILD_IDENTITY_FIELDS
            )
            if is_frozen_identity:
                frozen_identity_ids.add(identity_id)
            if manifest_path_value is None and not is_frozen_identity:
                raise ValueError(
                    "alternate validation build identity lacks execution manifest: "
                    f"{identity_id}"
                )
            if manifest_path_value is not None:
                if not isinstance(manifest_hash, str) or not re.fullmatch(
                    r"[0-9a-f]{64}", manifest_hash
                ):
                    raise ValueError(
                        f"validation execution manifest has invalid SHA256: {identity_id}"
                    )
                resolved_manifest = _validation_path(manifest_path_value, 2)
                try:
                    execution_manifest_bytes = resolved_manifest.read_bytes()
                except OSError as error:
                    raise ValueError(
                        f"validation execution manifest cannot be read: "
                        f"{identity_id}: {error}"
                    ) from error
                if hashlib.sha256(execution_manifest_bytes).hexdigest() != manifest_hash:
                    raise ValueError(
                        f"validation execution manifest hash mismatch: {identity_id}"
                    )
                try:
                    execution_manifest = json.loads(execution_manifest_bytes)
                except json.JSONDecodeError as error:
                    raise ValueError(
                        f"invalid execution manifest for {identity_id}: {error}"
                    ) from error
                if not isinstance(execution_manifest, dict):
                    raise ValueError("execution-build manifest is not an object")
                pending_execution_manifests.append(
                    (execution_manifest, identity_id, dict(value))
                )
            build_identities[identity_id] = dict(value)
        if sha256(local_campaign_manifest) != campaign_manifest_hash:
            raise ValueError(
                "schema-v2 trusted campaign manifest changed during validation"
            )
    if source_schema_version == 1:
        frozen_identity_ids = set()
    seen: set[tuple[str, str]] = set()
    coverage: dict[str, list[str]] = {}
    normalised_records: list[dict[str, object]] = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("validation index record is not an object")
        record = dict(record)
        mesh = str(record.get("mesh"))
        phase = str(record.get("phase"))
        if mesh not in MESH_BY_ID or phase not in PHASES:
            raise ValueError(f"unsupported validation record {mesh}/{phase}")
        if (mesh, phase) in seen:
            raise ValueError(f"duplicate validation record {mesh}/{phase}")
        seen.add((mesh, phase))
        if source_schema_version == 1:
            record["campaign_id"] = payload.get("campaign_id")
            record["gxtb_protocol_id"] = GXTB_PROTOCOL_ID
            record["build_id"] = identity_id
        if record.get("campaign_id") != payload.get("campaign_id"):
            raise ValueError(f"validation record campaign mismatch: {mesh}/{phase}")
        if record.get("gxtb_protocol_id") != GXTB_PROTOCOL_ID:
            raise ValueError(f"validation record protocol mismatch: {mesh}/{phase}")
        record_identity_id = str(record.get("build_id"))
        identity = build_identities.get(record_identity_id)
        if identity is None:
            raise ValueError(f"validation record has unknown build: {mesh}/{phase}")
        frozen_input = record.get("frozen_input")
        if frozen_input is None:
            if source_schema_version == 2:
                raise ValueError(
                    f"schema-v2 validation record lacks frozen_input: {mesh}/{phase}"
                )
            frozen_input = str(
                Path(str(record.get("output"))).parent
                / Path(str(record.get("input"))).name
            )
        frozen_input_hash = record.get(
            "frozen_input_sha256", record.get("input_sha256")
        )
        if source_schema_version == 2 and "frozen_input_sha256" not in record:
            raise ValueError(
                f"schema-v2 validation record lacks frozen_input_sha256: {mesh}/{phase}"
            )
        if frozen_input_hash != record.get("input_sha256"):
            raise ValueError(
                f"validation record {mesh}/{phase} source/frozen input hash mismatch"
            )
        files = {
            "input": (record.get("input"), record.get("input_sha256")),
            "frozen_input": (frozen_input, frozen_input_hash),
            "output": (record.get("output"), record.get("output_sha256")),
            "stamp": (record.get("stamp"), record.get("stamp_sha256")),
        }
        resolved: dict[str, Path] = {}
        artifact_bytes: dict[str, bytes] = {}
        for label, (file_name, expected_hash) in files.items():
            file_path, content = _read_validation_artifact_bytes(
                file_name,
                expected_hash,
                source_schema_version,
                f"{mesh}/{phase}/{label}",
            )
            resolved[label] = file_path
            artifact_bytes[label] = content
        expected_project = f"ice_{phase}_GXTB_{mesh}"
        expected_input_name = f"{expected_project}.inp"
        expected_output_name = f"{expected_project}.out"
        if source_schema_version in (1, 2):
            canonical_input = (
                ROOT / "gxtb_spglib_inputs" / mesh / expected_input_name
            ).resolve()
            canonical_output = (
                ROOT
                / "runs_gxtb_spglib"
                / mesh
                / phase
                / expected_output_name
            ).resolve()
            if resolved["input"] != canonical_input:
                raise ValueError(
                    f"validation record {mesh}/{phase} input path is not canonical"
                )
            if resolved["output"] != canonical_output:
                raise ValueError(
                    f"validation record {mesh}/{phase} output path is not canonical"
                )
            if resolved["frozen_input"] != canonical_output.parent / expected_input_name:
                raise ValueError(
                    f"validation record {mesh}/{phase} frozen input path mismatch"
                )
            if resolved["stamp"] != canonical_output.with_suffix(".run.json"):
                raise ValueError(
                    f"validation record {mesh}/{phase} stamp path mismatch"
                )
            if resolved["input"] == resolved["frozen_input"]:
                raise ValueError(
                    f"validation record {mesh}/{phase} source/frozen paths coincide"
                )
            if artifact_bytes["input"] != artifact_bytes["frozen_input"]:
                raise ValueError(
                    f"validation record {mesh}/{phase} source/frozen bytes differ"
                )
            for input_label in ("input", "frozen_input"):
                _validate_evidence_input(
                    resolved[input_label],
                    artifact_bytes[input_label],
                    mesh=mesh,
                    phase=phase,
                    label=f"record {mesh}/{phase}/{input_label}",
                )
        validated_energy = _evidence_total_energy(
            artifact_bytes["output"],
            f"record {mesh}/{phase}/output",
            expected_project=expected_project,
            expected_source_revision=str(identity["cp2k_source_revision"]),
            expected_tblite_source_revision=str(
                identity["tblite_source_revision"]
            ),
            allow_unknown_tblite_revision=(
                source_schema_version == 1
                or record_identity_id in frozen_identity_ids
            ),
        )
        try:
            stamp = json.loads(artifact_bytes["stamp"])
        except json.JSONDecodeError as error:
            raise ValueError(f"invalid stamp for {mesh}/{phase}: {error}") from error
        if not isinstance(stamp, dict):
            raise ValueError(f"invalid stamp for {mesh}/{phase}: not an object")
        expected_stamp_fields = {
            "campaign_id": payload.get("campaign_id"),
            **{field: identity.get(field) for field in BUILD_IDENTITY_FIELDS},
            "method": "GXTB",
            "mesh": mesh,
            "phase": phase,
            "gxtb_protocol_id": GXTB_PROTOCOL_ID,
            "input_sha256": record.get("input_sha256"),
            "output_sha256": record.get("output_sha256"),
            "input": expected_input_name,
            "output": expected_output_name,
            "input_contract_valid": True,
            "adopted_existing_output": False,
        }
        stamp_version = stamp.get("schema_version", 1)
        requires_v2_stamp = (
            source_schema_version == 2
            and record_identity_id not in frozen_identity_ids
        )
        if requires_v2_stamp and stamp_version != 2:
            raise ValueError(f"stamp mismatch for {mesh}/{phase}: schema_version")
        if stamp_version == 2:
            expected_stamp_fields.update(
                {
                    "schema_version": 2,
                    "build_id": record_identity_id,
                    "frozen_input": expected_input_name,
                    "frozen_input_sha256": frozen_input_hash,
                }
            )
        mismatches = [
            key
            for key, expected in expected_stamp_fields.items()
            if stamp.get(key) != expected
        ]
        if mismatches:
            raise ValueError(
                f"stamp mismatch for {mesh}/{phase}: {', '.join(mismatches)}"
            )
        if (
            "frozen_input" in stamp
            and stamp.get("frozen_input") != resolved["frozen_input"].name
        ):
            raise ValueError(f"stamp mismatch for {mesh}/{phase}: frozen_input")
        if (
            "frozen_input_sha256" in stamp
            and stamp.get("frozen_input_sha256") != frozen_input_hash
        ):
            raise ValueError(
                f"stamp mismatch for {mesh}/{phase}: frozen_input_sha256"
            )
        record["frozen_input"] = str(frozen_input)
        record["frozen_input_sha256"] = frozen_input_hash
        record["validated_energy_hartree"] = validated_energy
        normalised_records.append(record)
        coverage.setdefault(mesh, []).append(phase)
    reference_records = {
        (str(record["mesh"]), str(record["phase"])): record
        for record in normalised_records
        if str(record.get("build_id")) in frozen_identity_ids
    }
    for execution_manifest, execution_identity_id, execution_identity in (
        pending_execution_manifests
    ):
        _validate_execution_manifest_payload(
            execution_manifest,
            str(payload.get("campaign_id")),
            str(campaign_manifest_hash),
            execution_identity_id,
            execution_identity,
            frozen_identity,
            reference_records,
        )
    coverage = {
        mesh: [phase for phase in PHASES if phase in phases]
        for mesh, phases in coverage.items()
    }
    used_identity_ids = {str(record["build_id"]) for record in normalised_records}
    if source_schema_version == 2 and set(build_identities) != used_identity_ids:
        raise ValueError("schema-v2 validation index contains unused build identities")
    if source_schema_version == 2 and payload.get("validated_phase_coverage") != coverage:
        raise ValueError("schema-v2 validation coverage is not derived from records")
    if source_schema_version == 1:
        payload = {
            **payload,
            "source_identity": source_identity,
            "build_identities": build_identities,
        }
    payload["records"] = normalised_records
    payload["validated_phase_coverage"] = coverage
    payload["source_schema_version"] = source_schema_version
    payload["source_index_path"] = str(path.resolve())
    payload["source_index_sha256"] = hashlib.sha256(index_bytes).hexdigest()
    if (
        source_schema_version == 2
        and sha256(local_campaign_manifest) != campaign_manifest_hash
    ):
        raise ValueError("schema-v2 trusted campaign manifest changed during validation")
    return payload


def validation_coverage(
    validated_gxtb_meshes: set[str] | None,
    validation_index: dict[str, object] | None,
) -> dict[str, set[str]] | None:
    if validated_gxtb_meshes is None and validation_index is None:
        return None
    coverage = {
        mesh: set(PHASES) for mesh in (validated_gxtb_meshes or set())
    }
    if validation_index is not None:
        for record in validation_index.get("records", []):
            coverage.setdefault(str(record["mesh"]), set()).add(str(record["phase"]))
    return coverage


def build_convergence_report(
    results: dict[str, object],
    *,
    validation_index_path: Path | None = None,
    validation_index: dict[str, object] | None = None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Build a DMC-independent report from Ih-referenced g-xTB energies."""
    result_meshes = results.get("results", {})
    coverage = {
        mesh: set(phases)
        for mesh, phases in results.get("validated_gxtb_phases", {}).items()
    }
    comparisons: list[dict[str, object]] = []
    phase_rows: list[dict[str, object]] = []
    pair_by_key: dict[tuple[str, str], dict[str, object]] = {}
    nonreference = [phase for phase in PHASES if phase != "Ih"]

    for from_mesh, to_mesh in zip(CONVERGENCE_MESH_IDS, CONVERGENCE_MESH_IDS[1:]):
        from_method = result_meshes.get(from_mesh, {}).get("GXTB", {})
        to_method = result_meshes.get(to_mesh, {}).get("GXTB", {})
        from_relative = from_method.get("relative_kjmol", {})
        to_relative = to_method.get("relative_kjmol", {})
        common = [
            phase
            for phase in nonreference
            if phase in from_relative and phase in to_relative
        ]
        if not common:
            continue
        deltas = {
            phase: float(to_relative[phase]) - float(from_relative[phase])
            for phase in common
        }
        sum_squared = sum(value * value for value in deltas.values())
        observed_subset_rms = math.sqrt(sum_squared / len(deltas))
        full_set_rms_lower_bound = math.sqrt(sum_squared / len(nonreference))
        mad = sum(abs(value) for value in deltas.values()) / len(deltas)
        max_phase = max(common, key=lambda phase: abs(deltas[phase]))
        full_coverage = (
            coverage.get(from_mesh, set()) == set(PHASES)
            and coverage.get(to_mesh, set()) == set(PHASES)
            and len(common) == len(nonreference)
        )
        max_absolute_delta = abs(deltas[max_phase])
        passed = (
            max_absolute_delta <= CONVERGENCE_MAX_ABS_KJMOL
            and observed_subset_rms <= CONVERGENCE_RMS_KJMOL
            if full_coverage
            else None
        )
        pilot_definitively_rejects = (
            not full_coverage
            and (
                max_absolute_delta > CONVERGENCE_MAX_ABS_KJMOL
                or full_set_rms_lower_bound > CONVERGENCE_RMS_KJMOL
            )
        )
        comparison: dict[str, object] = {
            "from_mesh": from_mesh,
            "to_mesh": to_mesh,
            "coverage": "full" if full_coverage else "pilot",
            "eligible_for_stopping": full_coverage,
            "n_common_nonreference_phases": len(common),
            "expected_nonreference_phases": len(nonreference),
            "phases": common,
            "mean_absolute_delta_kjmol_per_h2o": mad,
            "rms_delta_kjmol_per_h2o": observed_subset_rms,
            "observed_subset_rms_delta_kjmol_per_h2o": observed_subset_rms,
            "full_set_rms_lower_bound_kjmol_per_h2o": full_set_rms_lower_bound,
            "max_absolute_delta_kjmol_per_h2o": max_absolute_delta,
            "max_absolute_delta_phase": max_phase,
            "passed_numeric_thresholds": passed,
            "passed_formal_stopping_pair": bool(passed) and full_coverage,
            "pilot_definitively_rejects_candidate": pilot_definitively_rejects,
            "pilot_inconclusive": (
                not full_coverage and not pilot_definitively_rejects
            ),
        }
        comparisons.append(comparison)
        pair_by_key[(from_mesh, to_mesh)] = comparison
        for phase in common:
            phase_rows.append(
                {
                    "from_mesh": from_mesh,
                    "to_mesh": to_mesh,
                    "coverage": comparison["coverage"],
                    "eligible_for_stopping": str(full_coverage).lower(),
                    "phase": phase,
                    "from_relative_kJmol_per_H2O": f"{float(from_relative[phase]):.12f}",
                    "to_relative_kJmol_per_H2O": f"{float(to_relative[phase]):.12f}",
                    "delta_kJmol_per_H2O": f"{deltas[phase]:.12f}",
                    "absolute_delta_kJmol_per_H2O": f"{abs(deltas[phase]):.12f}",
                }
            )

    complete_mesh_indices = [
        index
        for index, mesh in enumerate(CONVERGENCE_MESH_IDS)
        if coverage.get(mesh, set()) == set(PHASES)
        and set(
            result_meshes.get(mesh, {})
            .get("GXTB", {})
            .get("relative_kjmol", {})
        )
        >= set(PHASES)
    ]
    latest_complete_index = max(complete_mesh_indices, default=None)
    latest_pair: dict[str, object] | None = None
    trailing_pass_count = 0
    trailing_pass_start: int | None = None
    if latest_complete_index is not None and latest_complete_index > 0:
        latest_pair_index = latest_complete_index - 1
        latest_pair = pair_by_key.get(
            (
                CONVERGENCE_MESH_IDS[latest_pair_index],
                CONVERGENCE_MESH_IDS[latest_complete_index],
            )
        )
        pair_index = latest_pair_index
        while pair_index >= 0:
            pair = pair_by_key.get(
                (
                    CONVERGENCE_MESH_IDS[pair_index],
                    CONVERGENCE_MESH_IDS[pair_index + 1],
                )
            )
            if pair is None or not bool(pair["passed_formal_stopping_pair"]):
                break
            trailing_pass_count += 1
            trailing_pass_start = pair_index
            pair_index -= 1

    accepted: dict[str, object] | None = None
    if (
        trailing_pass_count >= CONVERGENCE_REQUIRED_CONSECUTIVE_PAIRS
        and trailing_pass_start is not None
    ):
        earlier = CONVERGENCE_MESH_IDS[trailing_pass_start]
        middle = CONVERGENCE_MESH_IDS[trailing_pass_start + 1]
        later = CONVERGENCE_MESH_IDS[trailing_pass_start + 2]
        accepted = {
            "earliest_reportable_mesh": middle,
            "validation_mesh": later,
            "supporting_pairs": [f"{earlier}->{middle}", f"{middle}->{later}"],
        }

    revoked_acceptance: dict[str, object] | None = None
    revoking_dense_pilots: list[dict[str, object]] = []
    if accepted is not None:
        validation_index = CONVERGENCE_MESH_IDS.index(
            str(accepted["validation_mesh"])
        )
        revoking_dense_pilots = [
            comparison
            for comparison in comparisons
            if comparison["coverage"] == "pilot"
            and bool(comparison["pilot_definitively_rejects_candidate"])
            and CONVERGENCE_MESH_IDS.index(str(comparison["to_mesh"]))
            > validation_index
        ]
        if revoking_dense_pilots:
            revoked_acceptance = accepted
            accepted = None

    if revoked_acceptance is not None:
        status = "not_converged"
        reason = (
            "a later dense pilot definitively rejects the previously accepted "
            "mesh sequence"
        )
    elif accepted is not None:
        status = "converged"
        reason = "two consecutive fully covered refinements pass both thresholds"
    elif latest_pair is not None and latest_pair["eligible_for_stopping"] and not bool(
        latest_pair["passed_formal_stopping_pair"]
    ):
        status = "not_converged"
        reason = "latest fully covered refinement fails at least one threshold"
    else:
        status = "insufficient_coverage"
        reason = "two consecutive passing fully covered refinements are unavailable"

    report: dict[str, object] = {
        "schema_version": 1,
        "benchmark": "DMC-ICE13",
        "method": "g-xTB",
        "quantity": "Ih-referenced relative energy",
        "unit": "kJ mol^-1 per H2O",
        "dmc_reference_used_for_convergence": False,
        "formula": (
            "delta(alpha;m,n)=2625.499638*[(E(alpha,m)/N(alpha)-"
            "E(Ih,m)/N(Ih))-(E(alpha,n)/N(alpha)-E(Ih,n)/N(Ih))]"
        ),
        "thresholds": {
            "max_absolute_delta_kjmol_per_h2o": CONVERGENCE_MAX_ABS_KJMOL,
            "rms_delta_kjmol_per_h2o": CONVERGENCE_RMS_KJMOL,
            "required_consecutive_fully_covered_pairs": (
                CONVERGENCE_REQUIRED_CONSECUTIVE_PAIRS
            ),
        },
        "core_meshes": [str(mesh["id"]) for mesh in MESHES],
        "supported_dense_extensions": [
            str(mesh["id"]) for mesh in DENSE_EXTENSION_MESHES
        ],
        "validated_phase_coverage": {
            mesh: [phase for phase in PHASES if phase in phases]
            for mesh, phases in coverage.items()
        },
        "comparisons": comparisons,
        "stopping_assessment": {
            "status": status,
            "reason": reason,
            "acceptance": accepted,
            "revoked_acceptance": revoked_acceptance,
            "revoking_dense_pilot_pairs": [
                f"{comparison['from_mesh']}->{comparison['to_mesh']}"
                for comparison in revoking_dense_pilots
            ],
            "latest_fully_validated_mesh": (
                CONVERGENCE_MESH_IDS[latest_complete_index]
                if latest_complete_index is not None
                else None
            ),
            "trailing_consecutive_passing_fully_covered_pairs": trailing_pass_count,
        },
    }
    if validation_index is not None and validation_index_path is not None:
        source_index_sha256 = validation_index.get("source_index_sha256")
        if (
            not isinstance(source_index_sha256, str)
            or sha256(validation_index_path) != source_index_sha256
        ):
            raise ValueError("validation index changed after verified read")
        report["validation_index"] = {
            "path": str(validation_index_path.resolve()),
            "sha256": source_index_sha256,
            "schema_version": validation_index.get("source_schema_version"),
            "campaign_id": validation_index.get("campaign_id"),
            "build_identities": validation_index.get("build_identities"),
        }
    return report, phase_rows


def mesh_parity_direction(from_mesh: str, to_mesh: str) -> str:
    """Return the odd/even direction for two consecutive cubic meshes."""
    from_size = int(str(MESH_BY_ID[from_mesh]["scheme"]).split()[1])
    to_size = int(str(MESH_BY_ID[to_mesh]["scheme"]).split()[1])
    from_parity = "odd" if from_size % 2 else "even"
    to_parity = "odd" if to_size % 2 else "even"
    if from_parity == to_parity:
        raise ValueError(
            f"adjacent mesh pair does not change parity: {from_mesh}->{to_mesh}"
        )
    return f"{from_parity}->{to_parity}"


def same_mesh_relative_energy(
    result_meshes: dict[str, object],
    mesh_id: str,
    method: str,
    phase: str,
) -> float | None:
    """Recompute a phase relative energy using Ih from the identical mesh."""
    method_result = result_meshes.get(mesh_id, {}).get(method, {})
    per_h2o = method_result.get("per_h2o_hartree", {})
    if not isinstance(per_h2o, dict):
        return None
    phase_energy = per_h2o.get(phase)
    ih_energy = per_h2o.get("Ih")
    if phase_energy is None or ih_energy is None:
        return None
    value = (float(phase_energy) - float(ih_energy)) * HARTREE_TO_KJMOL
    if not math.isfinite(value):
        return None
    stored_relative = method_result.get("relative_kjmol", {})
    if isinstance(stored_relative, dict) and phase in stored_relative:
        stored = float(stored_relative[phase])
        if not math.isclose(value, stored, rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(
                f"{method} {mesh_id}/{phase} relative energy is not bound to "
                "same-mesh Ih"
            )
    return value


def _phasewise_delta(
    from_mesh: str,
    to_mesh: str,
    from_value: float,
    to_value: float,
) -> dict[str, object]:
    value = to_value - from_value
    return {
        "from_mesh": from_mesh,
        "to_mesh": to_mesh,
        "value_kjmol_per_h2o": value,
        "absolute_value_kjmol_per_h2o": abs(value),
        "parity_direction": mesh_parity_direction(from_mesh, to_mesh),
    }


def _phasewise_delta_statistics(
    values: list[float],
    expected: int,
) -> dict[str, object]:
    """Summarize the chosen phase deltas without defining another gate."""
    return {
        "expected_delta_count": expected,
        "delta_count": len(values),
        "complete": len(values) == expected,
        "values_kjmol_per_h2o": values,
        "rms_delta_kjmol_per_h2o": (
            math.sqrt(sum(value * value for value in values) / len(values))
            if values
            else None
        ),
        "mean_absolute_delta_kjmol_per_h2o": (
            sum(abs(value) for value in values) / len(values)
            if values
            else None
        ),
        "max_absolute_delta_kjmol_per_h2o": (
            max((abs(value) for value in values), default=None)
        ),
        "diagnostic_only": True,
    }


def next_phasewise_mesh(values: dict[str, float]) -> str | None:
    """Return the earliest missing mesh that can create the next adjacent pair."""
    for left, right in zip(CONVERGENCE_MESH_IDS, CONVERGENCE_MESH_IDS[1:]):
        if left in values and right not in values:
            return right
    for left, right in reversed(
        list(zip(CONVERGENCE_MESH_IDS, CONVERGENCE_MESH_IDS[1:]))
    ):
        if right in values and left not in values:
            return left
    if not values:
        return CONVERGENCE_MESH_IDS[0]
    return None


def build_phasewise_kpoint_convergence_report(
    results: dict[str, object],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Build the phase-wise k-point-converged DMC13 result."""
    result_meshes = results.get("results", {})
    nonreference = [phase for phase in PHASES if phase != "Ih"]
    dmc_relative = {
        phase: DMC_ABS_KJMOL[phase] - DMC_ABS_KJMOL["Ih"]
        for phase in nonreference
    }
    method_reports: dict[str, object] = {}
    phase_rows: list[dict[str, object]] = []

    for method in METHODS:
        phase_convergence: dict[str, dict[str, object]] = {}
        diagnostics: dict[str, dict[str, object]] = {}
        unresolved: dict[str, dict[str, object]] = {}
        next_required_by_mesh: dict[str, list[str]] = {}
        exhausted_phases: list[str] = []

        for phase in nonreference:
            values = {
                mesh: value
                for mesh in CONVERGENCE_MESH_IDS
                if (
                    value := same_mesh_relative_energy(
                        result_meshes, mesh, method, phase
                    )
                )
                is not None
            }
            available_meshes = [
                mesh for mesh in CONVERGENCE_MESH_IDS if mesh in values
            ]
            stable_candidates: list[dict[str, object]] = []
            failed_pairs: list[dict[str, object]] = []
            revoked_candidates: list[dict[str, object]] = []

            for pair_index in range(len(CONVERGENCE_MESH_IDS) - 1):
                previous_mesh = CONVERGENCE_MESH_IDS[pair_index]
                converged_mesh = CONVERGENCE_MESH_IDS[pair_index + 1]
                if previous_mesh not in values or converged_mesh not in values:
                    continue
                delta = _phasewise_delta(
                    previous_mesh,
                    converged_mesh,
                    values[previous_mesh],
                    values[converged_mesh],
                )
                candidate: dict[str, object] = {
                    "smallest_required_mesh": converged_mesh,
                    "mesh_label": MESH_BY_ID[converged_mesh]["label"],
                    "mesh_n": round(
                        float(MESH_BY_ID[converged_mesh]["nk_total"])
                        ** (1.0 / 3.0)
                    ),
                    "nk_total": MESH_BY_ID[converged_mesh]["nk_total"],
                    "previous_mesh": previous_mesh,
                    "relative_energy_kjmol_per_h2o": values[converged_mesh],
                    "dmc_relative_kjmol_per_h2o": dmc_relative[phase],
                    "error_kjmol_per_h2o": (
                        values[converged_mesh] - dmc_relative[phase]
                    ),
                    "previous_relative_energy_kjmol_per_h2o": (
                        values[previous_mesh]
                    ),
                    "last_delta_kjmol_per_h2o": delta[
                        "value_kjmol_per_h2o"
                    ],
                    "absolute_last_delta_kjmol_per_h2o": delta[
                        "absolute_value_kjmol_per_h2o"
                    ],
                    "parity_direction": delta["parity_direction"],
                    "passes_phase_limit": (
                        float(delta["absolute_value_kjmol_per_h2o"])
                        <= PHASEWISE_KPOINT_MAX_ABS_KJMOL
                    ),
                    "reference_phase": "Ih",
                    "relative_energy_construction": (
                        "E(phase, mesh)/N(H2O, phase) - "
                        "E(Ih, same mesh)/N(H2O, Ih)"
                    ),
                    "chronological_pair_index": pair_index,
                }
                if not candidate["passes_phase_limit"]:
                    failed_pairs.append(candidate)
                    continue

                later_contradictions: list[dict[str, object]] = []
                for later_index in range(
                    pair_index + 1, len(CONVERGENCE_MESH_IDS) - 1
                ):
                    later_previous = CONVERGENCE_MESH_IDS[later_index]
                    later_mesh = CONVERGENCE_MESH_IDS[later_index + 1]
                    if (
                        later_previous not in values
                        or later_mesh not in values
                    ):
                        continue
                    later_delta = _phasewise_delta(
                        later_previous,
                        later_mesh,
                        values[later_previous],
                        values[later_mesh],
                    )
                    if (
                        float(later_delta["absolute_value_kjmol_per_h2o"])
                        > PHASEWISE_KPOINT_MAX_ABS_KJMOL
                    ):
                        later_contradictions.append(later_delta)
                candidate["later_available_contradictions"] = (
                    later_contradictions
                )
                candidate["revoked_by_later_evidence"] = bool(
                    later_contradictions
                )
                if later_contradictions:
                    revoked_candidates.append(candidate)
                    continue
                stable_candidates.append(candidate)

            diagnostics[phase] = {
                "available_meshes": available_meshes,
                "missing_meshes": [
                    mesh for mesh in CONVERGENCE_MESH_IDS if mesh not in values
                ],
                "failed_available_pairs": failed_pairs,
                "revoked_candidates": revoked_candidates,
                "stable_candidates": stable_candidates,
            }

            if stable_candidates:
                chosen = stable_candidates[0]
                phase_convergence[phase] = chosen
                phase_rows.append(
                    {
                        "method": method,
                        "method_label": method_label(method),
                        "phase": phase,
                        "status": "k-point-converged",
                        "smallest_required_mesh": chosen[
                            "smallest_required_mesh"
                        ],
                        "mesh_label": chosen["mesh_label"],
                        "mesh_n": chosen["mesh_n"],
                        "nk_total": chosen["nk_total"],
                        "previous_mesh": chosen["previous_mesh"],
                        "relative_energy_kJmol_per_H2O": chosen[
                            "relative_energy_kjmol_per_h2o"
                        ],
                        "DMC_relative_kJmol_per_H2O": chosen[
                            "dmc_relative_kjmol_per_h2o"
                        ],
                        "error_kJmol_per_H2O": chosen[
                            "error_kjmol_per_h2o"
                        ],
                        "previous_relative_energy_kJmol_per_H2O": chosen[
                            "previous_relative_energy_kjmol_per_h2o"
                        ],
                        "last_delta_kJmol_per_H2O": chosen[
                            "last_delta_kjmol_per_h2o"
                        ],
                        "absolute_last_delta_kJmol_per_H2O": chosen[
                            "absolute_last_delta_kjmol_per_h2o"
                        ],
                        "parity_direction": chosen["parity_direction"],
                        "next_required_mesh": None,
                    }
                )
                continue

            next_mesh = next_phasewise_mesh(values)
            unresolved[phase] = {
                **diagnostics[phase],
                "next_required_mesh": next_mesh,
            }
            if next_mesh is None:
                exhausted_phases.append(phase)
            else:
                next_required_by_mesh.setdefault(next_mesh, []).append(phase)
            phase_rows.append(
                {
                    "method": method,
                    "method_label": method_label(method),
                    "phase": phase,
                    "status": "unresolved",
                    "smallest_required_mesh": None,
                    "mesh_label": None,
                    "mesh_n": None,
                    "nk_total": None,
                    "previous_mesh": None,
                    "relative_energy_kJmol_per_H2O": None,
                    "DMC_relative_kJmol_per_H2O": None,
                    "error_kJmol_per_H2O": None,
                    "previous_relative_energy_kJmol_per_H2O": None,
                    "last_delta_kJmol_per_H2O": None,
                    "absolute_last_delta_kJmol_per_H2O": None,
                    "parity_direction": None,
                    "next_required_mesh": next_mesh,
                }
            )

        complete = len(phase_convergence) == len(nonreference)
        last_deltas = [
            float(phase_convergence[phase]["last_delta_kjmol_per_h2o"])
            for phase in nonreference
            if phase in phase_convergence
        ]
        delta_statistics = _phasewise_delta_statistics(
            last_deltas, len(nonreference)
        )

        converged_stats: dict[str, float] | None = None
        previous_stats: dict[str, float] | None = None
        mae_bound: dict[str, object] | None = None
        if complete:
            converged_errors = [
                float(
                    phase_convergence[phase][
                        "relative_energy_kjmol_per_h2o"
                    ]
                )
                - dmc_relative[phase]
                for phase in nonreference
            ]
            previous_errors = [
                float(
                    phase_convergence[phase][
                        "previous_relative_energy_kjmol_per_h2o"
                    ]
                )
                - dmc_relative[phase]
                for phase in nonreference
            ]
            converged_stats = stats(converged_errors)
            previous_stats = stats(previous_errors)
            mean_absolute_delta = float(
                delta_statistics[
                    "mean_absolute_delta_kjmol_per_h2o"
                ]
            )
            observed_mae_difference = abs(
                converged_stats["MAE"] - previous_stats["MAE"]
            )
            mae_bound = {
                "formula": (
                    "|MAE(converged)-MAE(previous mesh)| <= "
                    "mean_phase(|last delta|)"
                ),
                "observed_absolute_mae_difference_kjmol_per_h2o": (
                    observed_mae_difference
                ),
                "provable_absolute_mae_difference_upper_bound_kjmol_per_h2o": (
                    mean_absolute_delta
                ),
                "bound_satisfied": (
                    observed_mae_difference <= mean_absolute_delta + 1.0e-12
                ),
            }

        next_required_by_mesh = {
            mesh: next_required_by_mesh[mesh]
            for mesh in CONVERGENCE_MESH_IDS
            if mesh in next_required_by_mesh
        }
        unresolved_phases = [
            phase for phase in nonreference if phase in unresolved
        ]
        method_reports[method] = {
            "method_label": method_label(method),
            "result_label": "phase-wise k-point-converged MAE",
            "result_label_de": "phasenweise k-Punkt-konvergierter MAE",
            "status": (
                "phasewise_kpoint_converged"
                if complete
                else "unresolved_phases"
            ),
            "expected_phase_count": len(nonreference),
            "converged_phase_count": len(phase_convergence),
            "phasewise_kpoint_converged": complete,
            "phase_convergence": phase_convergence,
            "phase_diagnostics": diagnostics,
            "unresolved": unresolved,
            "unresolved_phases": unresolved_phases,
            "next_required_phases": unresolved_phases,
            "next_required_by_mesh": next_required_by_mesh,
            "required_same_mesh_reference_by_mesh": {
                mesh: "Ih" for mesh in next_required_by_mesh
            },
            "exhausted_mesh_sequence_phases": exhausted_phases,
            "last_delta_statistics_diagnostic": delta_statistics,
            "phasewise_kpoint_converged_stats_nonreference": (
                converged_stats
            ),
            "previous_mesh_comparison_stats_nonreference": previous_stats,
            "mae_mesh_difference_bound": mae_bound,
        }

    report: dict[str, object] = {
        "schema_version": 1,
        "benchmark": "DMC-ICE13",
        "protocol_id": "phasewise-kpoint-convergence-v1",
        "result_label": "phase-wise k-point-converged MAE",
        "result_label_de": "phasenweise k-Punkt-konvergierter MAE",
        "quantity": "same-mesh-Ih-referenced relative energy",
        "unit": "kJ mol^-1 per H2O",
        "mesh_sequence": CONVERGENCE_MESH_IDS,
        "reference_phase": "Ih",
        "nonreference_phases": nonreference,
        "dmc_reference_variant": "legacy_rounded_absolute_XI_0.16",
        "phase_convergence_rule": (
            "a phase is k-point converged at N^3 when "
            "|Erel(N^3)-Erel((N-1)^3)| <= 0.05 kJ mol^-1 per H2O"
        ),
        "dataset_completion_rule": (
            "the phase-wise result is complete when all 12 non-reference "
            "phases are k-point converged"
        ),
        "same_mesh_ih_required": True,
        "later_available_evidence_safety_check": True,
        "last_delta_statistics_are_diagnostic_only": True,
        "thresholds": {
            "per_phase_max_absolute_delta_kjmol_per_h2o": (
                PHASEWISE_KPOINT_MAX_ABS_KJMOL
            ),
        },
        "fixed_mesh_statistics_are_separate": True,
        "methods": method_reports,
    }
    return report, phase_rows


def write_phasewise_kpoint_convergence_artifacts(
    results: dict[str, object],
    output_prefix: str,
) -> dict[str, object]:
    report, phase_rows = build_phasewise_kpoint_convergence_report(results)
    json_path = data_output_path(
        "phasewise_kpoint_convergence.json", output_prefix
    )
    csv_path = data_output_path(
        "phasewise_kpoint_convergence.csv", output_prefix
    )
    atomic_write_bytes(
        json_path,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode(),
    )
    write_csv(
        csv_path,
        phase_rows,
        [
            "method",
            "method_label",
            "phase",
            "status",
            "smallest_required_mesh",
            "mesh_label",
            "mesh_n",
            "nk_total",
            "previous_mesh",
            "relative_energy_kJmol_per_H2O",
            "DMC_relative_kJmol_per_H2O",
            "error_kJmol_per_H2O",
            "previous_relative_energy_kJmol_per_H2O",
            "last_delta_kJmol_per_H2O",
            "absolute_last_delta_kJmol_per_H2O",
            "parity_direction",
            "next_required_mesh",
        ],
    )
    return report



def write_convergence_artifacts(
    results: dict[str, object],
    output_prefix: str,
    *,
    validation_index_path: Path | None = None,
    validation_index: dict[str, object] | None = None,
) -> dict[str, object]:
    report, phase_rows = build_convergence_report(
        results,
        validation_index_path=validation_index_path,
        validation_index=validation_index,
    )
    json_path = data_output_path("kpoint_convergence.json", output_prefix)
    csv_path = data_output_path("kpoint_convergence.csv", output_prefix)
    atomic_write_bytes(
        json_path,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode(),
    )
    write_csv(
        csv_path,
        phase_rows,
        [
            "from_mesh",
            "to_mesh",
            "coverage",
            "eligible_for_stopping",
            "phase",
            "from_relative_kJmol_per_H2O",
            "to_relative_kJmol_per_H2O",
            "delta_kJmol_per_H2O",
            "absolute_delta_kJmol_per_H2O",
        ],
    )
    return report


def validated_gxtb_output_paths(
    validation_index: dict[str, object],
    gxtb_run_root: Path,
) -> dict[tuple[str, str], Path]:
    schema_version = int(validation_index.get("source_schema_version", 1))
    paths: dict[tuple[str, str], Path] = {}
    for record in validation_index.get("records", []):
        mesh = str(record["mesh"])
        phase = str(record["phase"])
        indexed = _validation_path(record["output"], schema_version).resolve()
        expected = output_path(mesh, "GXTB", phase, gxtb_run_root).resolve()
        if indexed != expected:
            raise ValueError(
                f"validation output path mismatch for {mesh}/{phase}: "
                f"indexed {indexed}, expected {expected}"
            )
        paths[(mesh, phase)] = indexed
    return paths


def analyse(
    validated_gxtb_meshes: set[str] | None = None,
    *,
    gxtb_run_root: Path | None = None,
    output_prefix: str | None = None,
    validation_index_path: Path | None = None,
    validation_index_sha256: str | None = None,
) -> dict[str, object]:
    geometries = json.loads((DATA / "geometries.json").read_text())
    dmc_rel = {phase: DMC_ABS_KJMOL[phase] - DMC_ABS_KJMOL["Ih"] for phase in PHASES}
    validation_index = (
        read_validation_index(validation_index_path)
        if validation_index_path is not None
        else None
    )
    if validation_index is not None and validation_index_sha256 is not None:
        if validation_index.get("source_index_sha256") != validation_index_sha256:
            raise ValueError("validation index SHA256 pin mismatch")
    analysis_methods = METHODS if output_prefix else BASELINE_METHODS
    verified_gxtb_outputs: dict[tuple[str, str], Path] = {}
    verified_gxtb_energies: dict[tuple[str, str], float] = {}
    if validation_index is not None:
        if gxtb_run_root is None:
            raise ValueError("validation-index analysis requires --gxtb-run-root")
        verified_gxtb_outputs = validated_gxtb_output_paths(
            validation_index, gxtb_run_root
        )
        verified_gxtb_energies = {
            (str(record["mesh"]), str(record["phase"])): float(
                record["validated_energy_hartree"]
            )
            for record in validation_index.get("records", [])
        }
    coverage = validation_coverage(validated_gxtb_meshes, validation_index)
    analysis_meshes = [
        *MESHES,
        *[
            mesh
            for mesh in DENSE_EXTENSION_MESHES
            if output_prefix
            and coverage is not None
            and str(mesh["id"]) in coverage
        ],
    ]
    fully_validated_meshes = sorted(
        mesh for mesh, phases in (coverage or {}).items() if phases == set(PHASES)
    )
    prior_results: dict[str, object] = {}
    prior_path = DATA / "kpoint_results.json"
    if prior_path.is_file():
        try:
            prior_results = json.loads(prior_path.read_text()).get("results", {})
        except (json.JSONDecodeError, AttributeError):
            prior_results = {}
    results: dict[str, object] = {
        "meshes": MESHES,
        "supported_dense_extension_meshes": DENSE_EXTENSION_MESHES,
        "analysis_meshes": analysis_meshes,
        "methods": analysis_methods,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "gxtb_run_root": str(gxtb_run_root.resolve()) if gxtb_run_root else None,
        "validated_gxtb_meshes": fully_validated_meshes,
        "validated_gxtb_phases": {
            mesh: [phase for phase in PHASES if phase in phases]
            for mesh, phases in (coverage or {}).items()
        },
        "results": {},
    }
    relative_rows: list[dict[str, object]] = []
    stats_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []

    for mesh in analysis_meshes:
        mesh_id = str(mesh["id"])
        mesh_results: dict[str, object] = {}
        for method in analysis_methods:
            energies = {
                phase: (
                    verified_gxtb_energies.get((mesh_id, phase))
                    if method == "GXTB" and validation_index is not None
                    else parse_energy(
                        output_path(mesh_id, method, phase, gxtb_run_root)
                    )
                    if (
                        method != "GXTB"
                        or coverage is None
                        or phase in coverage.get(mesh_id, set())
                    )
                    else None
                )
                for phase in PHASES
            }
            complete = all(value is not None for value in energies.values())
            method_result: dict[str, object] = {"complete": complete, "energies_hartree": energies}
            available = {
                phase: float(value)
                for phase, value in energies.items()
                if value is not None
            }
            if "Ih" in available:
                per_h2o = {
                    phase: energy / geometries[phase]["counts"]["O"]
                    for phase, energy in available.items()
                }
                ih = per_h2o["Ih"]
                rel = {
                    phase: (energy - ih) * HARTREE_TO_KJMOL
                    for phase, energy in per_h2o.items()
                }
                method_result.update(
                    {
                        "per_h2o_hartree": per_h2o,
                        "relative_kjmol": rel,
                    }
                )
            if not complete:
                prior_method = (
                    prior_results.get(mesh_id, {}).get(method, {})
                    if isinstance(prior_results, dict)
                    else {}
                )
                if (
                    method != "GXTB"
                    and isinstance(prior_method, dict)
                    and prior_method.get("complete")
                ):
                    # Raw legacy runs are intentionally ignored by Git.  Keep a
                    # previously curated complete GFN1/GFN2 result when an
                    # additive GXTB-only analysis is performed in a fresh clone.
                    method_result = dict(prior_method)

            if method_result.get("complete"):
                rel = {
                    phase: float(method_result["relative_kjmol"][phase])
                    for phase in PHASES
                }
                err = {phase: rel[phase] - dmc_rel[phase] for phase in PHASES}
                err_nonref = [err[phase] for phase in PHASES if phase != "Ih"]
                s = stats(err_nonref)
                method_result.update(
                    {
                        "relative_errors_kjmol": err,
                        "stats_nonreference": s,
                    }
                )
                stats_rows.append(
                    {
                        "mesh": mesh_id,
                        "mesh_label": mesh["label"],
                        "nk_total": mesh["nk_total"],
                        "method": method_label(method),
                        "N": len(err_nonref),
                        **{key: f"{value:.6f}" for key, value in s.items()},
                    }
                )
                for phase in PHASES:
                    relative_rows.append(
                        {
                            "mesh": mesh_id,
                            "mesh_label": mesh["label"],
                            "method": method_label(method),
                            "phase": phase,
                            "DMC_relative_kJmol": f"{dmc_rel[phase]:.6f}",
                            "relative_kJmol": f"{rel[phase]:.6f}",
                            "error_kJmol": f"{err[phase]:.6f}",
                        }
                    )
                primary_errors = {
                    phase: rel[phase] - DMC_REL_PRIMARY_KJMOL[phase] for phase in PHASES
                }
                primary_stats = stats(
                    [primary_errors[phase] for phase in PHASES if phase != "Ih"]
                )
                for reference_variant, values in (
                    ("legacy_rounded_absolute_XI_0.16", s),
                    ("primary_explicit_relative_XI_0.15", primary_stats),
                ):
                    sensitivity_rows.append(
                        {
                            "mesh": mesh_id,
                            "mesh_label": mesh["label"],
                            "method": method_label(method),
                            "reference_variant": reference_variant,
                            "N": len(err_nonref),
                            **{key: f"{value:.6f}" for key, value in values.items()},
                        }
                    )
            mesh_results[method] = method_result
        results["results"][mesh_id] = mesh_results

    DATA.mkdir(exist_ok=True)
    data_output_path("kpoint_results.json", output_prefix).write_text(
        json.dumps(results, indent=2) + "\n"
    )
    write_csv(
        data_output_path("kpoint_relative_energies.csv", output_prefix),
        relative_rows,
        [
            "mesh",
            "mesh_label",
            "method",
            "phase",
            "DMC_relative_kJmol",
            "relative_kJmol",
            "error_kJmol",
        ],
    )
    write_csv(
        data_output_path("kpoint_stats.csv", output_prefix),
        stats_rows,
        [
            "mesh",
            "mesh_label",
            "nk_total",
            "method",
            "N",
            "ME",
            "MAE",
            "RMSE",
            "MaxAE",
        ],
    )
    write_csv(
        data_output_path("reference_sensitivity.csv", output_prefix),
        sensitivity_rows,
        [
            "mesh",
            "mesh_label",
            "method",
            "reference_variant",
            "N",
            "ME",
            "MAE",
            "RMSE",
            "MaxAE",
        ],
    )
    if output_prefix is not None and coverage is not None:
        write_convergence_artifacts(
            results,
            output_prefix,
            validation_index_path=validation_index_path,
            validation_index=validation_index,
        )
        write_phasewise_kpoint_convergence_artifacts(results, output_prefix)
    make_plots(stats_rows, output_prefix=output_prefix)
    return results


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_plots(
    stats_rows: list[dict[str, object]],
    *,
    output_prefix: str | None = None,
) -> None:
    if not stats_rows:
        return
    if shutil.which("gnuplot") is None or shutil.which("rsvg-convert") is None:
        return
    FIGURES.mkdir(exist_ok=True)
    dat = data_output_path("kpoint_stats_for_plot.dat", output_prefix)
    stats_by_mesh_method = {
        (str(row["mesh"]), str(row["method"])): row for row in stats_rows
    }
    complete_methods = [
        method
        for method in METHODS
        if all(
            (str(mesh["id"]), method_label(method)) in stats_by_mesh_method
            for mesh in MESHES
        )
    ]
    if not complete_methods:
        return
    plot_meshes = list(MESHES)
    if output_prefix is not None:
        completed_dense_indices = [
            index
            for index, mesh in enumerate(DENSE_EXTENSION_MESHES)
            if (
                str(mesh["id"]),
                method_label("GXTB"),
            )
            in stats_by_mesh_method
        ]
        if completed_dense_indices:
            plot_meshes.extend(
                DENSE_EXTENSION_MESHES[: max(completed_dense_indices) + 1]
            )
    with dat.open("w") as handle:
        columns = " ".join(f"{method}_MAE" for method in complete_methods)
        handle.write(f"# index mesh_label {columns}\n")
        for index, mesh in enumerate(plot_meshes, start=1):
            mesh_id = str(mesh["id"])
            values = [
                stats_by_mesh_method.get(
                    (mesh_id, method_label(method)), {"MAE": "NaN"}
                )["MAE"]
                for method in complete_methods
            ]
            handle.write(
                f"{index} \"{mesh['label']}\" "
                f"{' '.join(str(value) for value in values)}\n"
            )

    styles = []
    plots = []
    for index, method in enumerate(complete_methods, start=1):
        styles.append(
            f"set style line {index} lc rgb '{METHOD_COLORS[method]}' "
            f"lw 2.2 pt {5 + 2 * index} ps 0.9"
        )
        prefix = f"plot '{dat}'" if index == 1 else "     ''"
        plots.append(
            f"{prefix} using 1:{index + 2} with linespoints ls {index} "
            f"title '{method_label(method)} MAE'"
        )
    plot_command = ", \\\n".join(plots)

    figure_stem = (
        f"dmc_ice13_{output_prefix}_kpoint_mae"
        if output_prefix
        else "dmc_ice13_kpoint_mae"
    )
    svg = FIGURES / f"{figure_stem}.svg"
    if output_prefix is None:
        xrange = "[0.75:6.25]"
        xtics = (
            "('Gamma' 1, '1x1x1' 2, '2x2x2' 3, '3x3x3' 4, "
            "'4x4x4' 5, '5x5x5' 6)"
        )
    else:
        xrange = f"[0.75:{len(plot_meshes) + 0.25:.2f}]"
        xtics = "(" + ", ".join(
            f"'{mesh['label']}' {index}"
            for index, mesh in enumerate(plot_meshes, start=1)
        ) + ")"
    script = f"""
set terminal svg enhanced font 'Helvetica,12' size 840,500
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set output '{svg}'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set key top right spacing 1.2 samplen 2
set xlabel 'k-point mesh'
set ylabel 'Relative-energy error / kJ mol^{{-1}} per H_2O'
set xrange {xrange}
set yrange [0:*]
set xtics {xtics}
{chr(10).join(styles)}
{plot_command}
"""
    subprocess.run(["gnuplot"], input=script.encode(), check=True)
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n"
    )
    subprocess.run(
        ["rsvg-convert", str(svg), "-o", str(svg.with_suffix(".png"))],
        check=True,
    )
    subprocess.run(
        [
            "rsvg-convert",
            "-f",
            "pdf",
            str(svg),
            "-o",
            str(svg.with_suffix(".pdf")),
        ],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=["prepare", "analyse", "all"],
        nargs="?",
        default="all",
    )
    parser.add_argument(
        "--method",
        action="append",
        choices=METHODS,
        help="method(s) to prepare; analysis always merges every complete method",
    )
    parser.add_argument(
        "--mesh",
        action="append",
        choices=[str(mesh["id"]) for mesh in SUPPORTED_MESHES],
        help=(
            "mesh(es) to prepare; repeat as needed (default: frozen core meshes only)"
        ),
    )
    parser.add_argument(
        "--phase",
        action="append",
        choices=PHASES,
        help="phase(s) to prepare; repeat as needed (default: all phases)",
    )
    parser.add_argument(
        "--validated-gxtb-mesh",
        action="append",
        choices=[str(mesh["id"]) for mesh in SUPPORTED_MESHES],
        help=(
            "accept g-XTB outputs only for these externally hash-validated meshes; "
            "repeat as needed"
        ),
    )
    parser.add_argument(
        "--restrict-gxtb-to-validated",
        action="store_true",
        help="reject every g-XTB mesh not named by --validated-gxtb-mesh",
    )
    parser.add_argument(
        "--gxtb-input-root",
        type=Path,
        help="separate directory for generated g-xTB production inputs",
    )
    parser.add_argument(
        "--gxtb-run-root",
        type=Path,
        help="separate directory containing hash-validated g-xTB production runs",
    )
    parser.add_argument(
        "--output-prefix",
        help="write additive result/plot files instead of replacing the GFN1/GFN2 baseline",
    )
    parser.add_argument(
        "--validation-index",
        type=Path,
        help="runner-generated hash-valid phase index, including dense pilot coverage",
    )
    parser.add_argument(
        "--validation-index-sha256",
        help="required exact SHA256 pin for --validation-index",
    )
    args = parser.parse_args()
    if (args.validation_index is None) != (args.validation_index_sha256 is None):
        parser.error(
            "--validation-index and --validation-index-sha256 must be provided together"
        )
    if args.validation_index_sha256 is not None and not re.fullmatch(
        r"[0-9a-f]{64}", args.validation_index_sha256
    ):
        parser.error("--validation-index-sha256 must be 64 lowercase hex characters")
    if args.command in {"prepare", "all"}:
        prepare_inputs(args.method, args.gxtb_input_root, args.mesh, args.phase)
    if args.command in {"analyse", "all"}:
        validated_gxtb_meshes = (
            set(args.validated_gxtb_mesh or [])
            if args.restrict_gxtb_to_validated
            else None
        )
        analyse(
            validated_gxtb_meshes,
            gxtb_run_root=args.gxtb_run_root,
            output_prefix=args.output_prefix,
            validation_index_path=args.validation_index,
            validation_index_sha256=args.validation_index_sha256,
        )


if __name__ == "__main__":
    main()
