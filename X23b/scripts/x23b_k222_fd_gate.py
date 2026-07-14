#!/usr/bin/env python3
"""Finite-difference validation of X23b GXTB shifted-k222 derivatives.

The gate is deliberately additive and small: four frozen X23 reference
crystals, two normalized collective Cartesian directions, and two symmetric
strain directions.  Measurements are immutable and remain scientifically
unapproved until a separate approval artifact is written explicitly.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Mapping, Sequence

import x23b_common as common
import x23b_experimental_k222_preflight as preflight
import x23b_k222_force_stress_gate as derivative_parser
import x23b_kpoint_cellopt as cellopt


ROOT = Path(__file__).resolve().parents[1]
PHASE = "x23b_k222_finite_difference_gate_v1"
VARIANT = "frozen_reference_shifted_k222_spglib_fd"
SOURCE_POLICY = "experimental_reference"
MANIFEST_SCHEMA = 1
REPORT_SCHEMA = 1
APPROVAL_SCHEMA = 1
MANIFEST_NAME = "x23b_k222_fd_gate_manifest.json"
OUTPUT_NAME = "cp2k.out"
DEFAULT_SYSTEMS = (
    "ammonia",
    "14-cyclohexanedione",
    "acetic_acid",
    # Triclinic reference cell; this is the deliberately low-symmetry case.
    "ethylcarbamate",
)
DEFAULT_COORDINATE_STEP_BOHR = 1.0e-3
DEFAULT_STRAIN_STEP = 5.0e-4
DEFAULT_DIRECTION_COUNT = 2
BOHR_TO_ANGSTROM = 0.529177210903
# Exact SI definition: (1 GPa)*(1 Angstrom^3) expressed in Hartree.
GPA_ANGSTROM3_TO_HARTREE = 1.0e-21 / 4.3597447222071e-18
FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
STRAIN_DIRECTIONS = (
    {
        "id": "isotropic_linear",
        "description": "equal infinitesimal linear strain on x, y, and z",
        "generator": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    },
    {
        "id": "symmetric_xy_engineering_shear",
        "description": "symmetric xy shear; scalar parameter is engineering shear",
        "generator": [[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]],
    },
)
CSV_FIELDS = (
    "method",
    "system",
    "phase",
    "variant",
    "source_policy",
    "measurement_type",
    "direction_id",
    "scientific_status",
    "approved",
    "campaign_fingerprint_sha256",
    "manifest_sha256",
    "source_input_sha256",
    "structure_sha256",
    "baseline_input_sha256",
    "baseline_output_sha256",
    "minus_input_sha256",
    "minus_output_sha256",
    "plus_input_sha256",
    "plus_output_sha256",
    "step",
    "step_unit",
    "reference_volume_A3",
    "analytic_energy_derivative_hartree_per_parameter",
    "finite_difference_energy_derivative_hartree_per_parameter",
    "energy_derivative_error_hartree_per_parameter",
    "energy_derivative_relative_error",
    "force_projection_minus_F_dot_d_hartree_per_bohr",
    "stress_conjugation_sigma_colon_G_GPa",
    "finite_difference_stress_conjugation_GPa",
    "stress_conjugation_error_GPa",
    "gpa_A3_to_hartree",
)


def _fingerprint(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text)
    temporary.replace(path)


def _write_immutable(path: Path, text: str) -> None:
    if path.exists():
        if path.read_text() != text:
            raise ValueError(f"refusing to replace a different immutable artifact: {path}")
        return
    _atomic_text(path, text)


def _matrix_det(matrix: Sequence[Sequence[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _matrix_inverse(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    a, b, c = matrix
    det = _matrix_det(matrix)
    if abs(det) < 1.0e-14:
        raise ValueError("singular reference cell")
    # Inverse of a matrix whose rows are a, b, c.
    return [
        [
            (b[1] * c[2] - b[2] * c[1]) / det,
            (a[2] * c[1] - a[1] * c[2]) / det,
            (a[1] * b[2] - a[2] * b[1]) / det,
        ],
        [
            (b[2] * c[0] - b[0] * c[2]) / det,
            (a[0] * c[2] - a[2] * c[0]) / det,
            (a[2] * b[0] - a[0] * b[2]) / det,
        ],
        [
            (b[0] * c[1] - b[1] * c[0]) / det,
            (a[1] * c[0] - a[0] * c[1]) / det,
            (a[0] * b[1] - a[1] * b[0]) / det,
        ],
    ]


def _row_times_matrix(row: Sequence[float], matrix: Sequence[Sequence[float]]) -> list[float]:
    return [sum(row[k] * matrix[k][j] for k in range(3)) for j in range(3)]


def _cell(text: str) -> list[list[float]]:
    lines = text.splitlines()
    start, end = cellopt.section_bounds(lines, "CELL")
    vectors: dict[str, list[float]] = {}
    for line in lines[start + 1 : end]:
        fields = line.split()
        if len(fields) == 4 and fields[0].upper() in {"A", "B", "C"}:
            vectors[fields[0].upper()] = [float(value.replace("D", "E")) for value in fields[1:]]
    if set(vectors) != {"A", "B", "C"}:
        raise ValueError("input lacks exactly one explicit A/B/C cell")
    return [vectors[key] for key in "ABC"]


def _scaled_coordinates(text: str) -> tuple[list[str], list[list[float]]]:
    lines = text.splitlines()
    start, end = cellopt.section_bounds(lines, "COORD")
    if not any(line.strip().upper() == "SCALED" for line in lines[start + 1 : end]):
        raise ValueError("FD input requires SCALED coordinates")
    elements: list[str] = []
    values: list[list[float]] = []
    for line in lines[start + 1 : end]:
        fields = line.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Za-z]{1,3}", fields[0]):
            elements.append(fields[0].capitalize())
            values.append([float(value.replace("D", "E")) for value in fields[1:4]])
    if not values:
        raise ValueError("FD input contains no coordinates")
    return elements, values


def _replace_project_and_run_type(text: str, project: str, run_type: str) -> str:
    lines = text.splitlines()
    replacements = (
        (re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I), f'"{project}"', "PROJECT"),
        (re.compile(r"^(\s*RUN_TYPE\s+).*$", re.I), run_type, "RUN_TYPE"),
    )
    for pattern, value, label in replacements:
        matches = [index for index, line in enumerate(lines) if pattern.match(line)]
        if len(matches) != 1:
            raise ValueError(f"expected one {label}, found {len(matches)}")
        match = pattern.match(lines[matches[0]])
        assert match is not None
        lines[matches[0]] = f"{match.group(1)}{value}"
    return "\n".join(lines) + "\n"


def _replace_scaled_coordinates(text: str, coordinates: Sequence[Sequence[float]]) -> str:
    lines = text.splitlines()
    start, end = cellopt.section_bounds(lines, "COORD")
    rows = [
        index
        for index in range(start + 1, end)
        if len(lines[index].split()) >= 4
        and re.fullmatch(r"[A-Za-z]{1,3}", lines[index].split()[0])
    ]
    if len(rows) != len(coordinates):
        raise ValueError("coordinate replacement count differs")
    for index, values in zip(rows, coordinates):
        indent = re.match(r"\s*", lines[index]).group(0)
        element = lines[index].split()[0]
        lines[index] = (
            f"{indent}{element:<3s} {values[0]: .16f} {values[1]: .16f} {values[2]: .16f}"
        )
    return "\n".join(lines) + "\n"


def _replace_cell(text: str, matrix: Sequence[Sequence[float]]) -> str:
    lines = text.splitlines()
    start, end = cellopt.section_bounds(lines, "CELL")
    for label, values in zip("ABC", matrix):
        matches = [
            index
            for index in range(start + 1, end)
            if re.match(rf"^\s*{label}\s+", lines[index], flags=re.I)
        ]
        if len(matches) != 1:
            raise ValueError(f"expected one CELL/{label} vector")
        index = matches[0]
        indent = re.match(r"\s*", lines[index]).group(0)
        lines[index] = f"{indent}{label} {values[0]:.16f} {values[1]:.16f} {values[2]:.16f}"
    return "\n".join(lines) + "\n"


def coordinate_directions(n_atoms: int, count: int = DEFAULT_DIRECTION_COUNT) -> list[dict[str, object]]:
    """Return deterministic, translation-free, orthonormal 3N directions."""

    if n_atoms < 2 or count not in (1, 2):
        raise ValueError("coordinate FD supports one or two directions and at least two atoms")
    vectors: list[list[float]] = []
    for direction_index in range(count):
        flat = [
            math.sin((index + 1) * (math.sqrt(2.0) + 0.37 * direction_index))
            + math.cos((index + 1) * (math.sqrt(3.0) + 0.23 * direction_index))
            for index in range(3 * n_atoms)
        ]
        for axis in range(3):
            mean = sum(flat[3 * atom + axis] for atom in range(n_atoms)) / n_atoms
            for atom in range(n_atoms):
                flat[3 * atom + axis] -= mean
        for previous in vectors:
            projection = sum(left * right for left, right in zip(flat, previous))
            flat = [left - projection * right for left, right in zip(flat, previous)]
        norm = math.sqrt(sum(value * value for value in flat))
        if norm < 1.0e-12:
            raise ValueError("degenerate deterministic FD direction")
        vectors.append([value / norm for value in flat])
    records: list[dict[str, object]] = []
    for index, flat in enumerate(vectors, start=1):
        shaped = [[flat[3 * atom + axis] for axis in range(3)] for atom in range(n_atoms)]
        records.append(
            {
                "id": f"collective_cartesian_{index}",
                "description": "deterministic translation-free normalized Cartesian 3N direction",
                "normalization": "sum_i,alpha d_i,alpha^2 = 1",
                "vector_cartesian": shaped,
                "vector_sha256": _fingerprint(shaped),
            }
        )
    return records


def coordinate_displaced_input(
    baseline_text: str,
    project: str,
    direction: Sequence[Sequence[float]],
    signed_step_bohr: float,
) -> str:
    matrix = _cell(baseline_text)
    _, scaled = _scaled_coordinates(baseline_text)
    if len(direction) != len(scaled):
        raise ValueError("FD direction atom count differs")
    # For row cell vectors, Cartesian r = fractional @ cell.  Therefore
    # fractional displacement = Cartesian displacement @ inverse(cell).
    inverse = _matrix_inverse(matrix)
    displaced: list[list[float]] = []
    for fractional, cartesian_direction in zip(scaled, direction):
        delta_cart = [
            signed_step_bohr * BOHR_TO_ANGSTROM * float(value)
            for value in cartesian_direction
        ]
        delta_fractional = _row_times_matrix(delta_cart, inverse)
        displaced.append([left + right for left, right in zip(fractional, delta_fractional)])
    text = _replace_project_and_run_type(baseline_text, project, "ENERGY")
    return _replace_scaled_coordinates(text, displaced)


def strained_input(
    baseline_text: str,
    project: str,
    generator: Sequence[Sequence[float]],
    signed_step: float,
) -> str:
    matrix = _cell(baseline_text)
    deformation = [
        [float(i == j) + signed_step * float(generator[i][j]) for j in range(3)]
        for i in range(3)
    ]
    # CP2K stores cell vectors as columns internally.  With row vectors here,
    # H'=(I+hG)H becomes H'_row=H_row(I+hG)^T; G is explicitly symmetric.
    deformed = [
        _row_times_matrix(row, [list(column) for column in zip(*deformation)])
        for row in matrix
    ]
    text = _replace_project_and_run_type(baseline_text, project, "ENERGY")
    return _replace_cell(text, deformed)


def _validate_input(text: str, *, baseline: bool) -> None:
    preflight.validate_preflight_input(
        text if baseline else _replace_project_and_run_type(text, "validation", "ENERGY_FORCE"),
        # validate_preflight_input uses this only in diagnostics and metadata lookup.
        _system_from_project(text),
    )
    expected_run_type = "ENERGY_FORCE" if baseline else "ENERGY"
    if len(re.findall(rf"^\s*RUN_TYPE\s+{expected_run_type}\s*$", text, flags=re.I | re.M)) != 1:
        raise ValueError(f"FD input lacks unique RUN_TYPE {expected_run_type}")


def _system_from_project(text: str) -> str:
    project_match = re.search(r"^\s*PROJECT(?:_NAME)?\s+\"?([^\"\s]+)", text, flags=re.I | re.M)
    if project_match is None:
        raise ValueError("FD input has no project")
    project = project_match.group(1)
    for system in DEFAULT_SYSTEMS:
        if project.startswith(system.replace("-", "_")) or project.startswith(system):
            return system
    raise ValueError(f"FD project does not identify a selected system: {project}")


def _job(
    run_dir: Path,
    system: str,
    job_id: str,
    job_type: str,
    text: str,
    **metadata: object,
) -> dict[str, object]:
    job_dir = run_dir / job_id
    project = f"{system}_GXTB_fd_{job_id}".replace("-", "_")
    input_path = job_dir / f"{project}.inp"
    return {
        "job_id": job_id,
        "job_type": job_type,
        "run_dir": str(job_dir.resolve()),
        "input": str(input_path.resolve()),
        "input_sha256": _sha256_text(text),
        "output": str((job_dir / OUTPUT_NAME).resolve()),
        "input_text": text,
        **metadata,
    }


def _case_payload(
    output_root: Path,
    system: str,
    coordinate_step_bohr: float,
    strain_step: float,
    direction_count: int,
) -> dict[str, object]:
    metadata = next(row for row in cellopt.systems() if str(row["id"]) == system)
    source, structure = cellopt.experimental_reference_paths(system)
    cellopt.validate_experimental_reference_source(system, source, structure)
    baseline_project = f"{system}_GXTB_fd_baseline".replace("-", "_")
    baseline_text = preflight.preflight_input_text(source, system, baseline_project)
    elements, _ = _scaled_coordinates(baseline_text)
    directions = coordinate_directions(len(elements), direction_count)
    run_dir = output_root / "GXTB" / system / VARIANT
    jobs: list[dict[str, object]] = [
        _job(run_dir, system, "baseline", "baseline_energy_force", baseline_text)
    ]
    for direction in directions:
        vector = direction["vector_cartesian"]
        assert isinstance(vector, list)
        for sign, sign_name in ((-1, "minus"), (1, "plus")):
            job_id = f"coord_{direction['id']}_{sign_name}"
            project = f"{system}_GXTB_fd_{job_id}".replace("-", "_")
            text = coordinate_displaced_input(
                baseline_text,
                project,
                vector,
                sign * coordinate_step_bohr,
            )
            jobs.append(
                _job(
                    run_dir,
                    system,
                    job_id,
                    "coordinate_energy",
                    text,
                    direction_id=direction["id"],
                    direction_sha256=direction["vector_sha256"],
                    sign=sign,
                    step=coordinate_step_bohr,
                    step_unit="bohr",
                )
            )
    for strain in STRAIN_DIRECTIONS:
        for sign, sign_name in ((-1, "minus"), (1, "plus")):
            job_id = f"strain_{strain['id']}_{sign_name}"
            project = f"{system}_GXTB_fd_{job_id}".replace("-", "_")
            text = strained_input(
                baseline_text,
                project,
                strain["generator"],
                sign * strain_step,
            )
            jobs.append(
                _job(
                    run_dir,
                    system,
                    job_id,
                    "strain_energy",
                    text,
                    direction_id=strain["id"],
                    generator=strain["generator"],
                    sign=sign,
                    step=strain_step,
                    step_unit="dimensionless_strain_parameter",
                )
            )
    clean_jobs: list[dict[str, object]] = []
    for job in jobs:
        text = str(job.pop("input_text"))
        baseline = job["job_type"] == "baseline_energy_force"
        _validate_input(text, baseline=baseline)
        clean_jobs.append({**job, "_prepared_text": text})
    matrix = _cell(baseline_text)
    return {
        "method": "GXTB",
        "system": system,
        "source_policy": SOURCE_POLICY,
        "source_input": str(source.resolve(strict=True)),
        "source_input_sha256": common.sha256_file(source),
        "structure_path": str(structure.resolve(strict=True)),
        "structure_sha256": common.sha256_file(structure),
        "structure_source": str(metadata["structure_source"]),
        "atom_count": len(elements),
        "elements": elements,
        "reference_cell_A": matrix,
        "reference_volume_A3": abs(_matrix_det(matrix)),
        "coordinate_directions": directions,
        "jobs": clean_jobs,
    }


def _manifest_payload(
    output_root: Path,
    campaign_identity: Mapping[str, object],
    coordinate_step_bohr: float,
    strain_step: float,
    direction_count: int,
) -> dict[str, object]:
    if (
        not math.isfinite(coordinate_step_bohr)
        or not math.isfinite(strain_step)
        or not 0.0 < coordinate_step_bohr <= 5.0e-2
        or not 0.0 < strain_step <= 1.0e-2
    ):
        raise ValueError(
            "finite-difference steps must be finite, positive, and remain small "
            "(coordinate <= 0.05 bohr; strain <= 0.01)"
        )
    if direction_count not in (1, 2):
        raise ValueError("coordinate direction count must be one or two")
    common.validate_campaign_identity(campaign_identity)
    cases = [
        _case_payload(output_root, system, coordinate_step_bohr, strain_step, direction_count)
        for system in DEFAULT_SYSTEMS
    ]
    # Prepared text is intentionally omitted from the immutable JSON; its hash
    # and the on-disk input are sufficient, while the temporary value lets
    # prepare perform all validation before writing anything.
    payload: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "phase": PHASE,
        "variant": VARIANT,
        "source_policy": SOURCE_POLICY,
        "scientific_status": "prepared_not_measured",
        "campaign_identity": dict(campaign_identity),
        "protocol": {
            "systems": list(DEFAULT_SYSTEMS),
            "mesh": "MACDONALD 2 2 2 0.25 0.25 0.25",
            "symmetry": "SPGLIB reduced; complete eight-row mesh mapping required",
            "baseline_run_type": "ENERGY_FORCE with analytical stress in GPa",
            "perturbed_run_type": "ENERGY",
            "coordinate_step_bohr": coordinate_step_bohr,
            "coordinate_direction_count": direction_count,
            "coordinate_derivative_identity": "dE/dq = -sum_i F_i dot d_i",
            "strain_step": strain_step,
            "strain_directions": list(STRAIN_DIRECTIONS),
            "stress_derivative_identity": "dE/dh = -V*(sigma:G)*GPa*A^3_to_Eh",
            "gpa_A3_to_hartree": GPA_ANGSTROM3_TO_HARTREE,
            "bohr_to_angstrom": BOHR_TO_ANGSTROM,
        },
        "cases": cases,
    }
    serializable = json.loads(json.dumps(payload))
    for case in serializable["cases"]:
        for job in case["jobs"]:
            job.pop("_prepared_text", None)
    serializable["payload_sha256"] = _fingerprint(serializable)
    return {"serializable": serializable, "prepared": payload}


def manifest_path(output_root: Path) -> Path:
    return output_root.resolve() / MANIFEST_NAME


def prepare(
    output_root: Path,
    campaign_identity: Mapping[str, object],
    *,
    coordinate_step_bohr: float = DEFAULT_COORDINATE_STEP_BOHR,
    strain_step: float = DEFAULT_STRAIN_STEP,
    direction_count: int = DEFAULT_DIRECTION_COUNT,
) -> Path:
    output_root = output_root.resolve()
    path = manifest_path(output_root)
    if path.exists():
        load_manifest(output_root, campaign_identity)
        payload = json.loads(path.read_text())
        protocol = payload["protocol"]
        requested = (coordinate_step_bohr, strain_step, direction_count)
        frozen = (
            float(protocol["coordinate_step_bohr"]),
            float(protocol["strain_step"]),
            int(protocol["coordinate_direction_count"]),
        )
        if requested != frozen:
            raise ValueError(f"FD manifest already freezes steps/directions as {frozen}")
        return path
    bundle = _manifest_payload(
        output_root,
        campaign_identity,
        coordinate_step_bohr,
        strain_step,
        direction_count,
    )
    serializable = bundle["serializable"]
    prepared = bundle["prepared"]
    assert isinstance(serializable, dict) and isinstance(prepared, dict)
    for case in prepared["cases"]:
        for job in case["jobs"]:
            run_dir = Path(str(job["run_dir"]))
            input_path = Path(str(job["input"]))
            output = Path(str(job["output"]))
            if output.exists() or (run_dir / common.JOB_STAMP_NAME).exists():
                raise ValueError(f"refusing to prepare over pre-existing FD job state: {run_dir}")
            text = str(job["_prepared_text"])
            if input_path.exists() and input_path.read_text() != text:
                raise ValueError(f"stale FD input differs: {input_path}")
    # All cases have been validated before the first write.
    for case in prepared["cases"]:
        for job in case["jobs"]:
            input_path = Path(str(job["input"]))
            input_path.parent.mkdir(parents=True, exist_ok=True)
            _write_immutable(input_path, str(job["_prepared_text"]))
    _write_immutable(path, json.dumps(serializable, indent=2, sort_keys=True) + "\n")
    return path


def load_manifest(
    output_root: Path,
    campaign_identity: Mapping[str, object],
) -> dict[str, dict[str, object]]:
    path = manifest_path(output_root)
    payload = json.loads(path.read_text())
    digest = payload.pop("payload_sha256", None)
    if digest != _fingerprint(payload):
        raise ValueError(f"FD manifest payload fingerprint differs: {path}")
    payload["payload_sha256"] = digest
    if (
        payload.get("schema") != MANIFEST_SCHEMA
        or payload.get("phase") != PHASE
        or payload.get("variant") != VARIANT
        or payload.get("source_policy") != SOURCE_POLICY
        or payload.get("campaign_identity") != dict(campaign_identity)
    ):
        raise ValueError(f"invalid or foreign FD manifest: {path}")
    protocol = payload.get("protocol")
    if not isinstance(protocol, dict) or protocol.get("systems") != list(DEFAULT_SYSTEMS):
        raise ValueError("FD manifest does not contain the frozen four-system pilot")
    try:
        expected_payload = _manifest_payload(
            output_root.resolve(),
            campaign_identity,
            float(protocol["coordinate_step_bohr"]),
            float(protocol["strain_step"]),
            int(protocol["coordinate_direction_count"]),
        )["serializable"]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("FD manifest protocol is incomplete or invalid") from exc
    if payload != expected_payload:
        raise ValueError(
            "FD manifest is not the deterministic rendering of its frozen "
            "campaign/source/protocol"
        )
    cases: dict[str, dict[str, object]] = {}
    for case in payload.get("cases", []):
        if not isinstance(case, dict):
            raise ValueError("invalid FD case record")
        system = str(case.get("system", ""))
        if system in cases or system not in DEFAULT_SYSTEMS:
            raise ValueError(f"duplicate or unknown FD system: {system}")
        source, structure = cellopt.experimental_reference_paths(system)
        if Path(str(case["source_input"])).resolve() != source.resolve(strict=True):
            raise ValueError(f"noncanonical FD source input for {system}")
        if Path(str(case["structure_path"])).resolve() != structure.resolve(strict=True):
            raise ValueError(f"noncanonical FD structure for {system}")
        if common.sha256_file(source) != case.get("source_input_sha256"):
            raise ValueError(f"FD source input changed for {system}")
        if common.sha256_file(structure) != case.get("structure_sha256"):
            raise ValueError(f"FD structure changed for {system}")
        directions = case.get("coordinate_directions")
        if not isinstance(directions, list) or len(directions) != int(protocol["coordinate_direction_count"]):
            raise ValueError(f"FD direction count differs for {system}")
        expected_directions = coordinate_directions(int(case["atom_count"]), len(directions))
        if directions != expected_directions:
            raise ValueError(f"FD direction definition differs for {system}")
        jobs = case.get("jobs")
        if not isinstance(jobs, list) or len(jobs) != 1 + 2 * len(directions) + 2 * len(STRAIN_DIRECTIONS):
            raise ValueError(f"incomplete FD job set for {system}")
        seen: set[str] = set()
        for job in jobs:
            if not isinstance(job, dict):
                raise ValueError(f"invalid FD job for {system}")
            job_id = str(job.get("job_id", ""))
            if not job_id or job_id in seen:
                raise ValueError(f"duplicate FD job for {system}: {job_id}")
            seen.add(job_id)
            run_dir = Path(str(job["run_dir"])).resolve()
            expected_dir = output_root.resolve() / "GXTB" / system / VARIANT / job_id
            input_path = Path(str(job["input"])).resolve()
            output = Path(str(job["output"])).resolve()
            if run_dir != expected_dir or input_path.parent != run_dir or output != run_dir / OUTPUT_NAME:
                raise ValueError(f"noncanonical FD job paths for {system}/{job_id}")
            if common.sha256_file(input_path) != job.get("input_sha256"):
                raise ValueError(f"FD input changed for {system}/{job_id}")
            _validate_input(input_path.read_text(), baseline=job["job_type"] == "baseline_energy_force")
        cases[system] = case
    if tuple(cases) != DEFAULT_SYSTEMS:
        raise ValueError("FD manifest case order/coverage differs")
    return cases


def _manifest_protocol(output_root: Path) -> dict[str, object]:
    payload = json.loads(manifest_path(output_root).read_text())
    return dict(payload["protocol"])


def source_artifacts(output_root: Path, case: Mapping[str, object]) -> dict[str, Path]:
    return {
        "fd_manifest": manifest_path(output_root),
        "reference_input": Path(str(case["source_input"])),
        "reference_structure": Path(str(case["structure_path"])),
    }


def protocol_identity(output_root: Path, job: Mapping[str, object]) -> dict[str, object]:
    protocol = _manifest_protocol(output_root)
    return {
        "gate_phase": PHASE,
        "variant": VARIANT,
        "source_policy": SOURCE_POLICY,
        "manifest_payload_sha256": json.loads(manifest_path(output_root).read_text())["payload_sha256"],
        "mesh": protocol["mesh"],
        "coordinate_step_bohr": protocol["coordinate_step_bohr"],
        "strain_step": protocol["strain_step"],
        "job_id": job["job_id"],
        "job_type": job["job_type"],
        "direction_id": job.get("direction_id"),
        "direction_sha256": job.get("direction_sha256"),
        "generator": job.get("generator"),
        "sign": job.get("sign"),
    }


def _kpoint_summary(text: str) -> dict[str, object]:
    counts = [int(value) for value in re.findall(r"Number of Special K-points:\s*(\d+)", text, flags=re.I)]
    if not counts:
        raise ValueError("FD output has no evaluated k-point count")
    mesh_matches = re.findall(
        r"K-point Mesh:\s+(\d+)\s+(\d+)\s+(\d+)\s*$", text, flags=re.I | re.M
    )
    if not mesh_matches or [int(value) for value in mesh_matches[-1]] != [2, 2, 2]:
        raise ValueError("FD output has no shifted-k222 mesh report")
    lines = text.splitlines()
    header = max(
        index for index, line in enumerate(lines) if re.search(r"K-point Mesh:\s+2\s+2\s+2\s*$", line, re.I)
    )
    row_pattern = re.compile(
        rf"^\s*(\d+)\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}"
        r"\s+\d+\s+\d+\s+\d+\s*$"
    )
    rows: list[int] = []
    for line in lines[header + 1 :]:
        match = row_pattern.match(line)
        if match:
            rows.append(int(match.group(1)))
        elif rows:
            break
    if rows != list(range(1, 9)):
        raise ValueError("FD output has incomplete shifted-k222 full-mesh mapping")
    if not 1 <= counts[-1] <= 8:
        raise ValueError("FD output reports an invalid irreducible k-point count")
    return {"kpoint_count": counts[-1], "kpoint_mesh": [2, 2, 2], "kpoint_mesh_rows": rows}


def parse_energy_output(path: Path) -> dict[str, object]:
    text = path.read_text(errors="replace")
    if "PROGRAM ENDED" not in text:
        raise ValueError(f"CP2K did not end normally: {path}")
    if re.search(r"\*\*\*.*ABORT|SCF\s+(?:run\s+)?(?:did\s+)?NOT\s+converged", text, flags=re.I):
        raise ValueError(f"fatal/nonconverged marker in {path}")
    energies = re.findall(
        rf"ENERGY\| Total FORCE_EVAL[^\n]*?({FLOAT_PATTERN})\s*$", text, flags=re.I | re.M
    )
    if not energies:
        raise ValueError(f"total FORCE_EVAL energy missing from {path}")
    energy = float(energies[-1].replace("D", "E").replace("d", "e"))
    if not math.isfinite(energy):
        raise ValueError(f"nonfinite energy in {path}")
    return {"energy_hartree": energy, **_kpoint_summary(text)}


def parse_job_output(case: Mapping[str, object], job: Mapping[str, object]) -> dict[str, object]:
    output = Path(str(job["output"]))
    if job["job_type"] == "baseline_energy_force":
        parsed = derivative_parser.parse_cp2k_output(
            output,
            int(case["atom_count"]),
            [str(value) for value in case["elements"]],
        )
        if parsed["kpoint_mesh"] != [2, 2, 2] or parsed["kpoint_mesh_row_indices"] != list(range(1, 9)):
            raise ValueError(f"baseline output has incomplete shifted-k222 mapping: {output}")
        return parsed
    return parse_energy_output(output)


def _recorded_job_matches(
    output_root: Path,
    case: Mapping[str, object],
    job: Mapping[str, object],
    campaign_identity: Mapping[str, object],
) -> tuple[bool, str]:
    return common.recorded_job_stamp_matches(
        Path(str(job["run_dir"])),
        Path(str(job["input"])),
        "GXTB",
        PHASE,
        Path(str(job["output"])),
        campaign_identity=campaign_identity,
        accepted_status_prefixes=("converged",),
        protocol_identity=protocol_identity(output_root, job),
        source_artifacts=source_artifacts(output_root, case),
    )


def run_one(
    output_root: Path,
    case: Mapping[str, object],
    job: Mapping[str, object],
    cp2k: Path,
    threads: int,
    campaign_identity: Mapping[str, object],
) -> tuple[str, str, int, str]:
    system, job_id = str(case["system"]), str(job["job_id"])
    run_dir = Path(str(job["run_dir"]))
    input_path = Path(str(job["input"]))
    output = Path(str(job["output"]))
    with input_path.open() as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return system, job_id, common.BUSY_RETURN_CODE, "BUSY"
        stamp_path = run_dir / common.JOB_STAMP_NAME
        stamp_ok, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            "GXTB",
            PHASE,
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity(output_root, job),
            source_artifacts=source_artifacts(output_root, case),
        )
        if output.exists():
            recorded_ok, _ = _recorded_job_matches(output_root, case, job, campaign_identity)
            if not stamp_ok or not recorded_ok:
                return system, job_id, 1, "STALE_OUTPUT"
            try:
                parse_job_output(case, job)
            except ValueError:
                return system, job_id, 1, "STALE_OUTPUT"
            return system, job_id, 0, "SKIP"
        if stamp_path.exists():
            if not stamp_ok:
                return system, job_id, 1, "STALE_STAMP"
            try:
                prior_status = str(json.loads(stamp_path.read_text()).get("status", ""))
            except (json.JSONDecodeError, OSError):
                return system, job_id, 1, "STALE_STAMP"
            if prior_status.startswith("converged"):
                # A completed stamp with a missing output is evidence loss, not
                # an ordinary failed attempt which may be resumed in place.
                return system, job_id, 1, "STALE_STAMP"
        process = subprocess.run(
            [str(cp2k.resolve(strict=True)), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            env=common.thread_environment(threads),
            check=False,
        )
        (run_dir / "returncode.txt").write_text(f"{process.returncode}\n")
        details: dict[str, object] = {"returncode": process.returncode, "output": str(output)}
        status, action, code = "failed", "FAILED", process.returncode or 1
        if output.is_file():
            details["output_sha256"] = common.sha256_file(output)
        if process.returncode == 0 and output.is_file():
            try:
                details["parsed"] = parse_job_output(case, job)
            except ValueError as exc:
                details["parse_error"] = str(exc)
                status, action, code = "failed_parse", "INVALID_OUTPUT", 1
            else:
                status, action, code = "converged_measured_not_approved", "COMPLETED", 0
        common.write_job_stamp(
            run_dir,
            input_path,
            cp2k,
            "GXTB",
            PHASE,
            status,
            details=details,
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity(output_root, job),
            source_artifacts=source_artifacts(output_root, case),
        )
        return system, job_id, code, action


def _jobs_by_id(case: Mapping[str, object]) -> dict[str, dict[str, object]]:
    return {str(job["job_id"]): job for job in case["jobs"]}  # type: ignore[index]


def _completed_job(
    output_root: Path,
    case: Mapping[str, object],
    job: Mapping[str, object],
    campaign_identity: Mapping[str, object],
) -> dict[str, object]:
    valid, reason = _recorded_job_matches(output_root, case, job, campaign_identity)
    if not valid:
        raise ValueError(f"untrusted FD result for {case['system']}/{job['job_id']}: {reason}")
    parsed = parse_job_output(case, job)
    return {
        **parsed,
        "input_sha256": common.sha256_file(Path(str(job["input"]))),
        "output_sha256": common.sha256_file(Path(str(job["output"]))),
        "stamp_sha256": common.sha256_file(Path(str(job["run_dir"])) / common.JOB_STAMP_NAME),
    }


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-14)


def measured_rows(
    output_root: Path,
    campaign_identity: Mapping[str, object],
) -> list[dict[str, object]]:
    cases = load_manifest(output_root, campaign_identity)
    manifest_sha = common.sha256_file(manifest_path(output_root))
    rows: list[dict[str, object]] = []
    for system in DEFAULT_SYSTEMS:
        case = cases[system]
        jobs = _jobs_by_id(case)
        baseline_job = jobs["baseline"]
        baseline = _completed_job(output_root, case, baseline_job, campaign_identity)
        forces = [row["vector_hartree_per_bohr"] for row in baseline["forces"]]  # type: ignore[index]
        stress = baseline["stress_gpa"]
        common_fields: dict[str, object] = {
            "method": "GXTB",
            "system": system,
            "phase": PHASE,
            "variant": VARIANT,
            "source_policy": SOURCE_POLICY,
            "scientific_status": "measured_not_approved",
            "approved": False,
            "campaign_fingerprint_sha256": campaign_identity["fingerprint_sha256"],
            "manifest_sha256": manifest_sha,
            "source_input_sha256": case["source_input_sha256"],
            "structure_sha256": case["structure_sha256"],
            "baseline_input_sha256": baseline["input_sha256"],
            "baseline_output_sha256": baseline["output_sha256"],
            "reference_volume_A3": float(case["reference_volume_A3"]),
            "gpa_A3_to_hartree": GPA_ANGSTROM3_TO_HARTREE,
        }
        for direction in case["coordinate_directions"]:  # type: ignore[index]
            direction_id = str(direction["id"])
            minus_job = jobs[f"coord_{direction_id}_minus"]
            plus_job = jobs[f"coord_{direction_id}_plus"]
            minus = _completed_job(output_root, case, minus_job, campaign_identity)
            plus = _completed_job(output_root, case, plus_job, campaign_identity)
            step = float(minus_job["step"])
            fd = (float(plus["energy_hartree"]) - float(minus["energy_hartree"])) / (2.0 * step)
            vector = direction["vector_cartesian"]
            analytic = -sum(
                float(force_component) * float(direction_component)
                for force, direction_row in zip(forces, vector)
                for force_component, direction_component in zip(force, direction_row)
            )
            rows.append(
                {
                    **common_fields,
                    "measurement_type": "coordinate_directional_derivative",
                    "direction_id": direction_id,
                    "minus_input_sha256": minus["input_sha256"],
                    "minus_output_sha256": minus["output_sha256"],
                    "plus_input_sha256": plus["input_sha256"],
                    "plus_output_sha256": plus["output_sha256"],
                    "step": step,
                    "step_unit": "bohr",
                    "analytic_energy_derivative_hartree_per_parameter": analytic,
                    "finite_difference_energy_derivative_hartree_per_parameter": fd,
                    "energy_derivative_error_hartree_per_parameter": fd - analytic,
                    "energy_derivative_relative_error": _relative_error(fd, analytic),
                    "force_projection_minus_F_dot_d_hartree_per_bohr": analytic,
                    "stress_conjugation_sigma_colon_G_GPa": "",
                    "finite_difference_stress_conjugation_GPa": "",
                    "stress_conjugation_error_GPa": "",
                }
            )
        for strain in STRAIN_DIRECTIONS:
            direction_id = str(strain["id"])
            minus_job = jobs[f"strain_{direction_id}_minus"]
            plus_job = jobs[f"strain_{direction_id}_plus"]
            minus = _completed_job(output_root, case, minus_job, campaign_identity)
            plus = _completed_job(output_root, case, plus_job, campaign_identity)
            step = float(minus_job["step"])
            fd = (float(plus["energy_hartree"]) - float(minus["energy_hartree"])) / (2.0 * step)
            generator = strain["generator"]
            conjugation = sum(
                float(stress[i][j]) * float(generator[i][j])  # type: ignore[index]
                for i in range(3)
                for j in range(3)
            )
            volume = float(case["reference_volume_A3"])
            analytic = -volume * conjugation * GPA_ANGSTROM3_TO_HARTREE
            fd_conjugation = -fd / (volume * GPA_ANGSTROM3_TO_HARTREE)
            rows.append(
                {
                    **common_fields,
                    "measurement_type": "symmetric_strain_directional_derivative",
                    "direction_id": direction_id,
                    "minus_input_sha256": minus["input_sha256"],
                    "minus_output_sha256": minus["output_sha256"],
                    "plus_input_sha256": plus["input_sha256"],
                    "plus_output_sha256": plus["output_sha256"],
                    "step": step,
                    "step_unit": "dimensionless_strain_parameter",
                    "analytic_energy_derivative_hartree_per_parameter": analytic,
                    "finite_difference_energy_derivative_hartree_per_parameter": fd,
                    "energy_derivative_error_hartree_per_parameter": fd - analytic,
                    "energy_derivative_relative_error": _relative_error(fd, analytic),
                    "force_projection_minus_F_dot_d_hartree_per_bohr": "",
                    "stress_conjugation_sigma_colon_G_GPa": conjugation,
                    "finite_difference_stress_conjugation_GPa": fd_conjugation,
                    "stress_conjugation_error_GPa": fd_conjugation - conjugation,
                }
            )
    return rows


def _csv_text(rows: Sequence[Mapping[str, object]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return stream.getvalue()


def collect(
    output_root: Path,
    csv_path: Path,
    json_path: Path,
    campaign_identity: Mapping[str, object],
) -> dict[str, object]:
    rows = measured_rows(output_root, campaign_identity)
    expected_rows = len(DEFAULT_SYSTEMS) * (DEFAULT_DIRECTION_COUNT + len(STRAIN_DIRECTIONS))
    # A one-direction custom manifest is valid; derive its exact expected count.
    direction_count = int(_manifest_protocol(output_root)["coordinate_direction_count"])
    expected_rows = len(DEFAULT_SYSTEMS) * (direction_count + len(STRAIN_DIRECTIONS))
    if len(rows) != expected_rows:
        raise ValueError(f"incomplete FD report: {len(rows)}/{expected_rows} rows")
    csv_text = _csv_text(rows)
    _write_immutable(csv_path.resolve(), csv_text)
    report: dict[str, object] = {
        "schema": REPORT_SCHEMA,
        "phase": PHASE,
        "variant": VARIANT,
        "scientific_status": "measured_not_approved",
        "approved": False,
        "campaign_identity": dict(campaign_identity),
        "manifest": {
            "path": str(manifest_path(output_root)),
            "sha256": common.sha256_file(manifest_path(output_root)),
        },
        "measured_csv": {
            "path": str(csv_path.resolve()),
            "sha256": common.sha256_file(csv_path.resolve()),
        },
        "formulae": {
            "coordinate": "central dE/dq compared with -sum_i F_i dot d_i",
            "strain": "central dE/dh compared with -V*(sigma:G)*GPa*A^3_to_Eh",
            "stress_sign": "CP2K sigma = -(1/V) dE/d(strain)",
            "gpa_A3_to_hartree": GPA_ANGSTROM3_TO_HARTREE,
        },
        "systems": list(DEFAULT_SYSTEMS),
        "row_count": len(rows),
        "rows": rows,
    }
    _write_immutable(json_path.resolve(), json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def approve_report(
    report_json: Path,
    approval_json: Path,
    *,
    reviewer: str,
    coordinate_abs_tolerance_hartree_per_bohr: float,
    stress_abs_tolerance_gpa: float,
) -> dict[str, object]:
    if not reviewer.strip():
        raise ValueError("an explicit nonempty reviewer is required")
    if coordinate_abs_tolerance_hartree_per_bohr <= 0.0 or stress_abs_tolerance_gpa <= 0.0:
        raise ValueError("approval tolerances must be positive")
    report_path = report_json.resolve(strict=True)
    report = json.loads(report_path.read_text())
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("phase") != PHASE
        or report.get("scientific_status") != "measured_not_approved"
        or report.get("approved") is not False
    ):
        raise ValueError("approval source is not an unapproved FD measurement report")
    csv_record = report.get("measured_csv")
    manifest_record = report.get("manifest")
    if not isinstance(csv_record, dict) or not isinstance(manifest_record, dict):
        raise ValueError("FD report lacks measured artifact records")
    csv_path = Path(str(csv_record["path"])).resolve(strict=True)
    manifest = Path(str(manifest_record["path"])).resolve(strict=True)
    if common.sha256_file(csv_path) != csv_record.get("sha256"):
        raise ValueError("measured FD CSV changed after collection")
    if common.sha256_file(manifest) != manifest_record.get("sha256"):
        raise ValueError("FD manifest changed after collection")
    rows = report.get("rows")
    if not isinstance(rows, list) or len(rows) != int(report.get("row_count", -1)):
        raise ValueError("FD report rows are incomplete")
    if _csv_text(rows) != csv_path.read_text():
        raise ValueError("FD JSON rows do not reproduce the frozen measured CSV")
    checks: list[dict[str, object]] = []
    for row in rows:
        measurement_type = str(row["measurement_type"])
        if measurement_type == "coordinate_directional_derivative":
            error = abs(float(row["energy_derivative_error_hartree_per_parameter"]))
            tolerance = coordinate_abs_tolerance_hartree_per_bohr
            unit = "hartree_per_bohr"
        elif measurement_type == "symmetric_strain_directional_derivative":
            error = abs(float(row["stress_conjugation_error_GPa"]))
            tolerance = stress_abs_tolerance_gpa
            unit = "GPa"
        else:
            raise ValueError(f"unknown FD measurement type: {measurement_type}")
        checks.append(
            {
                "system": row["system"],
                "direction_id": row["direction_id"],
                "measurement_type": measurement_type,
                "absolute_error": error,
                "tolerance": tolerance,
                "unit": unit,
                "passed": error <= tolerance,
            }
        )
    passed = all(bool(check["passed"]) for check in checks)
    approval: dict[str, object] = {
        "schema": APPROVAL_SCHEMA,
        "phase": PHASE,
        "decision": "approved" if passed else "rejected",
        "reviewer": reviewer.strip(),
        "report_json": {"path": str(report_path), "sha256": common.sha256_file(report_path)},
        "measured_csv": csv_record,
        "manifest": manifest_record,
        "thresholds": {
            "coordinate_abs_tolerance_hartree_per_bohr": coordinate_abs_tolerance_hartree_per_bohr,
            "stress_abs_tolerance_gpa": stress_abs_tolerance_gpa,
        },
        "checks": checks,
        "passed_count": sum(bool(check["passed"]) for check in checks),
        "check_count": len(checks),
    }
    _write_immutable(approval_json.resolve(), json.dumps(approval, indent=2, sort_keys=True) + "\n")
    return approval


def _register_provenance_root(output_root: Path, **artifacts: Path | None) -> None:
    provenance = ROOT / "data" / common.GXTB_PROVENANCE_NAME
    if provenance.is_file():
        payload = json.loads(provenance.read_text())
        paths = payload.get("workflow_paths", {})
        if isinstance(paths, dict) and paths.get("k222_finite_difference_gate_root"):
            frozen = Path(str(paths["k222_finite_difference_gate_root"])).resolve()
            if frozen != output_root.resolve():
                raise ValueError(f"k222 finite-difference gate root is already frozen as {frozen}")
    common.update_gxtb_provenance(
        ROOT,
        cp2k=artifacts.get("cp2k"),
        cp2k_source=artifacts.get("cp2k_source"),
        save_tblite=artifacts.get("save_tblite"),
        save_tblite_source=artifacts.get("save_tblite_source"),
        campaign_manifest=artifacts.get("campaign_manifest"),
        workflow_paths={"k222_finite_difference_gate_root": output_root},
    )


def prepare_command(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / f".{output_root.name}.k222_fd_gate.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        campaign = common.load_campaign_identity(ROOT)
        path = prepare(
            output_root,
            campaign,
            coordinate_step_bohr=args.coordinate_step_bohr,
            strain_step=args.strain_step,
            direction_count=args.coordinate_directions,
        )
    print(f"Prepared immutable {len(DEFAULT_SYSTEMS)}-system k222 FD pilot: {path}")


def run_command(args: argparse.Namespace) -> None:
    common.require_gxtb_build_artifacts(
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    observed = common.validate_campaign_artifacts(
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )[0]
    campaign = common.load_campaign_identity(ROOT)
    if observed != campaign:
        raise ValueError("run artifacts differ from the frozen X23b campaign")
    cases = load_manifest(args.output_root, campaign)
    selected_systems = list(dict.fromkeys(args.system or list(DEFAULT_SYSTEMS)))
    selected: list[tuple[dict[str, object], dict[str, object]]] = []
    for system in selected_systems:
        case = cases[system]
        for job in case["jobs"]:  # type: ignore[index]
            if not args.job or str(job["job_id"]) in args.job:
                selected.append((case, job))
    if not selected:
        raise ValueError("no FD jobs selected")
    _register_provenance_root(
        args.output_root,
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    failed: list[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_one,
                args.output_root.resolve(),
                case,
                job,
                args.cp2k,
                args.threads_per_job,
                campaign,
            ): (case["system"], job["job_id"])
            for case, job in selected
        }
        for future in as_completed(futures):
            system, job_id, code, action = future.result()
            print(f"{action:14s} GXTB/{system}/{job_id} rc={code}", flush=True)
            if code != 0:
                failed.append(f"{system}/{job_id}")
    if failed:
        raise SystemExit(f"{len(failed)} FD job(s) failed or stale")


def collect_command(args: argparse.Namespace) -> None:
    campaign = common.load_campaign_identity(ROOT)
    report = collect(args.output_root, args.csv, args.json, campaign)
    _register_provenance_root(args.output_root)
    print(f"Collected {report['row_count']} measured, unapproved FD comparisons in {args.json}")


def approve_command(args: argparse.Namespace) -> None:
    approval = approve_report(
        args.report_json,
        args.approval_json,
        reviewer=args.reviewer,
        coordinate_abs_tolerance_hartree_per_bohr=(
            args.coordinate_abs_tolerance_hartree_per_bohr
        ),
        stress_abs_tolerance_gpa=args.stress_abs_tolerance_gpa,
    )
    print(f"Separate FD review decision: {approval['decision']} ({args.approval_json})")


def add_artifact_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--campaign-manifest", type=Path, default=common.DEFAULT_CAMPAIGN_MANIFEST)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.add_argument("--coordinate-step-bohr", type=float, default=DEFAULT_COORDINATE_STEP_BOHR)
    prepare_parser.add_argument("--strain-step", type=float, default=DEFAULT_STRAIN_STEP)
    prepare_parser.add_argument("--coordinate-directions", type=int, choices=(1, 2), default=DEFAULT_DIRECTION_COUNT)
    prepare_parser.set_defaults(function=prepare_command)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--system", action="append", choices=DEFAULT_SYSTEMS)
    run_parser.add_argument("--job", action="append", help="exact manifest job_id; repeat as needed")
    run_parser.add_argument("--jobs", type=int, default=2)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    add_artifact_arguments(run_parser)
    run_parser.set_defaults(function=run_command)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.add_argument("--json", type=Path, required=True)
    collect_parser.set_defaults(function=collect_command)

    approval_parser = subparsers.add_parser("approve")
    approval_parser.add_argument("--report-json", type=Path, required=True)
    approval_parser.add_argument("--approval-json", type=Path, required=True)
    approval_parser.add_argument("--reviewer", required=True)
    approval_parser.add_argument(
        "--coordinate-abs-tolerance-hartree-per-bohr", type=float, required=True
    )
    approval_parser.add_argument("--stress-abs-tolerance-gpa", type=float, required=True)
    approval_parser.set_defaults(function=approve_command)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        args.function(args)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
