#!/usr/bin/env python3
"""Fail-closed adaptive native-Bloch k-point convergence for X23b.

Two independent series are frozen for all 23 X23b crystals and for GFN1,
GFN2, and g-xTB on one production-ready CP2K/save_tblite build:

* fixed experimental crystal geometry single points; and
* independent full-cell optimizations started from that same geometry.

The workflow never constructs Born--von Karman supercells.  Every mesh is a
native CP2K ``MACDONALD N N N`` mesh with SPGLIB symmetry reduction and
``FULL_GRID F``.  One consecutive passing step is sufficient and the denser
result is selected.  Publication files are written atomically only after all
138 method/system/series tracks have converged.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable, Mapping

import x23b_common as common
import x23b_pipeline as pipeline


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
METHODS = tuple(common.METHODS)
SERIES = ("fixed_experimental_sp", "independent_cellopt")
MESH_NUMBERS = tuple(range(1, 9))
MESH_IDS = tuple(f"k{n}{n}{n}" for n in MESH_NUMBERS)
PROTOCOL_ID = "x23b-adaptive-native-bloch-kpoint-convergence-v1"
REQUIRED_CAMPAIGN_MARKER = "post5582"
WORKFLOW_SCHEMA = 1
OUTPUT_SCHEMA = 1
OUTPUT_STEM = "x23b_adaptive_kpoint_convergence"
HARTREE_TO_KJMOL = 2625.499638

ENERGY_TOLERANCE_KJMOL = 0.05
LENGTH_TOLERANCE_PERCENT = 0.05
VOLUME_TOLERANCE_PERCENT = 0.10
ANGLE_TOLERANCE_DEG = 0.05
COMPARISON_ABSOLUTE_TOLERANCE = 1.0e-10
MINIMUM_DENSE_MESH = {
    "fixed_experimental_sp": 2,
    "independent_cellopt": 3,
}

CSV_FIELDS = (
    "series",
    "method",
    "method_label",
    "system",
    "coarse_mesh",
    "dense_mesh",
    "eligible_for_stopping",
    "selected_step",
    "coarse_energy_hartree_per_molecule",
    "dense_energy_hartree_per_molecule",
    "energy_delta_kJmol",
    "energy_abs_delta_kJmol",
    "energy_pass",
    "coarse_a_A",
    "coarse_b_A",
    "coarse_c_A",
    "dense_a_A",
    "dense_b_A",
    "dense_c_A",
    "a_relative_delta_percent",
    "b_relative_delta_percent",
    "c_relative_delta_percent",
    "max_length_relative_delta_percent",
    "length_pass",
    "coarse_alpha_deg",
    "coarse_beta_deg",
    "coarse_gamma_deg",
    "dense_alpha_deg",
    "dense_beta_deg",
    "dense_gamma_deg",
    "alpha_abs_delta_deg",
    "beta_abs_delta_deg",
    "gamma_abs_delta_deg",
    "max_angle_abs_delta_deg",
    "angle_pass",
    "coarse_volume_A3",
    "dense_volume_A3",
    "volume_relative_delta_percent",
    "volume_pass",
    "all_required_criteria_pass",
    "coarse_output_sha256",
    "dense_output_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def method_label(method: str) -> str:
    return "g-xTB" if method == "GXTB" else f"{method}-xTB"


def require_post5582_campaign(binding: Mapping[str, object]) -> None:
    """Reject production manifests that predate the final post-5582 build."""
    identity = binding.get("campaign_identity")
    manifest = binding.get("campaign_manifest")
    if not isinstance(identity, Mapping) or not isinstance(manifest, Mapping):
        raise ValueError("workflow build lacks campaign identity records")
    campaign_ids = (
        str(identity.get("campaign_id", "")),
        str(manifest.get("campaign_id", "")),
    )
    normalized = tuple(re.sub(r"[^a-z0-9]", "", value.lower()) for value in campaign_ids)
    if any(REQUIRED_CAMPAIGN_MARKER not in value for value in normalized):
        raise ValueError(
            "X23b adaptive production is held for a production_ready post-5582 campaign"
        )
    if normalized[0] != normalized[1]:
        raise ValueError("campaign identity and manifest campaign IDs differ")


def systems() -> tuple[dict[str, object], ...]:
    metadata = json.loads((ROOT / "data" / "metadata.json").read_text())
    result = tuple(metadata["systems"])
    ids = tuple(str(row["id"]) for row in result)
    pipeline_ids = tuple(str(row["id"]) for row in pipeline.SYSTEMS)
    if len(result) != 23 or len(set(ids)) != 23 or ids != pipeline_ids:
        raise ValueError("X23b system order is not the exact canonical 23-system set")
    return result


def system_metadata(system: str) -> dict[str, object]:
    try:
        return next(row for row in systems() if str(row["id"]) == system)
    except StopIteration as error:
        raise ValueError(f"unknown X23b system: {system}") from error


def source_cif(system: str) -> Path:
    path = ROOT / "structures" / "cif" / f"{system}.cif"
    return path.resolve(strict=True)


def gamma_centered_shift(mesh: int) -> float:
    if mesh not in MESH_NUMBERS:
        raise ValueError(f"unsupported X23b convergence mesh: {mesh}")
    return 0.0 if mesh % 2 else (mesh - 1) / (2.0 * mesh)


def mesh_definition(mesh: int) -> dict[str, object]:
    shift = gamma_centered_shift(mesh)
    return {
        "id": f"k{mesh}{mesh}{mesh}",
        "label": f"{mesh}x{mesh}x{mesh}",
        "scheme": (
            f"MACDONALD {mesh} {mesh} {mesh} "
            f"{shift:.12g} {shift:.12g} {shift:.12g}"
        ),
    }


def replace_unique(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.I | re.M)
    if count != 1:
        raise ValueError(f"expected exactly one {label} while rendering X23b input")
    return updated


def render_input(series: str, method: str, system: str, mesh: int) -> str:
    if series not in SERIES or method not in METHODS:
        raise ValueError(f"invalid X23b convergence identity: {series}/{method}")
    metadata = system_metadata(system)
    geometry = pipeline.parse_cif(source_cif(system))
    text = pipeline.crystal_input(
        metadata, geometry, method, mesh_definition(mesh), "ENERGY"
    )
    project = f"x23b_kconv_{series}_{method}_{system}_k{mesh}{mesh}{mesh}".replace(
        "-", "_"
    )
    text = replace_unique(
        text,
        r"^(\s*PROJECT(?:_NAME)?\s+).*$",
        rf'\1"{project}"',
        "PROJECT",
    )
    if series == "independent_cellopt":
        text = replace_unique(
            text,
            r"^(\s*RUN_TYPE\s+)ENERGY\s*$",
            r"\1CELL_OPT",
            "RUN_TYPE ENERGY",
        )
        text += (
            "\n&MOTION\n"
            "  &CELL_OPT\n"
            "    OPTIMIZER CG\n"
            "    MAX_ITER 800\n"
            "    EXTERNAL_PRESSURE [bar] 0.0\n"
            "    KEEP_ANGLES F\n"
            "    &CG\n"
            "      &LINE_SEARCH\n"
            "        TYPE 2PNT\n"
            "      &END LINE_SEARCH\n"
            "    &END CG\n"
            "  &END CELL_OPT\n"
            "&END MOTION\n"
        )
    validate_input_contract(text, series, method, mesh)
    return text


def validate_input_contract(text: str, series: str, method: str, mesh: int) -> None:
    common.validate_method_input(text, method)
    shift = gamma_centered_shift(mesh)
    scheme = (
        rf"^\s*SCHEME\s+MACDONALD\s+{mesh}\s+{mesh}\s+{mesh}\s+"
        rf"{re.escape(f'{shift:.12g}')}\s+{re.escape(f'{shift:.12g}')}\s+"
        rf"{re.escape(f'{shift:.12g}')}\s*$"
    )
    required = {
        "native Bloch mesh": scheme,
        "SPGLIB symmetry": r"^\s*SYMMETRY\s+T\s*$",
        "reduced mesh": r"^\s*FULL_GRID\s+F\s*$",
        "SPGLIB backend": r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$",
        "SPGLIB reduction": r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$",
        "three-dimensional PBC": r"^\s*PERIODIC\s+XYZ\s*$",
    }
    for label, pattern in required.items():
        if len(re.findall(pattern, text, flags=re.I | re.M)) != 1:
            raise ValueError(f"input lacks exactly one {label}: {series}/{method}")
    forbidden = (
        "MULTIPLE_UNIT_CELL",
        "MULTIPLE_UNIT_CELL_REF",
        "BORN_VON_KARMAN",
        "SUPERCELL",
    )
    if any(token in text.upper() for token in forbidden):
        raise ValueError("Born--von Karman/supercell input is forbidden")
    if series == "fixed_experimental_sp":
        if not re.search(r"^\s*RUN_TYPE\s+ENERGY\s*$", text, re.I | re.M):
            raise ValueError("fixed experimental series is not a single point")
        if re.search(r"^\s*&MOTION\b", text, re.I | re.M):
            raise ValueError("fixed experimental series contains MOTION")
    elif series == "independent_cellopt":
        for pattern in (
            r"^\s*RUN_TYPE\s+CELL_OPT\s*$",
            r"^\s*&CELL_OPT\s*$",
            r"^\s*KEEP_ANGLES\s+F\s*$",
        ):
            if len(re.findall(pattern, text, re.I | re.M)) != 1:
                raise ValueError("independent cell series lacks its full-cell contract")
    else:
        raise ValueError(f"unknown convergence series: {series}")


def build_binding(
    *,
    cp2k: Path,
    cp2k_source: Path,
    save_tblite: Path,
    save_tblite_source: Path,
    campaign_manifest: Path,
) -> dict[str, object]:
    identity, cp2k_record, save_record, manifest_record = (
        common.validate_campaign_artifacts(
            cp2k=cp2k,
            cp2k_source=cp2k_source,
            save_tblite=save_tblite,
            save_tblite_source=save_tblite_source,
            campaign_manifest=campaign_manifest,
        )
    )
    binding = {
        "campaign_identity": identity,
        "campaign_manifest": manifest_record,
        "cp2k": {
            "path": str(cp2k.resolve(strict=True)),
            "sha256": str(cp2k_record["executable_sha256"]),
            "source_revision": str(cp2k_record["source_revision"]),
            "loaded_library_sha256": str(cp2k_record["loaded_library_sha256"]),
        },
        "save_tblite": {
            "path": str(save_tblite.resolve(strict=True)),
            "sha256": str(save_record["executable_sha256"]),
            "source_revision": str(save_record["source_revision"]),
            "static_library_sha256": str(save_record["static_library_sha256"]),
        },
    }
    require_post5582_campaign(binding)
    return binding


def job_relative_dir(series: str, method: str, system: str, mesh: str) -> Path:
    return Path("jobs") / series / method / system / mesh


def workflow_payload(output_root: Path, binding: Mapping[str, object]) -> dict[str, object]:
    system_ids = tuple(str(row["id"]) for row in systems())
    jobs: list[dict[str, object]] = []
    for series in SERIES:
        for method in METHODS:
            for system in system_ids:
                cif = source_cif(system)
                for mesh in MESH_NUMBERS:
                    mesh_id = f"k{mesh}{mesh}{mesh}"
                    relative_dir = job_relative_dir(series, method, system, mesh_id)
                    input_path = output_root / relative_dir / "input.inp"
                    jobs.append(
                        {
                            "series": series,
                            "method": method,
                            "system": system,
                            "mesh": mesh_id,
                            "mesh_number": mesh,
                            "run_dir": relative_dir.as_posix(),
                            "input": (relative_dir / "input.inp").as_posix(),
                            "input_sha256": sha256(input_path),
                            "output": (relative_dir / "cp2k.out").as_posix(),
                            "stamp": (
                                relative_dir / common.JOB_STAMP_NAME
                            ).as_posix(),
                            "source_cif": cif.relative_to(ROOT).as_posix(),
                            "source_cif_sha256": sha256(cif),
                        }
                    )
    return {
        "schema_version": WORKFLOW_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "status": "prepared_production_not_started",
        "benchmark": "X23b",
        "systems": list(system_ids),
        "methods": list(METHODS),
        "series": list(SERIES),
        "meshes": list(MESH_IDS),
        "native_kpoint_contract": {
            "scheme": "MACDONALD",
            "odd_mesh_shift": 0.0,
            "even_mesh_shift": "(N-1)/(2N)",
            "symmetry": True,
            "full_grid": False,
            "symmetry_backend": "SPGLIB",
            "symmetry_reduction_method": "SPGLIB",
            "born_von_karman_supercells": False,
        },
        "execution_contract": workflow_execution_contract(),
        "stopping": {
            "one_consecutive_step": True,
            "selected_value": "denser_mesh",
            "maximum_mesh": "k888",
            "fixed_experimental_sp": {
                "minimum_dense_mesh": "k222",
                "energy_abs_delta_kJmol_max": ENERGY_TOLERANCE_KJMOL,
            },
            "independent_cellopt": {
                "minimum_meshes": ["k111", "k222", "k333"],
                "minimum_dense_mesh": "k333",
                "max_length_relative_delta_percent": LENGTH_TOLERANCE_PERCENT,
                "volume_relative_delta_percent": VOLUME_TOLERANCE_PERCENT,
                "max_angle_abs_delta_deg": ANGLE_TOLERANCE_DEG,
                "all_criteria_required": True,
            },
        },
        "build": dict(binding),
        "generator": {
            "path": SCRIPT.relative_to(ROOT).as_posix(),
            "sha256": sha256(SCRIPT),
        },
        "jobs": jobs,
    }


def prepare(args: argparse.Namespace) -> tuple[Path, str]:
    output_root = args.output_root.resolve()
    if output_root.exists():
        raise ValueError(f"refusing to prepare over existing workflow root: {output_root}")
    binding = build_binding(
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    temporary = output_root.with_name(f".{output_root.name}.prepare.{os.getpid()}")
    if temporary.exists():
        shutil.rmtree(temporary)
    try:
        for series in SERIES:
            for method in METHODS:
                for system in (str(row["id"]) for row in systems()):
                    for mesh in MESH_NUMBERS:
                        mesh_id = f"k{mesh}{mesh}{mesh}"
                        run_dir = temporary / job_relative_dir(
                            series, method, system, mesh_id
                        )
                        run_dir.mkdir(parents=True, exist_ok=True)
                        (run_dir / "input.inp").write_text(
                            render_input(series, method, system, mesh)
                        )
        payload = workflow_payload(temporary, binding)
        # Paths are relative, so the exact same manifest remains valid after
        # the atomically prepared directory is renamed into place.
        manifest = temporary / "workflow_manifest.json"
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        output_root.parent.mkdir(parents=True, exist_ok=True)
        os.replace(temporary, output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    manifest = output_root / "workflow_manifest.json"
    return manifest, sha256(manifest)


def expected_job_keys() -> list[tuple[str, str, str, str]]:
    return [
        (series, method, str(system["id"]), mesh)
        for series in SERIES
        for method in METHODS
        for system in systems()
        for mesh in MESH_IDS
    ]


def load_workflow(path: Path, expected_sha256: str) -> dict[str, Any]:
    path = path.resolve(strict=True)
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("workflow manifest requires an exact SHA256 pin")
    if sha256(path) != expected_sha256:
        raise ValueError("workflow manifest SHA256 pin mismatch")
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError("workflow manifest is not a JSON object")
    expected_top = {
        "schema_version": WORKFLOW_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "benchmark": "X23b",
        "systems": [str(row["id"]) for row in systems()],
        "methods": list(METHODS),
        "series": list(SERIES),
        "meshes": list(MESH_IDS),
    }
    for key, expected in expected_top.items():
        if payload.get(key) != expected:
            raise ValueError(f"workflow manifest {key} differs from the protocol")
    stopping = payload.get("stopping")
    if not isinstance(stopping, dict) or stopping != workflow_stopping_contract():
        raise ValueError("workflow stopping contract differs from the protocol")
    generator = payload.get("generator")
    if (
        not isinstance(generator, dict)
        or generator.get("path") != SCRIPT.relative_to(ROOT).as_posix()
        or generator.get("sha256") != sha256(SCRIPT)
    ):
        raise ValueError("workflow generator fingerprint differs")
    build = payload.get("build")
    if not isinstance(build, dict):
        raise ValueError("workflow build binding is missing")
    identity = build.get("campaign_identity")
    if not isinstance(identity, dict):
        raise ValueError("workflow campaign identity is missing")
    common.validate_campaign_identity(identity)
    manifest_record = build.get("campaign_manifest")
    if not isinstance(manifest_record, dict):
        raise ValueError("workflow campaign-manifest record is missing")
    require_post5582_campaign(build)
    build_manifest = Path(str(manifest_record.get("path", ""))).resolve(strict=True)
    if sha256(build_manifest) != manifest_record.get("file_sha256"):
        raise ValueError("workflow campaign manifest has changed")
    if payload.get("execution_contract") != workflow_execution_contract():
        raise ValueError("workflow execution contract differs from the protocol")

    jobs = payload.get("jobs")
    keys = expected_job_keys()
    if not isinstance(jobs, list) or len(jobs) != len(keys):
        raise ValueError("workflow job matrix is incomplete or duplicated")
    root = path.parent
    for record, expected_key in zip(jobs, keys, strict=True):
        if not isinstance(record, dict):
            raise ValueError("workflow job matrix contains a non-object row")
        observed_key = tuple(str(record.get(key, "")) for key in (
            "series", "method", "system", "mesh"
        ))
        if observed_key != expected_key:
            raise ValueError("workflow job matrix is reordered, duplicated, or incomplete")
        series, method, system, mesh_id = expected_key
        mesh = int(record.get("mesh_number", 0))
        if mesh_id != f"k{mesh}{mesh}{mesh}" or mesh not in MESH_NUMBERS:
            raise ValueError(f"workflow mesh identity differs for {observed_key}")
        relative = job_relative_dir(series, method, system, mesh_id)
        expected_paths = {
            "run_dir": relative.as_posix(),
            "input": (relative / "input.inp").as_posix(),
            "output": (relative / "cp2k.out").as_posix(),
            "stamp": (relative / common.JOB_STAMP_NAME).as_posix(),
            "source_cif": source_cif(system).relative_to(ROOT).as_posix(),
        }
        for field, expected in expected_paths.items():
            if record.get(field) != expected:
                raise ValueError(f"workflow {field} differs for {observed_key}")
        input_path = root / str(record["input"])
        cif = ROOT / str(record["source_cif"])
        if sha256(input_path) != record.get("input_sha256"):
            raise ValueError(f"workflow input fingerprint differs for {observed_key}")
        if sha256(cif) != record.get("source_cif_sha256"):
            raise ValueError(f"workflow CIF fingerprint differs for {observed_key}")
        expected_text = render_input(series, method, system, mesh)
        if input_path.read_text() != expected_text:
            raise ValueError(f"workflow input is not reproducible for {observed_key}")
    return payload


def workflow_stopping_contract() -> dict[str, object]:
    return {
        "one_consecutive_step": True,
        "selected_value": "denser_mesh",
        "maximum_mesh": "k888",
        "fixed_experimental_sp": {
            "minimum_dense_mesh": "k222",
            "energy_abs_delta_kJmol_max": ENERGY_TOLERANCE_KJMOL,
        },
        "independent_cellopt": {
            "minimum_meshes": ["k111", "k222", "k333"],
            "minimum_dense_mesh": "k333",
            "max_length_relative_delta_percent": LENGTH_TOLERANCE_PERCENT,
            "volume_relative_delta_percent": VOLUME_TOLERANCE_PERCENT,
            "max_angle_abs_delta_deg": ANGLE_TOLERANCE_DEG,
            "all_criteria_required": True,
        },
    }


def workflow_execution_contract() -> dict[str, object]:
    return {
        "launcher": "direct frozen CP2K executable",
        "working_directory": "job run_dir",
        "arguments": ["-i", "input.inp", "-o", "cp2k.out"],
        "thread_count": "positive --threads-per-job recorded in every job stamp",
    }


def validate_runtime_binding(
    args: argparse.Namespace, workflow: Mapping[str, object]
) -> dict[str, object]:
    binding = build_binding(
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    if binding != workflow.get("build"):
        raise ValueError("runtime build differs from the frozen workflow build")
    return binding


def job_map(workflow: Mapping[str, object]) -> dict[tuple[str, str, str, str], dict[str, object]]:
    jobs = workflow["jobs"]
    assert isinstance(jobs, list)
    return {
        (
            str(row["series"]),
            str(row["method"]),
            str(row["system"]),
            str(row["mesh"]),
        ): row
        for row in jobs
    }


def resolved_job_paths(
    workflow_manifest: Path, record: Mapping[str, object]
) -> tuple[Path, Path, Path, Path]:
    root = workflow_manifest.resolve().parent
    return (
        root / str(record["run_dir"]),
        root / str(record["input"]),
        root / str(record["output"]),
        root / str(record["stamp"]),
    )


def protocol_identity(
    workflow_sha256: str,
    workflow: Mapping[str, object],
    record: Mapping[str, object],
) -> dict[str, object]:
    build = workflow["build"]
    assert isinstance(build, dict)
    identity = build["campaign_identity"]
    assert isinstance(identity, dict)
    return {
        "protocol_id": PROTOCOL_ID,
        "workflow_manifest_sha256": workflow_sha256,
        "campaign_fingerprint_sha256": identity["fingerprint_sha256"],
        "series": record["series"],
        "mesh": record["mesh"],
        "native_bloch": True,
        "spglib_reduced": True,
        "full_grid": False,
    }


def source_artifacts(
    workflow_manifest: Path, workflow: Mapping[str, object], record: Mapping[str, object]
) -> dict[str, Path]:
    build = workflow["build"]
    assert isinstance(build, dict)
    manifest = build["campaign_manifest"]
    assert isinstance(manifest, dict)
    return {
        "campaign_manifest": Path(str(manifest["path"])),
        "experimental_cif": ROOT / str(record["source_cif"]),
        "workflow_manifest": workflow_manifest,
    }


def phase_name(record: Mapping[str, object]) -> str:
    return f"x23b_kconv_{record['series']}_{record['mesh']}"


def job_command(cp2k: Path) -> list[str]:
    return [str(cp2k.resolve(strict=True)), "-i", "input.inp", "-o", "cp2k.out"]


def expected_stamp(
    workflow_manifest: Path,
    workflow_sha256: str,
    workflow: Mapping[str, object],
    record: Mapping[str, object],
    cp2k: Path,
) -> dict[str, object]:
    _, input_path, _, _ = resolved_job_paths(workflow_manifest, record)
    build = workflow["build"]
    assert isinstance(build, dict)
    identity = build["campaign_identity"]
    assert isinstance(identity, dict)
    return common.job_identity(
        input_path,
        cp2k,
        str(record["method"]),
        phase_name(record),
        campaign_identity=identity if record["method"] == "GXTB" else None,
        protocol_identity=protocol_identity(workflow_sha256, workflow, record),
        source_artifacts=source_artifacts(workflow_manifest, workflow, record),
    )


def validate_completed_stamp(
    workflow_manifest: Path,
    workflow_sha256: str,
    workflow: Mapping[str, object],
    record: Mapping[str, object],
    cp2k: Path,
) -> dict[str, object]:
    _, _, output, stamp_path = resolved_job_paths(workflow_manifest, record)
    if not output.is_file() or not stamp_path.is_file():
        raise ValueError(
            "partial output/stamp pair for "
            f"{record['method']}/{record['system']}/{record['mesh']}"
        )
    stamp = json.loads(stamp_path.read_text())
    expected = expected_stamp(
        workflow_manifest, workflow_sha256, workflow, record, cp2k
    )
    for key, value in expected.items():
        if stamp.get(key) != value:
            raise ValueError(
                f"stale {key} fingerprint for "
                f"{record['method']}/{record['system']}/{record['mesh']}"
            )
    details = stamp.get("details")
    if (
        not str(stamp.get("status", "")).startswith("converged")
        or not isinstance(details, dict)
        or details.get("returncode") != 0
        or details.get("output") != str(output.resolve())
        or details.get("output_sha256") != sha256(output)
    ):
        raise ValueError(
            "non-converged or stale output stamp for "
            f"{record['method']}/{record['system']}/{record['mesh']}"
        )
    threads = details.get("threads")
    if (
        details.get("command") != job_command(cp2k)
        or isinstance(threads, bool)
        or not isinstance(threads, int)
        or threads <= 0
    ):
        raise ValueError(
            "stale command/threading record for "
            f"{record['method']}/{record['system']}/{record['mesh']}"
        )
    return stamp


def output_is_fatal(text: str) -> bool:
    return bool(
        "PROGRAM ENDED" not in text
        or "ABORT" in text
        or re.search(r"SCF.*NOT|NOT.*SCF|DID NOT CONVERGE|convergence failure", text, re.I)
    )


def parse_energy(text: str) -> float:
    values = [
        float(match.group(1))
        for match in re.finditer(
            r"ENERGY\|\s+Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$", text, re.M
        )
    ]
    if not values:
        raise ValueError("CP2K output lacks a total energy")
    return values[-1]


def vector_length(vector: Iterable[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def vector_angle(left: list[float], right: list[float]) -> float:
    cosine = sum(a * b for a, b in zip(left, right, strict=True)) / (
        vector_length(left) * vector_length(right)
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def determinant(cell: list[list[float]]) -> float:
    a, b, c = cell
    return abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def parse_cell(text: str) -> dict[str, object]:
    vectors: dict[str, list[list[float]]] = {axis: [] for axis in "abc"}
    pattern = re.compile(
        r"CELL\|\s+Vector\s+([abc])\s+\[angstrom\]:\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
        re.I,
    )
    for match in pattern.finditer(text):
        vectors[match.group(1).lower()].append(
            [float(match.group(i)) for i in range(2, 5)]
        )
    if any(not vectors[axis] for axis in "abc"):
        raise ValueError("CP2K output lacks complete cell vectors")
    cell = [vectors[axis][-1] for axis in "abc"]
    lengths = [vector_length(vector) for vector in cell]
    angles = [
        vector_angle(cell[1], cell[2]),
        vector_angle(cell[0], cell[2]),
        vector_angle(cell[0], cell[1]),
    ]
    return {
        "vectors_A": cell,
        "lengths_A": lengths,
        "angles_deg": angles,
        "volume_A3": determinant(cell),
    }


def parse_observation(
    workflow_manifest: Path,
    workflow_sha256: str,
    workflow: Mapping[str, object],
    record: Mapping[str, object],
    cp2k: Path,
) -> dict[str, object]:
    _, _, output, stamp = resolved_job_paths(workflow_manifest, record)
    if output.exists() != stamp.exists():
        raise ValueError(
            "partial output/stamp pair for "
            f"{record['method']}/{record['system']}/{record['mesh']}"
        )
    if not output.exists():
        raise FileNotFoundError(output)
    validate_completed_stamp(
        workflow_manifest, workflow_sha256, workflow, record, cp2k
    )
    text = output.read_text(errors="ignore")
    if output_is_fatal(text):
        raise ValueError(f"fatal/incomplete CP2K output: {output}")
    series = str(record["series"])
    if series == "independent_cellopt":
        if (
            "GEOMETRY OPTIMIZATION COMPLETED" not in text
            and "CELL OPTIMIZATION COMPLETED" not in text
        ) or "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text:
            raise ValueError(f"cell optimization did not converge: {output}")
    observation: dict[str, object] = {
        "energy_hartree": parse_energy(text),
        "output": str(output.resolve()),
        "output_sha256": sha256(output),
        "stamp": str(stamp.resolve()),
        "stamp_sha256": sha256(stamp),
    }
    if series == "independent_cellopt":
        observation.update(parse_cell(text))
    return observation


def relative_percent(dense: float, coarse: float) -> float:
    if coarse == 0.0:
        raise ValueError("zero coarse value in relative convergence criterion")
    return 100.0 * abs(dense - coarse) / abs(coarse)


def within_threshold(value: float, threshold: float) -> bool:
    """Apply an inclusive threshold without rejecting round-off at its boundary."""
    return value <= threshold or math.isclose(
        value,
        threshold,
        rel_tol=0.0,
        abs_tol=COMPARISON_ABSOLUTE_TOLERANCE,
    )


def delta_row(
    series: str,
    method: str,
    system: str,
    coarse_mesh: int,
    dense_mesh: int,
    coarse: Mapping[str, object],
    dense: Mapping[str, object],
) -> dict[str, object]:
    row: dict[str, object] = {field: "" for field in CSV_FIELDS}
    row.update(
        {
            "series": series,
            "method": method,
            "method_label": method_label(method),
            "system": system,
            "coarse_mesh": f"k{coarse_mesh}{coarse_mesh}{coarse_mesh}",
            "dense_mesh": f"k{dense_mesh}{dense_mesh}{dense_mesh}",
            "eligible_for_stopping": dense_mesh >= MINIMUM_DENSE_MESH[series],
            "selected_step": False,
            "coarse_output_sha256": coarse["output_sha256"],
            "dense_output_sha256": dense["output_sha256"],
        }
    )
    if series == "fixed_experimental_sp":
        molecules = int(system_metadata(system)["molecules_per_cell"])
        coarse_energy = float(coarse["energy_hartree"]) / molecules
        dense_energy = float(dense["energy_hartree"]) / molecules
        delta = (dense_energy - coarse_energy) * HARTREE_TO_KJMOL
        passed = within_threshold(abs(delta), ENERGY_TOLERANCE_KJMOL)
        row.update(
            {
                "coarse_energy_hartree_per_molecule": coarse_energy,
                "dense_energy_hartree_per_molecule": dense_energy,
                "energy_delta_kJmol": delta,
                "energy_abs_delta_kJmol": abs(delta),
                "energy_pass": passed,
                "all_required_criteria_pass": passed,
            }
        )
        return row

    coarse_lengths = [float(value) for value in coarse["lengths_A"]]  # type: ignore[index]
    dense_lengths = [float(value) for value in dense["lengths_A"]]  # type: ignore[index]
    coarse_angles = [float(value) for value in coarse["angles_deg"]]  # type: ignore[index]
    dense_angles = [float(value) for value in dense["angles_deg"]]  # type: ignore[index]
    length_delta = [
        relative_percent(dense_value, coarse_value)
        for coarse_value, dense_value in zip(
            coarse_lengths, dense_lengths, strict=True
        )
    ]
    angle_delta = [
        abs(dense_value - coarse_value)
        for coarse_value, dense_value in zip(coarse_angles, dense_angles, strict=True)
    ]
    volume_delta = relative_percent(
        float(dense["volume_A3"]), float(coarse["volume_A3"])
    )
    length_pass = within_threshold(max(length_delta), LENGTH_TOLERANCE_PERCENT)
    volume_pass = within_threshold(volume_delta, VOLUME_TOLERANCE_PERCENT)
    angle_pass = within_threshold(max(angle_delta), ANGLE_TOLERANCE_DEG)
    row.update(
        {
            "coarse_a_A": coarse_lengths[0],
            "coarse_b_A": coarse_lengths[1],
            "coarse_c_A": coarse_lengths[2],
            "dense_a_A": dense_lengths[0],
            "dense_b_A": dense_lengths[1],
            "dense_c_A": dense_lengths[2],
            "a_relative_delta_percent": length_delta[0],
            "b_relative_delta_percent": length_delta[1],
            "c_relative_delta_percent": length_delta[2],
            "max_length_relative_delta_percent": max(length_delta),
            "length_pass": length_pass,
            "coarse_alpha_deg": coarse_angles[0],
            "coarse_beta_deg": coarse_angles[1],
            "coarse_gamma_deg": coarse_angles[2],
            "dense_alpha_deg": dense_angles[0],
            "dense_beta_deg": dense_angles[1],
            "dense_gamma_deg": dense_angles[2],
            "alpha_abs_delta_deg": angle_delta[0],
            "beta_abs_delta_deg": angle_delta[1],
            "gamma_abs_delta_deg": angle_delta[2],
            "max_angle_abs_delta_deg": max(angle_delta),
            "angle_pass": angle_pass,
            "coarse_volume_A3": float(coarse["volume_A3"]),
            "dense_volume_A3": float(dense["volume_A3"]),
            "volume_relative_delta_percent": volume_delta,
            "volume_pass": volume_pass,
            "all_required_criteria_pass": (
                length_pass and volume_pass and angle_pass
            ),
        }
    )
    return row


def assess_series(
    series: str,
    method: str,
    system: str,
    observations: Mapping[int, Mapping[str, object]],
) -> dict[str, object]:
    present = sorted(observations)
    if present and present != list(range(1, present[-1] + 1)):
        raise ValueError(f"non-contiguous mesh results for {series}/{method}/{system}")
    rows: list[dict[str, object]] = []
    selected: int | None = None
    for dense_mesh in range(2, (present[-1] if present else 0) + 1):
        row = delta_row(
            series,
            method,
            system,
            dense_mesh - 1,
            dense_mesh,
            observations[dense_mesh - 1],
            observations[dense_mesh],
        )
        eligible = bool(row["eligible_for_stopping"])
        passed = bool(row["all_required_criteria_pass"])
        if selected is None and eligible and passed:
            selected = dense_mesh
            row["selected_step"] = True
        rows.append(row)
    if selected is not None and present[-1] > selected:
        raise ValueError(
            f"completed mesh after first converged step for {series}/{method}/{system}"
        )
    minimum = MINIMUM_DENSE_MESH[series]
    if selected is not None:
        status = "converged"
    elif present and present[-1] == MESH_NUMBERS[-1]:
        status = "maximum_mesh_unconverged"
    else:
        status = "pending"
    required: list[int] = []
    if status == "pending":
        required_initial = list(range(1, minimum + 1))
        required = [mesh for mesh in required_initial if mesh not in observations]
        if not required and present:
            required = [present[-1] + 1]
    return {
        "status": status,
        "selected_mesh_number": selected,
        "selected_mesh": f"k{selected}{selected}{selected}" if selected else None,
        "selected_observation": observations.get(selected) if selected else None,
        "rows": rows,
        "required_meshes": required,
    }


def existing_observations(
    workflow_manifest: Path,
    workflow_sha256: str,
    workflow: Mapping[str, object],
    cp2k: Path,
) -> dict[tuple[str, str, str], dict[int, dict[str, object]]]:
    result: dict[tuple[str, str, str], dict[int, dict[str, object]]] = {}
    for key, record in job_map(workflow).items():
        series, method, system, _ = key
        _, _, output, stamp = resolved_job_paths(workflow_manifest, record)
        if output.exists() != stamp.exists():
            raise ValueError(f"partial output/stamp pair: {output}")
        if not output.exists():
            continue
        mesh = int(record["mesh_number"])
        result.setdefault((series, method, system), {})[mesh] = parse_observation(
            workflow_manifest,
            workflow_sha256,
            workflow,
            record,
            cp2k,
        )
    return result


def run_one(
    workflow_manifest: Path,
    workflow_sha256: str,
    workflow: Mapping[str, object],
    record: Mapping[str, object],
    cp2k: Path,
    threads: int,
) -> tuple[str, int, str]:
    run_dir, input_path, output, stamp_path = resolved_job_paths(
        workflow_manifest, record
    )
    identity = workflow["build"]["campaign_identity"]  # type: ignore[index]
    assert isinstance(identity, dict)
    with input_path.open() as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return str(input_path), common.BUSY_RETURN_CODE, "BUSY"
        if output.exists() or stamp_path.exists():
            try:
                parse_observation(
                    workflow_manifest,
                    workflow_sha256,
                    workflow,
                    record,
                    cp2k,
                )
            except (OSError, ValueError) as error:
                return str(input_path), 1, f"STALE:{error}"
            return str(input_path), 0, "SKIP_CONVERGED"
        command = job_command(cp2k)
        process = subprocess.run(
            command,
            cwd=run_dir,
            env=common.thread_environment(threads),
            check=False,
        )
        details: dict[str, object] = {
            "returncode": process.returncode,
            "output": str(output.resolve()),
            "command": command,
            "threads": threads,
        }
        status = "failed"
        if output.is_file():
            details["output_sha256"] = sha256(output)
        if process.returncode == 0 and output.is_file():
            try:
                text = output.read_text(errors="ignore")
                if output_is_fatal(text):
                    raise ValueError("fatal/incomplete output")
                parse_energy(text)
                if record["series"] == "independent_cellopt":
                    if (
                        "GEOMETRY OPTIMIZATION COMPLETED" not in text
                        and "CELL OPTIMIZATION COMPLETED" not in text
                    ) or "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text:
                        raise ValueError("cell optimization incomplete")
                    parse_cell(text)
                status = "converged"
            except ValueError as error:
                details["parse_error"] = str(error)
        common.write_job_stamp(
            run_dir,
            input_path,
            cp2k,
            str(record["method"]),
            phase_name(record),
            status,
            details=details,
            campaign_identity=identity if record["method"] == "GXTB" else None,
            protocol_identity=protocol_identity(
                workflow_sha256, workflow, record
            ),
            source_artifacts=source_artifacts(
                workflow_manifest, workflow, record
            ),
        )
        return str(input_path), process.returncode if status == "converged" else 1, status.upper()


def selected_tracks(args: argparse.Namespace) -> list[tuple[str, str, str]]:
    selected_series = tuple(args.series) if args.series else SERIES
    selected_methods = tuple(args.method) if args.method else METHODS
    selected_systems = (
        tuple(args.system)
        if args.system
        else tuple(str(row["id"]) for row in systems())
    )
    if len(set(selected_series)) != len(selected_series):
        raise ValueError("duplicate --series selection")
    if len(set(selected_methods)) != len(selected_methods):
        raise ValueError("duplicate --method selection")
    if len(set(selected_systems)) != len(selected_systems):
        raise ValueError("duplicate --system selection")
    unknown = set(selected_systems) - {str(row["id"]) for row in systems()}
    if unknown:
        raise ValueError("unknown X23b systems: " + ", ".join(sorted(unknown)))
    return [
        (series, method, system)
        for series in selected_series
        for method in selected_methods
        for system in selected_systems
    ]


def run(args: argparse.Namespace) -> None:
    workflow_manifest = args.workflow_manifest.resolve(strict=True)
    workflow = load_workflow(
        workflow_manifest, args.workflow_manifest_sha256
    )
    validate_runtime_binding(args, workflow)
    cp2k = args.cp2k.resolve(strict=True)
    records = job_map(workflow)
    tracks = selected_tracks(args)
    while True:
        observations = existing_observations(
            workflow_manifest,
            args.workflow_manifest_sha256,
            workflow,
            cp2k,
        )
        pending: list[dict[str, object]] = []
        unconverged: list[tuple[str, str, str]] = []
        for series, method, system in tracks:
            assessment = assess_series(
                series,
                method,
                system,
                observations.get((series, method, system), {}),
            )
            if assessment["status"] == "maximum_mesh_unconverged":
                unconverged.append((series, method, system))
            for mesh in assessment["required_meshes"]:  # type: ignore[union-attr]
                mesh_id = f"k{mesh}{mesh}{mesh}"
                pending.append(records[(series, method, system, mesh_id)])
        if unconverged:
            raise ValueError(
                "k888 reached without convergence for: "
                + ", ".join("/".join(key) for key in unconverged)
            )
        if not pending:
            return
        failures: list[str] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(
                    run_one,
                    workflow_manifest,
                    args.workflow_manifest_sha256,
                    workflow,
                    record,
                    cp2k,
                    args.threads_per_job,
                ): record
                for record in pending
            }
            for future in concurrent.futures.as_completed(futures):
                record = futures[future]
                _, returncode, action = future.result()
                label = f"{record['series']}/{record['method']}/{record['system']}/{record['mesh']}"
                print(f"{action:18s} {label}", flush=True)
                if returncode != 0:
                    failures.append(label)
        if failures:
            raise ValueError(f"{len(failures)} X23b convergence jobs failed")


def csv_text(rows: Iterable[Mapping[str, object]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})
    return stream.getvalue()


def tex_token(value: str) -> str:
    return "".join(
        piece[0].upper() + piece[1:]
        for piece in re.findall(r"[A-Za-z0-9]+", value)
    )


def tex_text(selections: Iterable[Mapping[str, object]]) -> str:
    lines = [
        "% AUTO-GENERATED by x23b_adaptive_kpoint_convergence.py",
        "% Do not edit numerical values by hand.",
        f"\\newcommand{{\\GXTBXKConvEnergyTolerance}}{{{ENERGY_TOLERANCE_KJMOL:.3f}}}",
        f"\\newcommand{{\\GXTBXKConvLengthTolerance}}{{{LENGTH_TOLERANCE_PERCENT:.3f}}}",
        f"\\newcommand{{\\GXTBXKConvVolumeTolerance}}{{{VOLUME_TOLERANCE_PERCENT:.3f}}}",
        f"\\newcommand{{\\GXTBXKConvAngleTolerance}}{{{ANGLE_TOLERANCE_DEG:.3f}}}",
    ]
    names: set[str] = set()
    for selection in selections:
        prefix = (
            "GXTBXKConv"
            + tex_token(str(selection["series"]))
            + tex_token(str(selection["method"]))
            + tex_token(str(selection["system"]))
        )
        values: dict[str, object] = {
            "Mesh": selection["selected_mesh_number"],
        }
        delta = selection["selected_delta"]
        assert isinstance(delta, dict)
        if selection["series"] == "fixed_experimental_sp":
            values["EnergyAbsDelta"] = delta["energy_abs_delta_kJmol"]
        else:
            values.update(
                {
                    "MaxLengthRelativeDelta": delta[
                        "max_length_relative_delta_percent"
                    ],
                    "VolumeRelativeDelta": delta[
                        "volume_relative_delta_percent"
                    ],
                    "MaxAngleAbsDelta": delta["max_angle_abs_delta_deg"],
                }
            )
        for suffix, raw in values.items():
            name = prefix + suffix
            if name in names:
                raise ValueError(f"duplicate TeX macro: {name}")
            names.add(name)
            value = str(raw) if suffix == "Mesh" else f"{float(raw):.9f}"
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    return "\n".join(lines) + "\n"


def finalize(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    workflow_manifest = args.workflow_manifest.resolve(strict=True)
    target = args.output_dir.resolve()
    outputs = tuple(target / f"{OUTPUT_STEM}.{suffix}" for suffix in ("csv", "json", "tex"))
    for output in outputs:
        output.unlink(missing_ok=True)
    temporaries = tuple(
        output.with_name(f".{output.name}.{os.getpid()}.tmp") for output in outputs
    )
    for temporary in temporaries:
        temporary.unlink(missing_ok=True)
    try:
        workflow = load_workflow(
            workflow_manifest, args.workflow_manifest_sha256
        )
        validate_runtime_binding(args, workflow)
        cp2k = args.cp2k.resolve(strict=True)
        observations = existing_observations(
            workflow_manifest,
            args.workflow_manifest_sha256,
            workflow,
            cp2k,
        )
        rows: list[dict[str, object]] = []
        selections: list[dict[str, object]] = []
        lineage: list[dict[str, object]] = []
        for series in SERIES:
            for method in METHODS:
                for system in (str(row["id"]) for row in systems()):
                    assessment = assess_series(
                        series,
                        method,
                        system,
                        observations.get((series, method, system), {}),
                    )
                    if assessment["status"] != "converged":
                        raise ValueError(
                            f"incomplete convergence track: {series}/{method}/{system}"
                        )
                    track_rows = assessment["rows"]
                    assert isinstance(track_rows, list)
                    rows.extend(track_rows)
                    selected_rows = [row for row in track_rows if row["selected_step"]]
                    if len(selected_rows) != 1:
                        raise ValueError(
                            f"track lacks exactly one selected step: {series}/{method}/{system}"
                        )
                    selected = int(assessment["selected_mesh_number"])
                    observation = assessment["selected_observation"]
                    assert isinstance(observation, dict)
                    selections.append(
                        {
                            "series": series,
                            "method": method,
                            "method_label": method_label(method),
                            "system": system,
                            "selected_mesh": assessment["selected_mesh"],
                            "selected_mesh_number": selected,
                            "selected_value": observation,
                            "selected_delta": selected_rows[0],
                        }
                    )
                    for mesh in range(1, selected + 1):
                        record = job_map(workflow)[
                            (series, method, system, f"k{mesh}{mesh}{mesh}")
                        ]
                        _, input_path, output, stamp = resolved_job_paths(
                            workflow_manifest, record
                        )
                        lineage.append(
                            {
                                "series": series,
                                "method": method,
                                "system": system,
                                "mesh": record["mesh"],
                                "input_sha256": sha256(input_path),
                                "output_sha256": sha256(output),
                                "stamp_sha256": sha256(stamp),
                            }
                        )
        if len(selections) != len(SERIES) * len(METHODS) * 23:
            raise ValueError("final X23b convergence selection coverage is not 138/138")
        csv_body = csv_text(rows)
        tex_body = tex_text(selections)
        payload = {
            "schema_version": OUTPUT_SCHEMA,
            "status": "publication_ready",
            "benchmark": "X23b",
            "protocol_id": PROTOCOL_ID,
            "coverage": {
                "systems": 23,
                "methods": list(METHODS),
                "series": list(SERIES),
                "tracks": len(selections),
                "exact_complete_coverage": True,
            },
            "stopping": workflow_stopping_contract(),
            "workflow_manifest": {
                "path": str(workflow_manifest),
                "sha256": args.workflow_manifest_sha256,
            },
            "build": workflow["build"],
            "selections": selections,
            "raw_consecutive_deltas": rows,
            "used_job_lineage": lineage,
            "generated_outputs": {
                "csv_sha256": hashlib.sha256(csv_body.encode()).hexdigest(),
                "tex_sha256": hashlib.sha256(tex_body.encode()).hexdigest(),
            },
        }
        target.mkdir(parents=True, exist_ok=True)
        temporaries[0].write_text(csv_body)
        temporaries[1].write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporaries[2].write_text(tex_body)
        for temporary, output in zip(temporaries, outputs, strict=True):
            os.replace(temporary, output)
    except BaseException:
        for path in (*temporaries, *outputs):
            path.unlink(missing_ok=True)
        raise
    return outputs  # type: ignore[return-value]


def add_build_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--campaign-manifest", type=Path, required=True)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    add_build_arguments(prepare_parser)

    run_parser = commands.add_parser("run")
    run_parser.add_argument("--workflow-manifest", type=Path, required=True)
    run_parser.add_argument("--workflow-manifest-sha256", required=True)
    add_build_arguments(run_parser)
    run_parser.add_argument("--jobs", type=int, default=1)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    run_parser.add_argument("--series", action="append", choices=SERIES)
    run_parser.add_argument("--method", action="append", choices=METHODS)
    run_parser.add_argument("--system", action="append")

    finalize_parser = commands.add_parser("finalize")
    finalize_parser.add_argument("--workflow-manifest", type=Path, required=True)
    finalize_parser.add_argument("--workflow-manifest-sha256", required=True)
    finalize_parser.add_argument("--output-dir", type=Path, required=True)
    add_build_arguments(finalize_parser)
    return result


def main() -> int:
    args = parser().parse_args()
    for field in (
        "cp2k",
        "cp2k_source",
        "save_tblite",
        "save_tblite_source",
        "campaign_manifest",
    ):
        setattr(args, field, getattr(args, field).expanduser().resolve())
    try:
        if args.command == "prepare":
            manifest, digest = prepare(args)
            print(manifest)
            print(digest)
        elif args.command == "run":
            if args.jobs <= 0 or args.threads_per_job <= 0:
                raise ValueError("--jobs and --threads-per-job must be positive")
            run(args)
        elif args.command == "finalize":
            print("\n".join(str(path) for path in finalize(args)))
        else:  # pragma: no cover
            raise AssertionError(args.command)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser().error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
