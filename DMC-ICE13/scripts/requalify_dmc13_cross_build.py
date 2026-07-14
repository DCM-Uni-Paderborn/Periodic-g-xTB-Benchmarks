#!/usr/bin/env python3
"""Requalify frozen DMC13 g-xTB energies against a new CP2K build."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import queue
import re
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any


HARTREE_TO_KJMOL = 2625.499638
TOTAL_ENERGY_TOLERANCE_HARTREE = 1.0e-10
RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O = 1.0e-3
PROTOCOL_ID = "dmc13-post-5582-cross-build-requalification-v1"
REQUIRED_CP2K_ANCESTOR = "c92cc08b45378b85150447011b5a4bb552f5b797"
SENTINEL_SCOPE = "sentinel"
FULL_PUBLICATION_SCOPE = "full-publication-matrix"
FIXED_COMPARISON_MESH = "k333"


def load_runner(repository: Path):
    path = repository / "scripts" / "run_dmc13_kpoint_jobs.py"
    spec = importlib.util.spec_from_file_location("dmc13_cross_build_runner", path)
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load DMC13 runner: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def finite_float(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def parse_selection(value: str, runner: Any) -> tuple[str, str]:
    mesh, separator, phase = value.partition(":")
    if not separator or mesh not in runner.SUPPORTED_MESHES or phase not in runner.PHASES:
        raise argparse.ArgumentTypeError(
            "selection must be MESH:PHASE with a supported DMC13 mesh and phase"
        )
    return mesh, phase


def validate_matrix(selections: list[tuple[str, str]]) -> None:
    if not selections:
        raise ValueError("the cross-build matrix must not be empty")
    if len(set(selections)) != len(selections):
        raise ValueError("the cross-build matrix contains duplicate selections")
    selected = set(selections)
    for mesh, phase in selections:
        if phase != "Ih" and (mesh, "Ih") not in selected:
            raise ValueError(f"{mesh}/{phase} lacks its same-mesh Ih reference")


def read_pinned_json(path: Path, expected_sha256: str, runner: Any) -> dict[str, Any]:
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError(f"invalid SHA256 pin for {path}")
    if runner.sha256(path) != expected_sha256:
        raise ValueError(f"SHA256 pin mismatch for {path}")
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read JSON artifact {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"JSON artifact is not an object: {path}")
    return payload


def validate_candidate_manifest(
    payload: dict[str, Any], identity: Any, args: argparse.Namespace, runner: Any
) -> None:
    manifest_path = args.candidate_build_manifest
    if payload.get("campaign_state") != "production_ready":
        raise ValueError("candidate build manifest is not production_ready")
    campaign_id = payload.get("campaign_id")
    if not isinstance(campaign_id, str) or not campaign_id:
        raise ValueError("candidate build manifest lacks /campaign_id")
    cp2k = payload.get("cp2k")
    save_tblite = payload.get("save_tblite")
    qualification = payload.get("qualification")
    if not all(isinstance(value, dict) for value in (cp2k, save_tblite, qualification)):
        raise ValueError(
            "candidate build manifest lacks /cp2k, /save_tblite, or /qualification"
        )

    def resolved_artifact(
        block: dict[str, Any], key: str, label: str
    ) -> Path:
        value = block.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(
                f"candidate build manifest lacks artifact path {label}/{key}"
            )
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = manifest_path.parent / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"candidate build manifest artifact is missing: {path}")
        expected = block.get(f"{key}_sha256")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError(
                f"candidate build manifest lacks valid {label}/{key}_sha256"
            )
        if runner.sha256(path) != expected:
            raise ValueError(
                f"candidate build manifest artifact hash mismatch: {label}/{key}"
            )
        return path

    def resolved_source(block: dict[str, Any], key: str, label: str) -> Path:
        value = block.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(f"candidate build manifest lacks /{label}/{key}")
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = manifest_path.parent / path
        path = path.resolve()
        if not path.is_dir():
            raise ValueError(f"candidate source path is missing: {path}")
        return path

    cp2k_binary = resolved_artifact(cp2k, "binary", "cp2k")
    cp2k_library = resolved_artifact(cp2k, "loaded_library", "cp2k")
    resolved_artifact(cp2k, "cmake_cache", "cp2k")
    tblite_cli = resolved_artifact(save_tblite, "cli", "save_tblite")
    tblite_library = resolved_artifact(
        save_tblite, "static_library", "save_tblite"
    )
    resolved_artifact(save_tblite, "cmake_cache", "save_tblite")
    expected_paths = {
        "cp2k/binary": (cp2k_binary, identity.cp2k),
        "cp2k/loaded_library": (cp2k_library, identity.cp2k_library),
        "save_tblite/cli": (tblite_cli, args.tblite),
        "save_tblite/static_library": (
            tblite_library,
            identity.tblite_static_library,
        ),
        "cp2k/source_path": (
            resolved_source(cp2k, "source_path", "cp2k"),
            args.cp2k_source,
        ),
        "save_tblite/source_path": (
            resolved_source(save_tblite, "source_path", "save_tblite"),
            args.tblite_source,
        ),
    }
    mismatched_paths = [
        label
        for label, (manifest_value, actual) in expected_paths.items()
        if manifest_value != actual.resolve()
    ]
    if mismatched_paths:
        raise ValueError(
            "candidate build manifest path mismatch: " + ", ".join(mismatched_paths)
        )
    expected_values = {
        "cp2k/revision": (cp2k.get("revision"), identity.cp2k_source_revision),
        "cp2k/required_upstream_ancestor": (
            cp2k.get("required_upstream_ancestor"),
            args.required_cp2k_ancestor,
        ),
        "cp2k/binary_sha256": (cp2k.get("binary_sha256"), identity.cp2k_sha256),
        "cp2k/loaded_library_sha256": (
            cp2k.get("loaded_library_sha256"),
            identity.cp2k_library_sha256,
        ),
        "save_tblite/revision": (
            save_tblite.get("revision"),
            identity.tblite_source_revision,
        ),
        "save_tblite/static_library_sha256": (
            save_tblite.get("static_library_sha256"),
            identity.tblite_static_library_sha256,
        ),
        "save_tblite/cli_sha256": (
            save_tblite.get("cli_sha256"),
            runner.sha256(args.tblite),
        ),
    }
    mismatched_values = [
        label
        for label, (manifest_value, actual) in expected_values.items()
        if manifest_value != actual
    ]
    if mismatched_values:
        raise ValueError(
            "candidate build manifest identity mismatch: "
            + ", ".join(mismatched_values)
        )
    for block, label in ((cp2k, "cp2k"), (save_tblite, "save_tblite")):
        for key in ("repository", "branch"):
            if not isinstance(block.get(key), str) or not block[key]:
                raise ValueError(f"candidate build manifest lacks /{label}/{key}")
    reported_revision = cp2k.get("reported_revision")
    if not isinstance(reported_revision, str) or not identity.cp2k_source_revision.startswith(
        reported_revision
    ):
        raise ValueError("candidate build manifest /cp2k/reported_revision mismatch")
    for block, source, label in (
        (cp2k, args.cp2k_source, "cp2k"),
        (save_tblite, args.tblite_source, "save_tblite"),
    ):
        process = subprocess.run(
            ["git", "diff", "--binary", "HEAD"],
            cwd=source,
            capture_output=True,
            check=False,
        )
        if process.returncode != 0:
            raise ValueError(f"cannot hash candidate {label} source diff")
        actual_diff_hash = hashlib.sha256(process.stdout).hexdigest()
        if block.get("source_diff_sha256") != actual_diff_hash:
            raise ValueError(
                f"candidate build manifest /{label}/source_diff_sha256 mismatch"
            )

    oracle = qualification.get("oracle")
    force_stress = qualification.get("force_stress_fd")
    if not isinstance(oracle, dict) or not isinstance(force_stress, dict):
        raise ValueError("candidate build manifest lacks qualification gates")
    if oracle.get("status") != "passed" or force_stress.get("status") != "passed":
        raise ValueError("candidate build qualification gates have not passed")
    if qualification.get("oracle_binary_is_production_binary") is not True:
        raise ValueError("candidate oracle binary is not the production binary")
    qualification_identity = {
        "cp2k_binary_sha256": identity.cp2k_sha256,
        "cp2k_loaded_library_sha256": identity.cp2k_library_sha256,
        "save_tblite_library_sha256": identity.tblite_static_library_sha256,
        "source_diff_sha256": cp2k.get("source_diff_sha256"),
    }
    for gate, gate_label in (
        (oracle, "oracle"),
        (force_stress, "force_stress_fd"),
    ):
        for key, expected in qualification_identity.items():
            if gate.get(key) != expected:
                raise ValueError(
                    f"candidate {gate_label} identity mismatch: {key}"
                )
        for key in ("input", "output", "affinity_record"):
            resolved_artifact(gate, key, f"qualification/{gate_label}")
    if oracle.get("launch_contract") is not None:
        resolved_artifact(
            oracle, "launch_contract", "qualification/oracle"
        )
    for key in ("exchange_duality_max_relative", "acp_duality_max_relative"):
        value = finite_float(oracle.get(key), f"candidate oracle {key}")
        if value < 0.0:
            raise ValueError(f"candidate oracle {key} is negative")
    for key in (
        "force_max_abs_error_hartree_per_bohr",
        "force_rms_error_hartree_per_bohr",
        "force_max_relative_error",
        "virial_max_abs_error_hartree",
        "virial_rms_error_hartree",
        "virial_max_relative_error",
        "debug_displacement_bohr",
        "debug_cell_displacement_bohr",
    ):
        value = finite_float(
            force_stress.get(key), f"candidate force_stress_fd {key}"
        )
        if value < 0.0:
            raise ValueError(f"candidate force_stress_fd {key} is negative")


def require_clean_source(source: Path, label: str) -> None:
    process = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=source,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0 or process.stdout.strip():
        raise ValueError(f"{label} source worktree is not clean")


def require_ancestor(source: Path, ancestor: str) -> None:
    process = subprocess.run(
        ["git", "merge-base", "--is-ancestor", ancestor, "HEAD"],
        cwd=source,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(
            f"candidate CP2K source does not contain required ancestor {ancestor}"
        )


def relative_path(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def artifact(path: Path, root: Path, runner: Any) -> dict[str, str]:
    return {"path": relative_path(path, root), "sha256": runner.sha256(path)}


def validation_record_map(index: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    records = index.get("records")
    if not isinstance(records, list):
        raise ValueError("reference validation index lacks records")
    mapped: dict[tuple[str, str], dict[str, Any]] = {}
    for value in records:
        if not isinstance(value, dict):
            raise ValueError("reference validation record is not an object")
        key = (str(value.get("mesh")), str(value.get("phase")))
        if key in mapped:
            raise ValueError(f"duplicate reference validation record {key}")
        mapped[key] = value
    return mapped


def derive_full_publication_matrix(
    summary: dict[str, Any],
    reference_index: dict[str, Any],
    reference_index_sha256: str,
    runner: Any,
) -> list[tuple[str, str]]:
    """Derive every paper-used g-xTB energy from pinned summary metadata."""
    if (
        summary.get("benchmark") != "DMC-ICE13"
        or summary.get("status") != "phasewise_kpoint_converged"
    ):
        raise ValueError("final summary is not a converged DMC-ICE13 summary")
    if summary.get("reference_phase") != "Ih":
        raise ValueError("final summary does not use ice Ih as reference")
    sources = summary.get("sources")
    methods = summary.get("methods")
    fixed = summary.get("fixed_k333_same_mesh_comparison")
    if not isinstance(sources, dict) or not isinstance(methods, dict):
        raise ValueError("final summary lacks sources or methods")
    validation_source = sources.get("validation_index")
    if (
        not isinstance(validation_source, dict)
        or validation_source.get("sha256") != reference_index_sha256
    ):
        raise ValueError("final summary validation-index SHA256 mismatch")
    if (
        not isinstance(fixed, dict)
        or fixed.get("mesh") != FIXED_COMPARISON_MESH
        or fixed.get("not_a_phasewise_converged_result") is not True
    ):
        raise ValueError("final summary lacks the pinned k333 comparator")
    gxtb = methods.get("GXTB")
    if (
        not isinstance(gxtb, dict)
        or gxtb.get("status") != "phasewise_kpoint_converged"
    ):
        raise ValueError("final summary lacks converged g-xTB results")
    phases = gxtb.get("phases")
    provenance = gxtb.get("provenance")
    if not isinstance(phases, dict) or not isinstance(provenance, dict):
        raise ValueError("final summary lacks g-xTB phase/provenance metadata")
    source_identity = reference_index.get("source_identity")
    if not isinstance(source_identity, dict):
        raise ValueError("reference validation index lacks source identity")
    if (
        provenance.get("cp2k_source_revision")
        != source_identity.get("cp2k_source_revision")
        or provenance.get("provider_source_revision")
        != source_identity.get("tblite_source_revision")
    ):
        raise ValueError("final summary and validation-index source identities differ")

    required: set[tuple[str, str]] = {
        (FIXED_COMPARISON_MESH, phase) for phase in runner.PHASES
    }
    expected_nonreference = set(runner.PHASES) - {"Ih"}
    if set(phases) != expected_nonreference:
        raise ValueError("final summary does not cover exactly 12 non-reference phases")
    for phase in runner.PHASES:
        if phase == "Ih":
            continue
        phase_payload = phases[phase]
        if not isinstance(phase_payload, dict):
            raise ValueError(f"final summary phase metadata is invalid: {phase}")
        for field in ("selected_mesh", "previous_mesh"):
            mesh = phase_payload.get(field)
            if mesh not in runner.SUPPORTED_MESHES:
                raise ValueError(f"final summary has invalid {phase}/{field}")
            required.add((str(mesh), phase))
            required.add((str(mesh), "Ih"))

    mesh_order = {mesh: index for index, mesh in enumerate(runner.SUPPORTED_MESHES)}
    phase_order = {phase: index for index, phase in enumerate(runner.PHASES)}
    return sorted(
        required,
        key=lambda value: (mesh_order[value[0]], phase_order[value[1]]),
    )


def select_requalification_matrix(
    args: argparse.Namespace,
    reference_index: dict[str, Any],
    runner: Any,
) -> tuple[list[tuple[str, str]], dict[str, Any] | None]:
    if args.scope == SENTINEL_SCOPE:
        if not args.selection:
            raise ValueError("sentinel scope requires at least one --selection")
        validate_matrix(args.selection)
        return list(args.selection), None
    if args.scope != FULL_PUBLICATION_SCOPE:
        raise ValueError(f"unsupported requalification scope: {args.scope}")
    if args.selection:
        raise ValueError(
            "full-publication-matrix scope is derived from the pinned final "
            "summary and does not accept --selection"
        )
    if args.final_summary is None or args.final_summary_sha256 is None:
        raise ValueError(
            "full-publication-matrix scope requires --final-summary and "
            "--final-summary-sha256"
        )
    summary = read_pinned_json(
        args.final_summary, args.final_summary_sha256, runner
    )
    selections = derive_full_publication_matrix(
        summary,
        reference_index,
        args.reference_validation_index_sha256,
        runner,
    )
    validate_matrix(selections)
    return selections, summary


def qualification_outcome(scope: str) -> dict[str, object]:
    if scope == SENTINEL_SCOPE:
        return {
            "status": "sentinel_passed",
            "old_results_reusable": False,
            "paper_freeze_authorized": False,
        }
    if scope == FULL_PUBLICATION_SCOPE:
        return {
            "status": "full_publication_matrix_passed",
            "old_results_reusable": True,
            "paper_freeze_authorized": True,
        }
    raise ValueError(f"unsupported requalification scope: {scope}")


def parse_reference_energy(
    root: Path,
    record: dict[str, Any],
    index: dict[str, Any],
    runner: Any,
) -> float:
    identities = index.get("build_identities")
    if not isinstance(identities, dict):
        raise ValueError("reference validation index lacks build identities")
    identity = identities.get(str(record.get("build_id")))
    if not isinstance(identity, dict):
        raise ValueError("reference record has unknown build identity")
    output = root / str(record.get("output"))
    return runner._evidence_total_energy(
        output.read_bytes(),
        f"reference {record['mesh']}/{record['phase']}",
        expected_project=f"ice_{record['phase']}_GXTB_{record['mesh']}",
        expected_source_revision=str(identity["cp2k_source_revision"]),
        expected_tblite_source_revision=str(identity["tblite_source_revision"]),
        allow_unknown_tblite_revision=True,
    )


def build_jobs(
    root: Path,
    output_root: Path,
    selections: list[tuple[str, str]],
    runner: Any,
) -> list[Any]:
    jobs = []
    for mesh, phase in selections:
        input_path = (
            root
            / runner.GXTB_INPUT_DIRECTORY
            / mesh
            / f"ice_{phase}_GXTB_{mesh}.inp"
        )
        jobs.append(
            runner.Job(
                mesh=mesh,
                method="GXTB",
                phase=phase,
                input_path=input_path,
                run_dir=output_root / mesh / phase,
                output_name=f"ice_{phase}_GXTB_{mesh}.out",
            )
        )
    return jobs


def resolve_identity(args: argparse.Namespace, runner: Any) -> Any:
    return runner.production_identity(
        PROTOCOL_ID,
        args.cp2k,
        args.cp2k_library,
        args.tblite_static_library,
        args.cp2k_source,
        args.tblite_source,
        args.tblite_source_revision,
        require_embedded_tblite_revision=True,
    )


def run_matrix(
    args: argparse.Namespace,
    runner: Any,
    identity: Any,
    selections: list[tuple[str, str]],
) -> list[Any]:
    cpu_sets = args.cpu_set or []
    runner.validate_cpu_sets(
        cpu_sets,
        args.jobs,
        args.mpi_ranks_per_job,
        args.threads_per_job,
    )
    output_root = args.output_root
    jobs = build_jobs(args.root, output_root, selections, runner)
    cpu_pool: queue.Queue[str] | None = None
    if cpu_sets:
        cpu_pool = queue.Queue()
        for value in cpu_sets:
            cpu_pool.put(value)
    stop_event = threading.Event()
    failures: list[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [
            executor.submit(
                runner.run_job_with_cpu_pool,
                identity,
                job,
                args.force,
                stop_event,
                args.threads_per_job,
                args.mpi_ranks_per_job,
                args.mpi_launcher,
                cpu_pool,
                args.taskset,
                args.mpi_launcher_arg,
            )
            for job in jobs
        ]
        for future in concurrent.futures.as_completed(futures):
            job, returncode = future.result()
            if returncode != 0 or not runner.stamp_valid(job, identity):
                failures.append(f"{job.mesh}/{job.phase} (exit {returncode})")
                stop_event.set()
    if failures:
        raise ValueError("candidate requalification jobs failed: " + ", ".join(failures))
    return jobs


def build_report(
    args: argparse.Namespace,
    runner: Any,
    identity: Any,
    jobs: list[Any],
    reference_index: dict[str, Any],
    selections: list[tuple[str, str]],
) -> dict[str, Any]:
    repository = args.root.parent
    records = validation_record_map(reference_index)
    job_map = {(job.mesh, job.phase): job for job in jobs}
    rows: list[dict[str, Any]] = []
    energies: dict[tuple[str, str], tuple[float, float]] = {}
    if set(job_map) != set(selections) or len(job_map) != len(selections):
        raise ValueError("executed matrix does not exactly match the requested matrix")
    for mesh, phase in selections:
        record = records.get((mesh, phase))
        if record is None:
            raise ValueError(f"reference validation index lacks {mesh}/{phase}")
        job = job_map[(mesh, phase)]
        source_input = job.input_path
        candidate_input = runner.frozen_input_path(job)
        candidate_output = job.run_dir / job.output_name
        candidate_stamp = runner.stamp_path(job)
        reference_output = args.root / str(record["output"])
        reference_stamp = args.root / str(record["stamp"])
        reference_energy = parse_reference_energy(
            args.root, record, reference_index, runner
        )
        candidate_energy = runner._evidence_total_energy(
            candidate_output.read_bytes(),
            f"candidate {mesh}/{phase}",
            expected_project=f"ice_{phase}_GXTB_{mesh}",
            expected_source_revision=identity.cp2k_source_revision,
            expected_tblite_source_revision=identity.tblite_source_revision,
            allow_unknown_tblite_revision=False,
        )
        delta = candidate_energy - reference_energy
        if abs(delta) > TOTAL_ENERGY_TOLERANCE_HARTREE:
            raise ValueError(
                f"{mesh}/{phase} total-energy delta {delta:.16g} Eh exceeds "
                f"{TOTAL_ENERGY_TOLERANCE_HARTREE:.1e} Eh"
            )
        energies[(mesh, phase)] = (reference_energy, candidate_energy)
        rows.append(
            {
                "mesh": mesh,
                "phase": phase,
                "reference_build_id": record["build_id"],
                "candidate_build_id": runner.build_id(
                    runner.execution_build_identity(identity)
                ),
                "source_input": artifact(source_input, repository, runner),
                "candidate_executed_input": artifact(
                    candidate_input, repository, runner
                ),
                "reference_output": artifact(reference_output, repository, runner),
                "reference_stamp": artifact(reference_stamp, repository, runner),
                "candidate_output": artifact(candidate_output, repository, runner),
                "candidate_stamp": artifact(candidate_stamp, repository, runner),
                "reference_total_energy_hartree": reference_energy,
                "candidate_total_energy_hartree": candidate_energy,
                "total_energy_delta_hartree": delta,
                "absolute_total_energy_delta_hartree": abs(delta),
            }
        )

    relative_rows: list[dict[str, Any]] = []
    for mesh, phase in selections:
        if phase == "Ih":
            continue
        reference_phase, candidate_phase = energies[(mesh, phase)]
        reference_ih, candidate_ih = energies[(mesh, "Ih")]
        phase_input = job_map[(mesh, phase)].input_path.read_text()
        ih_input = job_map[(mesh, "Ih")].input_path.read_text()
        phase_count = runner._cp2k_input_water_count(
            phase_input, f"{mesh}/{phase}"
        )
        ih_count = runner._cp2k_input_water_count(ih_input, f"{mesh}/Ih")
        reference_relative = (
            reference_phase / phase_count - reference_ih / ih_count
        ) * HARTREE_TO_KJMOL
        candidate_relative = (
            candidate_phase / phase_count - candidate_ih / ih_count
        ) * HARTREE_TO_KJMOL
        delta = candidate_relative - reference_relative
        if abs(delta) > RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O:
            raise ValueError(
                f"{mesh}/{phase} relative-energy delta {delta:.16g} "
                "kJ/mol/H2O exceeds tolerance"
            )
        relative_rows.append(
            {
                "mesh": mesh,
                "phase": phase,
                "reference_relative_energy_kjmol_per_h2o": reference_relative,
                "candidate_relative_energy_kjmol_per_h2o": candidate_relative,
                "relative_energy_delta_kjmol_per_h2o": delta,
                "absolute_relative_energy_delta_kjmol_per_h2o": abs(delta),
            }
        )

    max_total = max(row["absolute_total_energy_delta_hartree"] for row in rows)
    max_relative = max(
        (row["absolute_relative_energy_delta_kjmol_per_h2o"] for row in relative_rows),
        default=0.0,
    )
    full_publication_matrix = args.scope == FULL_PUBLICATION_SCOPE
    outcome = qualification_outcome(args.scope)
    reference: dict[str, Any] = {
        "validation_index": artifact(
            args.reference_validation_index, repository, runner
        ),
        "campaign_id": reference_index.get("campaign_id"),
        "source_identity": reference_index.get("source_identity"),
    }
    if full_publication_matrix:
        if args.final_summary is None:
            raise ValueError("full publication report lacks a final-summary path")
        reference["final_summary"] = artifact(
            args.final_summary, repository, runner
        )
    return {
        "schema_version": 1,
        "protocol_id": PROTOCOL_ID,
        "benchmark": "DMC-ICE13",
        "method": "g-xTB",
        "status": outcome["status"],
        "qualification_scope": args.scope,
        "acceptance": {
            "total_energy_tolerance_hartree": TOTAL_ENERGY_TOLERANCE_HARTREE,
            "relative_energy_tolerance_kjmol_per_h2o": (
                RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O
            ),
            "observed_max_abs_total_energy_delta_hartree": max_total,
            "observed_max_abs_relative_energy_delta_kjmol_per_h2o": max_relative,
            "old_results_reusable": outcome["old_results_reusable"],
            "paper_freeze_authorized": outcome["paper_freeze_authorized"],
        },
        "reference": reference,
        "candidate": {
            "build_manifest": artifact(
                args.candidate_build_manifest, repository, runner
            ),
            "build_id": runner.build_id(runner.execution_build_identity(identity)),
            "build_identity": runner.execution_build_identity(identity),
            "cp2k": artifact(identity.cp2k, repository, runner),
            "cp2k_library": artifact(identity.cp2k_library, repository, runner),
            "tblite_static_library": artifact(
                identity.tblite_static_library, repository, runner
            ),
            "required_cp2k_ancestor": args.required_cp2k_ancestor,
        },
        "matrix": {
            "derivation": (
                "pinned_final_summary_selected_and_previous_meshes_plus_all_k333"
                if full_publication_matrix
                else "explicit_sentinel_selections"
            ),
            "selection_count": len(selections),
            "selections": [f"{mesh}:{phase}" for mesh, phase in selections],
            "exact_coverage_passed": True,
            "total_energy_comparisons": rows,
            "same_mesh_ih_relative_energy_comparisons": relative_rows,
        },
        "parallelism": runner.execution_parallelism(
            args.jobs,
            args.threads_per_job,
            args.mpi_ranks_per_job,
            args.mpi_launcher,
            args.cpu_set,
            args.mpi_launcher_arg,
        ),
    }


def run(args: argparse.Namespace, runner: Any) -> Path:
    args.report.unlink(missing_ok=True)
    try:
        if args.required_cp2k_ancestor != REQUIRED_CP2K_ANCESTOR:
            raise ValueError(
                "the required CP2K ancestor is protocol-fixed at "
                f"{REQUIRED_CP2K_ANCESTOR}"
            )
        require_clean_source(args.cp2k_source, "CP2K")
        require_clean_source(args.tblite_source, "save_tblite")
        require_ancestor(args.cp2k_source, args.required_cp2k_ancestor)
        candidate_manifest = read_pinned_json(
            args.candidate_build_manifest,
            args.candidate_build_manifest_sha256,
            runner,
        )
        reference_index = runner.read_validation_index(
            args.reference_validation_index,
            args.root,
            expected_index_sha256=args.reference_validation_index_sha256,
        )
        selections, _ = select_requalification_matrix(
            args, reference_index, runner
        )
        records = validation_record_map(reference_index)
        for selection in selections:
            if selection not in records:
                raise ValueError(
                    f"reference validation index lacks "
                    f"{selection[0]}/{selection[1]}"
                )
        identity = resolve_identity(args, runner)
        validate_candidate_manifest(candidate_manifest, identity, args, runner)
        jobs = run_matrix(args, runner, identity, selections)
        report = build_report(
            args,
            runner,
            identity,
            jobs,
            reference_index,
            selections,
        )
        content = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
        runner.atomic_write_bytes(args.report, content)
        return args.report
    except Exception:
        args.report.unlink(missing_ok=True)
        raise


def main() -> None:
    repository = Path(__file__).resolve().parents[2]
    runner = load_runner(repository)
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--reference-validation-index", type=Path, required=True)
    parser.add_argument("--reference-validation-index-sha256", required=True)
    parser.add_argument(
        "--scope",
        choices=(SENTINEL_SCOPE, FULL_PUBLICATION_SCOPE),
        default=SENTINEL_SCOPE,
    )
    parser.add_argument("--final-summary", type=Path)
    parser.add_argument("--final-summary-sha256")
    parser.add_argument("--candidate-build-manifest", type=Path, required=True)
    parser.add_argument("--candidate-build-manifest-sha256", required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cp2k-library", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--tblite-static-library", type=Path, required=True)
    parser.add_argument("--tblite-source", type=Path, required=True)
    parser.add_argument("--tblite-source-revision", required=True)
    parser.add_argument(
        "--required-cp2k-ancestor", default=REQUIRED_CP2K_ANCESTOR
    )
    parser.add_argument(
        "--selection",
        action="append",
        type=lambda value: parse_selection(value, runner),
        default=[],
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--jobs", type=runner.positive_int, default=1)
    parser.add_argument("--threads-per-job", type=runner.positive_int, default=1)
    parser.add_argument("--mpi-ranks-per-job", type=runner.positive_int, default=1)
    parser.add_argument("--mpi-launcher", default="mpiexec")
    parser.add_argument("--mpi-launcher-arg", action="append", default=[])
    parser.add_argument("--cpu-set", action="append", default=[])
    parser.add_argument("--taskset", default="taskset")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for field in (
        "root",
        "reference_validation_index",
        "candidate_build_manifest",
        "cp2k",
        "cp2k_library",
        "cp2k_source",
        "tblite",
        "tblite_static_library",
        "tblite_source",
        "output_root",
        "report",
    ):
        setattr(args, field, getattr(args, field).expanduser().resolve())
    if args.final_summary is not None:
        args.final_summary = args.final_summary.expanduser().resolve()
    try:
        report = run(args, runner)
    except ValueError as error:
        parser.error(str(error))
    print(report)


if __name__ == "__main__":
    main()
