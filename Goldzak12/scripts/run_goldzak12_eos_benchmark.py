#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import benchmark_execution as execution  # noqa: E402
import run_goldzak12_benchmark as base  # noqa: E402


ROOT = base.ROOT
DEFAULT_SCALES = (0.82, 0.88, 0.94, 0.98, 1.00, 1.02, 1.06, 1.12, 1.20, 1.30, 1.45)
DEFAULT_RESULT_MESH = "k555"
MINIMUM_REDUCED_GXTB_FITS = 8
# Physical LC12 EOS curves vary by less than 2 Eh over the sampled interval.
# A much lower total energy is the numerical signature of the SCC charge collapse.
CHARGE_COLLAPSE_ENERGY_DROP_HARTREE = 25.0
CATASTROPHIC_CHARGE_COLLAPSE_ENERGY_HARTREE = -50.0
ADAPTIVE_SCALES = {
    ("MgO", "GFN2"): (0.90, 0.92, 0.926, 0.927, 0.928, 0.93),
    ("LiH", "GFN2"): (
        0.71,
        0.72,
        0.73,
        0.74,
        0.75,
        0.76,
        0.77,
        0.78,
        0.79,
        0.80,
        0.81,
        0.83,
        0.84,
        0.85,
        0.86,
        0.87,
        0.89,
        0.90,
        0.91,
        0.92,
        0.93,
    ),
}
BUILTIN_ADAPTIVE_SCALES = {key: tuple(values) for key, values in ADAPTIVE_SCALES.items()}
GXTB_BRANCH_DISCONTINUITY_HARTREE = 1.0
# A point that converges to a different SCC root can remain locally smooth
# enough to evade the single-mesh EOS test.  Across adjacent k meshes the
# physical total-energy shift is nevertheless a smooth function of volume.
# Flag an isolated departure from that shift before it can enter an EOS fit.
GXTB_CROSS_MESH_SHIFT_FLOOR_HARTREE = 1.0e-2
GXTB_CROSS_MESH_MAD_MULTIPLIER = 8.0
# The LC12 conventional cells all contain eight atoms. This tolerance only
# suppresses sub-meV/atom numerical noise in the GXTB EOS topology test; it is
# not an energy-discontinuity threshold.
GXTB_TOPOLOGY_TOLERANCE_HARTREE = 1.0e-4
GXTB_CLASSIFICATION_ACTIONS = {"exclude", "retain"}
FINAL_INPUT_LINEAGE_SCHEMA = 1


def gxtb_scale_manifest_path() -> Path:
    return ROOT / "data" / "gxtb_eos_scale_manifest.json"


def gxtb_classification_manifest_path() -> Path:
    return ROOT / "data" / "gxtb_eos_classifications.json"


def adaptive_scales_only(solid: str, method: str) -> tuple[float, ...]:
    built_in = set(BUILTIN_ADAPTIVE_SCALES.get((solid, method), ()))
    return tuple(value for value in ADAPTIVE_SCALES.get((solid, method), ()) if value not in built_in)


def restore_gxtb_scale_manifest(mesh: str, methods: tuple[str, ...]) -> None:
    """Restore previously requested adaptive GXTB points for a safe resume."""
    path = gxtb_scale_manifest_path()
    if "GXTB" not in methods or not path.exists():
        return
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return
    if payload.get("schema_version") != 1 or payload.get("eos_mesh") != mesh:
        return
    for record in payload.get("systems", []):
        if record.get("method") != "GXTB":
            continue
        solid = str(record.get("solid", ""))
        values = tuple(float(value) for value in record.get("adaptive_scales", []))
        if solid and values:
            key = (solid, "GXTB")
            ADAPTIVE_SCALES[key] = tuple(sorted(set(ADAPTIVE_SCALES.get(key, ())) | set(values)))


def write_gxtb_scale_manifest(
    mesh: str,
    base_scales: tuple[float, ...],
    methods: tuple[str, ...],
    solids: tuple[str, ...] = base.LC10_PAPER_SOLIDS,
) -> dict[str, object] | None:
    if "GXTB" not in methods:
        return None
    payload: dict[str, object] = {
        "schema_version": 1,
        "benchmark": "LC10 (fixed Goldzak12 subset)",
        "paper_systems": list(base.LC10_PAPER_SOLIDS),
        "diagnostic_only_systems": list(base.LC10_DIAGNOSTIC_ONLY_SOLIDS),
        "eos_mesh": mesh,
        "method": "GXTB",
        "base_scales": list(base_scales),
        "systems": [
            {
                "solid": ref.solid,
                "method": "GXTB",
                "requested_scales": list(scales_for(ref.solid, "GXTB", base_scales)),
                "adaptive_scales": list(adaptive_scales_only(ref.solid, "GXTB")),
            }
            for ref in base.REFERENCES
            if ref.solid in solids
        ],
    }
    path = gxtb_scale_manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def read_gxtb_scale_manifest() -> dict[str, object] | None:
    path = gxtb_scale_manifest_path()
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    return payload if payload.get("schema_version") == 1 else None


def load_gxtb_classifications(path: Path) -> dict[tuple[str, str, float], dict[str, str]]:
    """Load explicit branch exclusions/waivers; heuristic candidates never self-resolve."""
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    entries = payload.get("entries", [])
    refs = {ref.solid for ref in base.REFERENCES}
    result: dict[tuple[str, str, float], dict[str, str]] = {}
    for entry in entries:
        solid = str(entry.get("solid", ""))
        mesh = str(entry.get("mesh", ""))
        method = str(entry.get("method", "GXTB"))
        action = str(entry.get("action", ""))
        classification = str(entry.get("classification", "")).strip()
        rationale = str(entry.get("rationale", "")).strip()
        try:
            scale = round(float(entry["scale"]), 5)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid GXTB classification scale: {entry!r}") from exc
        if solid not in refs or method != "GXTB" or not mesh:
            raise ValueError(f"Invalid GXTB classification target: {entry!r}")
        if action not in GXTB_CLASSIFICATION_ACTIONS or not classification or not rationale:
            raise ValueError(
                "Every GXTB classification needs action=exclude|retain, classification, and rationale: "
                f"{entry!r}"
            )
        key = (solid, mesh, scale)
        if key in result:
            raise ValueError(f"Duplicate GXTB classification: {key}")
        result[key] = {
            "action": action,
            "classification": classification,
            "rationale": rationale,
        }
    return result


def gxtb_branch_candidates(
    points: list[tuple[float, float, float | None, bool]],
) -> dict[float, float]:
    """Flag, but never automatically exclude, large local EOS discontinuities."""
    valid = sorted((scale, float(energy)) for _, scale, energy, ok in points if ok and energy is not None)
    candidates: dict[float, float] = {}
    for index in range(1, len(valid) - 1):
        left_scale, left_energy = valid[index - 1]
        scale, energy = valid[index]
        right_scale, right_energy = valid[index + 1]
        fraction = (scale - left_scale) / (right_scale - left_scale)
        interpolated = left_energy + fraction * (right_energy - left_energy)
        residual = abs(energy - interpolated)
        if residual >= GXTB_BRANCH_DISCONTINUITY_HARTREE:
            candidates[round(scale, 5)] = residual
    # Local interpolation cannot diagnose an endpoint.  A jump far larger
    # than the documented physical LC12 energy span is therefore reported as
    # an endpoint candidate, still without automatically excluding it.
    if len(valid) >= 2:
        for endpoint, neighbor in ((valid[0], valid[1]), (valid[-1], valid[-2])):
            scale, energy = endpoint
            jump = abs(energy - neighbor[1])
            if jump >= CHARGE_COLLAPSE_ENERGY_DROP_HARTREE:
                candidates[round(scale, 5)] = jump
    return candidates


def previous_cubic_mesh(mesh: str) -> str | None:
    """Return the immediately preceding cubic mesh label, if one exists."""
    match = re.fullmatch(r"k([1-9][0-9]*)\1\1", mesh)
    if match is None:
        raise ValueError(f"non-cubic k mesh {mesh!r}")
    number = int(match.group(1))
    if number <= 1:
        return None
    previous = number - 1
    return f"k{previous}{previous}{previous}"


def gxtb_cross_mesh_branch_candidates(
    points: list[tuple[float, float, float | None, bool]],
    previous_points: list[tuple[float, float, float | None, bool]],
) -> dict[float, float]:
    """Flag isolated SCC-root changes from adjacent-mesh energy shifts.

    Absolute total energies are not compared directly.  For scales present in
    both meshes we form ``E_N(scale) - E_(N-1)(scale)`` and require a candidate
    to be anomalous both relative to the robust shift distribution and to the
    local interpolation of its neighbours.  This distinguishes a single
    alternate SCC solution from a smooth physical k-point correction.  The
    detector is fail-closed only: it reports candidates but never excludes a
    point without an explicit reviewed classification.
    """
    current = {
        round(float(scale), 5): float(energy)
        for _, scale, energy, ok in points
        if ok and energy is not None
    }
    previous = {
        round(float(scale), 5): float(energy)
        for _, scale, energy, ok in previous_points
        if ok and energy is not None
    }
    shifts = sorted(
        (scale, current[scale] - previous[scale])
        for scale in set(current) & set(previous)
    )
    if len(shifts) < 5:
        return {}

    values = np.array([shift for _, shift in shifts], dtype=float)
    centre = float(np.median(values))
    mad = float(np.median(np.abs(values - centre)))
    cutoff = max(
        GXTB_CROSS_MESH_SHIFT_FLOOR_HARTREE,
        GXTB_CROSS_MESH_MAD_MULTIPLIER * 1.4826 * mad,
    )
    candidates: dict[float, float] = {}
    for index in range(1, len(shifts) - 1):
        left_scale, left_shift = shifts[index - 1]
        scale, shift = shifts[index]
        right_scale, right_shift = shifts[index + 1]
        fraction = (scale - left_scale) / (right_scale - left_scale)
        interpolated = left_shift + fraction * (right_shift - left_shift)
        local_residual = abs(shift - interpolated)
        global_residual = abs(shift - centre)
        if local_residual >= cutoff and global_residual >= cutoff:
            candidates[scale] = local_residual

    # Endpoints have no two-sided interpolation.  Use a deliberately stricter
    # gate so that a smooth scale dependence is not mistaken for an SCC switch.
    endpoint_cutoff = 2.0 * cutoff
    for index, neighbor_index in ((0, 1), (len(shifts) - 1, len(shifts) - 2)):
        scale, shift = shifts[index]
        neighbor_shift = shifts[neighbor_index][1]
        residual = abs(shift - neighbor_shift)
        if residual >= endpoint_cutoff and abs(shift - centre) >= endpoint_cutoff:
            candidates[scale] = residual
    return candidates


def scale_tag(scale: float, method: str = "") -> str:
    # Preserve frozen GFN1/GFN2 paths while giving targeted GXTB points enough
    # precision to prevent adaptive-grid filename collisions.
    precision = 5 if method == "GXTB" else 3
    return f"s{scale:.{precision}f}".replace(".", "p")


def eos_project(solid: str, method: str, mesh: str, scale: float) -> str:
    return f"{solid}_{method}_eos_{mesh}_{scale_tag(scale, method)}"


def final_project(solid: str, method: str, mesh: str) -> str:
    return f"{solid}_{method}_eos_final_{mesh}"


def final_input_path(solid: str, method: str, mesh: str) -> Path:
    project = final_project(solid, method, mesh)
    return ROOT / "runs" / "eos_final_sp" / method / solid / mesh / f"{project}.inp"


def final_input_lineage_path(input_path: Path) -> Path:
    return input_path.with_suffix(input_path.suffix + ".eos.json")


def write_final_input_lineage(
    input_path: Path,
    fit: dict[str, object],
    mesh: str,
    *,
    valid: bool,
    reason: str,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": FINAL_INPUT_LINEAGE_SCHEMA,
        "benchmark": "LC12 (Goldzak12)",
        "valid": valid,
        "reason": reason,
        "solid": str(fit.get("solid", "")),
        "method": str(fit.get("method", "")),
        "eos_mesh": str(fit.get("eos_mesh", "")),
        "energy_mesh": mesh,
        "fit_status": str(fit.get("fit_status", "")),
        "a_eos_A": str(fit.get("a_eos_A", "")),
        "input": str(input_path.resolve()),
        "input_sha256": base.sha256(input_path) if input_path.is_file() else "",
        "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
    }
    path = final_input_lineage_path(input_path)
    base.write_file(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def final_input_lineage_issue(
    input_path: Path, fit: dict[str, object], mesh: str
) -> str | None:
    """Return why a GXTB final input is not derived from the current EOS fit."""
    if not input_path.is_file():
        return f"missing final input {input_path}"
    lineage_path = final_input_lineage_path(input_path)
    if not lineage_path.is_file():
        return f"missing EOS lineage for {input_path}"
    try:
        lineage = json.loads(lineage_path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid EOS lineage for {input_path}: {exc}"
    expected = {
        "schema_version": FINAL_INPUT_LINEAGE_SCHEMA,
        "valid": True,
        "solid": str(fit.get("solid", "")),
        "method": "GXTB",
        "eos_mesh": str(fit.get("eos_mesh", "")),
        "energy_mesh": mesh,
        "fit_status": "quadratic",
        "a_eos_A": str(fit.get("a_eos_A", "")),
        "input_sha256": base.sha256(input_path),
        "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
    }
    for key, value in expected.items():
        if lineage.get(key) != value:
            return f"EOS lineage {key} mismatch for {input_path}"
    refs = {ref.solid: ref for ref in base.REFERENCES}
    solid = str(fit.get("solid", ""))
    if solid not in refs or not fit.get("a_eos_A"):
        return f"invalid EOS fit identity for {input_path}"
    project = final_project(solid, "GXTB", mesh)
    expected_text = base.solid_input(
        refs[solid], "GXTB", "ENERGY", mesh, float(fit["a_eos_A"]), project
    )
    if input_path.read_text() != expected_text:
        return f"final input does not match the current EOS minimum for {input_path}"
    try:
        base.validate_method_input(input_path.read_text(), "GXTB")
    except ValueError as exc:
        return str(exc)
    return None


def invalidate_existing_gxtb_final_inputs(
    fits: list[dict[str, object]], meshes: list[str]
) -> None:
    """Mark pre-fit/reference-cell final inputs unusable until regenerated."""
    for fit in fits:
        if fit.get("method") != "GXTB":
            continue
        valid_fit = bool(fit.get("a_eos_A")) and fit.get("fit_status") == "quadratic"
        for mesh in meshes:
            input_path = final_input_path(str(fit["solid"]), "GXTB", mesh)
            if not input_path.is_file():
                continue
            if valid_fit and final_input_lineage_issue(input_path, fit, mesh) is None:
                continue
            reason = (
                "awaiting regeneration from the current quadratic EOS minimum"
                if valid_fit
                else "no valid quadratic EOS minimum; reference-cell/stale input is not runnable"
            )
            write_final_input_lineage(input_path, fit, mesh, valid=False, reason=reason)


def scales_for(solid: str, method: str, scales: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(sorted(set(scales) | set(ADAPTIVE_SCALES.get((solid, method), ()))))


def eos_job_specs(
    mesh: str,
    scales: tuple[float, ...],
    methods: tuple[str, ...] = base.METHODS,
    solids: tuple[str, ...] | None = None,
) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in base.REFERENCES:
        if solids is not None and ref.solid not in solids:
            continue
        for method in methods:
            for scale in scales_for(ref.solid, method, scales):
                project = eos_project(ref.solid, method, mesh, scale)
                a = ref.a_exp * scale
                run_dir = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale, method)
                inp = run_dir / f"{project}.inp"
                out = run_dir / f"{project}.out"
                text = base.solid_input(ref, method, "ENERGY", mesh, a, project)
                base.write_file(inp, text)
                specs.append((f"eos {method} {ref.solid} {mesh} {scale_tag(scale, method)}", inp, out, False))
    return specs


def strategy_path(output: Path) -> Path:
    return output.with_suffix(".strategy.json")


def write_strategy(output: Path, strategy: str, completed: bool) -> None:
    strategy_path(output).write_text(json.dumps({"strategy": strategy, "completed": completed}, indent=2) + "\n")


def read_strategy(output: Path) -> str:
    path = strategy_path(output)
    if not path.exists():
        return "unknown"
    try:
        return str(json.loads(path.read_text())["strategy"])
    except (json.JSONDecodeError, KeyError):
        return "unknown"


def output_charge_collapsed(output: Path, method: str) -> bool:
    if method != "GFN2":
        return False
    energy = base.parse_energy(output)
    return energy is not None and energy < CATASTROPHIC_CHARGE_COLLAPSE_ENERGY_HARTREE


def usable_output_ok(output: Path, method: str, require_opt: bool = False) -> bool:
    return base.output_ok(output, require_opt=require_opt) and not output_charge_collapsed(output, method)


def spec_method(spec: tuple[str, Path, Path, bool]) -> str:
    parts = spec[0].split()
    if len(parts) < 2 or parts[1] not in base.METHODS:
        raise ValueError(f"Cannot determine method from job label {spec[0]!r}")
    return parts[1]


def retries_exhausted(output: Path) -> bool:
    path = strategy_path(output)
    if not path.exists():
        return False
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return False
    return data.get("strategy") == "retry_m1_d001" and not bool(data.get("completed"))


def retry_input(inp: Path, iterations: int, memory: int, damping: float, label: str) -> Path:
    text = inp.read_text()
    marker = "        &END TBLITE\n      &END XTB"
    replacement = (
        "        &END TBLITE\n"
        "        SCC_MIXER TBLITE\n"
        "        &TBLITE_MIXER\n"
        f"          ITERATIONS {iterations}\n"
        f"          MEMORY {memory}\n"
        f"          DAMPING {damping:.6f}\n"
        "        &END TBLITE_MIXER\n"
        "      &END XTB"
    )
    if marker not in text:
        raise ValueError(f"Cannot insert TBLITE_MIXER into {inp}")
    text = text.replace(marker, replacement, 1)
    text = text.replace("      MAX_SCF 300", f"      MAX_SCF {iterations}", 1)
    path = inp.with_name(f"{inp.stem}_{label}.inp")
    base.write_file(path, text)
    return path


def cp2k_execution_command_contract(
    threads: int,
    execution_pool: execution.ExecutionPool | None,
) -> dict[str, object]:
    return {
        "driver": "cp2k",
        "omp_threads": threads,
        "mpi_ranks_per_job": (
            execution_pool.mpi_ranks_per_job if execution_pool is not None else 1
        ),
        "execution_mode": (
            "openmpi_ordered_pe_list" if execution_pool is not None else "direct"
        ),
        "execution_contract_sha256": (
            str(execution_pool.contract_sha256)
            if execution_pool is not None
            else None
        ),
    }


def run_jobs(
    specs: list[tuple[str, Path, Path, bool]],
    cp2k: Path,
    jobs: int,
    threads: int,
    force: bool,
    retry_scf: bool = True,
    campaign_fingerprint: dict[str, object] | None = None,
    execution_pool: execution.ExecutionPool | None = None,
    campaign_bind_all_methods: bool = False,
) -> None:
    if campaign_bind_all_methods and campaign_fingerprint is None:
        raise ValueError(
            "campaign_bind_all_methods requires a validated campaign fingerprint"
        )
    cp2k_identity = base.executable_fingerprint(cp2k)
    command_contract = cp2k_execution_command_contract(threads, execution_pool)

    def signature(spec: tuple[str, Path, Path, bool]) -> dict[str, object]:
        return base.job_signature(
            cp2k,
            spec[1],
            command_contract=command_contract,
            executable_identity=cp2k_identity,
            campaign_fingerprint=campaign_fingerprint,
        )

    pending: list[tuple[str, Path, Path, bool]] = []
    exhausted_before_run: list[tuple[str, Path, Path, bool]] = []
    execution_evidence_failures: list[tuple[str, str]] = []
    for spec in specs:
        method = spec_method(spec)
        if method == "GXTB" and campaign_fingerprint is None:
            raise ValueError("GXTB EOS jobs require a validated campaign fingerprint")
        base.validate_method_input(spec[1].read_text(), method)
        usable = usable_output_ok(spec[2], method, require_opt=spec[3])
        campaign_bound = method == "GXTB" or campaign_bind_all_methods
        stamped = not campaign_bound or base.job_stamp_matches(spec[2], signature(spec))
        if not force and usable and stamped:
            if execution_pool is None:
                continue
            issue = execution_pool.record_issue(spec[2], base.job_stamp_path(spec[2]))
            if issue is None:
                continue
            execution_evidence_failures.append((spec[0], issue))
            continue
        if not force and method != "GXTB" and retries_exhausted(spec[2]):
            exhausted_before_run.append(spec)
            continue
        pending.append(spec)
    if execution_evidence_failures:
        details = "; ".join(
            f"{label}: {issue}" for label, issue in execution_evidence_failures
        )
        raise RuntimeError(
            "scientifically completed outputs have missing/invalid separate execution "
            "evidence; refusing an implicit destructive rerun (use --force only after review): "
            + details
        )
    if not pending:
        if exhausted_before_run:
            raise RuntimeError(
                "Previously exhausted CP2K SCF retries remain unresolved: "
                + ", ".join(sorted(spec[0] for spec in exhausted_before_run))
            )
        print("No EOS jobs pending.")
        return

    def worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool, tuple[str, Path, Path, bool]]:
        label, inp, out, require_opt = spec
        method = spec_method(spec)
        observation: dict[str, object] | None = None
        if execution_pool is None:
            code = base.run_cp2k(cp2k, inp, out, threads)
        else:
            code, observation = execution_pool.run_cp2k(cp2k, inp, out)
        ok = usable_output_ok(out, method, require_opt=require_opt)
        if (
            execution_pool is not None
            and (
                observation is None
                or observation.get("runtime_affinity_gate") is not True
            )
        ):
            ok = False
            if code == 0:
                code = 97
        strategy = "native_gxtb_fdiis" if method == "GXTB" else "default_tblite_mixer"
        write_strategy(out, strategy, ok)
        if method == "GXTB" or campaign_bind_all_methods:
            stamp = base.job_stamp_path(out)
            if execution_pool is not None and not ok:
                stamp.unlink(missing_ok=True)
                return label, code, ok, spec
            base.write_job_stamp(out, signature(spec), completed=ok, return_code=code)
            if execution_pool is not None:
                assert observation is not None
                try:
                    execution_pool.write_record(out, observation, stamp)
                except Exception:
                    stamp.unlink(missing_ok=True)
                    raise
        return label, code, ok, spec

    parallel_label = (
        f", MPI ranks/job={execution_pool.mpi_ranks_per_job}, exact ordered PE-list pool"
        if execution_pool is not None
        else ""
    )
    print(
        f"Running {len(pending)} CP2K jobs with {jobs} worker(s), "
        f"OMP_NUM_THREADS={threads}{parallel_label}."
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in pending}
        done = 0
        failed: list[tuple[str, Path, Path, bool]] = []
        for future in concurrent.futures.as_completed(futures):
            label, code, ok, spec = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(pending):3d}] {status:14s} {label}", flush=True)
            if not ok:
                failed.append(spec)

    gxtb_failed = [spec for spec in failed if spec_method(spec) == "GXTB"]
    terminal_failed = exhausted_before_run + gxtb_failed
    if gxtb_failed:
        print(
            "Native g-XTB FDIIS did not converge for: "
            + ", ".join(spec[0] for spec in gxtb_failed)
            + "; no alternative mixer retry is permitted by the production protocol.",
            file=sys.stderr,
        )
    failed = [spec for spec in failed if spec_method(spec) != "GXTB"]
    if not retry_scf:
        terminal_failed.extend(failed)
        failed = []

    profiles = (
        ("retry_m1_d005", 1200, 1, 0.05),
        ("retry_m1_d001", 2400, 1, 0.01),
    )
    for profile, iterations, memory, damping in profiles if retry_scf else ():
        if not failed:
            break
        print(
            f"Retrying {len(failed)} failed job(s) with TBLITE_MIXER "
            f"MEMORY={memory}, DAMPING={damping}, ITERATIONS={iterations}."
        )

        def retry_worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool, tuple[str, Path, Path, bool]]:
            label, inp, out, require_opt = spec
            robust_inp = retry_input(inp, iterations, memory, damping, profile)
            code = base.run_cp2k(cp2k, robust_inp, out, threads)
            ok = usable_output_ok(out, spec_method(spec), require_opt=require_opt)
            write_strategy(out, profile, ok)
            return label, code, ok, spec

        next_failed: list[tuple[str, Path, Path, bool]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(jobs, 4)) as pool:
            futures = {pool.submit(retry_worker, spec): spec for spec in failed}
            done = 0
            for future in concurrent.futures.as_completed(futures):
                label, code, ok, spec = future.result()
                done += 1
                status = "ok" if ok else f"failed rc={code}"
                print(f"[retry {done:2d}/{len(failed):2d}] {status:14s} {label}", flush=True)
                if not ok:
                    next_failed.append(spec)
        failed = next_failed

    if failed:
        print("SCF retries exhausted for: " + ", ".join(spec[0] for spec in failed), file=sys.stderr)
        terminal_failed.extend(failed)
    if terminal_failed:
        raise RuntimeError(
            f"{len(terminal_failed)} CP2K job(s) failed after preserving completed jobs: "
            + ", ".join(sorted(spec[0] for spec in terminal_failed))
        )


def load_eos_points(
    ref: base.Reference,
    method: str,
    mesh: str,
    scales: tuple[float, ...],
    campaign_fingerprint: dict[str, object] | None = None,
) -> list[tuple[float, float, float | None, bool]]:
    points: list[tuple[float, float, float | None, bool]] = []
    for scale in scales:
        project = eos_project(ref.solid, method, mesh, scale)
        out = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale, method) / f"{project}.out"
        energy = base.parse_energy(out)
        ok = base.output_ok(out)
        if method == "GXTB":
            if campaign_fingerprint is None:
                raise ValueError("GXTB EOS collection requires a campaign fingerprint")
            issue = base.completed_stamp_campaign_issue(
                out,
                campaign_fingerprint,
                executable_role="cp2k",
                require_completed=ok,
            )
            if issue:
                raise RuntimeError(issue)
        points.append((ref.a_exp * scale, scale, energy, ok))
    return sorted(points)


def charge_collapsed_scales(points: list[tuple[float, float, float | None, bool]]) -> set[float]:
    converged = [energy for _, _, energy, ok in points if ok and energy is not None]
    if len(converged) < 3:
        return set()
    reference = float(np.median(converged))
    return {
        scale
        for _, scale, energy, ok in points
        if ok and energy is not None and energy < reference - CHARGE_COLLAPSE_ENERGY_DROP_HARTREE
    }


def fit_eos(points: list[tuple[float, float, float, bool]]) -> dict[str, object]:
    if len(points) < 3:
        return {"a_eos_A": "", "energy_fit_hartree": "", "fit_status": "insufficient_points", "n_points": len(points)}
    energies = np.array([p[2] for p in points], dtype=float)
    local_minima = [
        i
        for i in range(1, len(points) - 1)
        if points[i][2] < points[i - 1][2] and points[i][2] < points[i + 1][2]
    ]
    if not local_minima:
        return {
            "a_eos_A": "",
            "energy_fit_hartree": "",
            "fit_status": "no_local_minimum",
            "n_points": len(points),
            "grid_min_a_A": f"{points[int(np.argmin(energies))][0]:.10f}",
            "grid_min_scale": f"{points[int(np.argmin(energies))][1]:.5f}",
            "grid_min_energy_hartree": f"{points[int(np.argmin(energies))][2]:.12f}",
        }
    preferred = [i for i in local_minima if 0.88 <= points[i][1] <= 1.12]
    if preferred:
        idx = min(preferred, key=lambda i: points[i][2])
    else:
        idx = min(local_minima, key=lambda i: abs(points[i][1] - 1.0))
    lo = max(0, idx - 2)
    hi = min(len(points), idx + 3)
    if hi - lo < 3:
        lo = max(0, min(lo, len(points) - 3))
        hi = min(len(points), lo + 3)
    fit_points = points[lo:hi]
    x = np.array([p[0] for p in fit_points], dtype=float)
    y = np.array([p[2] for p in fit_points], dtype=float)
    coeff = np.polyfit(x, y, 2)
    fit_rmse = float(np.sqrt(np.mean((np.polyval(coeff, x) - y) ** 2)))
    status = "quadratic"
    if coeff[0] <= 0:
        a_min = points[idx][0]
        e_min = points[idx][2]
        status = "grid_min_negative_curvature"
    else:
        a_min = float(-coeff[1] / (2.0 * coeff[0]))
        e_min = float(np.polyval(coeff, a_min))
        if a_min < points[0][0] or a_min > points[-1][0]:
            a_min = points[idx][0]
            e_min = points[idx][2]
            status = "grid_min_outside_fit"
        elif fit_rmse > 2.0e-2 or e_min > points[idx][2] + 2.0e-2:
            return {
                "a_eos_A": "",
                "energy_fit_hartree": "",
                "fit_status": "poor_quadratic_fit",
                "fit_rmse_hartree": f"{fit_rmse:.12f}",
                "n_points": len(points),
                "grid_min_a_A": f"{points[idx][0]:.10f}",
                "grid_min_scale": f"{points[idx][1]:.5f}",
                "grid_min_energy_hartree": f"{points[idx][2]:.12f}",
            }
    return {
        "a_eos_A": f"{a_min:.10f}",
        "energy_fit_hartree": f"{e_min:.12f}",
        "fit_status": status,
        "fit_rmse_hartree": f"{fit_rmse:.12f}",
        "n_points": len(points),
        "grid_min_a_A": f"{points[idx][0]:.10f}",
        "grid_min_scale": f"{points[idx][1]:.5f}",
        "grid_min_energy_hartree": f"{points[idx][2]:.12f}",
    }


def gxtb_topology_reversals(
    points: list[tuple[float, float, float, bool]],
    tolerance_hartree: float = GXTB_TOPOLOGY_TOLERANCE_HARTREE,
) -> list[tuple[float, float, float]]:
    """Return energy-direction reversals relative to the sampled global minimum.

    A single-well EOS must decrease monotonically towards its sampled global
    minimum and increase monotonically away from it. This catches a lower
    compressed branch plus a higher local well even when a local quadratic fit
    happens to look acceptable. Adaptive points cannot remove an existing
    reversal; only an explicitly reviewed point exclusion can change the set
    passed to this gate.
    """
    ordered = sorted(
        (float(a), float(scale), float(energy))
        for a, scale, energy, ok in points
        if ok and energy is not None
    )
    if len(ordered) < 3:
        return []
    global_index = min(range(len(ordered)), key=lambda index: ordered[index][2])
    reversals: list[tuple[float, float, float]] = []
    for index in range(global_index):
        delta = ordered[index + 1][2] - ordered[index][2]
        if delta > tolerance_hartree:
            reversals.append((ordered[index][1], ordered[index + 1][1], delta))
    for index in range(global_index, len(ordered) - 1):
        delta = ordered[index + 1][2] - ordered[index][2]
        if delta < -tolerance_hartree:
            reversals.append((ordered[index][1], ordered[index + 1][1], delta))
    return reversals


def fit_gxtb_eos(points: list[tuple[float, float, float, bool]]) -> dict[str, object]:
    """Fit a GXTB EOS only after enforcing a sampled single-well topology."""
    reversals = gxtb_topology_reversals(points)
    if not reversals:
        return fit_eos(points)
    ordered = sorted(
        (float(a), float(scale), float(energy))
        for a, scale, energy, ok in points
        if ok and energy is not None
    )
    global_minimum = min(ordered, key=lambda point: point[2])
    return {
        "a_eos_A": "",
        "energy_fit_hartree": "",
        "fit_status": "nonmonotonic_branch",
        "n_points": len(ordered),
        "grid_min_a_A": f"{global_minimum[0]:.10f}",
        "grid_min_scale": f"{global_minimum[1]:.5f}",
        "grid_min_energy_hartree": f"{global_minimum[2]:.12f}",
        "topology_reversal_count": len(reversals),
        "topology_max_reversal_hartree": f"{max(abs(item[2]) for item in reversals):.12f}",
    }


def merge_keyed_rows(
    path: Path,
    rows: list[dict[str, object]],
    key_fields: tuple[str, ...],
    *,
    sort_key: Callable[[dict[str, object]], object] | None = None,
) -> list[dict[str, object]]:
    """Replace records by scientific identity without dropping other k meshes.

    The original single-EOS workflow replaced every row of a selected method.
    A genuine k-convergence campaign evaluates several independent EOS meshes,
    so method-only replacement would silently discard the preceding mesh.
    """
    new_keys = {
        tuple(str(row.get(field, "")) for field in key_fields) for row in rows
    }
    preserved = [
        dict(row)
        for row in base.read_csv(path)
        if tuple(str(row.get(field, "")) for field in key_fields) not in new_keys
    ]
    merged: list[dict[str, object]] = preserved + rows
    if sort_key is not None:
        merged.sort(key=sort_key)
    base.write_csv(path, merged)
    return merged


def make_eos_table(
    mesh: str,
    scales: tuple[float, ...],
    methods: tuple[str, ...] = base.METHODS,
    classifications: dict[tuple[str, str, float], dict[str, str]] | None = None,
    campaign_fingerprint: dict[str, object] | None = None,
    solids: tuple[str, ...] = base.LC10_PAPER_SOLIDS,
) -> list[dict[str, object]]:
    classifications = classifications or {}
    rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []
    branch_rows: list[dict[str, object]] = []
    for ref in base.REFERENCES:
        if ref.solid not in solids:
            continue
        for method in methods:
            requested_scales = scales_for(ref.solid, method, scales)
            all_points = load_eos_points(
                ref, method, mesh, requested_scales, campaign_fingerprint
            )
            legacy_collapsed_scales = (
                charge_collapsed_scales(all_points) if method == "GFN2" else set()
            )
            local_candidates = gxtb_branch_candidates(all_points) if method == "GXTB" else {}
            cross_mesh_candidates: dict[float, float] = {}
            if method == "GXTB":
                previous_mesh = previous_cubic_mesh(mesh)
                previous_root = (
                    ROOT / "runs" / "eos" / method / ref.solid / previous_mesh
                    if previous_mesh is not None
                    else None
                )
                if previous_root is not None and previous_root.is_dir():
                    previous_points = load_eos_points(
                        ref,
                        method,
                        previous_mesh,
                        requested_scales,
                        campaign_fingerprint,
                    )
                    cross_mesh_candidates = gxtb_cross_mesh_branch_candidates(
                        all_points, previous_points
                    )
            candidates = {
                scale: max(
                    local_candidates.get(scale, 0.0),
                    cross_mesh_candidates.get(scale, 0.0),
                )
                for scale in set(local_candidates) | set(cross_mesh_candidates)
            }
            explicit = {
                scale: classifications[(ref.solid, mesh, scale)]
                for scale in (round(value, 5) for value in requested_scales)
                if (ref.solid, mesh, scale) in classifications
            }
            explicit_excluded_scales = {
                scale for scale, entry in explicit.items() if entry["action"] == "exclude"
            }
            excluded_scales = legacy_collapsed_scales | explicit_excluded_scales
            points = [
                (a, scale, energy, ok)
                for a, scale, energy, ok in all_points
                if energy is not None and ok and round(scale, 5) not in excluded_scales
            ]
            fit = fit_gxtb_eos(points) if method == "GXTB" else fit_eos(points)
            rows.append(
                {
                    "solid": ref.solid,
                    "structure": ref.structure,
                    "method": method,
                    "eos_mesh": mesh,
                    "a_exp_A": ref.a_exp,
                    "n_requested": len(requested_scales),
                    "n_completed": len(points),
                    "n_converged_raw": sum(ok and energy is not None for _, _, energy, ok in all_points),
                    "n_charge_collapsed": len(legacy_collapsed_scales)
                    + sum(
                        entry["classification"] == "charge_collapse"
                        and entry["action"] == "exclude"
                        for entry in explicit.values()
                    ),
                    "n_explicit_excluded": len(explicit_excluded_scales),
                    "n_unresolved_branch_candidates": sum(
                        scale not in explicit for scale in candidates
                    ),
                    **fit,
                }
            )
            for a, scale, energy, ok in all_points:
                normalized_scale = round(scale, 5)
                project = eos_project(ref.solid, method, mesh, scale)
                output = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale, method) / f"{project}.out"
                entry = explicit.get(normalized_scale)
                candidate_residual = candidates.get(normalized_scale)
                local_candidate_residual = local_candidates.get(normalized_scale)
                cross_mesh_candidate_residual = cross_mesh_candidates.get(normalized_scale)
                excluded = normalized_scale in excluded_scales
                if not ok:
                    if entry is not None and entry["action"] == "exclude":
                        diagnostic = entry["classification"]
                        resolution = "explicit_failure_classification"
                    else:
                        diagnostic = "scf_failure"
                        resolution = "failed_job"
                elif energy is None:
                    diagnostic = "missing_energy"
                    resolution = "failed_job"
                elif entry is not None:
                    diagnostic = entry["classification"]
                    resolution = (
                        "explicit_exclusion" if entry["action"] == "exclude" else "explicit_waiver"
                    )
                elif candidate_residual is not None:
                    diagnostic = (
                        "cross_mesh_scc_root_candidate"
                        if cross_mesh_candidate_residual is not None
                        else "branch_discontinuity_candidate"
                    )
                    resolution = "unresolved_candidate"
                elif normalized_scale in legacy_collapsed_scales:
                    diagnostic = "charge_collapse"
                    resolution = "legacy_automatic_filter"
                else:
                    diagnostic = ""
                    resolution = ""
                point_rows.append(
                    {
                        "solid": ref.solid,
                        "method": method,
                        "mesh": mesh,
                        "scale": f"{scale:.5f}",
                        "a_A": f"{a:.10f}",
                        "energy_hartree": f"{energy:.12f}" if energy is not None else "",
                        "completed": ok,
                        "valid_for_eos": ok and energy is not None and not excluded,
                        "diagnostic": diagnostic,
                        "classification_resolution": resolution,
                        "classification_rationale": entry["rationale"] if entry is not None else "",
                        "scf_strategy": read_strategy(output),
                    }
                )
                if method == "GXTB" and (candidate_residual is not None or entry is not None):
                    branch_rows.append(
                        {
                            "solid": ref.solid,
                            "method": method,
                            "mesh": mesh,
                            "scale": f"{scale:.5f}",
                            "automatic_candidate": candidate_residual is not None,
                            "interpolation_residual_hartree": (
                                f"{candidate_residual:.12f}" if candidate_residual is not None else ""
                            ),
                            "single_mesh_residual_hartree": (
                                f"{local_candidate_residual:.12f}"
                                if local_candidate_residual is not None
                                else ""
                            ),
                            "cross_mesh_shift_residual_hartree": (
                                f"{cross_mesh_candidate_residual:.12f}"
                                if cross_mesh_candidate_residual is not None
                                else ""
                            ),
                            "detection_source": (
                                "single_mesh+cross_mesh"
                                if local_candidate_residual is not None
                                and cross_mesh_candidate_residual is not None
                                else (
                                    "cross_mesh"
                                    if cross_mesh_candidate_residual is not None
                                    else "single_mesh"
                                )
                            ),
                            "classification": entry["classification"] if entry is not None else "",
                            "action": entry["action"] if entry is not None else "",
                            "rationale": entry["rationale"] if entry is not None else "",
                            "resolution": resolution,
                            "interpretation": (
                                "numerical branch candidate; not a physical failure until explicitly classified"
                            ),
                        }
                    )
    merge_keyed_rows(
        ROOT / "data" / "eos_points.csv",
        point_rows,
        ("solid", "method", "mesh", "scale"),
        sort_key=lambda row: (
            [ref.solid for ref in base.REFERENCES].index(str(row["solid"])),
            base.METHODS.index(str(row["method"])),
            int(str(row["mesh"])[1:]),
            float(row["scale"]),
        ),
    )
    merge_keyed_rows(
        ROOT / "data" / "eos_fits.csv",
        rows,
        ("solid", "method", "eos_mesh"),
        sort_key=lambda row: (
            [ref.solid for ref in base.REFERENCES].index(str(row["solid"])),
            base.METHODS.index(str(row["method"])),
            int(str(row["eos_mesh"])[1:]),
        ),
    )
    if "GXTB" in methods:
        branch_path = ROOT / "data" / "gxtb_eos_branch_diagnostics.csv"
        selected_solid_set = set(solids)
        preserved_branch_rows = [
            dict(row)
            for row in base.read_csv(branch_path)
            if not (
                row.get("method") == "GXTB"
                and row.get("solid") in selected_solid_set
                and row.get("mesh") == mesh
            )
        ]
        combined_branch_rows: list[dict[str, object]] = (
            preserved_branch_rows + branch_rows
        )
        combined_branch_rows.sort(
            key=lambda row: (
                [ref.solid for ref in base.REFERENCES].index(str(row["solid"])),
                int(str(row["mesh"])[1:]),
                float(row["scale"]),
            )
        )
        base.write_csv(branch_path, combined_branch_rows)
        candidate_path = ROOT / "data" / "gxtb_eos_classification_candidates.json"
        candidate_path.parent.mkdir(parents=True, exist_ok=True)
        candidate_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "instructions": (
                        "Copy reviewed entries to gxtb_eos_classifications.json and set action "
                        "to exclude or retain plus a nonempty classification and rationale."
                    ),
                    "entries": [
                        {
                            "solid": row["solid"],
                            "method": "GXTB",
                            "mesh": row["mesh"],
                            "scale": float(row["scale"]),
                            "classification": "",
                            "action": "",
                            "rationale": "",
                        }
                        for row in combined_branch_rows
                        if row.get("resolution") == "unresolved_candidate"
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    return rows


def suggested_adaptive_scales(row: dict[str, object], requested: tuple[float, ...]) -> tuple[float, ...]:
    ordered = sorted(set(requested))
    grid_text = str(row.get("grid_min_scale", ""))
    grid = float(grid_text) if grid_text else 1.0
    if not ordered:
        return (0.98, 1.00, 1.02)
    index = min(range(len(ordered)), key=lambda idx: abs(ordered[idx] - grid))
    suggestions: set[float] = set()
    if index > 0:
        suggestions.add((ordered[index - 1] + ordered[index]) / 2.0)
    else:
        step = ordered[1] - ordered[0] if len(ordered) > 1 else 0.04
        suggestions.add(max(0.01, ordered[0] - step / 2.0))
    if index + 1 < len(ordered):
        suggestions.add((ordered[index] + ordered[index + 1]) / 2.0)
    else:
        step = ordered[-1] - ordered[-2] if len(ordered) > 1 else 0.04
        suggestions.add(ordered[-1] + step / 2.0)
    # A second pair close to the apparent minimum helps distinguish a narrow
    # smooth well from an SCC branch switch.
    local_step = min(
        (abs(value - grid) for value in ordered if abs(value - grid) > 1.0e-8),
        default=0.04,
    )
    suggestions.update({max(0.01, grid - local_step / 3.0), grid + local_step / 3.0})
    return tuple(
        sorted(round(value, 5) for value in suggestions if round(value, 5) not in ordered)
    )


def write_gxtb_adaptive_followup(
    fits: list[dict[str, object]], scales: tuple[float, ...]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for fit in fits:
        if fit.get("method") != "GXTB" or fit.get("fit_status") == "quadratic":
            continue
        solid = str(fit["solid"])
        requested = scales_for(solid, "GXTB", scales)
        suggestions = suggested_adaptive_scales(fit, requested)
        command_args = " ".join(f"--adaptive-scale {solid}={value:.5f}" for value in suggestions)
        rows.append(
            {
                "solid": solid,
                "method": "GXTB",
                "fit_status": fit.get("fit_status", ""),
                "n_requested": fit.get("n_requested", ""),
                "n_completed": fit.get("n_completed", ""),
                "grid_min_scale": fit.get("grid_min_scale", ""),
                "classification": "requires_adaptive_investigation",
                "interpretation": "numerical EOS/branch investigation required; not a physical failure",
                "suggested_scales": ";".join(f"{value:.5f}" for value in suggestions),
                "adaptive_investigated": bool(adaptive_scales_only(solid, "GXTB")),
                "rerun_arguments": command_args,
            }
        )
    base.write_csv(ROOT / "data" / "gxtb_adaptive_followup.csv", rows)
    lines = [
        "# GXTB adaptive EOS follow-up",
        "",
        "These entries are numerical EOS/branch diagnostics, not physical failures.",
        "Review `gxtb_eos_branch_diagnostics.csv`, run the suggested scales, and classify/waive",
        "every branch candidate explicitly in `gxtb_eos_classifications.json`.",
        "",
    ]
    if not rows:
        lines.append("No invalid GXTB EOS fits were detected.")
    else:
        lines += [
            "| solid | fit status | adaptive attempted | suggested scales |",
            "|---|---|---:|---|",
        ]
        for row in rows:
            lines.append(
                f"| {row['solid']} | {row['fit_status']} | {row['adaptive_investigated']} | "
                f"{row['suggested_scales']} |"
            )
    path = ROOT / "data" / "gxtb_adaptive_followup.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return rows


def enforce_gxtb_coverage(
    fits: list[dict[str, object]],
    followup: list[dict[str, object]],
    *,
    allow_reduced_coverage: bool,
    minimum_valid_fits: int = MINIMUM_REDUCED_GXTB_FITS,
) -> None:
    gxtb_fits = [row for row in fits if row.get("method") == "GXTB"]
    if not gxtb_fits:
        return
    branch_path = ROOT / "data" / "gxtb_eos_branch_diagnostics.csv"
    branch_rows = base.read_csv(branch_path)
    unresolved = [row for row in branch_rows if row.get("resolution") == "unresolved_candidate"]
    if unresolved:
        labels = ", ".join(f"{row['solid']}@{row['scale']}" for row in unresolved)
        raise RuntimeError(
            "Unresolved GXTB SCC-branch candidates require an explicit exclude/retain entry "
            f"with rationale in {gxtb_classification_manifest_path()}: {labels}"
        )
    valid = [row for row in gxtb_fits if row.get("fit_status") == "quadratic"]
    invalid = [row for row in gxtb_fits if row.get("fit_status") != "quadratic"]
    if not invalid:
        return
    labels = ", ".join(f"{row['solid']}={row['fit_status']}" for row in invalid)
    if not allow_reduced_coverage:
        raise RuntimeError(
            "Invalid GXTB EOS curves require adaptive follow-up before final single points: "
            + labels
            + ". See data/gxtb_adaptive_followup.csv; use --allow-reduced-coverage only after investigation."
        )
    if len(valid) < minimum_valid_fits:
        raise RuntimeError(
            f"Reduced GXTB coverage is not meaningful: {len(valid)} valid quadratic fits, "
            f"minimum {minimum_valid_fits}."
        )
    not_investigated = [row for row in followup if not truth(row["adaptive_investigated"])]
    if not_investigated:
        raise RuntimeError(
            "--allow-reduced-coverage requires adaptive scales for every invalid GXTB curve: "
            + ", ".join(str(row["solid"]) for row in not_investigated)
        )


def final_sp_specs(
    fits: list[dict[str, object]],
    meshes: list[str],
    solids: tuple[str, ...] | None = None,
) -> list[tuple[str, Path, Path, bool]]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    specs: list[tuple[str, Path, Path, bool]] = []
    for row in fits:
        if solids is not None and str(row.get("solid", "")) not in solids:
            continue
        a_text = row.get("a_eos_A", "")
        if a_text == "" or (row.get("method") == "GXTB" and row.get("fit_status") != "quadratic"):
            continue
        ref = refs[str(row["solid"])]
        method = str(row["method"])
        a = float(a_text)
        if method == "GXTB" and not row.get("eos_mesh"):
            raise ValueError(f"GXTB final input for {ref.solid} lacks its EOS mesh provenance")
        for mesh in meshes:
            project = final_project(ref.solid, method, mesh)
            run_dir = ROOT / "runs" / "eos_final_sp" / method / ref.solid / mesh
            inp = run_dir / f"{project}.inp"
            out = run_dir / f"{project}.out"
            base.write_file(inp, base.solid_input(ref, method, "ENERGY", mesh, a, project))
            if method == "GXTB":
                write_final_input_lineage(
                    inp,
                    row,
                    mesh,
                    valid=True,
                    reason="generated from the current quadratic EOS minimum",
                )
            specs.append((f"eos-final {method} {ref.solid} {mesh}", inp, out, False))
    return specs


GXTB_FIT_APPROVAL_FIELDS = (
    "solid",
    "method",
    "eos_mesh",
    "a_eos_A",
    "energy_fit_hartree",
    "fit_status",
    "fit_rmse_hartree",
    "n_requested",
    "n_completed",
    "n_converged_raw",
    "n_explicit_excluded",
    "n_unresolved_branch_candidates",
    "topology_reversal_count",
    "topology_max_reversal_hartree",
)


def gxtb_fit_approval_sha256(fits: list[dict[str, object]]) -> str:
    records = [
        {field: str(row.get(field, "")) for field in GXTB_FIT_APPROVAL_FIELDS}
        for row in fits
        if row.get("method") == "GXTB"
    ]
    records.sort(key=lambda row: row["solid"])
    payload = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def collect_results(
    fits: list[dict[str, object]],
    result_meshes: list[str],
    result_mesh: str,
    methods: tuple[str, ...] = base.METHODS,
    campaign_fingerprint: dict[str, object] | None = None,
    campaign_bind_all_methods: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    atom_e = base.atom_energies(
        methods,
        campaign_fingerprint,
        base.LC10_PAPER_ELEMENTS,
        campaign_bind_all_methods=campaign_bind_all_methods,
    )
    rows: list[dict[str, object]] = []
    for fit in fits:
        if not fit.get("a_eos_A") or (
            fit.get("method") == "GXTB" and fit.get("fit_status") != "quadratic"
        ):
            continue
        ref = refs[str(fit["solid"])]
        method = str(fit["method"])
        a_calc = float(fit["a_eos_A"])
        n_atoms = len(base.conventional_cell_atoms(ref))
        counts = base.atom_counts(ref)
        atom_sum = None
        if all((method, el) in atom_e for el in counts):
            atom_sum = sum(atom_e[(method, el)] * count for el, count in counts.items())
        for mesh in result_meshes:
            project = final_project(ref.solid, method, mesh)
            out = ROOT / "runs" / "eos_final_sp" / method / ref.solid / mesh / f"{project}.out"
            if method == "GXTB":
                if campaign_fingerprint is None:
                    raise ValueError("GXTB final collection requires a campaign fingerprint")
                issue = base.completed_stamp_campaign_issue(
                    out,
                    campaign_fingerprint,
                    executable_role="cp2k",
                    require_completed=base.output_ok(out),
                )
                if issue:
                    raise RuntimeError(issue)
            e_solid = base.parse_energy(out)
            ecoh = (atom_sum - e_solid) * base.HARTREE_TO_EV / n_atoms if atom_sum is not None and e_solid is not None else None
            rows.append(
                {
                    "solid": ref.solid,
                    "structure": ref.structure,
                    "method": method,
                    "eos_mesh": fit["eos_mesh"],
                    "energy_mesh": mesh,
                    "fit_status": fit["fit_status"],
                    "sp_completed": base.output_ok(out),
                    "sp_scf_strategy": read_strategy(out),
                    "a_calc_A": f"{a_calc:.8f}" if a_calc is not None else "",
                    "a_ref_exp_A": ref.a_exp,
                    "a_error_A": f"{(a_calc - ref.a_exp):.8f}" if a_calc is not None else "",
                    "a_abs_error_A": f"{abs(a_calc - ref.a_exp):.8f}" if a_calc is not None else "",
                    "ecoh_calc_eV_per_atom": f"{ecoh:.8f}" if ecoh is not None else "",
                    "ecoh_ref_exp_eV_per_atom": ref.ecoh_exp,
                    "ecoh_error_eV_per_atom": f"{(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                    "ecoh_abs_error_eV_per_atom": f"{abs(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                    "solid_energy_hartree": f"{e_solid:.12f}" if e_solid is not None else "",
                    "atom_reference_source": (
                        "save_tblite_cli"
                        if method == "GXTB" or campaign_bind_all_methods
                        else "tblite_cli"
                    ),
                    "a_HF_A": ref.a_hf,
                    "a_MP2_A": ref.a_mp2,
                    "a_SCS_MP2_A": ref.a_scs_mp2,
                    "a_SOS_MP2_A": ref.a_sos_mp2,
                    "ecoh_HF_eV_per_atom": ref.ecoh_hf,
                    "ecoh_MP2_eV_per_atom": ref.ecoh_mp2,
                    "ecoh_SCS_MP2_eV_per_atom": ref.ecoh_scs_mp2,
                    "ecoh_SOS_MP2_eV_per_atom": ref.ecoh_sos_mp2,
                }
            )
    rows = base.merge_method_rows(
        ROOT / "data" / "eos_results.csv",
        rows,
        methods,
        sort_key=lambda row: (
            [ref.solid for ref in base.REFERENCES].index(str(row["solid"])),
            base.METHODS.index(str(row["method"])),
            str(row["energy_mesh"]),
        ),
    )
    summary = summary_rows(rows, result_mesh, methods)
    lit_summary = literature_summary()
    base.write_csv(ROOT / "data" / "eos_summary.csv", summary + lit_summary)
    convergence = kpoint_convergence(rows)
    base.write_csv(ROOT / "data" / "eos_kpoint_convergence.csv", convergence)
    common_summary, common_solids = common_subset_summary(rows, result_mesh, methods)
    base.write_csv(ROOT / "data" / "eos_common_subset_summary.csv", common_summary)
    base.write_csv(
        ROOT / "data" / "eos_common_subset_systems.csv",
        [{"solid": solid} for solid in common_solids],
    )
    all_fits = [dict(row) for row in base.read_csv(ROOT / "data" / "eos_fits.csv")]
    write_markdown(rows, summary, lit_summary, convergence, result_mesh, all_fits, common_summary, common_solids)
    plot(rows, summary, lit_summary, result_mesh)
    plot_eos_diagnostics(all_fits)
    return rows, summary, convergence


def truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def summary_rows(
    rows: list[dict[str, object]],
    result_mesh: str,
    required_methods: tuple[str, ...] = (),
) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    available_methods = tuple(
        method
        for method in base.METHODS
        if method in required_methods or any(r["method"] == method for r in rows)
    )
    for method in available_methods:
        selected = [
            r
            for r in rows
            if r["method"] == method and r["energy_mesh"] == result_mesh and truth(r["sp_completed"])
        ]
        a_err = [float(r["a_error_A"]) for r in selected if r["a_error_A"] != ""]
        e_err = [float(r["ecoh_error_eV_per_atom"]) for r in selected if r["ecoh_error_eV_per_atom"] != ""]
        summary.append(
            {
                "source": "CP2K/save_tblite EOS" if method == "GXTB" else "CP2K/tblite EOS",
                "method": method,
                "n_complete": len(selected),
                "a_ME_A": mean(a_err),
                "a_MAE_A": mae(a_err),
                "a_RMSE_A": rmse(a_err),
                "ecoh_ME_eV_per_atom": mean(e_err),
                "ecoh_MAE_eV_per_atom": mae(e_err),
                "ecoh_RMSE_eV_per_atom": rmse(e_err),
            }
        )
    return summary


def common_subset_summary(
    rows: list[dict[str, object]],
    result_mesh: str,
    required_methods: tuple[str, ...] = (),
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    methods = tuple(
        method
        for method in base.METHODS
        if method in required_methods or any(r["method"] == method for r in rows)
    )
    valid_by_method: dict[str, set[str]] = {}
    by_key: dict[tuple[str, str], dict[str, object]] = {}
    for method in methods:
        selected = [
            row
            for row in rows
            if row["method"] == method
            and row["energy_mesh"] == result_mesh
            and truth(row["sp_completed"])
            and row.get("a_error_A", "") != ""
            and row.get("ecoh_error_eV_per_atom", "") != ""
        ]
        valid_by_method[method] = {str(row["solid"]) for row in selected}
        by_key.update({(str(row["solid"]), method): row for row in selected})
    common = set.intersection(*(valid_by_method[method] for method in methods)) if methods else set()
    ordered = tuple(ref.solid for ref in base.REFERENCES if ref.solid in common)
    summary: list[dict[str, object]] = []
    for method in methods:
        method_rows = [by_key[(solid, method)] for solid in ordered]
        a_err = [float(row["a_error_A"]) for row in method_rows]
        e_err = [float(row["ecoh_error_eV_per_atom"]) for row in method_rows]
        summary.append(
            {
                "method": method,
                "n_common": len(ordered),
                "systems": ";".join(ordered),
                "a_ME_A": mean(a_err),
                "a_MAE_A": mae(a_err),
                "a_RMSE_A": rmse(a_err),
                "ecoh_ME_eV_per_atom": mean(e_err),
                "ecoh_MAE_eV_per_atom": mae(e_err),
                "ecoh_RMSE_eV_per_atom": rmse(e_err),
            }
        )
    return summary, ordered


def literature_summary() -> list[dict[str, object]]:
    mapping = {
        "HF": ("a_hf", "ecoh_hf"),
        "MP2": ("a_mp2", "ecoh_mp2"),
        "SCS-MP2": ("a_scs_mp2", "ecoh_scs_mp2"),
        "SOS-MP2": ("a_sos_mp2", "ecoh_sos_mp2"),
    }
    rows: list[dict[str, object]] = []
    for name, (akey, ekey) in mapping.items():
        a_err = [getattr(ref, akey) - ref.a_exp for ref in base.REFERENCES]
        e_err = [getattr(ref, ekey) - ref.ecoh_exp for ref in base.REFERENCES]
        rows.append(
            {
                "source": "Goldzak2022",
                "method": name,
                "n_complete": len(base.REFERENCES),
                "a_ME_A": mean(a_err),
                "a_MAE_A": mae(a_err),
                "a_RMSE_A": rmse(a_err),
                "ecoh_ME_eV_per_atom": mean(e_err),
                "ecoh_MAE_eV_per_atom": mae(e_err),
                "ecoh_RMSE_eV_per_atom": rmse(e_err),
            }
        )
    return rows


def kpoint_convergence(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key = {(r["solid"], r["method"], r["energy_mesh"]): r for r in rows}
    conv: list[dict[str, object]] = []
    for ref in base.REFERENCES:
        for method in base.METHODS:
            dense = by_key.get((ref.solid, method, "k555"))
            if not dense or dense["ecoh_calc_eV_per_atom"] == "":
                continue
            e_dense = float(dense["ecoh_calc_eV_per_atom"])
            for mesh in ("k333", "k444"):
                row = by_key.get((ref.solid, method, mesh))
                if not row or row["ecoh_calc_eV_per_atom"] == "":
                    continue
                e = float(row["ecoh_calc_eV_per_atom"])
                conv.append(
                    {
                        "solid": ref.solid,
                        "method": method,
                        "mesh": mesh,
                        "reference_mesh": "k555",
                        "delta_ecoh_eV_per_atom": f"{(e - e_dense):.8f}",
                    }
                )
    return conv


def mean(values: list[float]) -> str:
    return f"{(sum(values) / len(values)):.8f}" if values else ""


def mae(values: list[float]) -> str:
    return f"{(sum(abs(v) for v in values) / len(values)):.8f}" if values else ""


def rmse(values: list[float]) -> str:
    return f"{math.sqrt(sum(v * v for v in values) / len(values)):.8f}" if values else ""


def write_markdown(
    rows: list[dict[str, object]],
    summary: list[dict[str, object]],
    lit_summary: list[dict[str, object]],
    convergence: list[dict[str, object]],
    result_mesh: str,
    fits: list[dict[str, object]],
    common_summary: list[dict[str, object]],
    common_solids: tuple[str, ...],
) -> None:
    selected = [r for r in rows if r["energy_mesh"] == result_mesh]
    by_key = {(r["solid"], r["method"]): r for r in selected}
    refs = {ref.solid: ref for ref in base.REFERENCES}
    lines = [
        f"# LC12 (Goldzak12) EOS results ({result_mesh} final energies)",
        "",
        "Solid energies use CP2K/tblite native Bloch k-points. Atomic references use the matching tblite CLI.",
        "",
        "## EOS fit coverage",
        "",
        "| method | valid fits | excluded EOS curves |",
        "|---|---:|---|",
    ]
    fit_labels = {
        "no_local_minimum": "no bracketed minimum",
        "poor_quadratic_fit": "discontinuous EOS",
        "nonmonotonic_branch": "nonmonotonic/multibranch EOS",
    }
    methods = tuple(method for method in base.METHODS if any(fit["method"] == method for fit in fits))
    for method in methods:
        method_fits = [fit for fit in fits if fit["method"] == method]
        excluded = [
            f"{fit['solid']} ({fit_labels.get(str(fit['fit_status']), str(fit['fit_status']))})"
            for fit in method_fits
            if fit.get("a_eos_A", "") == ""
            or (method == "GXTB" and fit.get("fit_status") != "quadratic")
        ]
        lines.append(f"| {method} | {len(method_fits) - len(excluded)}/{len(method_fits)} | {', '.join(excluded) or '-'} |")
    lines += [
        "",
        "## MAE comparison to experiment",
        "",
        "| source | method | n | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary + lit_summary:
        lines.append(
            f"| {row['source']} | {row['method']} | {row['n_complete']} | {row['a_ME_A']} | {row['a_MAE_A']} | "
            f"{row['ecoh_ME_eV_per_atom']} | {row['ecoh_MAE_eV_per_atom']} |"
        )
    header = ["solid", "a exp"]
    for method in methods:
        header += [f"a {method}", f"da {method}"]
    header.append("Ecoh exp")
    for method in methods:
        header += [f"Ecoh {method}", f"dE {method}"]
    lines += ["", "## Per-system GFN comparison", "", "| " + " | ".join(header) + " |"]
    lines.append("|---|" + "---:|" * (len(header) - 1))
    for ref in base.REFERENCES:
        cells = [ref.solid, f"{ref.a_exp:.3f}"]
        for method in methods:
            row = by_key.get((ref.solid, method), {})
            cells += [base.fmt(row.get("a_calc_A"), 4), base.fmt(row.get("a_error_A"), 4)]
        cells.append(f"{ref.ecoh_exp:.2f}")
        for method in methods:
            row = by_key.get((ref.solid, method), {})
            cells += [
                base.fmt(row.get("ecoh_calc_eV_per_atom"), 3),
                base.fmt(row.get("ecoh_error_eV_per_atom"), 3),
            ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += [
        "",
        "## Common valid GFN subset",
        "",
        "Systems: " + (", ".join(common_solids) if common_solids else "none"),
        "",
        "| method | n common | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in common_summary:
        lines.append(
            f"| {row['method']} | {row['n_common']} | {row['a_ME_A']} | {row['a_MAE_A']} | "
            f"{row['ecoh_ME_eV_per_atom']} | {row['ecoh_MAE_eV_per_atom']} |"
        )
    lines += [
        "",
        "## k-point convergence of cohesive energies",
        "",
        "| method | mesh vs k555 | mean abs delta (eV/atom) | max abs delta (eV/atom) |",
        "|---|---|---:|---:|",
    ]
    for method in methods:
        for mesh in ("k333", "k444"):
            vals = [abs(float(r["delta_ecoh_eV_per_atom"])) for r in convergence if r["method"] == method and r["mesh"] == mesh]
            lines.append(f"| {method} | {mesh} | {sum(vals) / len(vals):.6f} | {max(vals):.6f} |" if vals else f"| {method} | {mesh} |  |  |")
    (ROOT / "data" / "eos_results.md").write_text("\n".join(lines) + "\n")


def plot(rows: list[dict[str, object]], summary: list[dict[str, object]], lit_summary: list[dict[str, object]], result_mesh: str) -> None:
    selected = [r for r in rows if r["energy_mesh"] == result_mesh and truth(r["sp_completed"])]
    solids = [ref.solid for ref in base.REFERENCES]
    x = np.arange(len(solids))
    methods = tuple(method for method in base.METHODS if any(r["method"] == method for r in selected))
    width = min(0.8 / max(len(methods), 1), 0.36)
    colors = base.METHOD_COLORS
    if selected:
        for key, ylabel, name in [
            ("a_error_A", "lattice-constant error (A)", "goldzak12_eos_lattice_errors"),
            ("ecoh_error_eV_per_atom", "cohesive-energy error (eV/atom)", "goldzak12_eos_cohesive_errors"),
        ]:
            fig, ax = plt.subplots(figsize=(10.5, 4.6))
            for i, method in enumerate(methods):
                vals = []
                for solid in solids:
                    row = next((r for r in selected if r["solid"] == solid and r["method"] == method), None)
                    vals.append(float(row[key]) if row and row[key] != "" else np.nan)
                positions = x + (i - (len(methods) - 1) / 2.0) * width
                missing_positions = [position for position, value in zip(positions, vals) if np.isnan(value)]
                ax.bar(positions, vals, width, label=method, color=colors[method])
                for position in missing_positions:
                    ax.annotate("n/a", (position, 0.0), xytext=(0, 5), textcoords="offset points", ha="center", va="bottom", color="#666666", fontsize=8)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(solids, rotation=45, ha="right")
            ax.set_ylabel(ylabel)
            ax.set_title(f"LC12 (Goldzak12) EOS CP2K/tblite native-Bloch {result_mesh}")
            ax.legend(frameon=False)
            ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
            fig.tight_layout()
            out = ROOT / "figures" / name
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out.with_suffix(".png"), dpi=220)
            fig.savefig(out.with_suffix(".pdf"))
            plt.close(fig)

    plottable_summary = [
        row
        for row in summary
        if row["a_MAE_A"] != "" and row["ecoh_MAE_eV_per_atom"] != ""
    ]
    labels = [f"{r['method']}\n(n={r['n_complete']})" for r in lit_summary + plottable_summary]
    a_mae = [float(r["a_MAE_A"]) for r in lit_summary + plottable_summary]
    e_mae = [float(r["ecoh_MAE_eV_per_atom"]) for r in lit_summary + plottable_summary]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    method_colors = [base.METHOD_COLORS.get(str(row["method"]), "#999999") for row in plottable_summary]
    axes[0].bar(labels, a_mae, color=["#72B7B2"] * len(lit_summary) + method_colors)
    axes[1].bar(labels, e_mae, color=["#72B7B2"] * len(lit_summary) + method_colors)
    axes[0].set_ylabel("MAE a (A)")
    axes[1].set_ylabel("MAE Ecoh (eV/atom)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    fig.suptitle("LC12 (Goldzak12) comparison to zero-point corrected experiment")
    fig.tight_layout()
    out = ROOT / "figures" / "goldzak12_eos_mae_comparison"
    fig.savefig(out.with_suffix(".png"), dpi=220)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_eos_diagnostics(fits: list[dict[str, object]]) -> None:
    points_path = ROOT / "data" / "eos_points.csv"
    if not points_path.exists():
        return
    with points_path.open(newline="") as handle:
        points = list(csv.DictReader(handle))
    invalid = [
        fit
        for fit in fits
        if fit.get("a_eos_A", "") == ""
        or (fit.get("method") == "GXTB" and fit.get("fit_status") != "quadratic")
    ]
    if not invalid:
        return
    ncols = min(3, len(invalid))
    nrows = math.ceil(len(invalid) / ncols)
    fig, axes_array = plt.subplots(nrows, ncols, figsize=(5.4 * ncols, 4.2 * nrows), squeeze=False)
    axes = list(axes_array.flat)
    for ax, fit in zip(axes, invalid):
        solid = str(fit["solid"])
        method = str(fit["method"])
        selected = sorted(
            (row for row in points if row["solid"] == solid and row["method"] == method),
            key=lambda row: float(row["scale"]),
        )
        completed = [
            row
            for row in selected
            if row.get("valid_for_eos", row["completed"]) == "True" and row["energy_hartree"] != ""
        ]
        charge_collapsed = [row for row in selected if row.get("diagnostic") == "charge_collapse"]
        unstable = [
            row
            for row in selected
            if row["completed"] != "True" and row.get("diagnostic") != "charge_collapse"
        ]
        if not completed:
            ax.set_title(f"{method}/{solid}: no converged EOS points")
            ax.set_axis_off()
            continue
        energy_min = min(float(row["energy_hartree"]) for row in completed)
        scales = [float(row["scale"]) for row in completed]
        relative = [
            (float(row["energy_hartree"]) - energy_min) * base.HARTREE_TO_EV / 8.0 for row in completed
        ]
        color = base.METHOD_COLORS.get(method, "#666666")
        ax.plot(scales, relative, color=color, linewidth=1.2, alpha=0.75)
        ax.scatter(scales, relative, color=color, s=38, label="converged")
        marker_height = max(relative) * 1.08 if max(relative) > 0 else 1.0
        for row in unstable:
            scale = float(row["scale"])
            ax.axvline(scale, color="#888888", linewidth=0.8, linestyle="--", alpha=0.6)
            ax.scatter(scale, marker_height, color="#666666", marker="x", s=42, label="unstable SCC branch")
        for row in charge_collapsed:
            scale = float(row["scale"])
            ax.axvline(scale, color="#B91C1C", linewidth=0.8, linestyle=":", alpha=0.7)
            ax.scatter(scale, marker_height, color="#B91C1C", marker="x", s=48, label="charge-collapsed SCC solution")
        labels = {
            "no_local_minimum": "no bracketed minimum",
            "poor_quadratic_fit": "discontinuous EOS",
            "insufficient_points": "insufficient points",
            "nonmonotonic_branch": "nonmonotonic/multibranch EOS",
        }
        label = labels.get(str(fit["fit_status"]), str(fit["fit_status"]))
        ax.set_title(f"{method}/{solid}: {label}")
        ax.set_xlabel("lattice scale (a / experimental a)")
        ax.grid(color="#d0d0d0", linewidth=0.6, alpha=0.7)
    for ax in axes[len(invalid) :]:
        ax.set_axis_off()
    axes[0].set_ylabel("relative energy (eV/atom)")
    handles, labels = axes[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    axes[0].legend(unique.values(), unique.keys(), frameon=False)
    fig.suptitle("LC12 invalid-EOS diagnostics (native-Bloch k444)")
    fig.tight_layout()
    output = ROOT / "figures" / "goldzak12_eos_diagnostics"
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def add_cli_adaptive_scales(
    specs: list[str], methods: tuple[str, ...], parser: argparse.ArgumentParser
) -> None:
    if specs and len(methods) != 1:
        parser.error("--adaptive-scale requires exactly one selected --method")
    refs = {ref.solid for ref in base.REFERENCES}
    for spec in specs:
        try:
            solid, value = spec.split("=", 1)
            scale = float(value)
        except ValueError:
            parser.error(f"Bad --adaptive-scale {spec!r}; expected SOLID=SCALE")
        if solid not in refs:
            parser.error(f"Unknown --adaptive-scale solid {solid!r}")
        if scale <= 0.0:
            parser.error("Adaptive EOS scales must be positive")
        key = (solid, methods[0])
        ADAPTIVE_SCALES[key] = tuple(sorted(set(ADAPTIVE_SCALES.get(key, ())) | {scale}))


def final_stage_is_explicitly_approved(
    methods: tuple[str, ...],
    *,
    stop_after_eos: bool,
    fit_only: bool,
    approve_fits: bool,
) -> bool:
    if stop_after_eos or fit_only:
        return False
    return "GXTB" not in methods or approve_fits


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path)
    parser.add_argument("--tblite", type=Path, default=base.DEFAULT_TBLITE)
    parser.add_argument("--save-tblite", type=Path)
    parser.add_argument(
        "--campaign-manifest",
        type=Path,
        default=base.DEFAULT_GXTB_CAMPAIGN_MANIFEST,
        help="central frozen GXTB build manifest (source of executable/library paths and hashes)",
    )
    parser.add_argument(
        "--campaign-manifest-sha256",
        help="required external hash pin for MPI/affinity production execution",
    )
    parser.add_argument(
        "--cp2k-library",
        type=Path,
        help="optional exact libcp2k override; must match the campaign manifest",
    )
    parser.add_argument(
        "--save-tblite-library",
        type=Path,
        help="optional libtblite.a override; must match the campaign manifest",
    )
    parser.add_argument("--cp2k-source", type=Path, default=base.DEFAULT_CP2K_SOURCE)
    parser.add_argument("--tblite-source", type=Path, default=base.DEFAULT_TBLITE_SOURCE)
    parser.add_argument("--save-tblite-source", type=Path, default=base.DEFAULT_SAVE_TBLITE_SOURCE)
    parser.add_argument(
        "--method",
        action="append",
        choices=base.METHODS,
        help="method to run; repeat as needed (default: GFN1 and GFN2)",
    )
    parser.add_argument(
        "--solid",
        action="append",
        choices=tuple(ref.solid for ref in base.REFERENCES),
        help=(
            "solid to execute; repeat as needed (default: the fixed ten-system "
            "paper benchmark; LiH/MgO remain diagnostics only)"
        ),
    )
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mpi-ranks-per-job", type=int, default=1)
    parser.add_argument("--mpi-launcher", type=Path)
    parser.add_argument(
        "--mpi-launcher-arg",
        action="append",
        default=[],
        help=(
            "reserved compatibility option; production execution rejects every "
            "user-supplied MPI launcher argument"
        ),
    )
    parser.add_argument(
        "--pe-list",
        action="append",
        default=[],
        help=(
            "literal ordered CPU list (for example 96,97,98,99); repeat exactly "
            "once per --jobs worker"
        ),
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--stop-after-eos",
        action="store_true",
        help="run/collect EOS points and fits, but never launch final single points",
    )
    parser.add_argument(
        "--fit-only",
        action="store_true",
        help="collect already stamped EOS outputs and refresh fits without launching jobs",
    )
    parser.add_argument(
        "--approve-fits",
        action="store_true",
        help="explicitly approve the current GXTB fit fingerprint and launch final single points",
    )
    parser.add_argument("--eos-mesh", default="k444")
    parser.add_argument("--energy-mesh", action="append", default=[])
    parser.add_argument("--result-mesh", default=DEFAULT_RESULT_MESH)
    parser.add_argument("--scale", type=float, action="append", default=[])
    parser.add_argument(
        "--adaptive-scale",
        action="append",
        default=[],
        metavar="SOLID=SCALE",
        help="add a targeted EOS scale for the single selected method",
    )
    parser.add_argument(
        "--classification-manifest",
        type=Path,
        default=gxtb_classification_manifest_path(),
        help="reviewed per-point GXTB branch exclusions/waivers with rationale",
    )
    parser.add_argument(
        "--allow-reduced-coverage",
        action="store_true",
        help="deprecated and rejected: the paper benchmark requires exact LC10 coverage",
    )
    parser.add_argument(
        "--minimum-valid-fits",
        type=int,
        default=MINIMUM_REDUCED_GXTB_FITS,
        help="minimum quadratic GXTB fits when reduced coverage is explicitly allowed",
    )
    parser.add_argument(
        "--prune-transients",
        action="store_true",
        help="after validation, remove only large restart/WFN transients below GXTB run trees",
    )
    args = parser.parse_args()

    scales = tuple(args.scale) if args.scale else DEFAULT_SCALES
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]
    methods = base.selected_methods(args.method)
    if args.allow_reduced_coverage:
        parser.error(
            "the publication benchmark has fixed LC10 coverage; reduced coverage "
            "is not supported"
        )
    if args.jobs < 1 or args.threads < 1 or args.mpi_ranks_per_job < 1:
        parser.error("--jobs, --threads, and --mpi-ranks-per-job must be positive")
    if args.solid and len(args.solid) != len(set(args.solid)):
        parser.error("duplicate --solid selections are not allowed")
    selected_solids = tuple(args.solid or base.LC10_PAPER_SOLIDS)
    exact_lc10_scope = (
        len(selected_solids) == len(base.LC10_PAPER_SOLIDS)
        and set(selected_solids) == set(base.LC10_PAPER_SOLIDS)
    )
    if "GXTB" in methods and not exact_lc10_scope:
        parser.error(
            "the GXTB publication runner is restricted to the exact fixed LC10 set: "
            + ", ".join(base.LC10_PAPER_SOLIDS)
        )
    if args.fit_only and (args.force or args.approve_fits):
        parser.error("--fit-only cannot be combined with --force or --approve-fits")
    if args.stop_after_eos and args.approve_fits:
        parser.error("--stop-after-eos cannot be combined with --approve-fits")
    if not 1 <= args.minimum_valid_fits <= len(base.REFERENCES):
        parser.error(f"--minimum-valid-fits must be between 1 and {len(base.REFERENCES)}")
    campaign_fingerprint: dict[str, object] | None = None
    if "GXTB" in methods:
        try:
            campaign_fingerprint, campaign_paths = base.validated_gxtb_campaign_from_manifest(
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
            base.require_git_ancestor(
                args.cp2k_source,
                base.REQUIRED_CP2K_POST5582_ANCESTOR,
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
    else:
        args.cp2k = args.cp2k or base.DEFAULT_CP2K
        args.save_tblite = args.save_tblite or base.DEFAULT_SAVE_TBLITE
    execution_requested = bool(
        args.mpi_launcher
        or args.mpi_launcher_arg
        or args.pe_list
        or args.mpi_ranks_per_job != 1
    )
    execution_pool: execution.ExecutionPool | None = None
    if execution_requested:
        if methods != ("GXTB",):
            parser.error("MPI/affinity execution records are currently restricted to GXTB-only runs")
        if args.mpi_launcher is None:
            parser.error("--mpi-launcher is required with MPI/affinity execution")
        if not args.campaign_manifest_sha256 or not re.fullmatch(
            r"[0-9a-f]{64}", args.campaign_manifest_sha256
        ):
            parser.error(
                "MPI/affinity execution requires --campaign-manifest-sha256 as 64 lowercase hex digits"
            )
        observed_manifest_sha = base.sha256(args.campaign_manifest.resolve(strict=True))
        if observed_manifest_sha != args.campaign_manifest_sha256:
            parser.error(
                "campaign manifest hash pin mismatch: expected "
                f"{args.campaign_manifest_sha256}, observed {observed_manifest_sha}"
            )
        try:
            execution_pool = execution.ExecutionPool(
                concurrent_jobs=args.jobs,
                mpi_ranks_per_job=args.mpi_ranks_per_job,
                threads_per_rank=args.threads,
                mpi_launcher=args.mpi_launcher,
                mpi_launcher_args=args.mpi_launcher_arg,
                pe_lists=args.pe_list,
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
    restore_gxtb_scale_manifest(args.eos_mesh, methods)
    add_cli_adaptive_scales(args.adaptive_scale, methods, parser)
    if args.result_mesh not in energy_meshes:
        parser.error(f"--result-mesh {args.result_mesh} must also be supplied as --energy-mesh")
    try:
        classifications = load_gxtb_classifications(args.classification_manifest)
    except (ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    scale_manifest = write_gxtb_scale_manifest(
        args.eos_mesh, scales, methods, selected_solids
    )
    protocol = {
            "benchmark": "LC10 (fixed Goldzak12 subset)",
            "methods": methods,
            "selected_solids": selected_solids,
            "paper_systems": base.LC10_PAPER_SOLIDS,
            "diagnostic_only_systems": base.LC10_DIAGNOSTIC_ONLY_SOLIDS,
            "exact_lc10_scope": exact_lc10_scope,
            "cell_protocol": "cubic equation of state",
            "eos_mesh": args.eos_mesh,
            "energy_meshes": energy_meshes,
            "result_mesh": args.result_mesh,
            "scales": scales,
            "adaptive_scales": {f"{solid}/{method}": values for (solid, method), values in ADAPTIVE_SCALES.items()},
            "gxtb_scale_manifest": scale_manifest,
            "gxtb_scale_manifest_sha256": (
                base.sha256(gxtb_scale_manifest_path()) if scale_manifest is not None else None
            ),
            "gxtb_classification_manifest": str(args.classification_manifest.resolve()),
            "gxtb_classification_manifest_sha256": (
                base.sha256(args.classification_manifest)
                if args.classification_manifest.is_file()
                else None
            ),
            "allow_reduced_coverage": args.allow_reduced_coverage,
            "minimum_valid_gxtb_fits": args.minimum_valid_fits,
            "kpoint_scheme": "CP2K native Bloch MACDONALD with full SPGLIB symmetry reduction",
            "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
            "legacy_gxtb_full_grid_policy": base.LEGACY_GXTB_FULL_GRID_POLICY,
            "gxtb_energy_stress_policy": base.GXTB_ENERGY_STRESS_POLICY,
            "gxtb_atom_scf_policy": base.GXTB_ATOM_SCF_POLICY,
            "final_input_lineage_schema": FINAL_INPUT_LINEAGE_SCHEMA,
            "fit_approval_required": "GXTB" in methods,
            "fit_approved": False,
            "approved_gxtb_fit_sha256": None,
            "smearing_temperature_K": 300.0,
            "reported_energy": "Total energy extrapolated to T->0",
            "tblite_accuracy": 0.05,
            "legacy_gfn2_charge_collapse_filter_hartree_below_curve_median": CHARGE_COLLAPSE_ENERGY_DROP_HARTREE,
            "catastrophic_charge_collapse_energy_hartree": CATASTROPHIC_CHARGE_COLLAPSE_ENERGY_HARTREE,
            "gxtb_branch_candidate_interpolation_residual_hartree": GXTB_BRANCH_DISCONTINUITY_HARTREE,
            "gxtb_cross_mesh_shift_floor_hartree": GXTB_CROSS_MESH_SHIFT_FLOOR_HARTREE,
            "gxtb_cross_mesh_shift_mad_multiplier": GXTB_CROSS_MESH_MAD_MULTIPLIER,
            "gxtb_branch_policy": "automatic candidates require explicit per-point exclusion or retain waiver with rationale",
            "gxtb_topology_tolerance_hartree_per_eight_atom_cell": GXTB_TOPOLOGY_TOLERANCE_HARTREE,
            "gxtb_topology_policy": "sampled energy must decrease to the global minimum and increase away from it",
            "default_scf_strategy": "native g-XTB FDIIS" if methods == ("GXTB",) else "tblite modified-Broyden defaults",
            "scf_retry_strategies": [
                "TBLITE_MIXER ITERATIONS 1200 MEMORY 1 DAMPING 0.05",
                "TBLITE_MIXER ITERATIONS 2400 MEMORY 1 DAMPING 0.01",
            ] if any(method in base.LEGACY_METHODS for method in methods) else [],
            "gxtb_scc_mixer": "SCC_MIXER TBLITE (save_tblite native FDIIS)" if "GXTB" in methods else None,
            "cp2k_density_mixing": "METHOD DIRECT_P_MIXING; ALPHA 0.2",
            "atom_reference": "matching CLI, --method gxtb, explicit 2S spin" if "GXTB" in methods else "matching tblite CLI",
            "conventional_cell_atoms": 8,
            "eos_failure_policy": (
                "all ten fixed paper systems require valid quadratic fits; LiH/MgO "
                "branch studies are diagnostics outside the publication statistic"
            ),
            "execution_provenance": (
                {
                    "separate_from_scientific_job_stamp": True,
                    "record_schema": execution.SCHEMA_VERSION,
                    "contract": execution_pool.contract,
                    "contract_sha256": execution_pool.contract_sha256,
                }
                if execution_pool is not None
                else None
            ),
            "required_cp2k_ancestor": base.REQUIRED_CP2K_POST5582_ANCESTOR,
        }
    legacy = tuple(method for method in methods if method in base.LEGACY_METHODS)

    def write_provenance() -> None:
        if legacy:
            base.write_build_provenance(
                args.cp2k, args.tblite, args.cp2k_source, args.tblite_source, protocol
            )
        if "GXTB" in methods:
            assert campaign_fingerprint is not None
            base.write_gxtb_build_provenance(
                args.cp2k,
                args.save_tblite,
                args.cp2k_source,
                args.save_tblite_source,
                protocol,
                campaign_fingerprint,
                args.campaign_manifest,
            )

    write_provenance()
    if legacy:
        base.setup_inputs(args.eos_mesh, energy_meshes, legacy)
    if not args.fit_only:
        base.run_tblite_atom_jobs(
            args.tblite,
            args.jobs,
            args.force,
            methods,
            args.save_tblite,
            campaign_fingerprint,
            base.LC10_PAPER_ELEMENTS,
        )
        try:
            run_jobs(
                eos_job_specs(args.eos_mesh, scales, methods, selected_solids),
                args.cp2k,
                args.jobs,
                args.threads,
                args.force,
                campaign_fingerprint=campaign_fingerprint,
                execution_pool=execution_pool,
            )
        except RuntimeError:
            # Preserve a complete diagnostic/classification record even though the
            # production invocation must remain nonzero for every failed job.
            failed_fits = make_eos_table(
                args.eos_mesh,
                scales,
                methods,
                classifications,
                campaign_fingerprint,
                selected_solids,
            )
            if "GXTB" in methods:
                invalidate_existing_gxtb_final_inputs(failed_fits, energy_meshes)
                write_gxtb_adaptive_followup(failed_fits, scales)
            raise
    fits = make_eos_table(
        args.eos_mesh,
        scales,
        methods,
        classifications,
        campaign_fingerprint,
        selected_solids,
    )
    if "GXTB" in methods:
        invalidate_existing_gxtb_final_inputs(fits, energy_meshes)
    followup = write_gxtb_adaptive_followup(fits, scales) if "GXTB" in methods else []
    current_fit_sha = gxtb_fit_approval_sha256(fits) if "GXTB" in methods else None
    protocol["current_gxtb_fit_sha256"] = current_fit_sha
    write_provenance()
    if not final_stage_is_explicitly_approved(
        methods,
        stop_after_eos=args.stop_after_eos,
        fit_only=args.fit_only,
        approve_fits=args.approve_fits,
    ):
        print(
            "EOS staging complete; final single points were not launched. "
            "Review eos_fits.csv and branch diagnostics, then rerun with --approve-fits."
        )
        return 0
    enforce_gxtb_coverage(
        fits,
        followup,
        allow_reduced_coverage=args.allow_reduced_coverage,
        minimum_valid_fits=args.minimum_valid_fits,
    )
    if "GXTB" in methods:
        protocol["fit_approved"] = True
        protocol["approved_gxtb_fit_sha256"] = current_fit_sha
        write_provenance()
    try:
        run_jobs(
            final_sp_specs(fits, energy_meshes, selected_solids),
            args.cp2k,
            args.jobs,
            args.threads,
            args.force,
            campaign_fingerprint=campaign_fingerprint,
            execution_pool=execution_pool,
        )
    except RuntimeError:
        collect_results(
            fits, energy_meshes, args.result_mesh, methods, campaign_fingerprint
        )
        raise
    collect_results(fits, energy_meshes, args.result_mesh, methods, campaign_fingerprint)
    if args.prune_transients and "GXTB" in methods:
        count, size = base.prune_gxtb_transients()
        print(f"Pruned {count} validated GXTB transient file(s), {size} byte(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
