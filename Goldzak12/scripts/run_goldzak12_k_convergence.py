#!/usr/bin/env python3
"""Run the fail-closed, identical-set LC10 adaptive k-convergence campaign.

This is the publication runner.  It deliberately keeps two different studies
separate:

* independent EOS fits and an equilibrium single point at every k mesh drive
  the per-system selection of both a0 and Ecoh;
* the historical k333/k444/k555 single points on the fixed k444 EOS geometry
  remain a diagnostic series and never enter the convergence decision.

One adjacent n^3 -> (n+1)^3 step is sufficient when *both* raw changes pass.
The reported value is always the denser (n+1)^3 value.  There is no RMS gate
and no requirement for a second consecutive passing interval.  There is no
scientific maximum mesh: unresolved tracks keep advancing by one mesh.  The
optional ``--maximum-mesh`` is only a technical resource guard and reaching it
is an explicit error, never convergence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import benchmark_execution as execution  # noqa: E402
import run_goldzak12_benchmark as base  # noqa: E402
import run_goldzak12_eos_benchmark as eos  # noqa: E402


METHODS = ("GFN1", "GFN2", "GXTB")
PAPER_SYSTEMS = base.LC10_PAPER_SOLIDS
PAPER_ELEMENTS = base.LC10_PAPER_ELEMENTS
INITIAL_MESH_NUMBERS = (3, 4, 5)
LATTICE_THRESHOLD_A = 0.001
ECOH_THRESHOLD_KJMOL_PER_ATOM = 0.05
KJMOL_PER_EV = 96.48533212331002
ECOH_THRESHOLD_EV_PER_ATOM = ECOH_THRESHOLD_KJMOL_PER_ATOM / KJMOL_PER_EV
FIXED_GEOMETRY_EOS_MESH = "k444"
FIXED_GEOMETRY_ENERGY_MESHES = ("k333", "k444", "k555")
EQUILIBRIUM_LINEAGE_SCHEMA = 1
CONVERGENCE_SCHEMA = 2

VALUES_NAME = "lc10_independent_eos_k_values.csv"
STEPS_NAME = "lc10_adaptive_k_steps.csv"
SELECTION_NAME = "lc10_adaptive_k_selection.csv"
CONVERGENCE_NAME = "lc10_adaptive_k_convergence.json"
SCALE_MANIFEST_NAME = "lc10_k_convergence_scale_manifest.json"


def selected_methods(requested: Iterable[str] | None) -> tuple[str, ...]:
    """Return a unique selection in the canonical publication order."""
    selected = tuple(requested or ())
    if not selected:
        return METHODS
    if len(selected) != len(set(selected)):
        raise ValueError("--method may select each method at most once")
    unknown = sorted(set(selected) - set(METHODS))
    if unknown:
        raise ValueError(f"unknown method selection: {', '.join(unknown)}")
    return tuple(method for method in METHODS if method in selected)


def mesh_name(number: int) -> str:
    if number < 1:
        raise ValueError("k-mesh number must be positive")
    return f"k{number}{number}{number}"


def mesh_number(mesh: str) -> int:
    match = re.fullmatch(r"k([1-9][0-9]*)\1\1", mesh)
    if match is None:
        raise ValueError(f"non-cubic k mesh {mesh!r}")
    return int(match.group(1))


def equilibrium_project(solid: str, method: str, mesh: str) -> str:
    return f"{solid}_{method}_independent_eos_minimum_{mesh}"


def equilibrium_paths(root: Path, solid: str, method: str, mesh: str) -> tuple[Path, Path, Path]:
    run_dir = root / "runs" / "eos_k_convergence" / method / solid / mesh
    project = equilibrium_project(solid, method, mesh)
    input_path = run_dir / f"{project}.inp"
    output_path = run_dir / f"{project}.out"
    return input_path, output_path, input_path.with_suffix(".inp.eos.json")


def fit_fingerprint(
    fits: Iterable[dict[str, object]],
    methods: tuple[str, ...] = METHODS,
) -> str:
    fields = (
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
    records = [
        {field: str(row.get(field, "")) for field in fields}
        for row in fits
        if row.get("solid") in PAPER_SYSTEMS and row.get("method") in methods
    ]
    records.sort(
        key=lambda row: (
            methods.index(row["method"]),
            row["solid"],
            mesh_number(row["eos_mesh"]),
        )
    )
    encoded = json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def equilibrium_specs(
    root: Path,
    fits: Iterable[dict[str, object]],
    methods: tuple[str, ...] = METHODS,
) -> list[tuple[str, Path, Path, bool]]:
    refs = {ref.solid: ref for ref in base.LC10_PAPER_REFERENCES}
    specs: list[tuple[str, Path, Path, bool]] = []
    for fit in fits:
        solid = str(fit.get("solid", ""))
        method = str(fit.get("method", ""))
        mesh = str(fit.get("eos_mesh", ""))
        if solid not in refs or method not in methods:
            continue
        if fit.get("fit_status") != "quadratic" or not str(fit.get("a_eos_A", "")):
            raise RuntimeError(f"invalid independent EOS fit: {method}/{solid}/{mesh}")
        if int(fit.get("n_unresolved_branch_candidates", 0) or 0) != 0:
            raise RuntimeError(f"unresolved EOS branch candidate: {method}/{solid}/{mesh}")
        a0 = float(fit["a_eos_A"])
        input_path, output_path, lineage_path = equilibrium_paths(root, solid, method, mesh)
        text = base.solid_input(
            refs[solid], method, "ENERGY", mesh, a0, equilibrium_project(solid, method, mesh)
        )
        base.write_file(input_path, text)
        lineage = {
            "schema_version": EQUILIBRIUM_LINEAGE_SCHEMA,
            "benchmark": "LC10 independent-EOS k convergence",
            "solid": solid,
            "method": method,
            "eos_mesh": mesh,
            "energy_mesh": mesh,
            "a_eos_A": str(fit["a_eos_A"]),
            "fit_status": str(fit["fit_status"]),
            "fit_record_sha256": hashlib.sha256(
                json.dumps(fit, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
            "input_sha256": base.sha256(input_path),
            "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
        }
        lineage_path.write_text(json.dumps(lineage, indent=2, sort_keys=True) + "\n")
        specs.append((f"k-eq {method} {solid} {mesh}", input_path, output_path, False))
    return specs


def collect_independent_values(
    root: Path,
    fits: Iterable[dict[str, object]],
    campaign: dict[str, object],
    methods: tuple[str, ...] = METHODS,
) -> list[dict[str, object]]:
    refs = {ref.solid: ref for ref in base.LC10_PAPER_REFERENCES}
    atom_energies = base.atom_energies(
        methods,
        campaign,
        PAPER_ELEMENTS,
        campaign_bind_all_methods=True,
    )
    rows: list[dict[str, object]] = []
    for fit in fits:
        solid = str(fit.get("solid", ""))
        method = str(fit.get("method", ""))
        mesh = str(fit.get("eos_mesh", ""))
        if solid not in refs or method not in methods:
            continue
        if fit.get("fit_status") != "quadratic" or not str(fit.get("a_eos_A", "")):
            raise RuntimeError(f"invalid independent EOS fit: {method}/{solid}/{mesh}")
        input_path, output_path, lineage_path = equilibrium_paths(root, solid, method, mesh)
        if not base.output_ok(output_path):
            raise RuntimeError(f"incomplete equilibrium single point: {method}/{solid}/{mesh}")
        issue = base.completed_stamp_campaign_issue(
            output_path,
            campaign,
            executable_role="cp2k",
            require_completed=True,
        )
        if issue:
            raise RuntimeError(issue)
        lineage = json.loads(lineage_path.read_text())
        expected_lineage = {
            "schema_version": EQUILIBRIUM_LINEAGE_SCHEMA,
            "solid": solid,
            "method": method,
            "eos_mesh": mesh,
            "energy_mesh": mesh,
            "a_eos_A": str(fit["a_eos_A"]),
            "fit_status": "quadratic",
            "input_sha256": base.sha256(input_path),
            "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
        }
        for field, expected in expected_lineage.items():
            if lineage.get(field) != expected:
                raise RuntimeError(
                    f"equilibrium input lineage {field} mismatch: {method}/{solid}/{mesh}"
                )
        solid_energy = base.parse_energy(output_path)
        if solid_energy is None:
            raise RuntimeError(f"missing equilibrium energy: {method}/{solid}/{mesh}")
        ref = refs[solid]
        counts = base.atom_counts(ref)
        try:
            atom_sum = sum(atom_energies[(method, element)] * count for element, count in counts.items())
        except KeyError as error:
            raise RuntimeError(f"missing atom reference for {method}/{solid}: {error}") from error
        n_atoms = len(base.conventional_cell_atoms(ref))
        ecoh = (atom_sum - solid_energy) * base.HARTREE_TO_EV / n_atoms
        rows.append(
            {
                "solid": solid,
                "structure": ref.structure,
                "method": method,
                "mesh": mesh,
                "mesh_n": mesh_number(mesh),
                "a0_A": f"{float(fit['a_eos_A']):.10f}",
                "solid_energy_hartree": f"{solid_energy:.12f}",
                "ecoh_eV_per_atom": f"{ecoh:.12f}",
                "fit_energy_hartree": str(fit.get("energy_fit_hartree", "")),
                "fit_rmse_hartree": str(fit.get("fit_rmse_hartree", "")),
                "fit_status": "quadratic",
                "n_eos_points": str(fit.get("n_completed", "")),
                "input_sha256": base.sha256(input_path),
                "output_sha256": base.sha256(output_path),
                "lineage_sha256": base.sha256(lineage_path),
                "campaign_stamp_sha256": base.sha256(base.job_stamp_path(output_path)),
                "value_source": "single_point_at_own_independent_eos_minimum",
                "atom_reference_source": "save_tblite_cli",
            }
        )
    rows.sort(
        key=lambda row: (
            methods.index(str(row["method"])),
            PAPER_SYSTEMS.index(str(row["solid"])),
            int(row["mesh_n"]),
        )
    )
    return rows


def assess_convergence(
    values: Iterable[dict[str, object]],
    *,
    methods: tuple[str, ...] = METHODS,
    require_initial_meshes: bool = True,
    maximum_mesh: int | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[tuple[str, str, int]]]:
    if maximum_mesh is not None and maximum_mesh < max(INITIAL_MESH_NUMBERS):
        raise ValueError(
            f"technical maximum mesh must be at least {max(INITIAL_MESH_NUMBERS)}"
        )
    by_key: dict[tuple[str, str, int], dict[str, object]] = {}
    for row in values:
        method = str(row.get("method", ""))
        solid = str(row.get("solid", ""))
        number = int(row.get("mesh_n", mesh_number(str(row.get("mesh", "")))))
        key = (method, solid, number)
        if key in by_key:
            raise ValueError(f"duplicate independent-EOS value: {key}")
        by_key[key] = row

    steps: list[dict[str, object]] = []
    selections: list[dict[str, object]] = []
    pending: list[tuple[str, str, int]] = []
    for method in methods:
        for solid in PAPER_SYSTEMS:
            numbers = sorted(number for m, s, number in by_key if (m, s) == (method, solid))
            if require_initial_meshes and not set(INITIAL_MESH_NUMBERS) <= set(numbers):
                missing = sorted(set(INITIAL_MESH_NUMBERS) - set(numbers))
                raise ValueError(f"missing initial independent EOS meshes for {method}/{solid}: {missing}")
            if not numbers:
                pending.append((method, solid, INITIAL_MESH_NUMBERS[0]))
                continue
            if numbers != list(range(numbers[0], numbers[-1] + 1)):
                raise ValueError(f"non-consecutive independent EOS meshes for {method}/{solid}: {numbers}")
            selected: dict[str, object] | None = None
            for coarse_n, dense_n in zip(numbers, numbers[1:]):
                coarse = by_key[(method, solid, coarse_n)]
                dense = by_key[(method, solid, dense_n)]
                delta_a = float(dense["a0_A"]) - float(coarse["a0_A"])
                delta_ecoh = float(dense["ecoh_eV_per_atom"]) - float(coarse["ecoh_eV_per_atom"])
                lattice_passed = abs(delta_a) <= LATTICE_THRESHOLD_A or math.isclose(
                    abs(delta_a), LATTICE_THRESHOLD_A, rel_tol=0.0, abs_tol=5.0e-11
                )
                ecoh_passed = abs(delta_ecoh) <= ECOH_THRESHOLD_EV_PER_ATOM or math.isclose(
                    abs(delta_ecoh),
                    ECOH_THRESHOLD_EV_PER_ATOM,
                    rel_tol=0.0,
                    abs_tol=5.0e-13,
                )
                both_passed = lattice_passed and ecoh_passed
                step = {
                    "solid": solid,
                    "method": method,
                    "coarse_mesh": mesh_name(coarse_n),
                    "dense_mesh": mesh_name(dense_n),
                    "coarse_a0_A": f"{float(coarse['a0_A']):.10f}",
                    "dense_a0_A": f"{float(dense['a0_A']):.10f}",
                    "delta_a0_dense_minus_coarse_A": f"{delta_a:.10f}",
                    "abs_delta_a0_A": f"{abs(delta_a):.10f}",
                    "lattice_threshold_A": f"{LATTICE_THRESHOLD_A:.10f}",
                    "lattice_passed": lattice_passed,
                    "coarse_ecoh_eV_per_atom": f"{float(coarse['ecoh_eV_per_atom']):.12f}",
                    "dense_ecoh_eV_per_atom": f"{float(dense['ecoh_eV_per_atom']):.12f}",
                    "delta_ecoh_dense_minus_coarse_eV_per_atom": f"{delta_ecoh:.12f}",
                    "abs_delta_ecoh_eV_per_atom": f"{abs(delta_ecoh):.12f}",
                    "ecoh_threshold_eV_per_atom": f"{ECOH_THRESHOLD_EV_PER_ATOM:.12f}",
                    "ecoh_threshold_kJmol_per_atom": f"{ECOH_THRESHOLD_KJMOL_PER_ATOM:.8f}",
                    "ecoh_passed": ecoh_passed,
                    "both_passed": both_passed,
                    "decision_rule": "one_consecutive_step_AND",
                }
                steps.append(step)
                if selected is None and both_passed:
                    selected = {
                        **dense,
                        "converged_from_mesh": mesh_name(coarse_n),
                        "selected_mesh": mesh_name(dense_n),
                        "selection_status": "converged",
                        "selection_rule": "earliest_single_passing_adjacent_step_take_denser_value",
                        "selected_step_abs_delta_a0_A": step["abs_delta_a0_A"],
                        "selected_step_abs_delta_ecoh_eV_per_atom": step[
                            "abs_delta_ecoh_eV_per_atom"
                        ],
                    }
            if selected is not None:
                selections.append(selected)
            elif maximum_mesh is None or numbers[-1] < maximum_mesh:
                pending.append((method, solid, numbers[-1] + 1))
            else:
                selections.append(
                    {
                        "solid": solid,
                        "method": method,
                        "selection_status": "technical_resource_guard_reached",
                        "selected_mesh": "",
                        "last_evaluated_mesh": mesh_name(numbers[-1]),
                        "technical_resource_guard_mesh": mesh_name(maximum_mesh),
                        "selection_rule": (
                            "resource_guard_reached_without_scientific_convergence"
                        ),
                    }
                )
    steps.sort(
        key=lambda row: (
            methods.index(str(row["method"])),
            PAPER_SYSTEMS.index(str(row["solid"])),
            mesh_number(str(row["coarse_mesh"])),
        )
    )
    selections.sort(
        key=lambda row: (
            methods.index(str(row["method"])),
            PAPER_SYSTEMS.index(str(row["solid"])),
        )
    )
    return steps, selections, pending


def _csv_bytes(rows: list[dict[str, object]]) -> bytes:
    if not rows:
        return b""
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def write_convergence_artifacts(
    root: Path,
    values: list[dict[str, object]],
    steps: list[dict[str, object]],
    selections: list[dict[str, object]],
    pending: list[tuple[str, str, int]],
    *,
    campaign: dict[str, object],
    fits_sha256: str,
    methods: tuple[str, ...] = METHODS,
    maximum_mesh: int | None = None,
) -> dict[str, object]:
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    outputs = {
        VALUES_NAME: _csv_bytes(values),
        STEPS_NAME: _csv_bytes(steps),
        SELECTION_NAME: _csv_bytes(selections),
    }
    for name, content in outputs.items():
        temp = data / f".{name}.tmp.{os.getpid()}"
        temp.write_bytes(content)
        os.replace(temp, data / name)
    unconverged = [
        row for row in selections if row.get("selection_status") != "converged"
    ]
    resource_errors = [
        row
        for row in unconverged
        if row.get("selection_status") == "technical_resource_guard_reached"
    ]
    complete = (
        not pending
        and not unconverged
        and len(selections) == len(methods) * len(PAPER_SYSTEMS)
    )
    payload: dict[str, object] = {
        "schema_version": CONVERGENCE_SCHEMA,
        "benchmark": "LC10 (fixed Goldzak12 subset)",
        "status": (
            "converged"
            if complete
            else "technical_resource_limit_reached"
            if resource_errors
            else "incomplete"
        ),
        "methods": list(methods),
        "paper_systems": list(PAPER_SYSTEMS),
        "diagnostic_only_systems": list(base.LC10_DIAGNOSTIC_ONLY_SOLIDS),
        "algorithm": {
            "name": "one-step adaptive convergence",
            "initial_meshes": [mesh_name(number) for number in INITIAL_MESH_NUMBERS],
            "scientific_maximum_mesh": None,
            "technical_resource_guard_mesh": (
                mesh_name(maximum_mesh) if maximum_mesh is not None else None
            ),
            "technical_resource_guard_is_convergence": False,
            "required_consecutive_passing_steps": 1,
            "aggregate_rms_gate": False,
            "criteria_combination": "AND",
            "lattice_abs_delta_threshold_A": LATTICE_THRESHOLD_A,
            "cohesive_abs_delta_threshold_kJmol_per_atom": ECOH_THRESHOLD_KJMOL_PER_ATOM,
            "cohesive_abs_delta_threshold_eV_per_atom": ECOH_THRESHOLD_EV_PER_ATOM,
            "selected_value": "denser value from the earliest passing n->n+1 step",
        },
        "value_protocol": {
            "lattice": "quadratic a0 from an independent EOS at each mesh",
            "cohesive": "single point at that mesh's own independent-EOS minimum",
            "fixed_geometry_single_point_series_is_separate": True,
            "fixed_geometry_eos_mesh": FIXED_GEOMETRY_EOS_MESH,
            "fixed_geometry_energy_meshes": list(FIXED_GEOMETRY_ENERGY_MESHES),
        },
        "campaign_identity": campaign,
        "eos_fits_sha256": fits_sha256,
        "pending": [
            {"method": method, "solid": solid, "next_mesh": mesh_name(number)}
            for method, solid, number in pending
        ],
        "unconverged": unconverged,
        "resource_errors": resource_errors,
        "artifacts": {
            name: {
                "path": f"data/{name}",
                "sha256": hashlib.sha256(content).hexdigest(),
                "size_bytes": len(content),
            }
            for name, content in outputs.items()
        },
    }
    target = data / CONVERGENCE_NAME
    temp = data / f".{CONVERGENCE_NAME}.tmp.{os.getpid()}"
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temp, target)
    return payload


def selected_fits(
    fits: Iterable[dict[str, object]], selections: Iterable[dict[str, object]]
) -> list[dict[str, object]]:
    by_key = {
        (str(row["method"]), str(row["solid"]), str(row["eos_mesh"])): row
        for row in fits
    }
    result: list[dict[str, object]] = []
    for selected in selections:
        if selected.get("selection_status") != "converged":
            continue
        key = (
            str(selected["method"]),
            str(selected["solid"]),
            str(selected["selected_mesh"]),
        )
        result.append(by_key[key])
    return result


def scale_manifest(
    root: Path,
    fits: Iterable[dict[str, object]],
    scales: tuple[float, ...],
    methods: tuple[str, ...] = METHODS,
) -> dict[str, object]:
    records = []
    for fit in fits:
        method = str(fit.get("method", ""))
        solid = str(fit.get("solid", ""))
        mesh = str(fit.get("eos_mesh", ""))
        if method not in methods or solid not in PAPER_SYSTEMS:
            continue
        records.append(
            {
                "method": method,
                "solid": solid,
                "mesh": mesh,
                "requested_scales": list(eos.scales_for(solid, method, scales)),
            }
        )
    records.sort(
        key=lambda row: (
            methods.index(row["method"]),
            PAPER_SYSTEMS.index(row["solid"]),
            mesh_number(row["mesh"]),
        )
    )
    payload = {
        "schema_version": 1,
        "benchmark": "LC10 independent-EOS k convergence",
        "methods": list(methods),
        "paper_systems": list(PAPER_SYSTEMS),
        "records": records,
    }
    path = root / "data" / SCALE_MANIFEST_NAME
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def validate_campaign_outputs(
    specs: Iterable[tuple[str, Path, Path, bool]],
    campaign: dict[str, object],
) -> None:
    for label, _input_path, output_path, _require_opt in specs:
        issue = base.completed_stamp_campaign_issue(
            output_path,
            campaign,
            executable_role="cp2k",
            require_completed=True,
        )
        if issue:
            raise RuntimeError(f"{label}: {issue}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-manifest", type=Path, required=True)
    parser.add_argument("--campaign-manifest-sha256", required=True)
    parser.add_argument("--cp2k", type=Path)
    parser.add_argument("--cp2k-library", type=Path)
    parser.add_argument("--save-tblite", type=Path)
    parser.add_argument("--save-tblite-library", type=Path)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument(
        "--method",
        action="append",
        choices=METHODS,
        help="method to execute; repeat to select multiple methods (default: all)",
    )
    parser.add_argument(
        "--maximum-mesh",
        type=int,
        help=(
            "optional technical resource guard N for an NxNxN mesh; "
            "the default has no fixed mesh cap"
        ),
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mpi-ranks-per-job", type=int, default=1)
    parser.add_argument("--mpi-launcher", type=Path)
    parser.add_argument("--mpi-launcher-arg", action="append", default=[])
    parser.add_argument("--cpu-set", action="append", default=[])
    parser.add_argument("--taskset", default="taskset")
    parser.add_argument("--scale", type=float, action="append", default=[])
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--fit-only", action="store_true")
    parser.add_argument("--approve-fits", action="store_true")
    parser.add_argument("--stop-after-convergence", action="store_true")
    parser.add_argument(
        "--classification-manifest",
        type=Path,
        default=eos.gxtb_classification_manifest_path(),
    )
    args = parser.parse_args()

    try:
        run_methods = selected_methods(args.method)
    except ValueError as error:
        parser.error(str(error))

    if not re.fullmatch(r"[0-9a-f]{64}", args.campaign_manifest_sha256):
        parser.error("--campaign-manifest-sha256 must be 64 lowercase hexadecimal digits")
    observed_manifest_sha = base.sha256(args.campaign_manifest.resolve(strict=True))
    if observed_manifest_sha != args.campaign_manifest_sha256:
        parser.error(
            "campaign manifest hash pin mismatch: expected "
            f"{args.campaign_manifest_sha256}, observed {observed_manifest_sha}"
        )
    if args.jobs < 1 or args.threads < 1 or args.mpi_ranks_per_job < 1:
        parser.error("--jobs, --threads, and --mpi-ranks-per-job must be positive")
    if args.maximum_mesh is not None and args.maximum_mesh < max(INITIAL_MESH_NUMBERS):
        parser.error(
            f"--maximum-mesh must be at least {max(INITIAL_MESH_NUMBERS)}"
        )
    if args.fit_only and (args.force or args.approve_fits):
        parser.error("--fit-only cannot be combined with --force or --approve-fits")
    if args.stop_after_convergence and args.approve_fits:
        parser.error("--stop-after-convergence cannot be combined with --approve-fits")

    try:
        campaign, paths = base.validated_gxtb_campaign_from_manifest(
            args.campaign_manifest,
            args.cp2k_source,
            args.save_tblite_source,
            cp2k_override=args.cp2k,
            cp2k_library_override=args.cp2k_library,
            save_tblite_override=args.save_tblite,
            save_tblite_library_override=args.save_tblite_library,
        )
        base.require_git_ancestor(args.cp2k_source, base.REQUIRED_CP2K_POST5582_ANCESTOR)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    cp2k = paths["cp2k"]
    save_tblite = paths["save_tblite"]

    execution_requested = bool(
        args.mpi_launcher
        or args.mpi_launcher_arg
        or args.cpu_set
        or args.mpi_ranks_per_job != 1
    )
    pool: execution.ExecutionPool | None = None
    if execution_requested:
        if args.mpi_launcher is None:
            parser.error("--mpi-launcher is required with MPI/affinity execution")
        try:
            pool = execution.ExecutionPool(
                concurrent_jobs=args.jobs,
                mpi_ranks_per_job=args.mpi_ranks_per_job,
                threads_per_rank=args.threads,
                mpi_launcher=args.mpi_launcher,
                mpi_launcher_args=args.mpi_launcher_arg,
                cpu_sets=args.cpu_set,
                taskset=args.taskset,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))

    scales = tuple(args.scale) if args.scale else eos.DEFAULT_SCALES
    try:
        classifications = eos.load_gxtb_classifications(args.classification_manifest)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))

    protocol: dict[str, object] = {
        "benchmark": "LC10 (fixed Goldzak12 subset)",
        "methods": list(run_methods),
        "execution_methods": list(run_methods),
        "publication_methods": list(METHODS),
        "frozen_external_baseline_methods": [
            method for method in METHODS if method not in run_methods
        ],
        "selected_solids": list(PAPER_SYSTEMS),
        "paper_systems": list(PAPER_SYSTEMS),
        "diagnostic_only_systems": list(base.LC10_DIAGNOSTIC_ONLY_SOLIDS),
        "exact_lc10_scope": True,
        "single_cp2k_binary_for_selected_methods": True,
        "single_tblite_provider_for_selected_methods": "save_tblite",
        "campaign_manifest_sha256_external_pin": args.campaign_manifest_sha256,
        "adaptive_k_convergence": {
            "initial_meshes": [mesh_name(number) for number in INITIAL_MESH_NUMBERS],
            "scientific_maximum_mesh": None,
            "technical_resource_guard_mesh": (
                mesh_name(args.maximum_mesh)
                if args.maximum_mesh is not None
                else None
            ),
            "technical_resource_guard_is_convergence": False,
            "required_consecutive_passing_steps": 1,
            "aggregate_rms_gate": False,
            "criteria_combination": "AND",
            "lattice_abs_delta_threshold_A": LATTICE_THRESHOLD_A,
            "cohesive_abs_delta_threshold_kJmol_per_atom": ECOH_THRESHOLD_KJMOL_PER_ATOM,
            "cohesive_abs_delta_threshold_eV_per_atom": ECOH_THRESHOLD_EV_PER_ATOM,
            "selected_value_policy": "take denser n+1 value from earliest passing step",
            "equilibrium_energy_protocol": "single point at each mesh's own independent EOS minimum",
        },
        "fixed_geometry_single_point_diagnostic": {
            "separate_from_adaptive_selection": True,
            "eos_mesh": FIXED_GEOMETRY_EOS_MESH,
            "energy_meshes": list(FIXED_GEOMETRY_ENERGY_MESHES),
        },
        "fit_approval_required": True,
        "fit_approved": False,
        "approved_fit_sha256": None,
        "current_fit_sha256": None,
        "kpoint_scheme": "CP2K native Bloch MACDONALD with full SPGLIB symmetry reduction",
        "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
        "required_cp2k_ancestor": base.REQUIRED_CP2K_POST5582_ANCESTOR,
        "execution_provenance": (
            {
                "separate_from_scientific_job_stamp": True,
                "record_schema": execution.SCHEMA_VERSION,
                "contract": pool.contract,
                "contract_sha256": pool.contract_sha256,
            }
            if pool is not None
            else None
        ),
    }
    if run_methods == METHODS:
        protocol["single_cp2k_binary_for_all_methods"] = True
        protocol["single_tblite_provider_for_all_methods"] = "save_tblite"

    def write_provenance() -> None:
        if any(method in base.LEGACY_METHODS for method in run_methods):
            base.write_build_provenance(
                cp2k,
                save_tblite,
                args.cp2k_source,
                args.save_tblite_source,
                protocol,
            )
        if "GXTB" in run_methods:
            base.write_gxtb_build_provenance(
                cp2k,
                save_tblite,
                args.cp2k_source,
                args.save_tblite_source,
                protocol,
                campaign,
                args.campaign_manifest,
            )

    write_provenance()
    if not args.fit_only:
        base.run_tblite_atom_jobs(
            save_tblite,
            args.jobs,
            args.force,
            run_methods,
            save_tblite,
            campaign,
            PAPER_ELEMENTS,
            campaign_bind_all_methods=True,
        )

    targets: set[tuple[str, str, int]] = {
        (method, solid, number)
        for method in run_methods
        for solid in PAPER_SYSTEMS
        for number in INITIAL_MESH_NUMBERS
    }
    # Preserve and rebind any already archived adaptive meshes to this exact
    # campaign instead of silently omitting them from the final fingerprint.
    archived_maximum: dict[tuple[str, str], int] = {}
    for row in base.read_csv(base.ROOT / "data" / "eos_fits.csv"):
        method = str(row.get("method", ""))
        solid = str(row.get("solid", ""))
        if method not in run_methods or solid not in PAPER_SYSTEMS:
            continue
        try:
            number = mesh_number(str(row.get("eos_mesh", "")))
        except ValueError:
            continue
        if number >= min(INITIAL_MESH_NUMBERS):
            archived_maximum[(method, solid)] = max(
                number, archived_maximum.get((method, solid), 3)
            )
    for (method, solid), maximum in archived_maximum.items():
        targets.update((method, solid, number) for number in range(3, maximum + 1))
    completed_targets: set[tuple[str, str, int]] = set()
    final_values: list[dict[str, object]] = []
    final_steps: list[dict[str, object]] = []
    final_selections: list[dict[str, object]] = []
    final_pending: list[tuple[str, str, int]] = []

    while True:
        new_targets = sorted(
            targets - completed_targets,
            key=lambda item: (
                item[2],
                run_methods.index(item[0]),
                PAPER_SYSTEMS.index(item[1]),
            ),
        )
        if args.fit_only and new_targets:
            # fit-only consumes exactly what is archived; missing raw work fails below.
            new_targets = []
        for number in sorted({item[2] for item in new_targets}):
            mesh = mesh_name(number)
            for method in run_methods:
                solids = tuple(
                    solid
                    for solid in PAPER_SYSTEMS
                    if (method, solid, number) in new_targets
                )
                if not solids:
                    continue
                specs = eos.eos_job_specs(mesh, scales, (method,), solids)
                eos.run_jobs(
                    specs,
                    cp2k,
                    args.jobs,
                    args.threads,
                    args.force,
                    retry_scf=False,
                    campaign_fingerprint=campaign,
                    execution_pool=pool,
                    campaign_bind_all_methods=True,
                )
                validate_campaign_outputs(specs, campaign)
                eos.make_eos_table(
                    mesh,
                    scales,
                    (method,),
                    classifications,
                    campaign,
                    solids,
                )
                completed_targets.update((method, solid, number) for solid in solids)

        all_fits = [
            dict(row)
            for row in base.read_csv(base.ROOT / "data" / "eos_fits.csv")
            if row.get("method") in run_methods and row.get("solid") in PAPER_SYSTEMS
        ]
        expected_fit_keys = targets if not args.fit_only else {
            (str(row["method"]), str(row["solid"]), mesh_number(str(row["eos_mesh"])))
            for row in all_fits
        }
        fit_by_key = {
            (str(row["method"]), str(row["solid"]), mesh_number(str(row["eos_mesh"]))): row
            for row in all_fits
        }
        missing_fits = expected_fit_keys - set(fit_by_key)
        if missing_fits:
            raise RuntimeError(f"missing independent EOS fits: {sorted(missing_fits)}")
        invalid = [
            key
            for key in expected_fit_keys
            if fit_by_key[key].get("fit_status") != "quadratic"
            or not str(fit_by_key[key].get("a_eos_A", ""))
            or int(fit_by_key[key].get("n_unresolved_branch_candidates", 0) or 0) != 0
        ]
        if invalid:
            raise RuntimeError(
                "invalid independent EOS fits require scale/branch review before k adaptation: "
                + ", ".join("/".join((method, solid, mesh_name(number))) for method, solid, number in invalid)
            )
        active_fits = [
            fit_by_key[key]
            for key in sorted(
                expected_fit_keys,
                key=lambda item: (
                    run_methods.index(item[0]),
                    PAPER_SYSTEMS.index(item[1]),
                    item[2],
                ),
            )
        ]
        eq_specs = equilibrium_specs(base.ROOT, active_fits, run_methods)
        if not args.fit_only:
            eos.run_jobs(
                eq_specs,
                cp2k,
                args.jobs,
                args.threads,
                args.force,
                retry_scf=False,
                campaign_fingerprint=campaign,
                execution_pool=pool,
                campaign_bind_all_methods=True,
            )
        validate_campaign_outputs(eq_specs, campaign)
        final_values = collect_independent_values(
            base.ROOT, active_fits, campaign, run_methods
        )
        final_steps, final_selections, final_pending = assess_convergence(
            final_values,
            methods=run_methods,
            maximum_mesh=args.maximum_mesh,
        )
        scale_payload = scale_manifest(base.ROOT, active_fits, scales, run_methods)
        scale_path = base.ROOT / "data" / SCALE_MANIFEST_NAME
        protocol["k_convergence_scale_manifest_sha256"] = base.sha256(scale_path)
        protocol["k_convergence_scale_manifest"] = scale_payload
        protocol["current_fit_sha256"] = fit_fingerprint(active_fits, run_methods)
        convergence = write_convergence_artifacts(
            base.ROOT,
            final_values,
            final_steps,
            final_selections,
            final_pending,
            campaign=campaign,
            fits_sha256=base.sha256(base.ROOT / "data" / "eos_fits.csv"),
            methods=run_methods,
            maximum_mesh=args.maximum_mesh,
        )
        protocol["k_convergence_status"] = convergence["status"]
        protocol["k_convergence_artifact_sha256"] = base.sha256(
            base.ROOT / "data" / CONVERGENCE_NAME
        )
        write_provenance()
        resource_limited = [
            row
            for row in final_selections
            if row.get("selection_status") == "technical_resource_guard_reached"
        ]
        if resource_limited:
            assert args.maximum_mesh is not None
            labels = ", ".join(
                f"{row['method']}/{row['solid']}"
                for row in resource_limited
            )
            raise RuntimeError(
                "LC10 adaptive execution reached the technical "
                f"--maximum-mesh {mesh_name(args.maximum_mesh)} without "
                f"scientific convergence; no value was selected: {labels}"
            )
        if not final_pending:
            break
        if args.fit_only:
            raise RuntimeError(
                "archived independent EOS series is incomplete; next targets: "
                + ", ".join(
                    f"{method}/{solid}/{mesh_name(number)}"
                    for method, solid, number in final_pending
                )
            )
        targets.update(final_pending)

    expected_selections = len(run_methods) * len(PAPER_SYSTEMS)
    if len(final_selections) != expected_selections or any(
        row.get("selection_status") != "converged" for row in final_selections
    ):
        raise RuntimeError("LC10 adaptive k selection is incomplete")
    if args.stop_after_convergence or args.fit_only or not args.approve_fits:
        print(
            "Independent EOS k convergence is complete. Review the raw per-step "
            f"tables and fit fingerprint {protocol['current_fit_sha256']}, then rerun "
            "with --approve-fits for the separate fixed-geometry diagnostic series."
        )
        return 0

    protocol["fit_approved"] = True
    protocol["approved_fit_sha256"] = protocol["current_fit_sha256"]
    write_provenance()

    # Historical diagnostic only: all k meshes use the same k444 EOS geometry.
    fit_lookup = {
        (str(row["method"]), str(row["solid"]), str(row["eos_mesh"])): row
        for row in base.read_csv(base.ROOT / "data" / "eos_fits.csv")
    }
    fixed_fits = [
        fit_lookup[(method, solid, FIXED_GEOMETRY_EOS_MESH)]
        for method in run_methods
        for solid in PAPER_SYSTEMS
    ]
    fixed_specs = eos.final_sp_specs(
        fixed_fits,
        list(FIXED_GEOMETRY_ENERGY_MESHES),
        PAPER_SYSTEMS,
    )
    eos.run_jobs(
        fixed_specs,
        cp2k,
        args.jobs,
        args.threads,
        args.force,
        retry_scf=False,
        campaign_fingerprint=campaign,
        execution_pool=pool,
        campaign_bind_all_methods=True,
    )
    validate_campaign_outputs(fixed_specs, campaign)
    eos.collect_results(
        fixed_fits,
        list(FIXED_GEOMETRY_ENERGY_MESHES),
        "k555",
        run_methods,
        campaign,
        campaign_bind_all_methods=True,
    )
    protocol["fixed_geometry_single_point_diagnostic_complete"] = True
    write_provenance()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
