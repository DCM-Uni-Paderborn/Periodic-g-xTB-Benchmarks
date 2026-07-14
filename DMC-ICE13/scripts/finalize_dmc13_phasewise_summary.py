#!/usr/bin/env python3
"""Freeze the phase-wise DMC13 comparison into paper-ready CSV and JSON."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any


HARTREE_TO_KJMOL = 2625.499638
CONVERGENCE_THRESHOLD_KJMOL_PER_H2O = 0.05
METHODS = ("GFN1", "GFN2", "GXTB")
METHOD_LABELS = {
    "GFN1": "GFN1-xTB",
    "GFN2": "GFN2-xTB",
    "GXTB": "g-xTB",
}
NONREFERENCE_PHASES = (
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
SUMMARY_STEM = "dmc_ice13_gfn_gxtb_phasewise_summary"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return value


def finite_float(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def stats(errors: list[float]) -> dict[str, float]:
    if not errors:
        raise ValueError("cannot summarize an empty error set")
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(value) for value in errors) / len(errors),
        "RMSE": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "MaxAE": max(abs(value) for value in errors),
    }


def relative_path(path: Path, root: Path) -> str:
    resolved = path.resolve()
    root_resolved = root.resolve()
    try:
        return resolved.relative_to(root_resolved).as_posix()
    except ValueError:
        return str(resolved)


def artifact(path: Path, root: Path) -> dict[str, str]:
    if not path.is_file():
        raise ValueError(f"required artifact is missing: {path}")
    return {"path": relative_path(path, root), "sha256": sha256(path)}


def provenance_record(
    method: str,
    baseline: dict[str, Any],
    gxtb: dict[str, Any],
    baseline_path: Path,
    gxtb_path: Path,
    root: Path,
) -> dict[str, Any]:
    if method in {"GFN1", "GFN2"}:
        cp2k = baseline.get("cp2k", {})
        provider = baseline.get("tblite", {})
        if not isinstance(cp2k, dict) or not isinstance(provider, dict):
            raise ValueError("invalid frozen GFN1/GFN2 provenance")
        return {
            "artifact": artifact(baseline_path, root),
            "cp2k_source_revision": cp2k.get("source_revision"),
            "cp2k_executable_sha256": cp2k.get("executable_sha256"),
            "cp2k_library_sha256": cp2k.get("library_sha256"),
            "provider_name": "tblite",
            "provider_source_revision": provider.get("local_merge_revision"),
            "provider_main_revision": provider.get("main_revision"),
            "provider_executable_sha256": provider.get("executable_sha256"),
            "provider_library_sha256": provider.get("library_sha256"),
        }

    cp2k = gxtb.get("cp2k", {})
    provider = gxtb.get("save_tblite", {})
    campaign = gxtb.get("campaign", {})
    protocol = gxtb.get("protocol", {})
    if not all(
        isinstance(value, dict)
        for value in (cp2k, provider, campaign, protocol)
    ):
        raise ValueError("invalid g-xTB provenance")
    return {
        "artifact": artifact(gxtb_path, root),
        "campaign_id": campaign.get("id"),
        "campaign_manifest_sha256": campaign.get("manifest_sha256"),
        "execution_build_id": cp2k.get("execution_build_id"),
        "protocol_id": protocol.get("gxtb_protocol_id"),
        "cp2k_source_revision": cp2k.get("source_revision_validated"),
        "cp2k_executable_sha256": cp2k.get("sha256"),
        "cp2k_library_sha256": cp2k.get("loaded_library_sha256"),
        "provider_name": "save_tblite",
        "provider_source_revision": provider.get("source_revision_validated"),
        "provider_executable_sha256": provider.get("sha256"),
        "provider_library_sha256": provider.get("static_library_sha256"),
    }


def mesh_size(mesh: str) -> int:
    if not mesh.startswith("k"):
        raise ValueError(f"phase-wise result selected a non-k mesh: {mesh}")
    digits = mesh[1:]
    for size in range(1, 100):
        if digits == str(size) * 3:
            return size
    raise ValueError(f"cannot decode cubic mesh {mesh}")


def build_summary(root: Path) -> tuple[dict[str, Any], list[dict[str, object]]]:
    data = root / "data"
    phase_json_path = data / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json"
    phase_csv_path = data / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.csv"
    results_path = data / "dmc_ice13_gxtb_spglib_kpoint_results.json"
    validation_path = data / "dmc_ice13_gxtb_spglib_validation_index.json"
    geometries_path = data / "geometries.json"
    baseline_provenance_path = data / "build_provenance.json"
    gxtb_provenance_path = data / "build_provenance_gxtb_spglib.json"

    report = read_json(phase_json_path)
    results = read_json(results_path)
    geometries = read_json(geometries_path)
    baseline_provenance = read_json(baseline_provenance_path)
    gxtb_provenance = read_json(gxtb_provenance_path)
    sources = {
        "phasewise_json": artifact(phase_json_path, root),
        "phasewise_csv": artifact(phase_csv_path, root),
        "kpoint_results_json": artifact(results_path, root),
        "validation_index": artifact(validation_path, root),
        "baseline_provenance": artifact(baseline_provenance_path, root),
        "gxtb_provenance": artifact(gxtb_provenance_path, root),
        "geometries": artifact(geometries_path, root),
    }

    report_methods = report.get("methods")
    result_meshes = results.get("results")
    if not isinstance(report_methods, dict) or not isinstance(result_meshes, dict):
        raise ValueError("phase-wise report or k-point results has invalid structure")
    if report.get("dmc_reference_variant") != "legacy_rounded_absolute_XI_0.16":
        raise ValueError("unexpected DMC13 reference variant")
    if report.get("reference_phase") != "Ih" or not report.get("same_mesh_ih_required"):
        raise ValueError("phase-wise result is not bound to same-mesh ice Ih")

    method_payloads: dict[str, Any] = {}
    csv_rows: list[dict[str, object]] = []
    for method in METHODS:
        method_report = report_methods.get(method)
        if not isinstance(method_report, dict):
            raise ValueError(f"phase-wise report is missing {method}")
        if (
            method_report.get("status") != "phasewise_kpoint_converged"
            or method_report.get("phasewise_kpoint_converged") is not True
            or method_report.get("converged_phase_count") != len(NONREFERENCE_PHASES)
        ):
            raise ValueError(f"{method} is not phase-wise k-point converged")
        phase_convergence = method_report.get("phase_convergence")
        aggregate = method_report.get(
            "phasewise_kpoint_converged_stats_nonreference"
        )
        if not isinstance(phase_convergence, dict) or not isinstance(aggregate, dict):
            raise ValueError(f"{method} lacks complete phase-wise statistics")

        provenance = provenance_record(
            method,
            baseline_provenance,
            gxtb_provenance,
            baseline_provenance_path,
            gxtb_provenance_path,
            root,
        )
        phases: dict[str, Any] = {}
        selected_mesh_by_phase: dict[str, str] = {}
        previous_mesh_by_phase: dict[str, str] = {}
        absolute_last_delta_by_phase: dict[str, float] = {}
        errors: list[float] = []
        for phase in NONREFERENCE_PHASES:
            selected = phase_convergence.get(phase)
            if not isinstance(selected, dict):
                raise ValueError(f"{method}/{phase} lacks a selected mesh")
            mesh = str(selected.get("smallest_required_mesh", ""))
            previous_mesh = str(selected.get("previous_mesh", ""))
            mesh_result = result_meshes.get(mesh, {}).get(method, {})
            if not isinstance(mesh_result, dict):
                raise ValueError(f"{method}/{phase}/{mesh} result is missing")
            energies = mesh_result.get("energies_hartree")
            per_h2o = mesh_result.get("per_h2o_hartree")
            if not isinstance(energies, dict) or not isinstance(per_h2o, dict):
                raise ValueError(f"{method}/{phase}/{mesh} raw energies are missing")
            phase_total = finite_float(
                energies.get(phase), f"{method}/{phase}/{mesh} total energy"
            )
            ih_total = finite_float(
                energies.get("Ih"), f"{method}/Ih/{mesh} total energy"
            )
            phase_count = int(geometries.get(phase, {}).get("counts", {}).get("O", 0))
            ih_count = int(geometries.get("Ih", {}).get("counts", {}).get("O", 0))
            if phase_count <= 0 or ih_count <= 0:
                raise ValueError(f"invalid H2O count for {phase} or Ih")
            phase_per_h2o = phase_total / phase_count
            ih_per_h2o = ih_total / ih_count
            stored_phase_per_h2o = finite_float(
                per_h2o.get(phase), f"{method}/{phase}/{mesh} per-H2O energy"
            )
            stored_ih_per_h2o = finite_float(
                per_h2o.get("Ih"), f"{method}/Ih/{mesh} per-H2O energy"
            )
            if not math.isclose(
                phase_per_h2o, stored_phase_per_h2o, rel_tol=0.0, abs_tol=1.0e-12
            ) or not math.isclose(
                ih_per_h2o, stored_ih_per_h2o, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError(f"{method}/{phase}/{mesh} per-H2O energy mismatch")
            relative = (phase_per_h2o - ih_per_h2o) * HARTREE_TO_KJMOL
            reported_relative = finite_float(
                selected.get("relative_energy_kjmol_per_h2o"),
                f"{method}/{phase}/{mesh} relative energy",
            )
            if not math.isclose(relative, reported_relative, rel_tol=0.0, abs_tol=1.0e-9):
                raise ValueError(f"{method}/{phase}/{mesh} same-mesh-Ih mismatch")
            dmc_relative = finite_float(
                selected.get("dmc_relative_kjmol_per_h2o"),
                f"DMC/{phase} relative energy",
            )
            error = relative - dmc_relative
            reported_error = finite_float(
                selected.get("error_kjmol_per_h2o"),
                f"{method}/{phase} error",
            )
            if not math.isclose(error, reported_error, rel_tol=0.0, abs_tol=1.0e-9):
                raise ValueError(f"{method}/{phase} error mismatch")
            last_delta = finite_float(
                selected.get("last_delta_kjmol_per_h2o"),
                f"{method}/{phase} last delta",
            )
            absolute_last_delta = finite_float(
                selected.get("absolute_last_delta_kjmol_per_h2o"),
                f"{method}/{phase} absolute last delta",
            )
            if absolute_last_delta > CONVERGENCE_THRESHOLD_KJMOL_PER_H2O + 1.0e-12:
                raise ValueError(f"{method}/{phase} exceeds the phase convergence threshold")
            if not math.isclose(
                abs(last_delta), absolute_last_delta, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError(f"{method}/{phase} delta magnitude mismatch")

            selected_mesh_by_phase[phase] = mesh
            previous_mesh_by_phase[phase] = previous_mesh
            absolute_last_delta_by_phase[phase] = absolute_last_delta
            errors.append(error)
            phases[phase] = {
                "selected_mesh": mesh,
                "previous_mesh": previous_mesh,
                "mesh_n": mesh_size(mesh),
                "nk_total": int(selected.get("nk_total")),
                "phase_total_energy_hartree": phase_total,
                "phase_n_h2o": phase_count,
                "phase_energy_per_h2o_hartree": phase_per_h2o,
                "ih_total_energy_hartree": ih_total,
                "ih_n_h2o": ih_count,
                "ih_energy_per_h2o_hartree": ih_per_h2o,
                "relative_energy_kjmol_per_h2o": relative,
                "dmc_relative_energy_kjmol_per_h2o": dmc_relative,
                "error_kjmol_per_h2o": error,
                "last_delta_kjmol_per_h2o": last_delta,
                "absolute_last_delta_kjmol_per_h2o": absolute_last_delta,
                "parity_direction": selected.get("parity_direction"),
            }

        recomputed = stats(errors)
        metrics = {
            key: finite_float(aggregate.get(key), f"{method} {key}")
            for key in ("ME", "MAE", "RMSE", "MaxAE")
        }
        for key, value in metrics.items():
            if not math.isclose(value, recomputed[key], rel_tol=0.0, abs_tol=1.0e-9):
                raise ValueError(f"{method} aggregate {key} mismatch")

        primary_errors = list(errors)
        xi_index = NONREFERENCE_PHASES.index("XI")
        primary_errors[xi_index] += 0.01
        mesh_distribution = Counter(selected_mesh_by_phase.values())
        ordered_distribution = {
            mesh: mesh_distribution[mesh]
            for mesh in sorted(mesh_distribution, key=mesh_size)
        }
        mesh_sizes = [mesh_size(mesh) for mesh in selected_mesh_by_phase.values()]
        method_payloads[method] = {
            "method_label": METHOD_LABELS[method],
            "status": "phasewise_kpoint_converged",
            "n_nonreference_phases": len(NONREFERENCE_PHASES),
            "metrics_kjmol_per_h2o": metrics,
            "selected_mesh_by_phase": selected_mesh_by_phase,
            "previous_mesh_by_phase": previous_mesh_by_phase,
            "absolute_last_delta_kjmol_per_h2o_by_phase": (
                absolute_last_delta_by_phase
            ),
            "minimum_selected_mesh_n": min(mesh_sizes),
            "maximum_selected_mesh_n": max(mesh_sizes),
            "selected_mesh_distribution": ordered_distribution,
            "reference_sensitivity": {
                "legacy_rounded_absolute_XI_0.16": metrics,
                "primary_explicit_relative_XI_0.15": stats(primary_errors),
            },
            "phases": phases,
            "provenance": provenance,
        }
        csv_rows.append(
            {
                "method_id": method,
                "method_label": METHOD_LABELS[method],
                "status": "phasewise_kpoint_converged",
                "N_nonreference_phases": len(NONREFERENCE_PHASES),
                "reference_phase": "Ih",
                "ME_kJmol_per_H2O": f"{metrics['ME']:.9f}",
                "MAE_kJmol_per_H2O": f"{metrics['MAE']:.9f}",
                "RMSE_kJmol_per_H2O": f"{metrics['RMSE']:.9f}",
                "MaxAE_kJmol_per_H2O": f"{metrics['MaxAE']:.9f}",
                "convergence_threshold_kJmol_per_H2O": (
                    f"{CONVERGENCE_THRESHOLD_KJMOL_PER_H2O:.6f}"
                ),
                "minimum_selected_mesh_n": min(mesh_sizes),
                "maximum_selected_mesh_n": max(mesh_sizes),
                "selected_mesh_distribution": ";".join(
                    f"{mesh}:{count}" for mesh, count in ordered_distribution.items()
                ),
                "reference_variant": "legacy_rounded_absolute_XI_0.16",
                "phasewise_json_sha256": sources["phasewise_json"]["sha256"],
                "phasewise_csv_sha256": sources["phasewise_csv"]["sha256"],
                "validation_index_sha256": sources["validation_index"]["sha256"],
                "provenance_sha256": provenance["artifact"]["sha256"],
                "cp2k_source_revision": provenance.get("cp2k_source_revision"),
                "provider_name": provenance.get("provider_name"),
                "provider_source_revision": provenance.get(
                    "provider_source_revision"
                ),
                "provider_library_sha256": provenance.get(
                    "provider_library_sha256"
                ),
            }
        )

    summary = {
        "schema_version": 1,
        "benchmark": "DMC-ICE13",
        "status": "phasewise_kpoint_converged",
        "result_label": "phase-wise k-point-converged relative-energy errors",
        "quantity": "same-mesh-Ih-referenced relative lattice energy",
        "unit": "kJ mol^-1 per H2O",
        "reference_phase": "Ih",
        "reference_phase_in_statistics": False,
        "n_nonreference_phases": len(NONREFERENCE_PHASES),
        "nonreference_phases": list(NONREFERENCE_PHASES),
        "dmc_reference": {
            "citation": (
                "Della Pia, Zen, Alfe, and Michaelides, "
                "J. Chem. Phys. 157, 134701 (2022)"
            ),
            "doi": "10.1063/5.0102645",
            "reported_variant": "legacy_rounded_absolute_XI_0.16",
            "sensitivity_variant": "primary_explicit_relative_XI_0.15",
        },
        "conversion": {
            "hartree_to_kjmol": HARTREE_TO_KJMOL,
            "formula": (
                "Erel_phase(N) = 2625.499638 * "
                "[E_phase(N)/N_H2O_phase - E_Ih(N)/N_H2O_Ih]"
            ),
        },
        "convergence": {
            "threshold_kjmol_per_h2o": CONVERGENCE_THRESHOLD_KJMOL_PER_H2O,
            "rule": (
                "A phase is converged at N^3 when "
                "|Erel(N^3)-Erel((N-1)^3)| <= 0.05 kJ mol^-1 per H2O."
            ),
            "reported_value": "the higher-mesh N^3 value",
            "same_mesh_ih_required": True,
            "later_available_evidence_safety_check": True,
        },
        "sources": sources,
        "methods": method_payloads,
    }
    return summary, csv_rows


CSV_FIELDS = (
    "method_id",
    "method_label",
    "status",
    "N_nonreference_phases",
    "reference_phase",
    "ME_kJmol_per_H2O",
    "MAE_kJmol_per_H2O",
    "RMSE_kJmol_per_H2O",
    "MaxAE_kJmol_per_H2O",
    "convergence_threshold_kJmol_per_H2O",
    "minimum_selected_mesh_n",
    "maximum_selected_mesh_n",
    "selected_mesh_distribution",
    "reference_variant",
    "phasewise_json_sha256",
    "phasewise_csv_sha256",
    "validation_index_sha256",
    "provenance_sha256",
    "cp2k_source_revision",
    "provider_name",
    "provider_source_revision",
    "provider_library_sha256",
)


def finalize(root: Path) -> tuple[Path, Path]:
    data = root / "data"
    csv_path = data / f"{SUMMARY_STEM}.csv"
    json_path = data / f"{SUMMARY_STEM}.json"
    csv_path.unlink(missing_ok=True)
    json_path.unlink(missing_ok=True)

    summary, rows = build_summary(root)
    data.mkdir(parents=True, exist_ok=True)
    csv_tmp = csv_path.with_name(f".{csv_path.name}.{os.getpid()}.tmp")
    json_tmp = json_path.with_name(f".{json_path.name}.{os.getpid()}.tmp")
    try:
        with csv_tmp.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=CSV_FIELDS, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        json_tmp.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        os.replace(csv_tmp, csv_path)
        os.replace(json_tmp, json_path)
    finally:
        csv_tmp.unlink(missing_ok=True)
        json_tmp.unlink(missing_ok=True)
    return csv_path, json_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="DMC-ICE13 benchmark directory",
    )
    args = parser.parse_args()
    try:
        csv_path, json_path = finalize(args.root.resolve())
    except ValueError as error:
        parser.error(str(error))
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()
