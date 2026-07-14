#!/usr/bin/env python3
"""Create the fail-closed X23b GFN1/GFN2/g-XTB paper artifact.

The ordinary X23b collectors intentionally write method-owned working tables.
This script is the publication boundary.  It preserves the frozen GFN1/GFN2
tables, but refuses to add g-XTB until all 23 systems, the derivative gates,
the direct-ACP cross-build gate, and the k333--k444 convergence decision are
complete and cryptographically tied to the raw calculation artifacts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

import x23b_common as common


METHODS = ("GFN1", "GFN2", "GXTB")
METHOD_LABELS = {"GFN1": "GFN1-xTB", "GFN2": "GFN2-xTB", "GXTB": "g-xTB"}
FD_SYSTEMS = (
    "ammonia",
    "14-cyclohexanedione",
    "acetic_acid",
    "ethylcarbamate",
)
FD_PHASE = "x23b_k222_finite_difference_gate_v1"
PREFLIGHT_PHASE = "x23b_experimental_k222_preflight"
CELL_OPT_PHASE = "x23b_k222_cellopt"
CROSS_BUILD_PHASE = "x23b_gxtb_direct_acp_cross_build_v1"
CROSS_BUILD_FINAL_JOB_PHASE = "x23b_direct_acp_cross_build_final"
CROSS_BUILD_FROZEN_JOB_PHASE = "x23b_direct_acp_cross_build_frozen"
KPOINT_APPROVAL_PHASE = "x23b_k333_k444_convergence_v1"
SCHEMA_VERSION = 1
SUMMARY_STEM = "x23b_gfn_gxtb_paper_summary"
HARTREE_TO_KJMOL = 2625.499638
SUMMARY_FIELDS = (
    "method",
    "method_label",
    "quantity",
    "calculation",
    "mesh",
    "N",
    "ME",
    "MAE",
    "RMSE",
    "MaxAE",
)
PREFLIGHT_PROTOCOL = {
    "source_policy": "experimental_reference",
    "variant": "experimental_k222_preflight",
    "mesh": "MACDONALD 2 2 2 0.25 0.25 0.25",
    "symmetry": "SPGLIB reduced",
    "run_type": "ENERGY_FORCE",
    "stress": "ANALYTICAL GPa",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def relative_path(path: Path, root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def artifact(path: Path, root: Path) -> dict[str, object]:
    if not path.is_file():
        raise ValueError(f"required artifact is missing: {path}")
    return {
        "path": relative_path(path, root),
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"required CSV is missing or empty: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"required CSV has no records: {path}")
    return rows


def truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def finite_float(value: object, label: str) -> float:
    if value is None or isinstance(value, bool) or str(value).strip() == "":
        raise ValueError(f"{label} is missing")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is not numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def close(actual: float, expected: float, label: str, tolerance: float = 1.0e-8) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(
            f"{label} mismatch: found {actual:.15g}, expected {expected:.15g}"
        )


def unique_by(
    rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str
) -> dict[tuple[str, ...], dict[str, str]]:
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in fields)
        if key in result:
            raise ValueError(f"duplicate {label} record: {key}")
        result[key] = row
    return result


def stats(errors: list[float]) -> dict[str, float]:
    if not errors:
        raise ValueError("cannot summarize an empty X23b error set")
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(value) for value in errors) / len(errors),
        "RMSE": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "MaxAE": max(abs(value) for value in errors),
    }


def parse_energy(path: Path) -> float:
    text = path.read_text(errors="replace")
    if "PROGRAM ENDED" not in text:
        raise ValueError(f"CP2K output did not end normally: {path}")
    values = re.findall(
        r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.EeDd]+)\s*$",
        text,
        flags=re.M,
    )
    if not values:
        raise ValueError(f"CP2K total energy is missing: {path}")
    return finite_float(values[-1].replace("D", "E").replace("d", "e"), str(path))


def parse_volume(path: Path) -> float:
    values = re.findall(
        r"^\s*CELL\| Volume.*?([-+0-9.EeDd]+)\s*$",
        path.read_text(errors="replace"),
        flags=re.M,
    )
    if not values:
        raise ValueError(f"CP2K cell volume is missing: {path}")
    return finite_float(values[-1].replace("D", "E").replace("d", "e"), str(path))


def resolve_record_path(value: object, label: str) -> Path:
    if not value or not str(value).strip():
        raise ValueError(f"{label} path is missing")
    path = Path(str(value)).expanduser()
    if not path.is_file():
        raise ValueError(f"{label} artifact is missing: {path}")
    return path.resolve()


def validate_hash(path: Path, expected: object, label: str) -> None:
    expected_text = str(expected).lower()
    if not re.fullmatch(r"[0-9a-f]{64}", expected_text):
        raise ValueError(f"{label} has no valid SHA256")
    if sha256(path) != expected_text:
        raise ValueError(f"{label} fingerprint differs: {path}")


def metadata(root: Path) -> tuple[tuple[str, ...], dict[str, dict[str, object]], dict[str, object]]:
    path = root / "data" / "metadata.json"
    payload = read_json(path)
    raw = payload.get("systems")
    if not isinstance(raw, list) or len(raw) != 23:
        raise ValueError("X23b metadata must contain exactly 23 systems")
    records: dict[str, dict[str, object]] = {}
    for item in raw:
        if not isinstance(item, dict):
            raise ValueError("invalid X23b metadata record")
        system = str(item.get("id", ""))
        if not system or system in records:
            raise ValueError(f"duplicate or empty X23b metadata id: {system}")
        for field in ("ref_energy", "x23b_same_cell_ref_volume", "molecules_per_cell"):
            finite_float(item.get(field), f"metadata {system}/{field}")
        records[system] = item
    return tuple(sorted(records)), records, artifact(path, root)


def campaign_identity_from_manifest(path: Path) -> dict[str, object]:
    identity, state = common.declared_campaign_identity(path)
    if state != "production_ready":
        raise ValueError(f"campaign manifest is not production_ready: {path}")
    return identity


def locate_campaign_manifest(
    root: Path, record: Mapping[str, object], expected_identity: Mapping[str, object]
) -> tuple[Path, dict[str, object]]:
    expected_hash = str(record.get("file_sha256", "")).lower()
    if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
        raise ValueError("GXTB provenance lacks a valid campaign-manifest hash")
    campaign_id = str(expected_identity.get("campaign_id", ""))
    candidates = [
        Path(str(record.get("path", ""))),
        root.parent / "campaigns" / campaign_id / "build_manifest.json",
    ]
    campaigns = root.parent / "campaigns"
    if campaigns.is_dir():
        candidates.extend(sorted(campaigns.rglob("build_manifest.json")))
    path = next(
        (
            candidate
            for candidate in dict.fromkeys(candidates)
            if candidate.is_file() and sha256(candidate) == expected_hash
        ),
        None,
    )
    if path is None:
        raise ValueError("the exact GXTB campaign manifest is not locally available")
    identity = campaign_identity_from_manifest(path)
    expected = dict(expected_identity)
    if identity != expected:
        raise ValueError("campaign manifest identity differs from GXTB provenance")
    return path.resolve(), artifact(path, root)


def validate_build_provenance(
    root: Path,
) -> tuple[dict[str, object], dict[str, object], Path, dict[str, object], dict[str, object]]:
    legacy_path = root / "data" / "build_provenance.json"
    gxtb_path = root / "data" / "build_provenance_gxtb.json"
    legacy = read_json(legacy_path)
    gxtb = read_json(gxtb_path)
    if gxtb.get("method") not in (None, "GXTB"):
        raise ValueError("foreign method in GXTB build provenance")
    if gxtb.get("status") != "production_complete":
        raise ValueError(
            "GXTB build provenance is not production_complete "
            f"(found {gxtb.get('status', '<missing>')})"
        )
    validation = gxtb.get("validation")
    if not isinstance(validation, dict):
        raise ValueError("GXTB build provenance lacks validation coverage")
    required_counts = (
        "gas_optimizations",
        "experimental_k222_preflight",
        "k222_cell_optimizations",
        "k333_single_points",
        "k444_single_points",
    )
    for name in required_counts:
        if int(validation.get(f"{name}_expected", -1)) != 23 or int(
            validation.get(f"{name}_completed", -1)
        ) != 23:
            raise ValueError(f"GXTB provenance does not certify 23/23 {name}")
    workflow = gxtb.get("workflow_paths")
    if not isinstance(workflow, dict):
        raise ValueError("GXTB provenance lacks workflow paths")
    if workflow.get("k222_source_policy") != "experimental_reference":
        raise ValueError("GXTB k222 source policy is not experimental_reference")
    roots = workflow.get("final_single_point_roots")
    if not isinstance(roots, dict) or not roots.get("k333") or not roots.get("k444"):
        raise ValueError("GXTB provenance lacks distinct k333 and k444 roots")
    if Path(str(roots["k333"])).resolve() == Path(str(roots["k444"])).resolve():
        raise ValueError("GXTB k333 and k444 roots are not distinct")
    campaign = gxtb.get("campaign_identity")
    if not isinstance(campaign, dict):
        raise ValueError("GXTB provenance lacks campaign identity")
    common.validate_campaign_identity(campaign)
    manifest_record = gxtb.get("campaign_manifest")
    if not isinstance(manifest_record, dict):
        raise ValueError("GXTB provenance lacks campaign manifest record")
    manifest_path, manifest_artifact = locate_campaign_manifest(root, manifest_record, campaign)
    lineage = {
        "GFN1_GFN2": {
            "artifact": artifact(legacy_path, root),
            "cp2k": legacy.get("cp2k"),
            "provider_name": "tblite",
            "provider": legacy.get("tblite"),
            "repository_patches": legacy.get("repository_patches"),
            "raw_job_stamps": "not available for frozen legacy campaign; curated tables are hash-frozen",
        },
        "GXTB": {
            "artifact": artifact(gxtb_path, root),
            "campaign_identity": campaign,
            "campaign_manifest": manifest_artifact,
            "cp2k": gxtb.get("cp2k"),
            "provider_name": "save_tblite",
            "provider": gxtb.get("save_tblite"),
        },
    }
    return gxtb, campaign, manifest_path, manifest_artifact, lineage


def portable_job(
    root: Path,
    input_path: Path,
    output_path: Path,
    stamp_path: Path,
    campaign: Mapping[str, object],
    expected_phase: str | None,
    *,
    require_optimization: bool = False,
    expected_protocol: Mapping[str, object] | None = None,
) -> dict[str, object]:
    for path, label in (
        (input_path, "job input"),
        (output_path, "job output"),
        (stamp_path, "job stamp"),
    ):
        if not path.is_file():
            raise ValueError(f"{label} is missing: {path}")
    text = output_path.read_text(errors="replace")
    if "PROGRAM ENDED" not in text or "ENERGY| Total FORCE_EVAL" not in text:
        raise ValueError(f"incomplete CP2K output: {output_path}")
    if require_optimization and "GEOMETRY OPTIMIZATION COMPLETED" not in text:
        raise ValueError(f"unconverged optimization output: {output_path}")
    stamp = read_json(stamp_path)
    if stamp.get("schema") != common.JOB_STAMP_SCHEMA or stamp.get("method") != "GXTB":
        raise ValueError(f"invalid GXTB job stamp: {stamp_path}")
    if expected_phase is not None and stamp.get("phase") != expected_phase:
        raise ValueError(f"job phase differs in {stamp_path}")
    if stamp.get("campaign_identity") != dict(campaign):
        raise ValueError(f"job campaign identity differs in {stamp_path}")
    recorded_input = stamp.get("input")
    recorded_cp2k = stamp.get("cp2k")
    details = stamp.get("details")
    if not isinstance(recorded_input, dict) or recorded_input.get("sha256") != sha256(input_path):
        raise ValueError(f"job input fingerprint differs in {stamp_path}")
    if (
        not isinstance(recorded_cp2k, dict)
        or recorded_cp2k.get("sha256") != campaign.get("cp2k_executable_sha256")
    ):
        raise ValueError(f"job CP2K fingerprint differs in {stamp_path}")
    if not str(stamp.get("status", "")).startswith("converged"):
        raise ValueError(f"job stamp is not converged: {stamp_path}")
    if not isinstance(details, dict) or details.get("output_sha256") != sha256(output_path):
        raise ValueError(f"job output fingerprint differs in {stamp_path}")
    if expected_protocol is not None and stamp.get("protocol_identity") != dict(expected_protocol):
        raise ValueError(f"job protocol identity differs in {stamp_path}")
    source_lineage: dict[str, object] = {}
    raw_sources = stamp.get("source_artifacts")
    if raw_sources is not None:
        if not isinstance(raw_sources, dict):
            raise ValueError(f"invalid source-artifact record in {stamp_path}")
        for role, record in sorted(raw_sources.items()):
            if not isinstance(record, dict):
                raise ValueError(f"invalid {role} source record in {stamp_path}")
            source = resolve_record_path(record.get("path"), f"{role} source")
            validate_hash(source, record.get("sha256"), f"{role} source")
            source_lineage[str(role)] = artifact(source, root)
    return {
        "phase": stamp.get("phase"),
        "status": stamp.get("status"),
        "input": artifact(input_path, root),
        "output": artifact(output_path, root),
        "stamp": artifact(stamp_path, root),
        "source_artifacts": source_lineage,
    }


def require_lineage_hashes(
    job: Mapping[str, object], expected: Mapping[str, str], label: str
) -> None:
    raw = job.get("source_artifacts")
    if not isinstance(raw, dict) or set(raw) != set(expected):
        raise ValueError(f"{label} source-artifact roles differ")
    for role, expected_hash in expected.items():
        record = raw.get(role)
        if not isinstance(record, dict) or record.get("sha256") != expected_hash:
            raise ValueError(f"{label} {role} source fingerprint differs")


def validate_legacy(
    root: Path,
    systems: tuple[str, ...],
    meta: Mapping[str, Mapping[str, object]],
    rows_path: Path,
    volumes_path: Path,
) -> tuple[dict[str, dict[str, list[float]]], dict[str, object]]:
    rows = read_csv(rows_path)
    by_key = unique_by(rows, ("method", "system"), "legacy final-kpoint")
    expected = {(method, system) for method in ("GFN1", "GFN2") for system in systems}
    if set(by_key) != expected:
        raise ValueError(
            "legacy final-kpoint coverage differs: "
            f"missing={sorted(expected - set(by_key))}, unexpected={sorted(set(by_key) - expected)}"
        )
    volume_rows = read_csv(volumes_path)
    volume_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in volume_rows:
        label = row.get("method", "")
        method = {"GFN1-xTB": "GFN1", "GFN2-xTB": "GFN2"}.get(label, label)
        key = (method, row.get("system", ""))
        if key in volume_by_key:
            raise ValueError(f"duplicate legacy volume record: {key}")
        volume_by_key[key] = row
    if set(volume_by_key) != expected:
        raise ValueError("legacy volume coverage is not exactly 23/23 for GFN1 and GFN2")
    values: dict[str, dict[str, list[float]]] = {}
    convergence: dict[str, object] = {}
    for method in ("GFN1", "GFN2"):
        lattice_errors: list[float] = []
        volume_errors: list[float] = []
        deltas: list[float] = []
        for system in systems:
            row = by_key[(method, system)]
            reference = finite_float(meta[system]["ref_energy"], f"{system} reference energy")
            k333 = finite_float(row.get("k333_lattice_energy_kJmol"), f"{method}/{system} k333")
            k444 = finite_float(row.get("k444_lattice_energy_kJmol"), f"{method}/{system} k444")
            error = k333 - reference
            close(error, finite_float(row.get("k333_error_kJmol"), f"{method}/{system} k333 error"), f"{method}/{system} k333 error", 2.0e-9)
            delta = k444 - k333
            close(delta, finite_float(row.get("delta_k444_minus_k333_kJmol"), f"{method}/{system} k delta"), f"{method}/{system} k delta", 2.0e-9)
            vrow = volume_by_key[(method, system)]
            if not truth(vrow.get("complete")) or vrow.get("calculation") != "cell_opt" or vrow.get("mesh") != "k222":
                raise ValueError(f"incomplete legacy k222 volume for {method}/{system}")
            volume = finite_float(vrow.get("volume_A3"), f"{method}/{system} volume")
            volume_ref = finite_float(meta[system]["x23b_same_cell_ref_volume"], f"{system} volume reference")
            volume_error = 100.0 * (volume - volume_ref) / volume_ref
            close(volume_error, finite_float(vrow.get("volume_error_percent"), f"{method}/{system} volume error"), f"{method}/{system} volume error", 1.0e-5)
            lattice_errors.append(error)
            volume_errors.append(volume_error)
            deltas.append(delta)
        values[method] = {"lattice": lattice_errors, "volume": volume_errors}
        convergence[method] = {
            "N": 23,
            "mean_abs_change_kJmol": sum(abs(value) for value in deltas) / 23,
            "max_abs_change_kJmol": max(abs(value) for value in deltas),
            "per_system_delta_k444_minus_k333_kJmol": dict(zip(systems, deltas)),
        }
    return values, {
        "tables": {
            "final_kpoint_rows": artifact(rows_path, root),
            "cell_volumes": artifact(volumes_path, root),
        },
        "k333_to_k444": convergence,
    }


def validate_preflight(
    root: Path,
    systems: tuple[str, ...],
    campaign: Mapping[str, object],
    preflight_csv: Path,
    preflight_root: Path,
) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    rows = read_csv(preflight_csv)
    by_key = unique_by(rows, ("method", "system"), "preflight")
    expected = {("GXTB", system) for system in systems}
    if set(by_key) != expected:
        raise ValueError("experimental k222 preflight is not exactly 23/23")
    manifest_path = preflight_root / "experimental_k222_preflight_manifest.json"
    manifest = read_json(manifest_path)
    if (
        manifest.get("schema") != 1
        or manifest.get("phase") != PREFLIGHT_PHASE
        or manifest.get("source_policy") != "experimental_reference"
        or manifest.get("campaign_identity") != dict(campaign)
    ):
        raise ValueError("invalid experimental k222 preflight manifest")
    raw_records = manifest.get("systems")
    if not isinstance(raw_records, list):
        raise ValueError("preflight manifest has no system records")
    records = {str(item.get("system", "")): item for item in raw_records if isinstance(item, dict)}
    if set(records) != set(systems) or len(records) != 23:
        raise ValueError("preflight manifest is not exactly 23/23")
    lineage: dict[str, dict[str, object]] = {}
    for system in systems:
        row = by_key[("GXTB", system)]
        record = records[system]
        if (
            row.get("phase") != PREFLIGHT_PHASE
            or row.get("source_policy") != "experimental_reference"
            or not truth(row.get("program_ended"))
            or row.get("scientific_status") != "measured_not_approved"
            or truth(row.get("approved"))
            or row.get("campaign_fingerprint_sha256") != campaign.get("fingerprint_sha256")
        ):
            raise ValueError(f"invalid preflight status for GXTB/{system}")
        for field in (
            "energy_hartree",
            "max_force_hartree_per_bohr",
            "max_abs_stress_GPa",
            "pressure_GPa",
        ):
            finite_float(row.get(field), f"preflight {system}/{field}")
        input_path = resolve_record_path(row.get("input"), f"preflight {system} input")
        output_path = resolve_record_path(row.get("output"), f"preflight {system} output")
        validate_hash(input_path, row.get("input_sha256"), f"preflight {system} input")
        validate_hash(output_path, row.get("output_sha256"), f"preflight {system} output")
        for path_field, hash_field in (
            ("source_input", "source_input_sha256"),
            ("structure_path", "structure_sha256"),
        ):
            path = resolve_record_path(row.get(path_field), f"preflight {system}/{path_field}")
            validate_hash(path, row.get(hash_field), f"preflight {system}/{path_field}")
            if str(record.get(path_field)) != str(path) or record.get(hash_field) != row.get(hash_field):
                raise ValueError(f"preflight CSV/manifest lineage differs for {system}")
        if (
            str(record.get("input")) != str(input_path)
            or record.get("input_sha256") != row.get("input_sha256")
            or str(record.get("output")) != str(output_path)
        ):
            raise ValueError(f"preflight CSV/manifest raw paths differ for {system}")
        stamp_path = output_path.parent / common.JOB_STAMP_NAME
        job = portable_job(
            root,
            input_path,
            output_path,
            stamp_path,
            campaign,
            PREFLIGHT_PHASE,
            expected_protocol=PREFLIGHT_PROTOCOL,
        )
        require_lineage_hashes(
            job,
            {
                "reference_input": str(row["source_input_sha256"]),
                "reference_structure": str(row["structure_sha256"]),
            },
            f"preflight {system}",
        )
        lineage[system] = {
            **job,
            "source_input": artifact(resolve_record_path(row["source_input"], "source input"), root),
            "reference_structure": artifact(resolve_record_path(row["structure_path"], "reference structure"), root),
        }
    return lineage, {
        "phase": PREFLIGHT_PHASE,
        "decision": "accepted_by_exact_k222_cellopt_lineage",
        "coverage_expected": 23,
        "coverage_completed": 23,
        "coverage_accepted": 23,
        "scientific_status_before_acceptance": "measured_not_approved",
        "table": artifact(preflight_csv, root),
        "manifest": artifact(manifest_path, root),
    }


def validate_fd_gate(
    root: Path,
    campaign: Mapping[str, object],
    report_path: Path,
    approval_path: Path,
) -> dict[str, object]:
    approval = read_json(approval_path)
    report = read_json(report_path)
    if (
        approval.get("schema") != 1
        or approval.get("phase") != FD_PHASE
        or approval.get("decision") != "approved"
        or not str(approval.get("reviewer", "")).strip()
    ):
        raise ValueError("finite-difference pilot gate is not explicitly approved")
    if (
        report.get("schema") != 1
        or report.get("phase") != FD_PHASE
        or report.get("scientific_status") != "measured_not_approved"
        or report.get("approved") is not False
        or report.get("campaign_identity") != dict(campaign)
        or tuple(report.get("systems", [])) != FD_SYSTEMS
    ):
        raise ValueError("invalid finite-difference measurement report")
    report_record = approval.get("report_json")
    if not isinstance(report_record, dict):
        raise ValueError("FD approval lacks its report record")
    validate_hash(report_path, report_record.get("sha256"), "FD report")
    measured_record = report.get("measured_csv")
    manifest_record = report.get("manifest")
    if not isinstance(measured_record, dict) or not isinstance(manifest_record, dict):
        raise ValueError("FD report lacks measured CSV/manifest records")
    measured = resolve_record_path(measured_record.get("path"), "FD measured CSV")
    manifest_path = resolve_record_path(manifest_record.get("path"), "FD manifest")
    validate_hash(measured, measured_record.get("sha256"), "FD measured CSV")
    validate_hash(manifest_path, manifest_record.get("sha256"), "FD manifest")
    if approval.get("measured_csv") != measured_record or approval.get("manifest") != manifest_record:
        raise ValueError("FD approval/report artifact records differ")
    report_rows = report.get("rows")
    if not isinstance(report_rows, list) or len(report_rows) != int(report.get("row_count", -1)):
        raise ValueError("FD report row count is inconsistent")
    measured_rows = read_csv(measured)
    if len(measured_rows) != len(report_rows) or len(report_rows) not in (12, 16):
        raise ValueError("FD pilot is incomplete")
    with measured.open(newline="") as handle:
        measured_reader = csv.DictReader(handle)
        measured_fields = measured_reader.fieldnames
    if not measured_fields:
        raise ValueError("FD measured CSV has no header")
    import io

    rendered = io.StringIO(newline="")
    writer = csv.DictWriter(rendered, fieldnames=measured_fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(report_rows)
    if rendered.getvalue() != measured.read_text():
        raise ValueError("FD JSON rows do not exactly reproduce the measured CSV")
    checks = approval.get("checks")
    if (
        not isinstance(checks, list)
        or len(checks) != len(report_rows)
        or not all(isinstance(check, dict) and check.get("passed") is True for check in checks)
        or int(approval.get("passed_count", -1)) != len(checks)
        or int(approval.get("check_count", -1)) != len(checks)
    ):
        raise ValueError("FD approval checks are incomplete or failed")
    manifest = read_json(manifest_path)
    if (
        manifest.get("schema") != 1
        or manifest.get("phase") != FD_PHASE
        or manifest.get("campaign_identity") != dict(campaign)
        or manifest.get("scientific_status") != "prepared_not_measured"
    ):
        raise ValueError("invalid FD manifest")
    frozen_payload_hash = manifest.get("payload_sha256")
    unhashed = dict(manifest)
    unhashed.pop("payload_sha256", None)
    if frozen_payload_hash != fingerprint(unhashed):
        raise ValueError("FD manifest payload fingerprint is internally inconsistent")
    cases = manifest.get("cases")
    if not isinstance(cases, list) or tuple(str(case.get("system", "")) for case in cases if isinstance(case, dict)) != FD_SYSTEMS:
        raise ValueError("FD manifest does not contain the frozen four-system pilot")
    raw_jobs: dict[str, object] = {}
    protocol = manifest.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("FD manifest lacks protocol identity")
    for case in cases:
        assert isinstance(case, dict)
        system = str(case["system"])
        jobs = case.get("jobs")
        if not isinstance(jobs, list) or len(jobs) not in (7, 9):
            raise ValueError(f"incomplete FD job set for {system}")
        system_jobs: dict[str, object] = {}
        for job_record in jobs:
            if not isinstance(job_record, dict):
                raise ValueError(f"invalid FD job record for {system}")
            job_id = str(job_record.get("job_id", ""))
            input_path = resolve_record_path(job_record.get("input"), f"FD {system}/{job_id} input")
            output_path = resolve_record_path(job_record.get("output"), f"FD {system}/{job_id} output")
            validate_hash(input_path, job_record.get("input_sha256"), f"FD {system}/{job_id} input")
            job_protocol = {
                "gate_phase": FD_PHASE,
                "variant": "frozen_reference_shifted_k222_spglib_fd",
                "source_policy": "experimental_reference",
                "manifest_payload_sha256": manifest["payload_sha256"],
                "mesh": protocol["mesh"],
                "coordinate_step_bohr": protocol["coordinate_step_bohr"],
                "strain_step": protocol["strain_step"],
                "job_id": job_record["job_id"],
                "job_type": job_record["job_type"],
                "direction_id": job_record.get("direction_id"),
                "direction_sha256": job_record.get("direction_sha256"),
                "generator": job_record.get("generator"),
                "sign": job_record.get("sign"),
            }
            validated_job = portable_job(
                root,
                input_path,
                output_path,
                output_path.parent / common.JOB_STAMP_NAME,
                campaign,
                FD_PHASE,
                expected_protocol=job_protocol,
            )
            require_lineage_hashes(
                validated_job,
                {
                    "fd_manifest": sha256(manifest_path),
                    "reference_input": str(case["source_input_sha256"]),
                    "reference_structure": str(case["structure_sha256"]),
                },
                f"FD {system}/{job_id}",
            )
            system_jobs[job_id] = validated_job
        raw_jobs[system] = system_jobs
    return {
        "phase": FD_PHASE,
        "decision": "approved",
        "reviewer": approval["reviewer"],
        "pilot_systems": list(FD_SYSTEMS),
        "measurement_count": len(report_rows),
        "passed_count": len(checks),
        "thresholds": approval.get("thresholds"),
        "approval": artifact(approval_path, root),
        "report": artifact(report_path, root),
        "measured_csv": artifact(measured, root),
        "manifest": artifact(manifest_path, root),
        "raw_jobs": raw_jobs,
    }


def validate_gxtb_results(
    root: Path,
    systems: tuple[str, ...],
    meta: Mapping[str, Mapping[str, object]],
    campaign: Mapping[str, object],
    preflight: Mapping[str, Mapping[str, object]],
    provenance: Mapping[str, object],
    cellopt_csv: Path,
    k333_csv: Path,
    k444_csv: Path,
) -> tuple[dict[str, list[float]], dict[str, object], dict[str, object]]:
    workflow = provenance["workflow_paths"]
    assert isinstance(workflow, dict)
    cellopt_root = Path(str(workflow["k222_cellopt_root"])).resolve()
    final_roots = workflow["final_single_point_roots"]
    assert isinstance(final_roots, dict)
    roots = {mesh: Path(str(final_roots[mesh])).resolve() for mesh in ("k333", "k444")}
    cell_rows = unique_by(read_csv(cellopt_csv), ("method", "system"), "GXTB cellopt")
    expected = {("GXTB", system) for system in systems}
    if set(cell_rows) != expected:
        raise ValueError("GXTB k222 CELL_OPT coverage is not exactly 23/23")
    k_rows = {
        "k333": unique_by(read_csv(k333_csv), ("method", "system"), "GXTB k333"),
        "k444": unique_by(read_csv(k444_csv), ("method", "system"), "GXTB k444"),
    }
    if any(set(rows) != expected for rows in k_rows.values()):
        raise ValueError("GXTB k333/k444 coverage is not exactly 23/23")
    cell_manifest_path = cellopt_root.parent / "x23b_k222_cellopt_manifest.csv"
    cell_manifest = unique_by(read_csv(cell_manifest_path), ("method", "system"), "GXTB cellopt manifest")
    if set(cell_manifest) != expected:
        raise ValueError("GXTB cellopt manifest is not exactly 23/23")
    final_manifests: dict[str, dict[tuple[str, ...], dict[str, str]]] = {}
    for mesh in ("k333", "k444"):
        manifest = unique_by(read_csv(roots[mesh] / "manifest.csv"), ("method", "system"), f"GXTB {mesh} manifest")
        if set(manifest) != expected:
            raise ValueError(f"GXTB {mesh} manifest is not exactly 23/23")
        final_manifests[mesh] = manifest

    lattice_errors: list[float] = []
    volume_errors: list[float] = []
    deltas: list[float] = []
    lineage: dict[str, object] = {}
    for system in systems:
        c_row = cell_rows[("GXTB", system)]
        c_manifest = cell_manifest[("GXTB", system)]
        if (
            c_row.get("source_policy") != "experimental_reference"
            or c_manifest.get("source_policy") != "experimental_reference"
            or c_manifest.get("lineage")
            != "frozen_x23_reference_structure->experimental_k222_preflight->k222_cellopt_input"
            or not truth(c_row.get("program_ended"))
            or not truth(c_row.get("opt_completed"))
            or truth(c_row.get("max_iter_reached"))
        ):
            raise ValueError(f"invalid GXTB k222 CELL_OPT status/lineage for {system}")
        pf = preflight[system]
        for manifest_field, lineage_field in (
            ("preflight_input_sha256", "input"),
            ("preflight_output_sha256", "output"),
            ("preflight_stamp_sha256", "stamp"),
        ):
            record = pf[lineage_field]
            assert isinstance(record, dict)
            if c_manifest.get(manifest_field) != record.get("sha256") or c_row.get(manifest_field) != record.get("sha256"):
                raise ValueError(f"k222 cellopt does not bind accepted preflight {lineage_field} for {system}")
        c_input = resolve_record_path(c_manifest.get("input"), f"GXTB/{system} k222 input")
        c_output = resolve_record_path(c_row.get("output"), f"GXTB/{system} k222 output")
        validate_hash(c_input, c_manifest.get("input_sha256"), f"GXTB/{system} k222 input")
        c_job = portable_job(
            root,
            c_input,
            c_output,
            c_output.parent / common.JOB_STAMP_NAME,
            campaign,
            CELL_OPT_PHASE,
            require_optimization=True,
            expected_protocol={
                "source_policy": "experimental_reference",
                "variant": c_manifest["variant"],
                "lineage": c_manifest["lineage"],
            },
        )
        require_lineage_hashes(
            c_job,
            {
                "reference_input": str(c_manifest["source_sha256"]),
                "reference_structure": str(c_manifest["structure_sha256"]),
                "preflight_input": str(c_manifest["preflight_input_sha256"]),
                "preflight_output": str(c_manifest["preflight_output_sha256"]),
                "preflight_stamp": str(c_manifest["preflight_stamp_sha256"]),
            },
            f"GXTB/{system} k222 cellopt",
        )
        restarts = list(c_output.parent.glob("*-1.restart")) + list(c_output.parent.glob("*-1_*.restart"))
        if not restarts:
            raise ValueError(f"GXTB/{system} final k222 restart is missing")
        c_restart = max(restarts, key=lambda path: path.stat().st_mtime_ns)
        gas_stem = f"{system}_GXTB_mol_geoopt"
        gas_dir = root / "runs" / "molecule_geoopt" / "GXTB" / gas_stem
        gas_input = gas_dir / f"{gas_stem}.inp"
        gas_output = gas_dir / f"{gas_stem}.out"
        gas_job = portable_job(
            root,
            gas_input,
            gas_output,
            gas_dir / common.JOB_STAMP_NAME,
            campaign,
            "x23b_molecule_geoopt",
            require_optimization=True,
        )
        gas_restarts = list(gas_dir.glob("*-1.restart")) + list(gas_dir.glob("*-1_*.restart"))
        if not gas_restarts:
            raise ValueError(f"GXTB/{system} gas optimization restart is missing")
        gas_restart = max(gas_restarts, key=lambda path: path.stat().st_mtime_ns)
        gas_energy = parse_energy(gas_output)
        close(gas_energy, finite_float(c_row.get("gas_energy_hartree"), f"GXTB/{system} gas energy"), f"GXTB/{system} gas energy", 1.0e-10)
        cell_energy = parse_energy(c_output)
        close(
            cell_energy,
            finite_float(c_row.get("energy_hartree"), f"GXTB/{system} k222 energy"),
            f"GXTB/{system} k222 raw/table energy",
            1.0e-10,
        )
        raw_volume = parse_volume(c_output)
        close(
            raw_volume,
            finite_float(c_row.get("volume_A3"), f"GXTB/{system} volume"),
            f"GXTB/{system} raw/table volume",
            1.0e-6,
        )
        n_molecules = int(meta[system]["molecules_per_cell"])
        source_lattice = (gas_energy - cell_energy / n_molecules) * HARTREE_TO_KJMOL
        close(
            source_lattice,
            finite_float(c_row.get("lattice_energy_kJmol"), f"GXTB/{system} k222 lattice"),
            f"GXTB/{system} k222 lattice",
            2.0e-6,
        )

        mesh_lineage: dict[str, object] = {}
        lattice: dict[str, float] = {}
        for mesh in ("k333", "k444"):
            result = k_rows[mesh][("GXTB", system)]
            manifest = final_manifests[mesh][("GXTB", system)]
            target = int(mesh[1])
            if (
                result.get("target_mesh") != mesh
                or not truth(result.get("program_ended"))
                or int(manifest.get("mesh", -1)) != target
                or manifest.get("source_policy") != "experimental_reference"
            ):
                raise ValueError(f"invalid GXTB {mesh} status/manifest for {system}")
            for field, actual in (
                ("source_input_sha256", sha256(c_input)),
                ("source_output_sha256", sha256(c_output)),
                ("source_restart_sha256", sha256(c_restart)),
            ):
                if manifest.get(field) != actual:
                    raise ValueError(f"GXTB {mesh} geometry/restart lineage differs for {system}: {field}")
            target_input = resolve_record_path(manifest.get("input"), f"GXTB/{system} {mesh} input")
            target_output = resolve_record_path(result.get("output"), f"GXTB/{system} {mesh} output")
            validate_hash(target_input, manifest.get("input_sha256"), f"GXTB/{system} {mesh} input")
            try:
                source_protocol = json.loads(manifest["source_protocol_identity"])
                source_hashes = json.loads(manifest["source_artifact_hashes"])
            except (KeyError, TypeError, json.JSONDecodeError) as error:
                raise ValueError(f"invalid GXTB {mesh} source protocol for {system}") from error
            if not isinstance(source_protocol, dict) or not isinstance(source_hashes, dict):
                raise ValueError(f"invalid GXTB {mesh} source lineage for {system}")
            expected_protocol = {
                "source_policy": "experimental_reference",
                "source_variant": manifest["source_variant"],
                "source_protocol_identity": source_protocol,
                "target_mesh": mesh,
            }
            target_job = portable_job(
                root,
                target_input,
                target_output,
                target_output.parent / common.JOB_STAMP_NAME,
                campaign,
                f"x23b_final_k{target}{target}{target}_on_k222",
                expected_protocol=expected_protocol,
            )
            require_lineage_hashes(
                target_job,
                {
                    **{str(name): str(value) for name, value in source_hashes.items()},
                    "cellopt_input": sha256(c_input),
                    "cellopt_output": sha256(c_output),
                    "cellopt_restart": sha256(c_restart),
                },
                f"GXTB/{system} {mesh}",
            )
            target_energy = parse_energy(target_output)
            table_energy = finite_float(result.get("target_energy_hartree"), f"GXTB/{system} {mesh} energy")
            close(target_energy, table_energy, f"GXTB/{system} {mesh} raw/table energy", 1.0e-10)
            lattice_value = (gas_energy - target_energy / n_molecules) * HARTREE_TO_KJMOL
            close(lattice_value, finite_float(result.get("target_lattice_energy_kJmol"), f"GXTB/{system} {mesh} lattice"), f"GXTB/{system} {mesh} lattice", 2.0e-6)
            lattice[mesh] = lattice_value
            close(
                cell_energy,
                finite_float(result.get("source_energy_hartree"), f"GXTB/{system} {mesh} source energy"),
                f"GXTB/{system} {mesh} source energy",
                1.0e-10,
            )
            close(
                source_lattice,
                finite_float(result.get("source_lattice_energy_kJmol"), f"GXTB/{system} {mesh} source lattice"),
                f"GXTB/{system} {mesh} source lattice",
                2.0e-6,
            )
            close(
                lattice_value - source_lattice,
                finite_float(result.get("delta_target_minus_source_kJmol"), f"GXTB/{system} {mesh} target-source delta"),
                f"GXTB/{system} {mesh} target-source delta",
                2.0e-6,
            )
            close(
                lattice_value - finite_float(meta[system]["ref_energy"], f"{system} reference energy"),
                finite_float(result.get("target_error_kJmol"), f"GXTB/{system} {mesh} error"),
                f"GXTB/{system} {mesh} error",
                2.0e-6,
            )
            mesh_lineage[mesh] = target_job
        reference = finite_float(meta[system]["ref_energy"], f"{system} lattice reference")
        lattice_error = lattice["k333"] - reference
        close(lattice_error, finite_float(k_rows["k333"][("GXTB", system)].get("target_error_kJmol"), f"GXTB/{system} k333 error"), f"GXTB/{system} k333 error", 2.0e-6)
        volume = raw_volume
        volume_ref = finite_float(meta[system]["x23b_same_cell_ref_volume"], f"{system} volume reference")
        volume_error = 100.0 * (volume - volume_ref) / volume_ref
        close(volume_error, finite_float(c_row.get("volume_error_percent"), f"GXTB/{system} volume error"), f"GXTB/{system} volume error", 2.0e-5)
        delta = lattice["k444"] - lattice["k333"]
        lattice_errors.append(lattice_error)
        volume_errors.append(volume_error)
        deltas.append(delta)
        lineage[system] = {
            "gas_optimization": {**gas_job, "final_restart": artifact(gas_restart, root)},
            "preflight": preflight[system],
            "k222_cell_optimization": {**c_job, "final_restart": artifact(c_restart, root)},
            "final_single_points": mesh_lineage,
            "reported_lattice_energy_kJmol": lattice["k333"],
            "volume_error_percent": volume_error,
            "delta_k444_minus_k333_kJmol": delta,
        }
    convergence = {
        "N": 23,
        "mean_abs_change_kJmol": sum(abs(value) for value in deltas) / 23,
        "max_abs_change_kJmol": max(abs(value) for value in deltas),
        "per_system_delta_k444_minus_k333_kJmol": dict(zip(systems, deltas)),
    }
    tables = {
        "cellopt": artifact(cellopt_csv, root),
        "k333": artifact(k333_csv, root),
        "k444": artifact(k444_csv, root),
        "cellopt_manifest": artifact(cell_manifest_path, root),
        "k333_manifest": artifact(roots["k333"] / "manifest.csv", root),
        "k444_manifest": artifact(roots["k444"] / "manifest.csv", root),
    }
    return {"lattice": lattice_errors, "volume": volume_errors}, convergence, {
        "tables": tables,
        "systems": lineage,
    }


def validate_cross_build_gate(
    root: Path,
    systems: tuple[str, ...],
    final_manifest_path: Path,
    final_campaign: Mapping[str, object],
    approval_path: Path,
) -> dict[str, object]:
    approval = read_json(approval_path)
    if (
        approval.get("schema") != 1
        or approval.get("phase") != CROSS_BUILD_PHASE
        or approval.get("decision") != "approved"
        or not str(approval.get("reviewer", "")).strip()
    ):
        raise ValueError("direct-ACP cross-build gate is not explicitly approved")
    comparison_record = approval.get("comparison_csv")
    final_record = approval.get("final_campaign_manifest")
    frozen_record = approval.get("frozen_campaign_manifest")
    if not all(isinstance(record, dict) for record in (comparison_record, final_record, frozen_record)):
        raise ValueError("cross-build approval lacks comparison/build records")
    assert isinstance(comparison_record, dict)
    assert isinstance(final_record, dict)
    assert isinstance(frozen_record, dict)
    comparison_path = resolve_record_path(comparison_record.get("path"), "cross-build comparison CSV")
    frozen_manifest_path = resolve_record_path(frozen_record.get("path"), "frozen cross-build manifest")
    validate_hash(comparison_path, comparison_record.get("sha256"), "cross-build comparison CSV")
    validate_hash(final_manifest_path, final_record.get("sha256"), "final cross-build manifest")
    validate_hash(frozen_manifest_path, frozen_record.get("sha256"), "frozen cross-build manifest")
    if sha256(final_manifest_path) == sha256(frozen_manifest_path):
        raise ValueError("cross-build gate compares the final build with itself")
    if campaign_identity_from_manifest(final_manifest_path) != dict(final_campaign):
        raise ValueError("cross-build final manifest is not the production campaign")
    frozen_campaign = campaign_identity_from_manifest(frozen_manifest_path)
    if frozen_campaign == dict(final_campaign):
        raise ValueError("cross-build final and frozen campaign identities are identical")
    rows = read_csv(comparison_path)
    by_key = unique_by(rows, ("method", "system"), "cross-build comparison")
    expected = {("GXTB", system) for system in systems}
    if set(by_key) != expected:
        raise ValueError("direct-ACP cross-build comparison is not exactly 23/23")
    threshold = finite_float(approval.get("absolute_energy_tolerance_hartree"), "cross-build tolerance")
    if threshold <= 0.0:
        raise ValueError("cross-build tolerance must be positive")
    deltas: list[float] = []
    lineage: dict[str, object] = {}
    for system in systems:
        row = by_key[("GXTB", system)]
        if not truth(row.get("passed")):
            raise ValueError(f"cross-build comparison failed for {system}")
        jobs: dict[str, object] = {}
        energies: dict[str, float] = {}
        for role, campaign, job_phase in (
            ("final", final_campaign, CROSS_BUILD_FINAL_JOB_PHASE),
            ("frozen", frozen_campaign, CROSS_BUILD_FROZEN_JOB_PHASE),
        ):
            input_path = resolve_record_path(row.get(f"{role}_input"), f"cross-build {system} {role} input")
            output_path = resolve_record_path(row.get(f"{role}_output"), f"cross-build {system} {role} output")
            stamp_path = resolve_record_path(row.get(f"{role}_stamp"), f"cross-build {system} {role} stamp")
            for path, suffix in ((input_path, "input"), (output_path, "output"), (stamp_path, "stamp")):
                validate_hash(path, row.get(f"{role}_{suffix}_sha256"), f"cross-build {system} {role} {suffix}")
            jobs[role] = portable_job(
                root, input_path, output_path, stamp_path, campaign, job_phase
            )
            energies[role] = parse_energy(output_path)
        delta = energies["final"] - energies["frozen"]
        close(delta, finite_float(row.get("delta_final_minus_frozen_hartree"), f"cross-build {system} delta"), f"cross-build {system} delta", 1.0e-12)
        if abs(delta) > threshold:
            raise ValueError(f"cross-build energy delta exceeds approval tolerance for {system}")
        deltas.append(delta)
        lineage[system] = {**jobs, "delta_final_minus_frozen_hartree": delta}
    maximum = max(abs(value) for value in deltas)
    close(maximum, finite_float(approval.get("max_abs_delta_hartree"), "approved cross-build MaxAE"), "approved cross-build MaxAE", 1.0e-12)
    checks = approval.get("checks")
    if not isinstance(checks, list) or len(checks) != 23 or not all(
        isinstance(check, dict) and check.get("passed") is True for check in checks
    ):
        raise ValueError("cross-build approval checks are incomplete or failed")
    return {
        "phase": CROSS_BUILD_PHASE,
        "decision": "approved",
        "reviewer": approval["reviewer"],
        "coverage_expected": 23,
        "coverage_completed": 23,
        "coverage_passed": 23,
        "absolute_energy_tolerance_hartree": threshold,
        "max_abs_delta_hartree": maximum,
        "approval": artifact(approval_path, root),
        "comparison_csv": artifact(comparison_path, root),
        "final_campaign_manifest": artifact(final_manifest_path, root),
        "frozen_campaign_manifest": artifact(frozen_manifest_path, root),
        "frozen_campaign_identity": frozen_campaign,
        "systems": lineage,
    }


def validate_kpoint_approval(
    root: Path,
    convergence: Mapping[str, Mapping[str, object]],
    approval_path: Path,
) -> dict[str, object]:
    approval = read_json(approval_path)
    if (
        approval.get("schema") != 1
        or approval.get("phase") != KPOINT_APPROVAL_PHASE
        or approval.get("decision") != "approved"
        or not str(approval.get("reviewer", "")).strip()
    ):
        raise ValueError("k333--k444 convergence is not explicitly approved")
    checks = approval.get("checks")
    if not isinstance(checks, list) or len(checks) != 3:
        raise ValueError("k333--k444 approval must contain three method checks")
    by_method = {str(check.get("method", "")): check for check in checks if isinstance(check, dict)}
    if set(by_method) != set(METHODS):
        raise ValueError("k333--k444 approval does not cover GFN1/GFN2/GXTB")
    validated: dict[str, object] = {}
    for method in METHODS:
        check = by_method[method]
        observed = convergence[method]
        if check.get("passed") is not True or int(check.get("N", -1)) != 23:
            raise ValueError(f"k333--k444 convergence check failed for {method}")
        mean_tolerance = finite_float(check.get("mean_abs_tolerance_kJmol"), f"{method} k-point mean tolerance")
        max_tolerance = finite_float(check.get("max_abs_tolerance_kJmol"), f"{method} k-point max tolerance")
        actual_mean = finite_float(observed.get("mean_abs_change_kJmol"), f"{method} k-point mean")
        actual_max = finite_float(observed.get("max_abs_change_kJmol"), f"{method} k-point max")
        close(actual_mean, finite_float(check.get("mean_abs_change_kJmol"), f"approved {method} k-point mean"), f"approved {method} k-point mean", 1.0e-10)
        close(actual_max, finite_float(check.get("max_abs_change_kJmol"), f"approved {method} k-point max"), f"approved {method} k-point max", 1.0e-10)
        if actual_mean > mean_tolerance or actual_max > max_tolerance:
            raise ValueError(f"k333--k444 convergence exceeds approved thresholds for {method}")
        validated[method] = {
            **observed,
            "mean_abs_tolerance_kJmol": mean_tolerance,
            "max_abs_tolerance_kJmol": max_tolerance,
            "passed": True,
        }
    return {
        "phase": KPOINT_APPROVAL_PHASE,
        "decision": "approved",
        "reviewer": approval["reviewer"],
        "approval": artifact(approval_path, root),
        "methods": validated,
    }


def summary_rows(values: Mapping[str, Mapping[str, list[float]]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for method in METHODS:
        for quantity, calculation, mesh in (
            ("lattice_energy_kJmol", "cell_opt_single_point", "k333"),
            ("volume_error_percent", "cell_opt", "k222"),
        ):
            key = "lattice" if quantity == "lattice_energy_kJmol" else "volume"
            errors = values[method][key]
            metrics = stats(errors)
            rows.append(
                {
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "quantity": quantity,
                    "calculation": calculation,
                    "mesh": mesh,
                    "N": len(errors),
                    **{name: f"{value:.12f}" for name, value in metrics.items()},
                }
            )
    return rows


def csv_text(rows: list[dict[str, object]]) -> str:
    import io

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def output_paths(root: Path, csv_path: Path | None, json_path: Path | None) -> tuple[Path, Path]:
    data = root / "data"
    return (
        csv_path.resolve() if csv_path is not None else data / f"{SUMMARY_STEM}.csv",
        json_path.resolve() if json_path is not None else data / f"{SUMMARY_STEM}.json",
    )


def finalize(args: argparse.Namespace) -> tuple[Path, Path, dict[str, object]]:
    root = args.root.resolve()
    paths = {
        "legacy_rows": args.legacy_rows or root / "data" / "x23b_final_geometry_kpoint_rows.csv",
        "legacy_volumes": args.legacy_volumes or root / "data" / "x23b_cell_volumes.csv",
        "preflight_csv": args.preflight_csv or root / "data" / "gxtb_staging" / "x23b_experimental_k222_preflight.csv",
        "fd_report": args.fd_report or root / "data" / "gxtb_staging" / "x23b_k222_fd_measured.json",
        "fd_approval": args.fd_approval or root / "data" / "gxtb_staging" / "x23b_k222_fd_approval.json",
        "cellopt_csv": args.cellopt_csv or root / "data" / "gxtb_staging" / "x23b_k222_cellopt_results.csv",
        "k333_csv": args.k333_csv or root / "data" / "gxtb_staging" / "x23b_k333_results.csv",
        "k444_csv": args.k444_csv or root / "data" / "gxtb_staging" / "x23b_k444_results.csv",
        "cross_build_approval": args.cross_build_approval or root / "data" / "gxtb_staging" / "x23b_direct_acp_cross_build_approval.json",
        "kpoint_approval": args.kpoint_approval or root / "data" / "gxtb_staging" / "x23b_k333_k444_convergence_approval.json",
    }
    paths = {name: Path(path).resolve() for name, path in paths.items()}
    out_csv, out_json = output_paths(root, args.output_csv, args.output_json)
    for final in (out_csv, out_json):
        final.unlink(missing_ok=True)
    temporary_csv = out_csv.with_name(f".{out_csv.name}.{os.getpid()}.tmp")
    temporary_json = out_json.with_name(f".{out_json.name}.{os.getpid()}.tmp")
    for temporary in (temporary_csv, temporary_json):
        temporary.unlink(missing_ok=True)
    try:
        systems, meta, metadata_artifact = metadata(root)
        provenance, campaign, manifest_path, manifest_artifact, build_lineage = validate_build_provenance(root)
        legacy_values, legacy_lineage = validate_legacy(
            root, systems, meta, paths["legacy_rows"], paths["legacy_volumes"]
        )
        workflow = provenance["workflow_paths"]
        assert isinstance(workflow, dict)
        preflight_root = Path(str(workflow["experimental_k222_preflight_root"])).resolve()
        preflight, preflight_gate = validate_preflight(
            root, systems, campaign, paths["preflight_csv"], preflight_root
        )
        fd_gate = validate_fd_gate(root, campaign, paths["fd_report"], paths["fd_approval"])
        gxtb_values, gxtb_convergence, gxtb_lineage = validate_gxtb_results(
            root,
            systems,
            meta,
            campaign,
            preflight,
            provenance,
            paths["cellopt_csv"],
            paths["k333_csv"],
            paths["k444_csv"],
        )
        cross_build = validate_cross_build_gate(
            root,
            systems,
            manifest_path,
            campaign,
            paths["cross_build_approval"],
        )
        convergence = {
            "GFN1": legacy_lineage["k333_to_k444"]["GFN1"],
            "GFN2": legacy_lineage["k333_to_k444"]["GFN2"],
            "GXTB": gxtb_convergence,
        }
        kpoint_gate = validate_kpoint_approval(root, convergence, paths["kpoint_approval"])
        values = {**legacy_values, "GXTB": gxtb_values}
        rows = summary_rows(values)
        payload: dict[str, object] = {
            "schema_version": SCHEMA_VERSION,
            "benchmark": "X23b",
            "publication_status": "publication_ready",
            "methods": list(METHODS),
            "method_labels": METHOD_LABELS,
            "coverage": {
                "required": 23,
                "common": 23,
                "systems": list(systems),
                "exact_common_coverage": True,
            },
            "protocol": {
                "geometry": "method-owned k222 KEEP_ANGLES CELL_OPT",
                "reported_lattice_energy": "k333 single point on the method-owned final k222 geometry",
                "volume": "final k222 CELL_OPT volume",
                "convergence_check": "k444 minus k333 on the identical final k222 geometry",
                "gxtb_source_policy": "experimental_reference accepted only through exact 23/23 preflight-to-cellopt hash lineage",
            },
            "summary": rows,
            "gates": {
                "direct_acp_cross_build": cross_build,
                "finite_difference_pilot": fd_gate,
                "experimental_k222_preflight": preflight_gate,
                "k333_to_k444_convergence": kpoint_gate,
            },
            "provenance": {
                "metadata": metadata_artifact,
                "builds": build_lineage,
                "production_campaign_manifest": manifest_artifact,
                "GFN1_GFN2": legacy_lineage,
                "GXTB": gxtb_lineage,
            },
        }
        csv_body = csv_text(rows)
        payload["publication_csv_sha256"] = hashlib.sha256(csv_body.encode()).hexdigest()
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        temporary_csv.write_text(csv_body)
        temporary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary_csv, out_csv)
        os.replace(temporary_json, out_json)
        return out_csv, out_json, payload
    except BaseException:
        for path in (temporary_csv, temporary_json, out_csv, out_json):
            path.unlink(missing_ok=True)
        raise


def parser() -> argparse.ArgumentParser:
    default_root = Path(__file__).resolve().parents[1]
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--root", type=Path, default=default_root)
    result.add_argument("--legacy-rows", type=Path)
    result.add_argument("--legacy-volumes", type=Path)
    result.add_argument("--preflight-csv", type=Path)
    result.add_argument("--fd-report", type=Path)
    result.add_argument("--fd-approval", type=Path)
    result.add_argument("--cellopt-csv", type=Path)
    result.add_argument("--k333-csv", type=Path)
    result.add_argument("--k444-csv", type=Path)
    result.add_argument("--cross-build-approval", type=Path)
    result.add_argument("--kpoint-approval", type=Path)
    result.add_argument("--output-csv", type=Path)
    result.add_argument("--output-json", type=Path)
    return result


def main() -> None:
    args = parser().parse_args()
    csv_path, json_path, payload = finalize(args)
    print(
        f"Frozen X23b paper artifact: {csv_path} and {json_path}; "
        f"common coverage={payload['coverage']['common']}/23"
    )


if __name__ == "__main__":
    main()
