#!/usr/bin/env python3
"""Freeze the validated LC12 GFN1/GFN2/g-XTB comparison for publication.

The normal collectors deliberately keep method-selective working tables.  This
finalizer is the stricter publication boundary: it recomputes every reported
quantity from the archived raw outputs, requires the approved g-XTB EOS
fingerprint and all k333/k444/k555 final single points, and writes one compact
CSV plus a complete JSON lineage manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import run_goldzak12_benchmark as base
import run_goldzak12_eos_benchmark as eos


METHODS = ("GFN1", "GFN2", "GXTB")
METHOD_LABELS = {
    "GFN1": "GFN1-xTB",
    "GFN2": "GFN2-xTB",
    "GXTB": "g-xTB",
}
EOS_MESH = "k444"
ENERGY_MESHES = ("k333", "k444", "k555")
RESULT_MESH = "k555"
SUMMARY_STEM = "lc12_gfn_gxtb_paper_summary"
SCHEMA_VERSION = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def optional_artifact(path: Path, root: Path) -> dict[str, object] | None:
    return artifact(path, root) if path.is_file() else None


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"required CSV is missing or empty: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"required CSV has no records: {path}")
    return rows


def read_optional_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return value


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


def close(actual: float, expected: float, label: str, tolerance: float) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(
            f"{label} mismatch: found {actual:.15g}, expected {expected:.15g}"
        )


def unique_by(
    rows: Iterable[dict[str, str]],
    fields: tuple[str, ...],
    label: str,
) -> dict[tuple[str, ...], dict[str, str]]:
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in fields)
        if key in result:
            raise ValueError(f"duplicate {label} record: {key}")
        result[key] = row
    return result


def stats(errors: list[float]) -> dict[str, float]:
    if not errors:
        raise ValueError("cannot summarize an empty LC12 error set")
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(value) for value in errors) / len(errors),
        "RMSE": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "MaxAE": max(abs(value) for value in errors),
    }


def fit_is_valid(row: Mapping[str, object], method: str) -> bool:
    if not str(row.get("a_eos_A", "")).strip():
        return False
    return method != "GXTB" or row.get("fit_status") == "quadratic"


def portable_gxtb_stamp(
    result: Path,
    local_input: Path,
    campaign: dict[str, object],
    executable_role: str,
    root: Path,
) -> dict[str, object]:
    """Validate a stamp without depending on its host-specific absolute path."""
    base.validate_campaign_identity(campaign)
    stamp_path = base.job_stamp_path(result)
    stamp = read_json(stamp_path)
    if stamp.get("schema_version") != 1 or stamp.get("completed") is not True:
        raise ValueError(f"incomplete campaign stamp for {result}")
    if int(stamp.get("return_code", -1)) != 0:
        raise ValueError(f"nonzero campaign return code for {result}")
    if stamp.get("campaign_identity") != campaign:
        raise ValueError(f"campaign identity mismatch for {result}")
    if not local_input.is_file() or stamp.get("input_sha256") != sha256(local_input):
        raise ValueError(f"input fingerprint mismatch for {result}")
    expected_field = {
        "cp2k": "cp2k_executable_sha256",
        "save_tblite": "save_tblite_executable_sha256",
    }.get(executable_role)
    if expected_field is None:
        raise ValueError(f"unknown executable role {executable_role}")
    if stamp.get("executable_sha256") != campaign.get(expected_field):
        raise ValueError(f"{executable_role} executable fingerprint mismatch for {result}")
    return artifact(stamp_path, root)


def validate_manifest(
    root: Path,
    provenance: dict[str, Any],
    campaign: dict[str, object],
) -> tuple[Path, dict[str, object]]:
    record = provenance.get("campaign_manifest")
    if not isinstance(record, dict):
        raise ValueError("g-XTB provenance lacks its campaign manifest record")
    campaign_id = str(campaign.get("campaign_id", ""))
    expected_hash = str(record.get("file_sha256", ""))
    if not expected_hash:
        raise ValueError("g-XTB provenance lacks its campaign-manifest hash")
    candidates = [
        Path(str(record.get("path", ""))),
        root.parent / "campaigns" / campaign_id / "build_manifest.json",
    ]
    campaign_root = root.parent / "campaigns"
    if campaign_root.is_dir():
        candidates.extend(sorted(campaign_root.rglob("build_manifest.json")))
    manifest_path = next(
        (
            path
            for path in dict.fromkeys(candidates)
            if path.is_file() and sha256(path) == expected_hash
        ),
        None,
    )
    if manifest_path is None:
        raise ValueError(
            "the exact frozen g-XTB campaign manifest is not locally available"
        )
    manifest = read_json(manifest_path)
    if base.campaign_identity_from_manifest(manifest, manifest_path) != campaign:
        raise ValueError("g-XTB campaign manifest identity mismatch")
    return manifest_path, artifact(manifest_path, root)


def validate_build_provenance(
    root: Path,
    fits: list[dict[str, str]],
) -> tuple[dict[str, Any], dict[str, object], dict[str, object]]:
    data = root / "data"
    legacy_path = data / "build_provenance.json"
    gxtb_path = data / "build_provenance_gxtb.json"
    legacy = read_json(legacy_path)
    gxtb = read_json(gxtb_path)
    legacy_protocol = legacy.get("protocol")
    protocol = gxtb.get("protocol")
    if not isinstance(legacy_protocol, dict) or not isinstance(protocol, dict):
        raise ValueError("LC12 build provenance lacks its protocol")
    if legacy_protocol.get("result_mesh") != RESULT_MESH:
        raise ValueError("GFN1/GFN2 provenance has the wrong result mesh")
    if protocol.get("eos_mesh") != EOS_MESH or protocol.get("result_mesh") != RESULT_MESH:
        raise ValueError("g-XTB provenance has the wrong EOS/result mesh")
    if set(protocol.get("energy_meshes", [])) != set(ENERGY_MESHES):
        raise ValueError("g-XTB provenance does not contain k333/k444/k555")
    if protocol.get("fit_approval_required") is not True or protocol.get("fit_approved") is not True:
        raise ValueError("g-XTB EOS fits are not explicitly approved")
    gxtb_fits = [row for row in fits if row.get("method") == "GXTB"]
    current_fit_hash = eos.gxtb_fit_approval_sha256(gxtb_fits)
    if protocol.get("approved_gxtb_fit_sha256") != current_fit_hash:
        raise ValueError("approved g-XTB EOS fingerprint differs from eos_fits.csv")
    if protocol.get("current_gxtb_fit_sha256") not in (None, current_fit_hash):
        raise ValueError("current g-XTB EOS fingerprint differs from eos_fits.csv")
    campaign = gxtb.get("campaign_identity")
    if not isinstance(campaign, dict):
        raise ValueError("g-XTB provenance lacks the complete campaign identity")
    base.validate_campaign_identity(campaign)
    manifest_path, manifest_record = validate_manifest(root, gxtb, campaign)
    scale_path = data / "gxtb_eos_scale_manifest.json"
    if not scale_path.is_file():
        raise ValueError("g-XTB scale manifest is missing")
    if protocol.get("gxtb_scale_manifest_sha256") != sha256(scale_path):
        raise ValueError("g-XTB scale-manifest hash differs from provenance")
    classification_path = data / "gxtb_eos_classifications.json"
    recorded_classification_hash = protocol.get(
        "gxtb_classification_manifest_sha256"
    )
    if classification_path.is_file():
        if recorded_classification_hash != sha256(classification_path):
            raise ValueError(
                "g-XTB classification-manifest hash differs from provenance"
            )
    elif recorded_classification_hash is not None:
        raise ValueError("the approved g-XTB classification manifest is missing")
    return (
        gxtb,
        campaign,
        {
            "GFN1_GFN2": {
                "artifact": artifact(legacy_path, root),
                "cp2k": legacy.get("cp2k"),
                "provider_name": "tblite",
                "provider": legacy.get("tblite"),
                "repository_patches": legacy.get("repository_patches"),
            },
            "GXTB": {
                "artifact": artifact(gxtb_path, root),
                "campaign_manifest": manifest_record,
                "campaign_manifest_path_used": relative_path(manifest_path, root),
                "campaign_identity": campaign,
                "cp2k": gxtb.get("cp2k"),
                "provider_name": "save_tblite",
                "provider": gxtb.get("save_tblite"),
            },
        },
    )


def validate_scale_manifest(
    root: Path,
    point_rows: list[dict[str, str]],
) -> None:
    manifest = read_json(root / "data" / "gxtb_eos_scale_manifest.json")
    if manifest.get("eos_mesh") != EOS_MESH:
        raise ValueError("g-XTB scale manifest has the wrong EOS mesh")
    expected: set[tuple[str, str]] = set()
    systems = manifest.get("systems")
    if not isinstance(systems, list):
        raise ValueError("g-XTB scale manifest has no system list")
    for item in systems:
        if not isinstance(item, dict) or item.get("method") != "GXTB":
            raise ValueError("invalid g-XTB scale-manifest record")
        solid = str(item.get("solid", ""))
        for scale in item.get("requested_scales", []):
            expected.add((solid, f"{float(scale):.5f}"))
    actual = {
        (row.get("solid", ""), f"{float(row.get('scale', 'nan')):.5f}")
        for row in point_rows
        if row.get("method") == "GXTB" and row.get("mesh") == EOS_MESH
    }
    if expected != actual:
        raise ValueError(
            "g-XTB EOS point coverage differs from its scale manifest: "
            f"missing {sorted(expected - actual)}, unexpected {sorted(actual - expected)}"
        )


def validate_atom_references(
    root: Path,
    campaign: dict[str, object],
) -> tuple[dict[tuple[str, str], float], dict[str, object], dict[str, object]]:
    data = root / "data"
    legacy_path = data / "atom_energies_tblite_cli.csv"
    gxtb_path = data / "atom_energies_save_tblite_cli_gxtb.csv"
    check_path = data / "atom_reference_cp2k_vs_save_tblite_gxtb.csv"
    rows = read_csv(legacy_path) + read_csv(gxtb_path)
    by_key = unique_by(rows, ("method", "element"), "atom reference")
    elements = tuple(sorted(base.ELEMENT_MULTIPLICITY))
    expected = {(method, element) for method in METHODS for element in elements}
    if set(by_key) != expected:
        raise ValueError(
            "atom-reference coverage differs: "
            f"missing {sorted(expected - set(by_key))}, unexpected {sorted(set(by_key) - expected)}"
        )
    energies: dict[tuple[str, str], float] = {}
    lineage: dict[str, object] = {method: {} for method in METHODS}
    for method, element in sorted(expected):
        row = by_key[(method, element)]
        expected_source = "save_tblite_cli" if method == "GXTB" else "tblite_cli"
        if row.get("source") != expected_source:
            raise ValueError(f"{method}/{element} atom-reference source mismatch")
        json_path = root / "runs" / "atoms_cli" / method / element / f"atom_{element}_{method}.json"
        out_path = json_path.with_suffix(".out")
        xyz_path = json_path.parent / f"atom_{element}.xyz"
        raw_energy = base.parse_tblite_json_energy(json_path)
        if raw_energy is None:
            raise ValueError(f"cannot parse {method}/{element} atom energy")
        table_energy = finite_float(row.get("energy_hartree"), f"{method}/{element} atom energy")
        close(raw_energy, table_energy, f"{method}/{element} atom energy", 1.0e-10)
        energies[(method, element)] = raw_energy
        record: dict[str, object] = {
            "energy_hartree": raw_energy,
            "json": artifact(json_path, root),
            "stdout": artifact(out_path, root),
            "geometry": artifact(xyz_path, root),
        }
        if method == "GXTB":
            record["campaign_stamp"] = portable_gxtb_stamp(
                json_path, xyz_path, campaign, "save_tblite", root
            )
        cast_lineage = lineage[method]
        assert isinstance(cast_lineage, dict)
        cast_lineage[element] = record

    checks = read_csv(check_path)
    check_by_element = unique_by(checks, ("element",), "g-XTB atom check")
    if set(key[0] for key in check_by_element) != set(elements):
        raise ValueError("g-XTB CP2K/save_tblite atom check is incomplete")
    check_lineage: dict[str, object] = {}
    for element in elements:
        row = check_by_element[(element,)]
        if row.get("method") != "GXTB" or not truth(row.get("passed")):
            raise ValueError(f"g-XTB atom check failed for {element}")
        if row.get("campaign_stamp_issue", ""):
            raise ValueError(f"g-XTB atom check has a campaign-stamp issue for {element}")
        cli_energy = finite_float(row.get("cli_energy_hartree"), f"GXTB/{element} checked CLI energy")
        close(cli_energy, energies[("GXTB", element)], f"GXTB/{element} checked CLI energy", 1.0e-10)
        cp2k_output = root / "runs" / "atoms" / "GXTB" / element / f"atom_{element}_GXTB.out"
        cp2k_input = cp2k_output.with_suffix(".inp")
        cp2k_energy = base.parse_energy(cp2k_output)
        if cp2k_energy is None or not base.output_ok(cp2k_output):
            raise ValueError(f"g-XTB CP2K atom check output is incomplete for {element}")
        close(
            cp2k_energy,
            finite_float(row.get("cp2k_energy_hartree"), f"GXTB/{element} checked CP2K energy"),
            f"GXTB/{element} checked CP2K energy",
            1.0e-10,
        )
        delta = cp2k_energy - cli_energy
        close(
            delta,
            finite_float(row.get("delta_cp2k_minus_cli_hartree"), f"GXTB/{element} atom delta"),
            f"GXTB/{element} atom delta",
            1.0e-10,
        )
        tolerance = finite_float(row.get("tolerance_hartree"), f"GXTB/{element} atom tolerance")
        if abs(delta) > tolerance:
            raise ValueError(f"g-XTB atom check exceeds tolerance for {element}")
        check_lineage[element] = {
            "cp2k_energy_hartree": cp2k_energy,
            "delta_cp2k_minus_cli_hartree": delta,
            "tolerance_hartree": tolerance,
            "input": artifact(cp2k_input, root),
            "output": artifact(cp2k_output, root),
            "campaign_stamp": portable_gxtb_stamp(
                cp2k_output, cp2k_input, campaign, "cp2k", root
            ),
        }
    return energies, lineage, {
        "table": artifact(check_path, root),
        "elements": check_lineage,
    }


def validate_eos_and_collect_lineage(
    root: Path,
    fits: list[dict[str, str]],
    points: list[dict[str, str]],
    campaign: dict[str, object],
    protocol: dict[str, Any],
) -> tuple[
    dict[tuple[str, str], dict[str, str]],
    dict[str, dict[str, object]],
]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    fit_by_key = unique_by(fits, ("solid", "method"), "EOS fit")
    expected_fit_keys = {(solid, method) for solid in refs for method in METHODS}
    if set(fit_by_key) != expected_fit_keys:
        raise ValueError(
            "EOS-fit coverage differs: "
            f"missing {sorted(expected_fit_keys - set(fit_by_key))}, "
            f"unexpected {sorted(set(fit_by_key) - expected_fit_keys)}"
        )
    for row in fits:
        if row.get("eos_mesh") != EOS_MESH:
            raise ValueError(f"{row.get('method')}/{row.get('solid')} has the wrong EOS mesh")
    point_by_key = unique_by(points, ("solid", "method", "mesh", "scale"), "EOS point")
    del point_by_key  # duplicate detection is the purpose of this index
    validate_scale_manifest(root, points)
    branch_rows = read_optional_csv(root / "data" / "gxtb_eos_branch_diagnostics.csv")
    unresolved_branches = [
        row for row in branch_rows if row.get("resolution") == "unresolved_candidate"
    ]
    if unresolved_branches:
        raise ValueError("g-XTB branch diagnostics still contain unresolved candidates")

    points_by_fit: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in points:
        key = (row.get("solid", ""), row.get("method", ""))
        if key not in expected_fit_keys or row.get("mesh") != EOS_MESH:
            raise ValueError(f"unexpected EOS point identity: {key}/{row.get('mesh')}")
        points_by_fit.setdefault(key, []).append(row)

    valid_gxtb = 0
    invalid_gxtb: list[str] = []
    failed_gxtb_points: list[str] = []
    followup_by_solid = {
        row.get("solid", ""): row
        for row in read_optional_csv(root / "data" / "gxtb_adaptive_followup.csv")
    }
    lineage: dict[str, dict[str, object]] = {method: {} for method in METHODS}
    for key in sorted(expected_fit_keys, key=lambda item: (METHODS.index(item[1]), item[0])):
        solid, method = key
        fit = fit_by_key[key]
        fit_points = points_by_fit.get(key, [])
        if int(fit.get("n_requested", -1)) != len(fit_points):
            raise ValueError(f"{method}/{solid} n_requested differs from eos_points.csv")
        raw_converged = sum(
            truth(row.get("completed")) and bool(row.get("energy_hartree", ""))
            for row in fit_points
        )
        if int(fit.get("n_converged_raw", -1)) != raw_converged:
            raise ValueError(f"{method}/{solid} n_converged_raw mismatch")
        if method == "GXTB" and int(
            fit.get("n_unresolved_branch_candidates", -1)
        ) != 0:
            raise ValueError(f"g-XTB/{solid} fit retains unresolved branch candidates")
        accepted = [row for row in fit_points if truth(row.get("valid_for_eos"))]
        if int(fit.get("n_completed", -1)) != len(accepted):
            raise ValueError(f"{method}/{solid} accepted-point count mismatch")
        numeric_points: list[tuple[float, float, float, bool]] = []
        point_lineage: list[dict[str, object]] = []
        for row in sorted(fit_points, key=lambda item: float(item["scale"])):
            scale = finite_float(row.get("scale"), f"{method}/{solid} EOS scale")
            completed = truth(row.get("completed"))
            accepted_point = truth(row.get("valid_for_eos"))
            project = eos.eos_project(solid, method, EOS_MESH, scale)
            run_dir = (
                root
                / "runs"
                / "eos"
                / method
                / solid
                / EOS_MESH
                / eos.scale_tag(scale, method)
            )
            input_path = run_dir / f"{project}.inp"
            output_path = run_dir / f"{project}.out"
            record: dict[str, object] = {
                "scale": scale,
                "a_A": finite_float(row.get("a_A"), f"{method}/{solid}@{scale} lattice value"),
                "completed": completed,
                "valid_for_eos": accepted_point,
                "diagnostic": row.get("diagnostic", ""),
                "classification_resolution": row.get("classification_resolution", ""),
                "classification_rationale": row.get("classification_rationale", ""),
            }
            if completed:
                if not base.output_ok(output_path):
                    raise ValueError(f"completed EOS point has an invalid output: {output_path}")
                raw_energy = base.parse_energy(output_path)
                if raw_energy is None:
                    raise ValueError(f"completed EOS point has no energy: {output_path}")
                table_energy = finite_float(
                    row.get("energy_hartree"), f"{method}/{solid}@{scale} EOS energy"
                )
                close(raw_energy, table_energy, f"{method}/{solid}@{scale} EOS energy", 1.0e-10)
                record.update(
                    {
                        "energy_hartree": raw_energy,
                        "input": artifact(input_path, root),
                        "output": artifact(output_path, root),
                    }
                )
                if method == "GXTB":
                    if row.get("scf_strategy") not in {
                        "native_gxtb_fdiis",
                        "native_gxtb_fdiis_adaptive",
                    }:
                        raise ValueError(
                            f"g-XTB EOS point used a non-production mixer: {solid}@{scale}"
                        )
                    record["campaign_stamp"] = portable_gxtb_stamp(
                        output_path, input_path, campaign, "cp2k", root
                    )
                    if row.get("classification_resolution") == "unresolved_candidate":
                        raise ValueError(
                            f"g-XTB EOS point has an unresolved branch candidate: {solid}@{scale}"
                        )
                    if not accepted_point and (
                        row.get("classification_resolution") != "explicit_exclusion"
                        or not row.get("classification_rationale")
                    ):
                        raise ValueError(
                            f"g-XTB EOS point was excluded without review: {solid}@{scale}"
                        )
                if accepted_point:
                    numeric_points.append((float(row["a_A"]), scale, raw_energy, True))
            elif accepted_point:
                raise ValueError(f"incomplete EOS point is marked valid: {method}/{solid}@{scale}")
            elif method == "GXTB":
                failed_gxtb_points.append(f"{solid}@{scale:.5f}")
                if (
                    row.get("classification_resolution") != "explicit_failure_classification"
                    or not row.get("classification_rationale")
                ):
                    raise ValueError(
                        f"g-XTB failed EOS point is not explicitly classified: {solid}@{scale}"
                    )
                record["input"] = artifact(input_path, root)
                record["output"] = optional_artifact(output_path, root)
                record["campaign_stamp"] = optional_artifact(
                    base.job_stamp_path(output_path), root
                )
            point_lineage.append(record)

        recomputed = eos.fit_gxtb_eos(numeric_points) if method == "GXTB" else eos.fit_eos(numeric_points)
        if recomputed.get("fit_status") != fit.get("fit_status"):
            raise ValueError(f"{method}/{solid} EOS fit status is not reproducible")
        for field, tolerance in (
            ("a_eos_A", 5.0e-8),
            ("energy_fit_hartree", 5.0e-10),
            ("fit_rmse_hartree", 5.0e-10),
            ("grid_min_a_A", 5.0e-8),
            ("grid_min_scale", 5.0e-8),
            ("grid_min_energy_hartree", 5.0e-10),
        ):
            stored = str(fit.get(field, "")).strip()
            regenerated = str(recomputed.get(field, "")).strip()
            if bool(stored) != bool(regenerated):
                raise ValueError(f"{method}/{solid} EOS {field} presence mismatch")
            if stored:
                close(float(stored), float(regenerated), f"{method}/{solid} EOS {field}", tolerance)
        if fit_is_valid(fit, method):
            if method == "GXTB":
                valid_gxtb += 1
        elif method == "GXTB":
            invalid_gxtb.append(solid)
        lineage[method][solid] = {
            "fit": dict(fit),
            "eos_points": point_lineage,
        }

    allow_reduced = bool(protocol.get("allow_reduced_coverage"))
    minimum = int(protocol.get("minimum_valid_gxtb_fits", -1))
    if failed_gxtb_points and not allow_reduced:
        raise ValueError(
            "g-XTB has explicitly classified failed EOS points without reduced-coverage approval: "
            + ", ".join(failed_gxtb_points)
        )
    if invalid_gxtb and not allow_reduced:
        raise ValueError(
            "g-XTB EOS coverage is reduced without explicit approval: "
            + ", ".join(invalid_gxtb)
        )
    if allow_reduced and valid_gxtb < minimum:
        raise ValueError(
            f"g-XTB has {valid_gxtb} valid fits, below the approved minimum {minimum}"
        )
    for solid in invalid_gxtb:
        row = followup_by_solid.get(solid)
        if (
            row is None
            or not truth(row.get("adaptive_investigated"))
            or not row.get("classification")
            or not row.get("interpretation")
        ):
            raise ValueError(f"g-XTB invalid EOS fit lacks adaptive investigation: {solid}")
        system = lineage["GXTB"][solid]
        system["reporting_status"] = "excluded_no_valid_eos_minimum"
        system["adaptive_followup"] = dict(row)
        discarded: dict[str, object] = {}
        fit = fit_by_key[(solid, "GXTB")]
        for mesh in ENERGY_MESHES:
            input_path = (
                root
                / "runs"
                / "eos_final_sp"
                / "GXTB"
                / solid
                / mesh
                / f"{eos.final_project(solid, 'GXTB', mesh)}.inp"
            )
            if not input_path.is_file():
                continue
            lineage_path = eos.final_input_lineage_path(input_path)
            invalid_lineage = read_json(lineage_path)
            if (
                invalid_lineage.get("schema_version") != eos.FINAL_INPUT_LINEAGE_SCHEMA
                or invalid_lineage.get("valid") is not False
                or invalid_lineage.get("solid") != solid
                or invalid_lineage.get("method") != "GXTB"
                or invalid_lineage.get("energy_mesh") != mesh
                or invalid_lineage.get("fit_status") != fit.get("fit_status")
                or invalid_lineage.get("input_sha256") != sha256(input_path)
            ):
                raise ValueError(
                    f"invalid-fit g-XTB final input is not explicitly invalidated: {solid}/{mesh}"
                )
            output_path = input_path.with_suffix(".out")
            discarded[mesh] = {
                "status": "not_reported_invalid_eos_fit",
                "input": artifact(input_path, root),
                "invalidation_lineage": artifact(lineage_path, root),
                "preexisting_output": optional_artifact(output_path, root),
                "preexisting_campaign_stamp": optional_artifact(
                    base.job_stamp_path(output_path), root
                ),
            }
        if discarded:
            system["discarded_final_artifacts"] = discarded
    for solid in lineage["GXTB"]:
        lineage["GXTB"][solid].setdefault(
            "reporting_status", "reported_at_approved_eos_minimum"
        )
    return fit_by_key, lineage


def validate_final_results(
    root: Path,
    results: list[dict[str, str]],
    fit_by_key: dict[tuple[str, str], dict[str, str]],
    eos_lineage: dict[str, dict[str, object]],
    atom_energies: dict[tuple[str, str], float],
    campaign: dict[str, object],
) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    result_by_key = unique_by(results, ("solid", "method", "energy_mesh"), "final result")
    valid_fit_keys = {
        key for key, fit in fit_by_key.items() if fit_is_valid(fit, key[1])
    }
    expected_results = {
        (solid, method, mesh)
        for solid, method in valid_fit_keys
        for mesh in ENERGY_MESHES
    }
    if set(result_by_key) != expected_results:
        raise ValueError(
            "final single-point coverage differs: "
            f"missing {sorted(expected_results - set(result_by_key))}, "
            f"unexpected {sorted(set(result_by_key) - expected_results)}"
        )

    selected: dict[str, list[dict[str, object]]] = {method: [] for method in METHODS}
    for solid, method in sorted(valid_fit_keys, key=lambda item: (METHODS.index(item[1]), item[0])):
        ref = refs[solid]
        fit = fit_by_key[(solid, method)]
        a_calc = finite_float(fit.get("a_eos_A"), f"{method}/{solid} fitted lattice constant")
        atom_sum = sum(
            atom_energies[(method, element)] * count
            for element, count in base.atom_counts(ref).items()
        )
        n_atoms = len(base.conventional_cell_atoms(ref))
        final_lineage: dict[str, object] = {}
        selected_row: dict[str, object] | None = None
        for mesh in ENERGY_MESHES:
            row = result_by_key[(solid, method, mesh)]
            if not truth(row.get("sp_completed")):
                raise ValueError(f"incomplete final single point: {method}/{solid}/{mesh}")
            if row.get("eos_mesh") != EOS_MESH or row.get("fit_status") != fit.get("fit_status"):
                raise ValueError(f"final result fit identity mismatch: {method}/{solid}/{mesh}")
            project = eos.final_project(solid, method, mesh)
            run_dir = root / "runs" / "eos_final_sp" / method / solid / mesh
            input_path = run_dir / f"{project}.inp"
            output_path = run_dir / f"{project}.out"
            if not base.output_ok(output_path):
                raise ValueError(f"invalid final single-point output: {output_path}")
            raw_energy = base.parse_energy(output_path)
            if raw_energy is None:
                raise ValueError(f"missing final single-point energy: {output_path}")
            close(
                raw_energy,
                finite_float(row.get("solid_energy_hartree"), f"{method}/{solid}/{mesh} solid energy"),
                f"{method}/{solid}/{mesh} solid energy",
                1.0e-10,
            )
            close(
                a_calc,
                finite_float(row.get("a_calc_A"), f"{method}/{solid}/{mesh} lattice constant"),
                f"{method}/{solid}/{mesh} lattice constant",
                5.0e-8,
            )
            close(
                float(ref.a_exp),
                finite_float(row.get("a_ref_exp_A"), f"{solid} lattice reference"),
                f"{solid} lattice reference",
                1.0e-10,
            )
            a_error = a_calc - ref.a_exp
            close(
                a_error,
                finite_float(row.get("a_error_A"), f"{method}/{solid}/{mesh} lattice error"),
                f"{method}/{solid}/{mesh} lattice error",
                5.0e-8,
            )
            ecoh = (atom_sum - raw_energy) * base.HARTREE_TO_EV / n_atoms
            close(
                ecoh,
                finite_float(row.get("ecoh_calc_eV_per_atom"), f"{method}/{solid}/{mesh} cohesive energy"),
                f"{method}/{solid}/{mesh} cohesive energy",
                5.0e-8,
            )
            close(
                float(ref.ecoh_exp),
                finite_float(row.get("ecoh_ref_exp_eV_per_atom"), f"{solid} cohesive reference"),
                f"{solid} cohesive reference",
                1.0e-10,
            )
            ecoh_error = ecoh - ref.ecoh_exp
            close(
                ecoh_error,
                finite_float(row.get("ecoh_error_eV_per_atom"), f"{method}/{solid}/{mesh} cohesive error"),
                f"{method}/{solid}/{mesh} cohesive error",
                5.0e-8,
            )
            expected_atom_source = "save_tblite_cli" if method == "GXTB" else "tblite_cli"
            if row.get("atom_reference_source") != expected_atom_source:
                raise ValueError(f"{method}/{solid}/{mesh} atom-reference source mismatch")
            record: dict[str, object] = {
                "solid_energy_hartree": raw_energy,
                "cohesive_energy_eV_per_atom": ecoh,
                "input": artifact(input_path, root),
                "output": artifact(output_path, root),
            }
            if method == "GXTB":
                lineage_path = eos.final_input_lineage_path(input_path)
                input_lineage = read_json(lineage_path)
                expected_lineage = {
                    "schema_version": eos.FINAL_INPUT_LINEAGE_SCHEMA,
                    "valid": True,
                    "solid": solid,
                    "method": "GXTB",
                    "eos_mesh": EOS_MESH,
                    "energy_mesh": mesh,
                    "fit_status": "quadratic",
                    "a_eos_A": str(fit.get("a_eos_A", "")),
                    "input_sha256": sha256(input_path),
                    "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
                }
                for field, expected in expected_lineage.items():
                    if input_lineage.get(field) != expected:
                        raise ValueError(
                            f"g-XTB final-input EOS lineage {field} mismatch: {solid}/{mesh}"
                        )
                record["eos_input_lineage"] = artifact(lineage_path, root)
                record["campaign_stamp"] = portable_gxtb_stamp(
                    output_path, input_path, campaign, "cp2k", root
                )
            final_lineage[mesh] = record
            if mesh == RESULT_MESH:
                selected_row = {
                    "solid": solid,
                    "structure": ref.structure,
                    "method": method,
                    "a_calc_A": a_calc,
                    "a_ref_A": ref.a_exp,
                    "a_error_A": a_error,
                    "ecoh_calc_eV_per_atom": ecoh,
                    "ecoh_ref_eV_per_atom": ref.ecoh_exp,
                    "ecoh_error_eV_per_atom": ecoh_error,
                }
        assert selected_row is not None
        selected[method].append(selected_row)
        system_lineage = eos_lineage[method][solid]
        assert isinstance(system_lineage, dict)
        system_lineage["final_single_points"] = final_lineage
        system_lineage["reported_result"] = selected_row
    return selected, eos_lineage


def scope_record(
    method: str,
    scope: str,
    rows: list[dict[str, object]],
    ordered_systems: tuple[str, ...],
) -> dict[str, object]:
    by_solid = {str(row["solid"]): row for row in rows}
    chosen = [by_solid[solid] for solid in ordered_systems]
    a_stats = stats([float(row["a_error_A"]) for row in chosen])
    e_stats = stats([float(row["ecoh_error_eV_per_atom"]) for row in chosen])
    return {
        "method_id": method,
        "method_label": METHOD_LABELS[method],
        "scope": scope,
        "n_systems": len(chosen),
        "coverage_denominator": len(base.REFERENCES),
        "systems": ";".join(ordered_systems),
        "eos_mesh": EOS_MESH,
        "result_mesh": RESULT_MESH,
        "lattice_ME_A": a_stats["ME"],
        "lattice_MAE_A": a_stats["MAE"],
        "lattice_RMSE_A": a_stats["RMSE"],
        "lattice_MaxAE_A": a_stats["MaxAE"],
        "cohesive_ME_eV_per_atom": e_stats["ME"],
        "cohesive_MAE_eV_per_atom": e_stats["MAE"],
        "cohesive_RMSE_eV_per_atom": e_stats["RMSE"],
        "cohesive_MaxAE_eV_per_atom": e_stats["MaxAE"],
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.9f}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )


def build_summary(root: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    data = root / "data"
    fit_path = data / "eos_fits.csv"
    point_path = data / "eos_points.csv"
    result_path = data / "eos_results.csv"
    fits = read_csv(fit_path)
    points = read_csv(point_path)
    results = read_csv(result_path)
    gxtb_provenance, campaign, build_provenance = validate_build_provenance(root, fits)
    protocol = gxtb_provenance["protocol"]
    assert isinstance(protocol, dict)
    atom_energies, atom_lineage, atom_check = validate_atom_references(root, campaign)
    fit_by_key, eos_lineage = validate_eos_and_collect_lineage(
        root, fits, points, campaign, protocol
    )
    selected, system_lineage = validate_final_results(
        root,
        results,
        fit_by_key,
        eos_lineage,
        atom_energies,
        campaign,
    )

    reference_order = tuple(ref.solid for ref in base.REFERENCES)
    available_systems = {
        method: tuple(solid for solid in reference_order if any(row["solid"] == solid for row in selected[method]))
        for method in METHODS
    }
    common_systems = tuple(
        solid
        for solid in reference_order
        if all(solid in available_systems[method] for method in METHODS)
    )
    if not common_systems:
        raise ValueError("the three-method LC12 common subset is empty")
    reduced_gxtb_coverage = (
        len(available_systems["GXTB"]) < len(base.REFERENCES)
    )
    summary_rows: list[dict[str, object]] = []
    method_payload: dict[str, object] = {}
    for method in METHODS:
        available = scope_record(
            method, "method_available_coverage", selected[method], available_systems[method]
        )
        common = scope_record(
            method, "three_method_common_subset", selected[method], common_systems
        )
        summary_rows.extend((available, common))
        method_payload[method] = {
            "method_label": METHOD_LABELS[method],
            "available_coverage": available,
            "three_method_common_subset": common,
            "atom_references": atom_lineage[method],
            "systems": system_lineage[method],
        }

    source_paths = {
        "reference": data / "reference_goldzak2022.csv",
        "eos_fits": fit_path,
        "eos_points": point_path,
        "eos_results": result_path,
        "legacy_atom_references": data / "atom_energies_tblite_cli.csv",
        "gxtb_atom_references": data / "atom_energies_save_tblite_cli_gxtb.csv",
        "gxtb_atom_check": data / "atom_reference_cp2k_vs_save_tblite_gxtb.csv",
        "gxtb_scale_manifest": data / "gxtb_eos_scale_manifest.json",
        "legacy_build_provenance": data / "build_provenance.json",
        "gxtb_build_provenance": data / "build_provenance_gxtb.json",
    }
    sources = {name: artifact(path, root) for name, path in source_paths.items()}
    for name, path in (
        ("gxtb_branch_diagnostics", data / "gxtb_eos_branch_diagnostics.csv"),
        ("gxtb_adaptive_followup", data / "gxtb_adaptive_followup.csv"),
        ("gxtb_classifications", data / "gxtb_eos_classifications.json"),
    ):
        item = optional_artifact(path, root)
        if item is not None:
            sources[name] = item
    payload: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "benchmark": "LC12 (Goldzak12)",
        "status": (
            "publication_ready_reduced_coverage"
            if reduced_gxtb_coverage
            else "publication_ready"
        ),
        "methods": method_payload,
        "protocol": {
            "eos_mesh": EOS_MESH,
            "energy_meshes": list(ENERGY_MESHES),
            "result_mesh": RESULT_MESH,
            "lattice_unit": "angstrom",
            "cohesive_energy_unit": "eV/atom",
            "common_subset_systems": list(common_systems),
            "common_subset_count": len(common_systems),
            "gxtb_fit_approval_sha256": protocol["approved_gxtb_fit_sha256"],
            "gxtb_allow_reduced_coverage": bool(protocol.get("allow_reduced_coverage")),
            "gxtb_reduced_coverage_reported": reduced_gxtb_coverage,
            "gxtb_minimum_valid_fits": int(protocol.get("minimum_valid_gxtb_fits", -1)),
        },
        "summary_rows": summary_rows,
        "atom_reference_acceptance": atom_check,
        "build_provenance": build_provenance,
        "sources": sources,
    }
    return payload, summary_rows


def finalize(root: Path) -> tuple[Path, Path]:
    data = root / "data"
    csv_path = data / f"{SUMMARY_STEM}.csv"
    json_path = data / f"{SUMMARY_STEM}.json"
    csv_temp = data / f".{SUMMARY_STEM}.csv.tmp.{os.getpid()}"
    json_temp = data / f".{SUMMARY_STEM}.json.tmp.{os.getpid()}"
    csv_path.unlink(missing_ok=True)
    json_path.unlink(missing_ok=True)
    csv_temp.unlink(missing_ok=True)
    json_temp.unlink(missing_ok=True)
    try:
        payload, rows = build_summary(root)
        write_csv(csv_temp, rows)
        payload["paper_summary_csv"] = {
            "path": relative_path(csv_path, root),
            "sha256": sha256(csv_temp),
            "size_bytes": csv_temp.stat().st_size,
        }
        json_temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(csv_temp, csv_path)
        os.replace(json_temp, json_path)
    except BaseException:
        # A two-file publication bundle has no useful partially committed state.
        csv_path.unlink(missing_ok=True)
        json_path.unlink(missing_ok=True)
        raise
    finally:
        csv_temp.unlink(missing_ok=True)
        json_temp.unlink(missing_ok=True)
    return csv_path, json_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Goldzak12 benchmark directory",
    )
    args = parser.parse_args()
    try:
        csv_path, json_path = finalize(args.root.resolve())
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"LC12 publication summary: {csv_path}")
    print(f"LC12 publication lineage: {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
