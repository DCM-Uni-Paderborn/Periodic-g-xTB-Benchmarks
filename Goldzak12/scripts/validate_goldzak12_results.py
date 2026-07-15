#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import run_goldzak12_benchmark as base
import run_goldzak12_eos_benchmark as eos


ROOT = base.ROOT


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def validate_job_stamp(
    result: Path,
    campaign_fingerprint: dict[str, object],
    *,
    executable_role: str,
    require_completed: bool = True,
) -> str | None:
    return base.completed_stamp_campaign_issue(
        result,
        campaign_fingerprint,
        executable_role=executable_role,
        require_completed=require_completed,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eos-mesh", default="k444")
    parser.add_argument("--energy-mesh", action="append", default=[])
    parser.add_argument("--result-mesh", default=eos.DEFAULT_RESULT_MESH)
    parser.add_argument("--method", action="append", choices=base.METHODS)
    parser.add_argument(
        "--allow-reduced-coverage",
        action="store_true",
        help="deprecated and rejected: the paper benchmark requires exact LC10 coverage",
    )
    parser.add_argument(
        "--minimum-valid-fits",
        type=int,
        default=eos.MINIMUM_REDUCED_GXTB_FITS,
    )
    args = parser.parse_args()
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]
    methods = base.selected_methods(args.method)
    if args.allow_reduced_coverage:
        parser.error(
            "the paper benchmark has fixed LC10 coverage; reduced coverage is not allowed"
        )
    if not 1 <= args.minimum_valid_fits <= len(base.REFERENCES):
        parser.error(f"--minimum-valid-fits must be between 1 and {len(base.REFERENCES)}")

    expected_pairs = {
        (solid, method) for solid in base.LC10_PAPER_SOLIDS for method in methods
    }
    problems: list[str] = []

    if "GXTB" in methods:
        atom_rows = [
            row
            for row in read_csv(ROOT / "data" / "atom_energies_save_tblite_cli_gxtb.csv")
            if row["method"] == "GXTB"
        ]
        elements = list(base.LC10_PAPER_ELEMENTS)
        if sorted(row["element"] for row in atom_rows) != elements:
            problems.append(f"GXTB atom references: expected {len(elements)} elements, found {len(atom_rows)}")
        for row in atom_rows:
            element = row["element"]
            if row.get("source") != "save_tblite_cli":
                problems.append(f"GXTB/{element} atom reference is not labelled save_tblite_cli")
            if row.get("spin_2S") != str(base.ELEMENT_MULTIPLICITY[element] - 1):
                problems.append(f"GXTB/{element} atom spin metadata is inconsistent")

    points = [
        row
        for row in read_csv(ROOT / "data" / "eos_points.csv")
        if row["method"] in methods and row["solid"] in base.LC10_PAPER_SOLIDS
    ]
    expected_point_keys: set[tuple[str, str, str]] = set()
    scale_manifest = eos.read_gxtb_scale_manifest() if "GXTB" in methods else None
    if "GXTB" in methods and scale_manifest is None:
        problems.append("Missing or invalid gxtb_eos_scale_manifest.json")
    elif scale_manifest is not None and scale_manifest.get("eos_mesh") != args.eos_mesh:
        problems.append("GXTB scale manifest EOS mesh does not match validation request")
    for solid, method in expected_pairs:
        if method != "GXTB":
            expected_scales = eos.scales_for(solid, method, eos.DEFAULT_SCALES)
        else:
            if scale_manifest is None:
                expected_scales = eos.scales_for(solid, method, eos.DEFAULT_SCALES)
            elif scale_manifest.get("eos_mesh") != args.eos_mesh:
                expected_scales = ()
            else:
                record = next(
                    (
                        item
                        for item in scale_manifest.get("systems", [])
                        if item.get("solid") == solid and item.get("method") == "GXTB"
                    ),
                    None,
                )
                if record is None:
                    problems.append(f"GXTB scale manifest has no record for {solid}")
                    expected_scales = ()
                else:
                    expected_scales = tuple(float(value) for value in record.get("requested_scales", []))
        expected_point_keys.update(
            (solid, method, f"{scale:.5f}") for scale in expected_scales
        )
    point_keys = {(row["solid"], row["method"], f"{float(row['scale']):.5f}") for row in points}
    if point_keys != expected_point_keys or len(points) != len(expected_point_keys):
        problems.append(
            "EOS point coverage differs from the persisted scale manifest/protocol: "
            f"missing {sorted(expected_point_keys - point_keys)}, "
            f"unexpected {sorted(point_keys - expected_point_keys)}"
        )
    for row in (point for point in points if point["method"] == "GXTB"):
        scale = float(row["scale"])
        project = eos.eos_project(row["solid"], "GXTB", row["mesh"], scale)
        input_path = (
            ROOT
            / "runs"
            / "eos"
            / "GXTB"
            / row["solid"]
            / row["mesh"]
            / eos.scale_tag(scale, "GXTB")
            / f"{project}.inp"
        )
        if not input_path.is_file():
            problems.append(f"Missing GXTB EOS input: {input_path}")
            continue
        try:
            base.validate_method_input(input_path.read_text(), "GXTB")
        except ValueError as exc:
            problems.append(f"Invalid GXTB EOS input {input_path}: {exc}")
    failed_points = [row for row in points if not truth(row["completed"])]
    gxtb_failed_points = [row for row in failed_points if row["method"] == "GXTB"]
    unclassified_failed_points = [
        row
        for row in gxtb_failed_points
        if row.get("classification_resolution") != "explicit_failure_classification"
        or not row.get("classification_rationale")
    ]
    if unclassified_failed_points:
        labels = ", ".join(
            f"{row['method']}/{row['solid']}@{row['scale']}:{row.get('diagnostic', '')}"
            for row in unclassified_failed_points
        )
        problems.append(
            f"EOS points neither completed nor explicitly classified "
            f"({len(unclassified_failed_points)}): {labels}"
        )
    if gxtb_failed_points:
        problems.append("The fixed LC10 benchmark does not allow failed g-XTB EOS points")
    unclassified_exclusions = [
        row
        for row in points
        if row["method"] == "GXTB"
        and truth(row["completed"])
        and not truth(row.get("valid_for_eos", row["completed"]))
        and (
            not row.get("diagnostic")
            or row.get("classification_resolution")
            not in {
                "explicit_exclusion",
                "explicit_failure_classification",
                "legacy_automatic_filter",
            }
        )
    ]
    if unclassified_exclusions:
        problems.append("EOS points excluded without an explicit classification")
    unresolved_points = [
        row for row in points if row.get("classification_resolution") == "unresolved_candidate"
    ]
    if unresolved_points:
        problems.append(
            "Unresolved GXTB SCC-branch candidates: "
            + ", ".join(f"{row['solid']}@{row['scale']}" for row in unresolved_points)
        )

    fits = [
        row
        for row in read_csv(ROOT / "data" / "eos_fits.csv")
        if row["method"] in methods and row["solid"] in base.LC10_PAPER_SOLIDS
    ]
    fit_pairs = {(row["solid"], row["method"]) for row in fits if row["eos_mesh"] == args.eos_mesh}
    if fit_pairs != expected_pairs:
        problems.append(f"EOS fit coverage differs: missing {sorted(expected_pairs - fit_pairs)}")
    legacy_bad_fits = [
        row
        for row in fits
        if row["method"] in base.LEGACY_METHODS
        and (row["a_eos_A"] == "" or row["fit_status"] != "quadratic")
    ]
    if legacy_bad_fits:
        labels = ", ".join(
            f"{row['method']}/{row['solid']}={row['fit_status']}"
            for row in legacy_bad_fits
        )
        problems.append(f"Invalid fixed-LC10 legacy EOS fits: {labels}")

    gxtb_fits = [row for row in fits if row["method"] == "GXTB"]
    valid_gxtb_fits = [row for row in gxtb_fits if row["fit_status"] == "quadratic" and row["a_eos_A"]]
    invalid_gxtb_fits = [row for row in gxtb_fits if row not in valid_gxtb_fits]
    if "GXTB" in methods and not valid_gxtb_fits:
        problems.append("GXTB has zero valid quadratic EOS fits")
    if invalid_gxtb_fits:
        problems.append(
            "GXTB requires exact 10/10 quadratic fit coverage: "
            + ", ".join(f"{row['solid']}={row['fit_status']}" for row in invalid_gxtb_fits)
        )
    for fit in invalid_gxtb_fits:
        for mesh in energy_meshes:
            input_path = eos.final_input_path(fit["solid"], "GXTB", mesh)
            if not input_path.is_file():
                continue
            lineage_path = eos.final_input_lineage_path(input_path)
            try:
                lineage = json.loads(lineage_path.read_text())
            except (FileNotFoundError, json.JSONDecodeError, OSError):
                problems.append(f"Stale GXTB final input lacks invalid EOS lineage: {input_path}")
                continue
            if (
                lineage.get("schema_version") != eos.FINAL_INPUT_LINEAGE_SCHEMA
                or lineage.get("valid") is not False
                or lineage.get("input_sha256") != base.sha256(input_path)
            ):
                problems.append(f"Stale GXTB final input is not explicitly invalidated: {input_path}")

    results = [
        row
        for row in read_csv(ROOT / "data" / "eos_results.csv")
        if row["method"] in methods and row["solid"] in base.LC10_PAPER_SOLIDS
    ]
    valid_pairs = {
        (row["solid"], row["method"])
        for row in fits
        if row["a_eos_A"] != ""
        and (row["method"] != "GXTB" or row["fit_status"] == "quadratic")
    }
    expected_results = {(solid, method, mesh) for solid, method in valid_pairs for mesh in energy_meshes}
    result_keys = {(row["solid"], row["method"], row["energy_mesh"]) for row in results}
    if result_keys != expected_results:
        problems.append(
            "Final result coverage differs: "
            f"missing {sorted(expected_results - result_keys)}, "
            f"unexpected {sorted(result_keys - expected_results)}"
        )
    fit_by_pair = {(row["solid"], row["method"]): row for row in fits}
    for row in (result for result in results if result["method"] == "GXTB"):
        fit = fit_by_pair.get((row["solid"], "GXTB"))
        if fit is None:
            problems.append(f"GXTB final result has no EOS fit: {row['solid']}/{row['energy_mesh']}")
            continue
        input_path = eos.final_input_path(row["solid"], "GXTB", row["energy_mesh"])
        issue = eos.final_input_lineage_issue(
            input_path,
            fit,
            row["energy_mesh"],
        )
        if issue:
            problems.append(issue)
    failed_sp = [row for row in results if not truth(row["sp_completed"])]
    if failed_sp:
        labels = ", ".join(f"{row['method']}/{row['solid']}/{row['energy_mesh']}" for row in failed_sp)
        problems.append(f"Incomplete final single points ({len(failed_sp)}): {labels}")

    provenance_paths = []
    if any(method in base.LEGACY_METHODS for method in methods):
        provenance_paths.append(ROOT / "data" / "build_provenance.json")
    if "GXTB" in methods:
        provenance_paths.append(ROOT / "data" / "build_provenance_gxtb.json")
        bad_strategies = [
            row
            for row in points
            if row["method"] == "GXTB"
            and row.get("scf_strategy", "")
            not in {"native_gxtb_fdiis", "native_gxtb_fdiis_adaptive"}
        ]
        if bad_strategies:
            problems.append("GXTB EOS points contain a non-production mixer strategy")
        branch_rows = read_csv(ROOT / "data" / "gxtb_eos_branch_diagnostics.csv")
        unresolved_branches = [
            row
            for row in branch_rows
            if row.get("resolution") == "unresolved_candidate"
            or (
                row.get("automatic_candidate") == "True"
                and (
                    not row.get("classification")
                    or not row.get("rationale")
                    or row.get("action") not in {"exclude", "retain"}
                )
            )
        ]
        if unresolved_branches:
            problems.append("GXTB branch diagnostics contain candidates without explicit classification/waiver")
    gxtb_provenance: dict[str, object] | None = None
    for provenance_path in provenance_paths:
        if not provenance_path.exists():
            problems.append(f"Missing provenance: {provenance_path.name}")
            continue
        provenance = json.loads(provenance_path.read_text())
        if provenance["protocol"]["result_mesh"] != args.result_mesh:
            problems.append(f"{provenance_path.name} result mesh does not match validation request")
        if provenance_path.name == "build_provenance_gxtb.json":
            gxtb_provenance = provenance

    if "GXTB" in methods and gxtb_provenance is not None:
        protocol = gxtb_provenance.get("protocol", {})
        if protocol.get("kpoint_mesh_contract") != base.KPOINT_MESH_CONTRACT:
            problems.append("GXTB provenance does not record the LC12 SPGLIB mesh contract")
        if protocol.get("legacy_gxtb_full_grid_policy") != base.LEGACY_GXTB_FULL_GRID_POLICY:
            problems.append("GXTB provenance does not exclude legacy unreduced full-grid results")
        if protocol.get("gxtb_energy_stress_policy") != base.GXTB_ENERGY_STRESS_POLICY:
            problems.append("GXTB provenance does not record the LC12 energy-only stress policy")
        if protocol.get("final_input_lineage_schema") != eos.FINAL_INPUT_LINEAGE_SCHEMA:
            problems.append("GXTB provenance final-input lineage schema is missing or incompatible")
        if protocol.get("fit_approval_required") is not True or protocol.get("fit_approved") is not True:
            problems.append("GXTB final results lack explicit EOS-fit approval")
        current_fit_sha = eos.gxtb_fit_approval_sha256(
            [dict(row) for row in fits if row["method"] == "GXTB"]
        )
        if protocol.get("approved_gxtb_fit_sha256") != current_fit_sha:
            problems.append("GXTB approved fit fingerprint differs from the current EOS fits")
        if protocol.get("allow_reduced_coverage") not in (None, False):
            problems.append("GXTB provenance must disable reduced coverage for fixed LC10")
        selected_solids = tuple(str(value) for value in protocol.get("selected_solids", ()))
        if (
            protocol.get("exact_lc10_scope") is not True
            or len(selected_solids) != len(base.LC10_PAPER_SOLIDS)
            or set(selected_solids) != set(base.LC10_PAPER_SOLIDS)
            or tuple(protocol.get("paper_systems", ())) != base.LC10_PAPER_SOLIDS
            or tuple(protocol.get("diagnostic_only_systems", ()))
            != base.LC10_DIAGNOSTIC_ONLY_SOLIDS
        ):
            problems.append("GXTB provenance does not certify the exact fixed LC10 scope")
        scale_path = eos.gxtb_scale_manifest_path()
        if not scale_path.is_file() or protocol.get("gxtb_scale_manifest_sha256") != base.sha256(scale_path):
            problems.append("GXTB scale manifest hash differs from build provenance")
        classification_path = Path(
            str(protocol.get("gxtb_classification_manifest", eos.gxtb_classification_manifest_path()))
        )
        classification_hash = (
            base.sha256(classification_path) if classification_path.is_file() else None
        )
        if protocol.get("gxtb_classification_manifest_sha256") != classification_hash:
            problems.append("GXTB classification manifest hash differs from build provenance")
        campaign = gxtb_provenance.get("campaign_identity", {})
        if not isinstance(campaign, dict) or not campaign:
            problems.append("GXTB provenance lacks the complete campaign identity")
            campaign = {}
        else:
            try:
                base.validate_campaign_identity(campaign)
                manifest_record = gxtb_provenance.get("campaign_manifest")
                if not isinstance(manifest_record, dict) or not manifest_record.get("path"):
                    raise ValueError("campaign manifest record is missing")
                manifest_path = Path(str(manifest_record["path"])).resolve(strict=True)
                manifest = json.loads(manifest_path.read_text())
                if base.campaign_identity_from_manifest(manifest, manifest_path) != campaign:
                    raise ValueError(
                        "current campaign manifest build identity differs from LC12 provenance"
                    )
            except (json.JSONDecodeError, OSError, TypeError, ValueError) as exc:
                problems.append(f"GXTB campaign provenance is invalid: {exc}")
        stamp_problems: list[str] = []
        for element in base.LC10_PAPER_ELEMENTS:
            atom_json = (
                ROOT
                / "runs"
                / "atoms_cli"
                / "GXTB"
                / element
                / f"atom_{element}_GXTB.json"
            )
            issue = validate_job_stamp(
                atom_json, campaign, executable_role="save_tblite"
            )
            if issue:
                stamp_problems.append(issue)
        for row in (point for point in points if point["method"] == "GXTB"):
            scale = float(row["scale"])
            project = eos.eos_project(row["solid"], "GXTB", row["mesh"], scale)
            output = (
                ROOT
                / "runs"
                / "eos"
                / "GXTB"
                / row["solid"]
                / row["mesh"]
                / eos.scale_tag(scale, "GXTB")
                / f"{project}.out"
            )
            issue = validate_job_stamp(
                output,
                campaign,
                executable_role="cp2k",
                require_completed=truth(row["completed"]),
            )
            if issue:
                stamp_problems.append(issue)
        for row in (result for result in results if result["method"] == "GXTB"):
            project = eos.final_project(row["solid"], "GXTB", row["energy_mesh"])
            output = (
                ROOT
                / "runs"
                / "eos_final_sp"
                / "GXTB"
                / row["solid"]
                / row["energy_mesh"]
                / f"{project}.out"
            )
            issue = validate_job_stamp(output, campaign, executable_role="cp2k")
            if issue:
                stamp_problems.append(issue)
        if stamp_problems:
            preview = "; ".join(stamp_problems[:5])
            problems.append(
                f"GXTB per-job provenance stamps failed ({len(stamp_problems)}): {preview}"
            )

    if problems:
        print("LC10 validation FAILED")
        for problem in problems:
            print(f"- {problem}")
        return 1

    print(
        f"LC10 validation passed: {len(points)} EOS points "
        f"({len(failed_points)} documented nonessential/classified failures), "
        f"{len(valid_pairs)}/{len(fits)} valid fits, {len(results)} final single points."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
