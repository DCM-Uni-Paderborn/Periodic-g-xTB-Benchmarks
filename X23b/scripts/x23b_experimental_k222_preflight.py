#!/usr/bin/env python3
"""Prepare, run, and collect the frozen-structure X23b GXTB k222 preflight.

This is an additive derivative-measurement phase.  It never edits or consumes
the quarantined Gamma CELL_OPT directories and it deliberately separates a
technically complete calculation from a later scientific approval decision.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Mapping

import x23b_common as common
import x23b_k222_force_stress_gate as derivative_gate
import x23b_kpoint_cellopt as cellopt


ROOT = Path(__file__).resolve().parents[1]
PHASE = "x23b_experimental_k222_preflight"
VARIANT = "experimental_k222_preflight"
SOURCE_POLICY = "experimental_reference"
MANIFEST_SCHEMA = 1
MANIFEST_NAME = "experimental_k222_preflight_manifest.json"
OUTPUT_NAME = "cp2k.out"
PROTOCOL_IDENTITY = {
    "source_policy": SOURCE_POLICY,
    "variant": VARIANT,
    "mesh": "MACDONALD 2 2 2 0.25 0.25 0.25",
    "symmetry": "SPGLIB reduced",
    "run_type": "ENERGY_FORCE",
    "stress": "ANALYTICAL GPa",
}
CSV_FIELDS = (
    "method",
    "system",
    "phase",
    "variant",
    "source_policy",
    "campaign_fingerprint_sha256",
    "program_ended",
    "scientific_status",
    "approved",
    "energy_hartree",
    "atom_count",
    "kpoint_count",
    "kpoint_count_source",
    "start_volume_A3",
    "max_force_hartree_per_bohr",
    "rms_force_hartree_per_bohr",
    "rms_force_component_hartree_per_bohr",
    "net_force_hartree_per_bohr",
    "stress_xx_GPa",
    "stress_xy_GPa",
    "stress_xz_GPa",
    "stress_yx_GPa",
    "stress_yy_GPa",
    "stress_yz_GPa",
    "stress_zx_GPa",
    "stress_zy_GPa",
    "stress_zz_GPa",
    "max_abs_stress_GPa",
    "pressure_GPa",
    "pressure_bar",
    "source_input",
    "source_input_sha256",
    "structure_path",
    "structure_sha256",
    "structure_source",
    "input",
    "input_sha256",
    "output",
    "output_sha256",
)


def manifest_path(output_root: Path) -> Path:
    return output_root.resolve() / MANIFEST_NAME


def systems() -> list[dict[str, object]]:
    return cellopt.systems()


def _atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _replace_unique(lines: list[str], pattern: re.Pattern[str], value: str, label: str) -> None:
    matches = [index for index, line in enumerate(lines) if pattern.match(line)]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {label}, found {len(matches)}")
    index = matches[0]
    match = pattern.match(lines[index])
    assert match is not None
    lines[index] = f"{match.group(1)}{value}"


def preflight_input_text(source: Path, system: str, project: str) -> str:
    source_input, structure = cellopt.experimental_reference_paths(system)
    if source.resolve(strict=True) != source_input.resolve(strict=True):
        raise ValueError(f"noncanonical preflight source for {system}: {source}")
    cellopt.validate_experimental_reference_source(system, source, structure)
    lines = source.read_text().splitlines()
    _replace_unique(
        lines,
        re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I),
        f'"{project}"',
        "PROJECT",
    )
    _replace_unique(
        lines,
        re.compile(r"^(\s*RUN_TYPE\s+).*$", re.I),
        "ENERGY_FORCE",
        "RUN_TYPE",
    )

    k_start, k_end = cellopt.section_bounds(lines, "KPOINTS")
    if not any(re.match(r"^\s*VERBOSE\s+T\s*$", line, flags=re.I) for line in lines[k_start:k_end]):
        indent = re.match(r"\s*", lines[k_start]).group(0) + "  "
        lines.insert(k_end, f"{indent}VERBOSE T")

    force_start, force_end = cellopt.section_bounds(lines, "FORCE_EVAL")
    if any(line.strip().upper() == "&PRINT" for line in lines[force_start:force_end]):
        raise ValueError(f"frozen source unexpectedly contains FORCE_EVAL/PRINT: {source}")
    indent = re.match(r"\s*", lines[force_start]).group(0) + "  "
    lines[force_end:force_end] = [
        f"{indent}&PRINT",
        f"{indent}  &FORCES ON",
        f"{indent}    NDIGITS 12",
        f"{indent}  &END FORCES",
        f"{indent}  &STRESS_TENSOR ON",
        f"{indent}    STRESS_UNIT GPa",
        f"{indent}  &END STRESS_TENSOR",
        f"{indent}&END PRINT",
    ]
    text = "\n".join(lines) + "\n"
    validate_preflight_input(text, system)
    return text


def validate_preflight_input(text: str, system: str) -> None:
    common.validate_method_input(text, "GXTB")
    required = {
        "RUN_TYPE ENERGY_FORCE": r"^\s*RUN_TYPE\s+ENERGY_FORCE\s*$",
        "analytical stress": r"^\s*STRESS_TENSOR\s+ANALYTICAL\s*$",
        "shifted k222": (
            r"^\s*SCHEME\s+MACDONALD\s+2\s+2\s+2\s+"
            r"0\.25\s+0\.25\s+0\.25\s*$"
        ),
        "verbose k-point mapping": r"^\s*VERBOSE\s+T\s*$",
        "SPGLIB symmetry": r"^\s*SYMMETRY\s+T\s*$",
        "reduced grid": r"^\s*FULL_GRID\s+F\s*$",
        "SPGLIB backend": r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$",
        "SPGLIB reduction": r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$",
        "force print": r"^\s*&FORCES\s+ON\s*$",
        "force precision": r"^\s*NDIGITS\s+12\s*$",
        "stress print": r"^\s*&STRESS_TENSOR\s+ON\s*$",
        "GPa stress unit": r"^\s*STRESS_UNIT\s+GPa\s*$",
    }
    for description, pattern in required.items():
        if len(re.findall(pattern, text, flags=re.I | re.M)) != 1:
            raise ValueError(f"preflight input lacks unique {description} for {system}")
    if re.search(r"^\s*&MOTION\b", text, flags=re.I | re.M):
        raise ValueError("experimental k222 preflight must not contain MOTION")
    metadata = next(row for row in systems() if str(row["id"]) == system)
    if len(cellopt._coord_elements(text)) != int(metadata["n_atoms_crystal"]):
        raise ValueError(f"preflight atom count differs for {system}")


def _cell_volume(text: str) -> float:
    lines = text.splitlines()
    start, end = cellopt.section_bounds(lines, "CELL")
    vectors: dict[str, list[float]] = {}
    for line in lines[start + 1 : end]:
        fields = line.split()
        if len(fields) == 4 and fields[0].upper() in {"A", "B", "C"}:
            vectors[fields[0].upper()] = [float(value.replace("D", "E")) for value in fields[1:]]
    if set(vectors) != {"A", "B", "C"}:
        raise ValueError("preflight input lacks explicit A/B/C cell vectors")
    a, b, c = (vectors[name] for name in ("A", "B", "C"))
    volume = abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("invalid preflight start volume")
    return volume


def source_artifacts(record: Mapping[str, object]) -> dict[str, Path]:
    return {
        "reference_input": Path(str(record["source_input"])),
        "reference_structure": Path(str(record["structure_path"])),
    }


def _case_record(output_root: Path, metadata: Mapping[str, object]) -> dict[str, object]:
    system = str(metadata["id"])
    source, structure = cellopt.experimental_reference_paths(system)
    cellopt.validate_experimental_reference_source(system, source, structure)
    run_dir = output_root / "GXTB" / system / VARIANT
    project = f"{system}_GXTB_{VARIANT}".replace("-", "_")
    input_path = run_dir / f"{project}.inp"
    output = run_dir / OUTPUT_NAME
    text = preflight_input_text(source, system, project)
    if abs(_cell_volume(text) - float(metadata["input_volume"])) > 1.0e-5:
        raise ValueError(f"preflight start volume differs from metadata for {system}")
    if output.exists():
        raise ValueError(f"refusing to prepare over existing preflight output: {output}")
    run_dir.mkdir(parents=True, exist_ok=True)
    if input_path.exists() and input_path.read_text() != text:
        raise ValueError(f"stale preflight input differs: {input_path}")
    input_path.write_text(text)
    return {
        "method": "GXTB",
        "system": system,
        "phase": PHASE,
        "variant": VARIANT,
        "source_policy": SOURCE_POLICY,
        "source_input": str(source.resolve(strict=True)),
        "source_input_sha256": common.sha256_file(source),
        "structure_path": str(structure.resolve(strict=True)),
        "structure_sha256": common.sha256_file(structure),
        "structure_source": str(metadata["structure_source"]),
        "input": str(input_path.resolve()),
        "input_sha256": common.sha256_file(input_path),
        "run_dir": str(run_dir.resolve()),
        "output": str(output.resolve()),
        "start_volume_A3": _cell_volume(text),
        "atom_count": int(metadata["n_atoms_crystal"]),
    }


def prepare(output_root: Path, campaign_identity: Mapping[str, object]) -> Path:
    output_root = output_root.resolve()
    common.validate_campaign_identity(campaign_identity)
    path = manifest_path(output_root)
    records = [_case_record(output_root, row) for row in systems()]
    payload: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "phase": PHASE,
        "variant": VARIANT,
        "source_policy": SOURCE_POLICY,
        "campaign_identity": dict(campaign_identity),
        "systems": records,
    }
    if path.exists():
        observed = json.loads(path.read_text())
        if observed != payload:
            raise ValueError(f"existing preflight manifest differs: {path}")
    else:
        _atomic_json(path, payload)
    return path


def load_manifest(
    output_root: Path,
    campaign_identity: Mapping[str, object],
) -> dict[str, dict[str, object]]:
    output_root = output_root.resolve()
    path = manifest_path(output_root)
    payload = json.loads(path.read_text())
    if (
        payload.get("schema") != MANIFEST_SCHEMA
        or payload.get("phase") != PHASE
        or payload.get("variant") != VARIANT
        or payload.get("source_policy") != SOURCE_POLICY
    ):
        raise ValueError(f"invalid experimental preflight manifest contract: {path}")
    if payload.get("campaign_identity") != dict(campaign_identity):
        raise ValueError("preflight manifest campaign fingerprint differs")
    raw_records = payload.get("systems")
    if not isinstance(raw_records, list):
        raise ValueError("preflight manifest has no system records")
    records: dict[str, dict[str, object]] = {}
    metadata = {str(row["id"]): row for row in systems()}
    for record in raw_records:
        if not isinstance(record, dict):
            raise ValueError("invalid preflight system record")
        system = str(record.get("system", ""))
        if system in records or system not in metadata:
            raise ValueError(f"duplicate or unknown preflight system: {system}")
        if (
            record.get("method") != "GXTB"
            or record.get("phase") != PHASE
            or record.get("variant") != VARIANT
            or record.get("source_policy") != SOURCE_POLICY
        ):
            raise ValueError(f"preflight identity differs for {system}")
        expected_source, expected_structure = cellopt.experimental_reference_paths(system)
        source = Path(str(record["source_input"])).resolve()
        structure = Path(str(record["structure_path"])).resolve()
        if source != expected_source.resolve(strict=True) or structure != expected_structure.resolve(strict=True):
            raise ValueError(f"noncanonical preflight lineage for {system}")
        if common.sha256_file(source) != record.get("source_input_sha256"):
            raise ValueError(f"preflight source input fingerprint differs for {system}")
        if common.sha256_file(structure) != record.get("structure_sha256"):
            raise ValueError(f"preflight structure fingerprint differs for {system}")
        if record.get("structure_source") != metadata[system]["structure_source"]:
            raise ValueError(f"preflight structure source label differs for {system}")
        cellopt.validate_experimental_reference_source(system, source, structure)
        input_path = Path(str(record["input"])).resolve()
        run_dir = Path(str(record["run_dir"])).resolve()
        output = Path(str(record["output"])).resolve()
        if input_path.parent != run_dir or output.parent != run_dir or output.name != OUTPUT_NAME:
            raise ValueError(f"preflight path contract differs for {system}")
        try:
            run_dir.relative_to(output_root)
        except ValueError as exc:
            raise ValueError(f"preflight run directory escapes output root: {run_dir}") from exc
        expected_run_dir = output_root / "GXTB" / system / VARIANT
        expected_project = f"{system}_GXTB_{VARIANT}".replace("-", "_")
        if run_dir != expected_run_dir or input_path != run_dir / f"{expected_project}.inp":
            raise ValueError(f"noncanonical preflight run/input path for {system}")
        if common.sha256_file(input_path) != record.get("input_sha256"):
            raise ValueError(f"preflight generated input fingerprint differs for {system}")
        expected_text = preflight_input_text(source, system, expected_project)
        if input_path.read_text() != expected_text:
            raise ValueError(f"preflight input is not reproducible for {system}")
        if abs(_cell_volume(expected_text) - float(record["start_volume_A3"])) > 1.0e-9:
            raise ValueError(f"preflight start-volume manifest differs for {system}")
        if int(record["atom_count"]) != int(metadata[system]["n_atoms_crystal"]):
            raise ValueError(f"preflight atom-count manifest differs for {system}")
        records[system] = record
    expected = set(metadata)
    if set(records) != expected:
        raise ValueError(
            f"complete 23-system preflight manifest required: "
            f"missing={sorted(expected - set(records))}, unexpected={sorted(set(records) - expected)}"
        )
    return records


def derivative_summary(output: Path, record: Mapping[str, object]) -> dict[str, object]:
    input_path = Path(str(record["input"]))
    elements = cellopt._coord_elements(input_path.read_text())
    parsed = derivative_gate.parse_cp2k_output(
        output,
        int(record["atom_count"]),
        elements,
    )
    if parsed["kpoint_mesh"] != [2, 2, 2]:
        raise ValueError(f"preflight output reports the wrong k-point mesh: {output}")
    count = int(parsed["kpoint_count"])
    if not 1 <= count <= 8:
        raise ValueError(f"invalid shifted-k222 irreducible count in {output}: {count}")
    if parsed["kpoint_mesh_row_indices"] != list(range(1, 9)):
        raise ValueError(f"incomplete shifted-k222 full-mesh mapping in {output}")
    forces = [list(row["vector_hartree_per_bohr"]) for row in parsed["forces"]]
    stress = [list(row) for row in parsed["stress_gpa"]]
    flat = [float(value) for vector in forces for value in vector]
    flat_stress = [float(value) for row in stress for value in row]
    values = [float(parsed["energy_hartree"]), *flat, *flat_stress]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"nonfinite preflight energy/force/stress in {output}")
    norms = [math.sqrt(sum(float(value) ** 2 for value in vector)) for vector in forces]
    net = [sum(float(vector[axis]) for vector in forces) for axis in range(3)]
    pressure = -(stress[0][0] + stress[1][1] + stress[2][2]) / 3.0
    return {
        "energy_hartree": float(parsed["energy_hartree"]),
        "atom_count": len(forces),
        "kpoint_count": count,
        "kpoint_count_source": parsed["kpoint_count_source"],
        "max_force_hartree_per_bohr": max(norms),
        "rms_force_hartree_per_bohr": math.sqrt(sum(value * value for value in norms) / len(norms)),
        "rms_force_component_hartree_per_bohr": math.sqrt(
            sum(value * value for value in flat) / len(flat)
        ),
        "net_force_hartree_per_bohr": math.sqrt(sum(value * value for value in net)),
        "stress_gpa": stress,
        "max_abs_stress_GPa": max(abs(value) for value in flat_stress),
        "pressure_GPa": pressure,
        "pressure_bar": pressure * 10000.0,
    }


def validate_completed_case(
    output_root: Path,
    system: str,
    campaign_identity: Mapping[str, object],
) -> dict[str, object]:
    records = load_manifest(output_root, campaign_identity)
    if system not in records:
        raise ValueError(f"preflight manifest has no {system}")
    record = records[system]
    input_path = Path(str(record["input"]))
    output = Path(str(record["output"]))
    valid, reason = common.recorded_job_stamp_matches(
        output.parent,
        input_path,
        "GXTB",
        PHASE,
        output,
        campaign_identity=campaign_identity,
        accepted_status_prefixes=("converged",),
        protocol_identity=PROTOCOL_IDENTITY,
        source_artifacts=source_artifacts(record),
    )
    if not valid:
        raise ValueError(f"untrusted experimental k222 preflight for {system}: {reason}")
    summary = derivative_summary(output, record)
    return {**record, **summary}


def run_one(
    record: Mapping[str, object],
    cp2k: Path,
    threads: int,
    campaign_identity: Mapping[str, object],
) -> tuple[str, int, str]:
    system = str(record["system"])
    input_path = Path(str(record["input"]))
    run_dir = Path(str(record["run_dir"]))
    output = Path(str(record["output"]))
    common.validate_method_input(input_path.read_text(), "GXTB")
    with input_path.open() as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return system, common.BUSY_RETURN_CODE, "BUSY"
        stamp_ok, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            "GXTB",
            PHASE,
            campaign_identity=campaign_identity,
            protocol_identity=PROTOCOL_IDENTITY,
            source_artifacts=source_artifacts(record),
        )
        if output.exists():
            recorded_ok, _ = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                "GXTB",
                PHASE,
                output,
                campaign_identity=campaign_identity,
                accepted_status_prefixes=("converged",),
                protocol_identity=PROTOCOL_IDENTITY,
                source_artifacts=source_artifacts(record),
            )
            if not stamp_ok or not recorded_ok:
                return system, 1, "STALE_OUTPUT"
            try:
                derivative_summary(output, record)
            except ValueError:
                return system, 1, "STALE_OUTPUT"
            return system, 0, "SKIP"

        process = subprocess.run(
            [str(cp2k.resolve(strict=True)), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            env=common.thread_environment(threads),
            check=False,
        )
        (run_dir / "returncode.txt").write_text(f"{process.returncode}\n")
        details: dict[str, object] = {"returncode": process.returncode, "output": str(output)}
        action = "FAILED"
        status = "failed"
        code = process.returncode
        if output.is_file():
            details["output_sha256"] = common.sha256_file(output)
        if process.returncode == 0 and output.is_file():
            try:
                summary = derivative_summary(output, record)
            except ValueError as exc:
                details["parse_error"] = str(exc)
                action = "INVALID_OUTPUT"
                status = "failed_parse"
                code = 1
            else:
                details["derivative_summary"] = summary
                action = "COMPLETED"
                status = "converged_measured_not_approved"
                code = 0
        common.write_job_stamp(
            run_dir,
            input_path,
            cp2k,
            "GXTB",
            PHASE,
            status,
            details=details,
            campaign_identity=campaign_identity,
            protocol_identity=PROTOCOL_IDENTITY,
            source_artifacts=source_artifacts(record),
        )
        return system, code, action


def _register_provenance_root(output_root: Path, **artifacts: Path | None) -> None:
    provenance = ROOT / "data" / common.GXTB_PROVENANCE_NAME
    if provenance.is_file():
        payload = json.loads(provenance.read_text())
        paths = payload.get("workflow_paths", {})
        if isinstance(paths, dict) and paths.get("experimental_k222_preflight_root"):
            frozen = Path(str(paths["experimental_k222_preflight_root"])).resolve()
            if frozen != output_root.resolve():
                raise ValueError(f"experimental preflight root is already frozen as {frozen}")
    common.update_gxtb_provenance(
        ROOT,
        cp2k=artifacts.get("cp2k"),
        cp2k_source=artifacts.get("cp2k_source"),
        save_tblite=artifacts.get("save_tblite"),
        save_tblite_source=artifacts.get("save_tblite_source"),
        campaign_manifest=artifacts.get("campaign_manifest"),
        workflow_paths={"experimental_k222_preflight_root": output_root},
    )


def prepare_command(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / f".{output_root.name}.experimental_k222_preflight.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        campaign = common.load_campaign_identity(ROOT)
        path = prepare(output_root, campaign)
    print(f"Prepared 23 additive experimental-k222 preflight inputs: {path}")


def run_command(args: argparse.Namespace) -> None:
    common.require_gxtb_build_artifacts(
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    _register_provenance_root(
        args.output_root,
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    campaign = common.load_campaign_identity(ROOT)
    records = load_manifest(args.output_root, campaign)
    selected = sorted(set(args.system)) if args.system else sorted(records)
    failed: list[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_one,
                records[system],
                args.cp2k,
                args.threads_per_job,
                campaign,
            ): system
            for system in selected
        }
        for future in as_completed(futures):
            system, code, action = future.result()
            print(f"{action:14s} GXTB/{system} rc={code}", flush=True)
            if code != 0:
                failed.append(system)
    _register_provenance_root(args.output_root)
    if failed:
        raise SystemExit(f"{len(failed)} experimental k222 preflight job(s) failed")


def _format(value: float, digits: int = 12) -> str:
    return f"{value:.{digits}f}"


def collect_command(args: argparse.Namespace) -> None:
    campaign = common.load_campaign_identity(ROOT)
    records = load_manifest(args.output_root, campaign)
    rows: list[dict[str, object]] = []
    for system in sorted(records):
        result = validate_completed_case(args.output_root, system, campaign)
        stress = result["stress_gpa"]
        assert isinstance(stress, list)
        output = Path(str(result["output"]))
        row: dict[str, object] = {
            "method": "GXTB",
            "system": system,
            "phase": PHASE,
            "variant": VARIANT,
            "source_policy": SOURCE_POLICY,
            "campaign_fingerprint_sha256": campaign["fingerprint_sha256"],
            "program_ended": True,
            "scientific_status": "measured_not_approved",
            "approved": False,
            "energy_hartree": _format(float(result["energy_hartree"])),
            "atom_count": result["atom_count"],
            "kpoint_count": result["kpoint_count"],
            "kpoint_count_source": result["kpoint_count_source"],
            "start_volume_A3": _format(float(result["start_volume_A3"]), 8),
            "max_force_hartree_per_bohr": _format(float(result["max_force_hartree_per_bohr"])),
            "rms_force_hartree_per_bohr": _format(float(result["rms_force_hartree_per_bohr"])),
            "rms_force_component_hartree_per_bohr": _format(
                float(result["rms_force_component_hartree_per_bohr"])
            ),
            "net_force_hartree_per_bohr": _format(float(result["net_force_hartree_per_bohr"])),
            "max_abs_stress_GPa": _format(float(result["max_abs_stress_GPa"])),
            "pressure_GPa": _format(float(result["pressure_GPa"])),
            "pressure_bar": _format(float(result["pressure_bar"]), 6),
            "source_input": result["source_input"],
            "source_input_sha256": result["source_input_sha256"],
            "structure_path": result["structure_path"],
            "structure_sha256": result["structure_sha256"],
            "structure_source": result["structure_source"],
            "input": result["input"],
            "input_sha256": result["input_sha256"],
            "output": result["output"],
            "output_sha256": common.sha256_file(output),
        }
        for axis, name in enumerate("xyz"):
            for component, other in enumerate("xyz"):
                row[f"stress_{name}{other}_GPa"] = _format(float(stress[axis][component]))
        rows.append(row)
    _register_provenance_root(args.output_root)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Collected 23 parsed preflight measurements in {args.csv}")


def add_artifact_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument(
        "--campaign-manifest",
        type=Path,
        default=common.DEFAULT_CAMPAIGN_MANIFEST,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.set_defaults(function=prepare_command)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--system", action="append", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--jobs", type=int, default=4)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    add_artifact_arguments(run_parser)
    run_parser.set_defaults(function=run_command)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.set_defaults(function=collect_command)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
