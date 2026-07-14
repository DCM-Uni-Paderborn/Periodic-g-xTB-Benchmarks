#!/usr/bin/env python3
"""Run and evaluate the X23b GXTB shifted-k222 SPGLIB derivative gate.

The gate is intentionally separate from production.  It admits the central
campaign while its state is ``validation_in_progress``, but it never changes
that state.  The normal X23b production entry points continue to require
``production_ready``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

import x23b_common as common


REPOSITORY = Path(__file__).resolve().parents[2]
GATE_ROOT = REPOSITORY / "X23b" / "validation" / "gxtb_k222_force_stress"
DEFAULT_SPEC = GATE_ROOT / "gate_spec.json"
DEFAULT_RUN_ROOT = (
    REPOSITORY / "X23b" / "runs" / "validation" / "gxtb_k222_force_stress_v1_final"
)
JOB_STAMP = "gate_job_provenance.json"
OUTPUT_NAME = "cp2k.out"
LAUNCHER_LOG = "launcher.log"
FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def repo_path(value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else REPOSITORY / path


def _required_line(text: str, pattern: str, description: str) -> None:
    if re.search(pattern, text, flags=re.I | re.M) is None:
        raise ValueError(f"gate input is missing {description}")


def _coord_elements(text: str) -> list[str]:
    lines = text.splitlines()
    start = next(
        (index for index, line in enumerate(lines) if line.strip().upper() == "&COORD"),
        None,
    )
    if start is None:
        raise ValueError("gate input has no COORD section")
    elements: list[str] = []
    for line in lines[start + 1 :]:
        stripped = line.strip()
        if stripped.upper() == "&END COORD":
            return elements
        if not stripped or stripped.upper() == "SCALED" or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Za-z]{1,3}", fields[0]):
            elements.append(fields[0].capitalize())
    raise ValueError("unterminated COORD section in gate input")


def _normalized_pair_input(text: str) -> str:
    differing_keywords = (
        "SYMMETRY",
        "FULL_GRID",
        "SYMMETRY_BACKEND",
        "SYMMETRY_REDUCTION_METHOD",
    )
    normalized: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if re.match(r"^PROJECT(?:_NAME)?\s+", stripped, flags=re.I):
            normalized.append("  PROJECT <VARIANT>")
        elif any(re.match(rf"^{keyword}\s+", stripped, flags=re.I) for keyword in differing_keywords):
            continue
        else:
            normalized.append(line.rstrip())
    return "\n".join(normalized) + "\n"


def validate_gate_inputs(spec: Mapping[str, object]) -> dict[str, Path]:
    inputs = spec.get("inputs")
    structure = spec.get("structure")
    if not isinstance(inputs, Mapping) or not isinstance(structure, Mapping):
        raise ValueError("gate specification has no inputs/structure records")
    structure_path = repo_path(structure.get("path"))
    if common.sha256_file(structure_path) != structure.get("sha256"):
        raise ValueError("frozen ammonia structure differs from the gate specification")

    paths: dict[str, Path] = {}
    texts: dict[str, str] = {}
    for variant in ("full", "spglib"):
        record = inputs.get(variant)
        if not isinstance(record, Mapping):
            raise ValueError(f"gate specification has no {variant} input")
        path = repo_path(record.get("path"))
        if common.sha256_file(path) != record.get("sha256"):
            raise ValueError(f"versioned {variant} gate input differs from gate_spec.json")
        paths[variant] = path
        texts[variant] = path.read_text()

    for variant, text in texts.items():
        _required_line(text, r"^\s*RUN_TYPE\s+ENERGY_FORCE\s*$", "RUN_TYPE ENERGY_FORCE")
        _required_line(text, r"^\s*STRESS_TENSOR\s+ANALYTICAL\s*$", "analytical stress")
        _required_line(text, r"^\s*METHOD\s+GXTB\s*$", "METHOD GXTB")
        _required_line(text, r"^\s*SCC_MIXER\s+TBLITE\s*$", "native TBLITE mixer")
        _required_line(text, r"^\s*EPS_DEFAULT\s+1\.0E-12\s*$", "EPS_DEFAULT 1e-12")
        _required_line(text, r"^\s*EPS_SCF\s+1\.0E-9\s*$", "EPS_SCF 1e-9")
        _required_line(text, r"^\s*MAX_SCF\s+300\s*$", "MAX_SCF 300")
        _required_line(text, r"^\s*SCF_GUESS\s+MOPAC\s*$", "SCF_GUESS MOPAC")
        _required_line(
            text,
            r"^\s*SCHEME\s+MACDONALD\s+2\s+2\s+2\s+0\.25\s+0\.25\s+0\.25\s*$",
            "shifted MACDONALD 2x2x2 mesh",
        )
        _required_line(text, r"^\s*VERBOSE\s+T\s*$", "auditable k-point-count output")
        _required_line(text, r"^\s*&FORCES\s+ON\s*$", "explicit force output")
        _required_line(text, r"^\s*NDIGITS\s+12\s*$", "12-digit force output")
        _required_line(text, r"^\s*&STRESS_TENSOR\s+ON\s*$", "explicit stress output")
        _required_line(text, r"^\s*STRESS_UNIT\s+GPa\s*$", "GPa stress output")
        if re.search(r"^\s*&MOTION\b", text, flags=re.I | re.M):
            raise ValueError(f"{variant} gate input must be a single ENERGY_FORCE evaluation")

    _required_line(texts["full"], r"^\s*SYMMETRY\s+F\s*$", "full-grid SYMMETRY F")
    _required_line(texts["full"], r"^\s*FULL_GRID\s+T\s*$", "FULL_GRID T")
    if "SPGLIB" in texts["full"].upper():
        raise ValueError("full-grid gate input must not request SPGLIB reduction")
    common.validate_method_input(texts["spglib"], "GXTB")

    if _normalized_pair_input(texts["full"]) != _normalized_pair_input(texts["spglib"]):
        raise ValueError(
            "full and SPGLIB inputs differ outside PROJECT and the k-point reduction contract"
        )
    elements = {variant: _coord_elements(text) for variant, text in texts.items()}
    expected_atoms = int(structure.get("atom_count", 0))
    if elements["full"] != elements["spglib"] or len(elements["full"]) != expected_atoms:
        raise ValueError("gate input atom ordering/count is inconsistent")
    return paths


def load_spec(spec_path: Path) -> dict[str, object]:
    path = spec_path.resolve(strict=True)
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict) or payload.get("schema") != 1:
        raise ValueError(f"unsupported gate specification: {path}")
    campaign = payload.get("campaign_identity")
    if not isinstance(campaign, dict):
        raise ValueError("gate specification has no campaign identity")
    common.validate_campaign_identity(campaign)
    tolerances = payload.get("tolerances")
    if not isinstance(tolerances, dict) or any(float(value) <= 0.0 for value in tolerances.values()):
        raise ValueError("gate specification has invalid tolerances")
    validate_gate_inputs(payload)
    return payload


def allowed_states(spec: Mapping[str, object]) -> tuple[str, ...]:
    values = spec.get("allowed_campaign_states")
    if not isinstance(values, list) or not values or not all(isinstance(value, str) for value in values):
        raise ValueError("gate specification has no allowed campaign states")
    return tuple(values)


def validate_manifest_choice(spec: Mapping[str, object], campaign_manifest: Path) -> Path:
    expected = repo_path(spec.get("campaign_manifest")).resolve(strict=True)
    observed = campaign_manifest.resolve(strict=True)
    if observed != expected:
        raise ValueError(f"gate is pinned to campaign manifest {expected}, not {observed}")
    identity, _ = common.declared_campaign_identity(
        observed,
        allowed_campaign_states=allowed_states(spec),
    )
    if identity != spec.get("campaign_identity"):
        raise ValueError("central campaign build identity differs from the gate specification")
    return observed


def validate_artifacts(
    spec: Mapping[str, object],
    *,
    cp2k: Path,
    cp2k_source: Path,
    save_tblite: Path,
    save_tblite_source: Path,
    campaign_manifest: Path,
) -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    manifest_path = validate_manifest_choice(spec, campaign_manifest)
    records = common.validate_campaign_artifacts(
        cp2k=cp2k,
        cp2k_source=cp2k_source,
        save_tblite=save_tblite,
        save_tblite_source=save_tblite_source,
        campaign_manifest=manifest_path,
        allowed_campaign_states=allowed_states(spec),
    )
    if records[0] != spec.get("campaign_identity"):
        raise ValueError("observed artifacts differ from the gate campaign identity")
    return records


def command_for(cp2k: Path, input_path: Path, output: Path) -> list[str]:
    return [str(cp2k.resolve(strict=True)), "-i", input_path.name, "-o", output.name]


def immutable_stamp_fields(payload: Mapping[str, object]) -> dict[str, object]:
    return {
        key: payload.get(key)
        for key in (
            "schema",
            "gate_id",
            "variant",
            "campaign_identity",
            "specification",
            "structure",
            "source_input",
            "run_input",
            "command",
            "thread_environment",
        )
    }


def build_stamp(
    spec: Mapping[str, object],
    spec_path: Path,
    variant: str,
    run_root: Path,
    records: tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]],
) -> tuple[Path, dict[str, object]]:
    identity, cp2k_record, save_record, manifest_record = records
    input_records = spec["inputs"]
    structure = spec["structure"]
    assert isinstance(input_records, Mapping) and isinstance(structure, Mapping)
    variant_record = input_records[variant]
    assert isinstance(variant_record, Mapping)
    source_input = repo_path(variant_record["path"]).resolve(strict=True)
    run_dir = run_root.resolve() / variant
    run_input = run_dir / source_input.name
    output = run_dir / OUTPUT_NAME
    command = command_for(Path(str(cp2k_record["path"])), run_input, output)
    env = common.thread_environment(1)
    selected_env = {
        key: env[key]
        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "OMP_WAIT_POLICY",
        )
    }
    payload: dict[str, object] = {
        "schema": 1,
        "gate_id": spec["gate_id"],
        "variant": variant,
        "status": "prepared",
        "prepared_at_utc": now_utc(),
        "campaign_identity": identity,
        "campaign_manifest": manifest_record,
        "specification": {
            "path": str(spec_path.resolve(strict=True)),
            "sha256": common.sha256_file(spec_path),
        },
        "structure": {
            "path": str(repo_path(structure["path"]).resolve(strict=True)),
            "sha256": structure["sha256"],
        },
        "source_input": {
            "path": str(source_input),
            "sha256": variant_record["sha256"],
        },
        "run_input": {
            "path": str(run_input),
            "sha256": variant_record["sha256"],
        },
        "output": {"path": str(output)},
        "command": command,
        "thread_environment": selected_env,
        "cp2k": cp2k_record,
        "save_tblite": save_record,
    }
    return run_dir, payload


def prepare_variant(
    spec: Mapping[str, object],
    spec_path: Path,
    variant: str,
    run_root: Path,
    records: tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]],
) -> Path:
    run_dir, expected = build_stamp(spec, spec_path, variant, run_root, records)
    run_dir.mkdir(parents=True, exist_ok=True)
    source_input = Path(str(expected["source_input"]["path"]))  # type: ignore[index]
    run_input = Path(str(expected["run_input"]["path"]))  # type: ignore[index]
    output = Path(str(expected["output"]["path"]))  # type: ignore[index]
    if output.exists():
        raise ValueError(
            f"refusing to prepare over existing {variant} output {output}; check it or use a new run root"
        )
    if run_input.exists():
        if common.sha256_file(run_input) != expected["run_input"]["sha256"]:  # type: ignore[index]
            raise ValueError(f"stale {variant} run input differs from the gate specification")
    else:
        shutil.copyfile(source_input, run_input)
    stamp_path = run_dir / JOB_STAMP
    if stamp_path.exists():
        observed = json.loads(stamp_path.read_text())
        if observed.get("status") != "prepared" or immutable_stamp_fields(observed) != immutable_stamp_fields(expected):
            raise ValueError(f"stale or foreign {variant} gate stamp: {stamp_path}")
    else:
        atomic_json(stamp_path, expected)
    return stamp_path


def prepare_command(args: argparse.Namespace) -> int:
    spec_path = args.spec.resolve(strict=True)
    spec = load_spec(spec_path)
    records = validate_artifacts(
        spec,
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    for variant in ("full", "spglib"):
        stamp = prepare_variant(spec, spec_path, variant, args.run_root, records)
        print(f"prepared {variant}: {stamp}")
    print("Run the two stamped jobs with:")
    for variant in ("full", "spglib"):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "run",
            "--variant",
            variant,
            "--spec",
            str(spec_path),
            "--run-root",
            str(args.run_root.resolve()),
            "--cp2k",
            str(args.cp2k.resolve()),
            "--cp2k-source",
            str(args.cp2k_source.resolve()),
            "--save-tblite",
            str(args.save_tblite.resolve()),
            "--save-tblite-source",
            str(args.save_tblite_source.resolve()),
            "--campaign-manifest",
            str(args.campaign_manifest.resolve()),
        ]
        print(shlex.join(command))
    print("Then evaluate the gate with:")
    print(
        shlex.join(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "check",
                "--spec",
                str(spec_path),
                "--run-root",
                str(args.run_root.resolve()),
                "--campaign-manifest",
                str(args.campaign_manifest.resolve()),
            ]
        )
    )
    return 0


def run_command(args: argparse.Namespace) -> int:
    spec_path = args.spec.resolve(strict=True)
    spec = load_spec(spec_path)
    records = validate_artifacts(
        spec,
        cp2k=args.cp2k,
        cp2k_source=args.cp2k_source,
        save_tblite=args.save_tblite,
        save_tblite_source=args.save_tblite_source,
        campaign_manifest=args.campaign_manifest,
    )
    stamp_path = prepare_variant(spec, spec_path, args.variant, args.run_root, records)
    stamp = json.loads(stamp_path.read_text())
    output = Path(str(stamp["output"]["path"]))
    if output.exists():
        raise ValueError(f"refusing to overwrite existing gate output: {output}")
    stamp["status"] = "running"
    stamp["started_at_utc"] = now_utc()
    atomic_json(stamp_path, stamp)

    run_dir = stamp_path.parent
    launcher_log = run_dir / LAUNCHER_LOG
    with launcher_log.open("w") as handle:
        process = subprocess.run(
            list(stamp["command"]),
            cwd=run_dir,
            env=common.thread_environment(1),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    stamp["ended_at_utc"] = now_utc()
    stamp["exit_code"] = process.returncode
    stamp["status"] = "completed" if process.returncode == 0 and output.is_file() else "failed"
    stamp["launcher_log"] = {
        "path": str(launcher_log),
        "sha256": common.sha256_file(launcher_log),
    }
    if output.is_file():
        stamp["output"]["sha256"] = common.sha256_file(output)
    atomic_json(stamp_path, stamp)
    print(f"{args.variant}: {stamp['status']} (exit {process.returncode}); {output}")
    return process.returncode if process.returncode != 0 else (0 if output.is_file() else 1)


def parse_cp2k_output(
    path: Path,
    expected_atoms: int,
    expected_elements: list[str] | None = None,
    fallback_kpoint_count: int | None = None,
) -> dict[str, object]:
    text = path.read_text(errors="replace")
    if "PROGRAM ENDED" not in text:
        raise ValueError(f"CP2K did not end normally: {path}")
    fatal_patterns = (
        r"\*\*\*.*ABORT",
        r"SCF\s+(?:run\s+)?(?:did\s+)?NOT\s+converged",
    )
    for pattern in fatal_patterns:
        if re.search(pattern, text, flags=re.I):
            raise ValueError(f"fatal/nonconverged marker in {path}: {pattern}")

    energy: float | None = None
    lines = text.splitlines()
    for line in lines:
        if "ENERGY| Total FORCE_EVAL" in line:
            try:
                energy = float(line.split()[-1].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
    if energy is None:
        raise ValueError(f"total FORCE_EVAL energy missing from {path}")

    legacy_force_headers = [
        index
        for index, line in enumerate(lines)
        if re.match(r"^\s*ATOMIC FORCES in \[a\.u\.\]\s*$", line, flags=re.I)
    ]
    compact_force_headers = [
        index
        for index, line in enumerate(lines)
        if re.match(
            r"^\s*FORCES\|\s*Atomic forces \[hartree/bohr\]\s*$",
            line,
            flags=re.I,
        )
    ]
    legacy_force_row = re.compile(
        rf"^\s*(\d+)\s+(\d+)\s+([A-Za-z]{{1,3}})\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
    )
    compact_force_row = re.compile(
        rf"^\s*FORCES\|\s*(\d+)\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})(?:\s+{FLOAT_PATTERN})?\s*$",
        flags=re.I,
    )
    forces: list[dict[str, object]] = []
    saw_force_sum = False
    if compact_force_headers:
        for line in lines[compact_force_headers[-1] + 1 :]:
            if re.match(r"^\s*FORCES\|\s*Sum\b", line, flags=re.I):
                saw_force_sum = True
                break
            match = compact_force_row.match(line)
            if not match:
                continue
            atom = int(match.group(1))
            element = expected_elements[atom - 1] if expected_elements else None
            forces.append(
                {
                    "atom": atom,
                    "kind": None,
                    "element": element,
                    "vector_hartree_per_bohr": [
                        float(match.group(index).replace("D", "E").replace("d", "e"))
                        for index in (2, 3, 4)
                    ],
                }
            )
    elif legacy_force_headers:
        for line in lines[legacy_force_headers[-1] + 1 :]:
            if re.match(r"^\s*SUM OF ATOMIC FORCES", line, flags=re.I):
                saw_force_sum = True
                break
            match = legacy_force_row.match(line)
            if not match:
                continue
            forces.append(
                {
                    "atom": int(match.group(1)),
                    "kind": int(match.group(2)),
                    "element": match.group(3).capitalize(),
                    "vector_hartree_per_bohr": [
                        float(match.group(index).replace("D", "E").replace("d", "e"))
                        for index in (4, 5, 6)
                    ],
                }
            )
    else:
        raise ValueError(f"atomic force block missing from {path}")
    if not saw_force_sum or len(forces) != expected_atoms:
        raise ValueError(
            f"incomplete atomic force block in {path}: {len(forces)}/{expected_atoms} atoms"
        )
    if [row["atom"] for row in forces] != list(range(1, expected_atoms + 1)):
        raise ValueError(f"noncanonical atom ordering in force block: {path}")
    if expected_elements is not None:
        if len(expected_elements) != expected_atoms:
            raise ValueError("expected force element list has the wrong atom count")
        observed_elements = [row["element"] for row in forces]
        if observed_elements != expected_elements:
            raise ValueError(f"force block element ordering differs from the gate input: {path}")

    stress_headers = [
        index
        for index, line in enumerate(lines)
        if re.match(r"^\s*STRESS\| Analytical stress tensor \[GPa\]\s*$", line, flags=re.I)
    ]
    if not stress_headers:
        raise ValueError(f"analytical GPa stress tensor missing from {path}")
    stress_row = re.compile(
        rf"^\s*STRESS\|\s+([xyz])\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})\s*$",
        flags=re.I,
    )
    stress: list[list[float]] = []
    axes: list[str] = []
    for line in lines[stress_headers[-1] + 1 :]:
        match = stress_row.match(line)
        if not match:
            continue
        axes.append(match.group(1).lower())
        stress.append(
            [
                float(match.group(index).replace("D", "E").replace("d", "e"))
                for index in (2, 3, 4)
            ]
        )
        if len(stress) == 3:
            break
    if axes != ["x", "y", "z"]:
        raise ValueError(f"incomplete analytical stress tensor in {path}")

    counts = [
        int(value)
        for value in re.findall(r"Number of Special K-points:\s*(\d+)", text, flags=re.I)
    ]
    if not counts:
        counts = [
            int(value)
            for value in re.findall(
                r"BRILLOUIN\| List of Kpoints[^\n]*?\s(\d+)\s*$",
                text,
                flags=re.I | re.M,
            )
        ]
    if not counts:
        counts = [
            int(value)
            for value in re.findall(
                r"KPOINTS\| Number of kpoints per group\s+(\d+)\s*$",
                text,
                flags=re.I | re.M,
            )
        ]
    if counts:
        kpoint_count = counts[-1]
        kpoint_count_source = "cp2k_output"
    elif fallback_kpoint_count is not None:
        kpoint_count = fallback_kpoint_count
        kpoint_count_source = "hashed_input_contract"
    else:
        raise ValueError(f"number of evaluated k-points missing from {path}")

    mesh_matches = re.findall(
        r"K-point Mesh:\s+(\d+)\s+(\d+)\s+(\d+)\s*$",
        text,
        flags=re.I | re.M,
    )
    kpoint_mesh = [int(value) for value in mesh_matches[-1]] if mesh_matches else None
    mesh_rows: list[int] = []
    if mesh_matches:
        mesh_header = max(
            index
            for index, line in enumerate(lines)
            if re.search(r"K-point Mesh:\s+\d+\s+\d+\s+\d+\s*$", line, flags=re.I)
        )
        mesh_row = re.compile(
            rf"^\s*(\d+)\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}"
            r"\s+\d+\s+\d+\s+\d+\s*$"
        )
        for line in lines[mesh_header + 1 :]:
            match = mesh_row.match(line)
            if match:
                mesh_rows.append(int(match.group(1)))
            elif mesh_rows:
                break

    return {
        "energy_hartree": energy,
        "forces": forces,
        "stress_gpa": stress,
        "kpoint_count": kpoint_count,
        "kpoint_count_source": kpoint_count_source,
        "kpoint_mesh": kpoint_mesh,
        "kpoint_mesh_rows": len(mesh_rows),
        "kpoint_mesh_row_indices": mesh_rows,
    }


def validate_completed_stamp(
    stamp_path: Path,
    spec: Mapping[str, object],
    spec_path: Path,
    variant: str,
) -> tuple[Path, dict[str, object]]:
    if not stamp_path.is_file():
        raise ValueError(f"missing {variant} gate stamp: {stamp_path}")
    stamp = json.loads(stamp_path.read_text())
    if (
        stamp.get("schema") != 1
        or stamp.get("gate_id") != spec.get("gate_id")
        or stamp.get("variant") != variant
        or stamp.get("status") != "completed"
        or stamp.get("exit_code") != 0
    ):
        raise ValueError(f"{variant} gate stamp does not record a completed successful job")
    if stamp.get("campaign_identity") != spec.get("campaign_identity"):
        raise ValueError(f"{variant} gate stamp has a foreign campaign identity")
    manifest_record = stamp.get("campaign_manifest")
    expected_manifest = repo_path(spec["campaign_manifest"]).resolve(strict=True)
    if (
        not isinstance(manifest_record, dict)
        or Path(str(manifest_record.get("path", ""))).resolve() != expected_manifest
        or manifest_record.get("campaign_id") != spec["campaign_identity"]["campaign_id"]
        or manifest_record.get("campaign_state") not in allowed_states(spec)
        or not re.fullmatch(r"[0-9a-f]{64}", str(manifest_record.get("file_sha256", "")))
    ):
        raise ValueError(f"{variant} gate stamp lacks the pinned campaign-manifest record")
    specification = stamp.get("specification")
    if not isinstance(specification, dict) or specification.get("sha256") != common.sha256_file(spec_path):
        raise ValueError(f"{variant} gate stamp has a stale specification")
    input_record = spec["inputs"][variant]  # type: ignore[index]
    run_input = stamp.get("run_input")
    if not isinstance(run_input, dict) or run_input.get("sha256") != input_record["sha256"]:
        raise ValueError(f"{variant} gate stamp has a stale input fingerprint")
    run_input_path = Path(str(run_input.get("path", ""))).resolve(strict=True)
    if common.sha256_file(run_input_path) != input_record["sha256"]:
        raise ValueError(f"{variant} run input changed after execution")

    campaign = spec["campaign_identity"]
    cp2k = stamp.get("cp2k")
    save = stamp.get("save_tblite")
    if not isinstance(cp2k, dict) or not isinstance(save, dict):
        raise ValueError(f"{variant} gate stamp lacks build records")
    exact_fields = (
        (cp2k.get("executable_sha256"), campaign["cp2k_executable_sha256"]),
        (cp2k.get("loaded_library_sha256"), campaign["cp2k_loaded_library_sha256"]),
        (cp2k.get("cmake_cache_sha256"), campaign["cp2k_cmake_cache_sha256"]),
        (cp2k.get("embedded_source_revision"), campaign["cp2k_embedded_source_revision"]),
        (cp2k.get("source_revision"), campaign["cp2k_source_revision"]),
        (save.get("executable_sha256"), campaign["save_tblite_executable_sha256"]),
        (save.get("static_library_sha256"), campaign["save_tblite_library_sha256"]),
        (save.get("cmake_cache_sha256"), campaign["save_tblite_cmake_cache_sha256"]),
        (save.get("source_revision"), campaign["save_tblite_source_revision"]),
    )
    if any(observed != expected for observed, expected in exact_fields):
        raise ValueError(f"{variant} gate stamp build records differ from the campaign")

    output = stamp.get("output")
    if not isinstance(output, dict) or not output.get("sha256"):
        raise ValueError(f"{variant} gate stamp has no output fingerprint")
    output_path = Path(str(output.get("path", ""))).resolve(strict=True)
    if common.sha256_file(output_path) != output["sha256"]:
        raise ValueError(f"{variant} output changed after execution")
    return output_path, stamp


def rms(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def check_command(args: argparse.Namespace) -> int:
    spec_path = args.spec.resolve(strict=True)
    spec = load_spec(spec_path)
    manifest_path = validate_manifest_choice(spec, args.campaign_manifest)
    identity, campaign_state = common.declared_campaign_identity(
        manifest_path,
        allowed_campaign_states=allowed_states(spec),
    )
    expected_atoms = int(spec["structure"]["atom_count"])  # type: ignore[index]
    expected_elements = _coord_elements(
        repo_path(spec["inputs"]["full"]["path"]).read_text()  # type: ignore[index]
    )
    parsed: dict[str, dict[str, object]] = {}
    stamps: dict[str, dict[str, object]] = {}
    for variant in ("full", "spglib"):
        output, stamp = validate_completed_stamp(
            args.run_root.resolve() / variant / JOB_STAMP,
            spec,
            spec_path,
            variant,
        )
        fallback = (
            int(spec["inputs"]["full"]["expected_kpoint_count"])  # type: ignore[index]
            if variant == "full"
            else None
        )
        parsed[variant] = parse_cp2k_output(
            output,
            expected_atoms,
            expected_elements,
            fallback_kpoint_count=fallback,
        )
        stamps[variant] = stamp

    full_forces = parsed["full"]["forces"]
    reduced_forces = parsed["spglib"]["forces"]
    assert isinstance(full_forces, list) and isinstance(reduced_forces, list)
    full_elements = [row["element"] for row in full_forces]
    reduced_elements = [row["element"] for row in reduced_forces]
    if full_elements != reduced_elements:
        raise ValueError("full and SPGLIB force blocks use different atom ordering")
    force_deltas = [
        float(left) - float(right)
        for full_row, reduced_row in zip(full_forces, reduced_forces)
        for left, right in zip(
            full_row["vector_hartree_per_bohr"],
            reduced_row["vector_hartree_per_bohr"],
        )
    ]
    full_stress = parsed["full"]["stress_gpa"]
    reduced_stress = parsed["spglib"]["stress_gpa"]
    assert isinstance(full_stress, list) and isinstance(reduced_stress, list)
    stress_deltas = [
        float(left) - float(right)
        for full_row, reduced_row in zip(full_stress, reduced_stress)
        for left, right in zip(full_row, reduced_row)
    ]
    metrics = {
        "energy_abs_hartree": abs(
            float(parsed["full"]["energy_hartree"])
            - float(parsed["spglib"]["energy_hartree"])
        ),
        "force_max_abs_hartree_per_bohr": max(abs(value) for value in force_deltas),
        "force_rms_hartree_per_bohr": rms(force_deltas),
        "stress_max_abs_gpa": max(abs(value) for value in stress_deltas),
        "stress_rms_gpa": rms(stress_deltas),
    }
    tolerances = spec["tolerances"]
    full_record = spec["inputs"]["full"]  # type: ignore[index]
    reduced_record = spec["inputs"]["spglib"]  # type: ignore[index]
    checks = {
        "energy": metrics["energy_abs_hartree"] <= float(tolerances["energy_max_abs_hartree"]),
        "forces": metrics["force_max_abs_hartree_per_bohr"]
        <= float(tolerances["force_max_abs_hartree_per_bohr"]),
        "stress": metrics["stress_max_abs_gpa"] <= float(tolerances["stress_max_abs_gpa"]),
        "full_kpoint_count": parsed["full"]["kpoint_count"]
        == int(full_record["expected_kpoint_count"]),
        "spglib_reduces_mesh": int(parsed["spglib"]["kpoint_count"])
        <= int(reduced_record["maximum_kpoint_count"])
        and int(parsed["spglib"]["kpoint_count"]) < int(parsed["full"]["kpoint_count"]),
        "spglib_reports_complete_mesh_mapping": parsed["spglib"]["kpoint_mesh"]
        == reduced_record["expected_full_mesh"]
        and int(parsed["spglib"]["kpoint_mesh_rows"])
        == int(reduced_record["expected_full_mesh_rows"])
        and parsed["spglib"]["kpoint_mesh_row_indices"]
        == list(range(1, int(reduced_record["expected_full_mesh_rows"]) + 1)),
    }
    passed = all(checks.values())
    jobs = {
        variant: {
            "input_sha256": stamps[variant]["run_input"]["sha256"],
            "output_path": stamps[variant]["output"]["path"],
            "output_sha256": stamps[variant]["output"]["sha256"],
            "kpoint_count": parsed[variant]["kpoint_count"],
            "kpoint_count_source": parsed[variant]["kpoint_count_source"],
            "kpoint_mesh": parsed[variant]["kpoint_mesh"],
            "kpoint_mesh_rows": parsed[variant]["kpoint_mesh_rows"],
            "energy_hartree": parsed[variant]["energy_hartree"],
            "forces": parsed[variant]["forces"],
            "stress_gpa": parsed[variant]["stress_gpa"],
        }
        for variant in ("full", "spglib")
    }
    result = {
        "schema": 1,
        "gate_id": spec["gate_id"],
        "checked_at_utc": now_utc(),
        "passed": passed,
        "campaign_state_at_check": campaign_state,
        "campaign_identity": identity,
        "campaign_manifest": {
            "path": str(manifest_path),
            "sha256": common.sha256_file(manifest_path),
        },
        "specification": {
            "path": str(spec_path),
            "sha256": common.sha256_file(spec_path),
        },
        "tolerances": tolerances,
        "metrics": metrics,
        "checks": checks,
        "jobs": jobs,
    }
    result_path = args.result or (args.run_root.resolve() / "gate_result.json")
    atomic_json(result_path, result)
    print(json.dumps({"passed": passed, "checks": checks, "metrics": metrics}, indent=2))
    print(f"result: {result_path}")
    return 0 if passed else 1


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


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="verify artifacts and stamp two clean run dirs")
    prepare.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    prepare.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    add_artifact_arguments(prepare)
    prepare.set_defaults(func=prepare_command)

    run = subparsers.add_parser("run", help="run one exactly stamped CP2K variant")
    run.add_argument("--variant", choices=("full", "spglib"), required=True)
    run.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    run.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    add_artifact_arguments(run)
    run.set_defaults(func=run_command)

    check = subparsers.add_parser("check", help="parse both outputs and apply the frozen gate")
    check.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    check.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    check.add_argument("--result", type=Path)
    check.add_argument(
        "--campaign-manifest",
        type=Path,
        default=common.DEFAULT_CAMPAIGN_MANIFEST,
    )
    check.set_defaults(func=check_command)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        return int(args.func(args))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
