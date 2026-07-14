#!/usr/bin/env python3
"""Prepare, run, and collect native-Bloch X23b k-point cell optimizations."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Mapping

import x23b_common as common
import x23b_pipeline as pipeline


ROOT = Path(__file__).resolve().parents[1]
METHODS = common.METHODS
HARTREE_TO_KJMOL = 2625.499638
VARIANT = "k222_cellopt_keep_angles"
EXPERIMENTAL_VARIANT = "k222_cellopt_keep_angles_from_experimental_reference"
SOURCE_POLICIES = ("gamma_cellopt_restart", "experimental_reference")
LINEAGE_BY_POLICY = {
    "gamma_cellopt_restart": "gamma_cellopt_restart->k222_cellopt_input",
    "experimental_reference": (
        "frozen_x23_reference_structure->experimental_k222_preflight->k222_cellopt_input"
    ),
}
MANIFEST_SCHEMA = 2
MANIFEST_FIELDS = (
    "schema",
    "method",
    "system",
    "variant",
    "source_policy",
    "source_kind",
    "lineage",
    "source_path",
    "source_sha256",
    "structure_path",
    "structure_sha256",
    "structure_source",
    "preflight_input",
    "preflight_input_sha256",
    "preflight_output",
    "preflight_output_sha256",
    "preflight_stamp",
    "preflight_stamp_sha256",
    "source_restart",
    "input",
    "input_sha256",
    "run_dir",
)
FIELDS = (
    "method",
    "system",
    "variant",
    "mesh",
    "source_variant",
    "source",
    "source_policy",
    "source_path",
    "source_sha256",
    "structure_path",
    "structure_sha256",
    "structure_source",
    "lineage",
    "preflight_input",
    "preflight_input_sha256",
    "preflight_output",
    "preflight_output_sha256",
    "preflight_stamp",
    "preflight_stamp_sha256",
    "returncode",
    "program_ended",
    "opt_completed",
    "max_iter_reached",
    "last_step",
    "energy_hartree",
    "gas_energy_hartree",
    "lattice_energy_kJmol",
    "x23b_ref_lattice_energy_kJmol",
    "error_kJmol",
    "volume_A3",
    "x23b_same_cell_ref_volume_A3",
    "volume_error_percent",
    "last_pressure_bar",
    "last_max_step",
    "last_rms_step",
    "last_max_gradient",
    "last_rms_gradient",
    "source_restart",
    "run_dir",
    "output",
)

CELL_OPT_LIMITS = {
    "pressure": 100.0,
    "max_step": 0.003,
    "rms_step": 0.0015,
    "max_gradient": 0.00045,
    "rms_gradient": 0.0003,
}

CELL_OPT_PATTERNS = {
    "step": re.compile(r"^\s*OPT\| Step number\s+(\d+)\s*$"),
    "pressure": re.compile(
        r"^\s*OPT\| (?:Internal pressure|Pressure deviation) \[?bar\]?\s+([-+0-9.Ee]+)\s*$"
    ),
    "max_step": re.compile(r"^\s*OPT\| Maximum step size\s+([-+0-9.Ee]+)\s*$"),
    "rms_step": re.compile(r"^\s*OPT\| RMS step size\s+([-+0-9.Ee]+)\s*$"),
    "max_gradient": re.compile(r"^\s*OPT\| Maximum gradient\s+([-+0-9.Ee]+)\s*$"),
    "rms_gradient": re.compile(r"^\s*OPT\| RMS gradient\s+([-+0-9.Ee]+)\s*$"),
}


def manifest_path(output_root: Path) -> Path:
    return output_root.resolve().parent / "x23b_k222_cellopt_manifest.csv"


def variant_for_policy(policy: str) -> str:
    if policy == "gamma_cellopt_restart":
        return VARIANT
    if policy == "experimental_reference":
        return EXPERIMENTAL_VARIANT
    raise ValueError(f"unsupported k222 source policy: {policy}")


def experimental_reference_paths(system: str) -> tuple[Path, Path]:
    return (
        ROOT / "inputs" / "crystal_sp" / "k222" / "GXTB" / f"{system}_GXTB_k222_sp.inp",
        ROOT / "structures" / "cif" / f"{system}.cif",
    )


def _normalized_manifest_row(row: Mapping[str, str]) -> dict[str, str]:
    """Normalize legacy Gamma rows without weakening the GXTB V1 contract."""

    normalized = {field: str(row.get(field, "") or "") for field in MANIFEST_FIELDS}
    normalized["schema"] = normalized["schema"] or "1"
    if normalized["schema"] == "1":
        normalized["variant"] = normalized["variant"] or VARIANT
        normalized["source_policy"] = normalized["source_policy"] or "gamma_cellopt_restart"
        normalized["source_kind"] = normalized["source_kind"] or "gamma_cellopt_restart"
        normalized["source_path"] = normalized["source_path"] or normalized["source_restart"]
        normalized["lineage"] = normalized["lineage"] or LINEAGE_BY_POLICY[
            "gamma_cellopt_restart"
        ]
    return normalized


def _manifest_metadata(system: str) -> dict[str, object]:
    try:
        return next(row for row in systems() if str(row["id"]) == system)
    except StopIteration as exc:
        raise ValueError(f"unknown X23b system in manifest: {system}") from exc


def validate_manifest_row(row: Mapping[str, str], output_root: Path) -> dict[str, str]:
    normalized = _normalized_manifest_row(row)
    method, system = normalized["method"], normalized["system"]
    if method not in METHODS or not system:
        raise ValueError(f"invalid k222 manifest identity: {method}/{system}")
    _manifest_metadata(system)
    try:
        schema = int(normalized["schema"])
    except ValueError as exc:
        raise ValueError(f"invalid k222 manifest schema for {method}/{system}") from exc
    if schema not in (1, MANIFEST_SCHEMA):
        raise ValueError(f"unsupported k222 manifest schema {schema} for {method}/{system}")

    input_path = Path(normalized["input"]).resolve()
    run_dir = Path(normalized["run_dir"]).resolve()
    if input_path.parent != run_dir:
        raise ValueError(f"manifest input/run_dir mismatch for {method}/{system}")
    try:
        input_path.relative_to(output_root.resolve())
    except ValueError as exc:
        raise ValueError(f"manifest input escapes output root: {input_path}") from exc
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    # Schema 1 is retained solely for the frozen GFN1/GFN2 and historical
    # Gamma manifests.  A GXTB experimental-reference row must always carry
    # the complete schema-2 lineage and hashes.
    if schema == 1:
        if method == "GXTB":
            raise ValueError(f"policyless legacy GXTB manifest row rejected for {system}")
        return normalized

    expected_variant = variant_for_policy(normalized["source_policy"])
    if normalized["variant"] != expected_variant:
        raise ValueError(f"source policy/variant mismatch for {method}/{system}")
    if normalized["source_kind"] != normalized["source_policy"]:
        raise ValueError(f"source kind/policy mismatch for {method}/{system}")
    if normalized["lineage"] != LINEAGE_BY_POLICY[normalized["source_policy"]]:
        raise ValueError(f"noncanonical source lineage for {method}/{system}")
    expected_run_dir = output_root.resolve() / method / system / expected_variant
    if run_dir != expected_run_dir:
        raise ValueError(f"noncanonical k222 run directory for {method}/{system}: {run_dir}")
    if not normalized["input_sha256"] or common.sha256_file(input_path) != normalized["input_sha256"]:
        raise ValueError(f"generated input fingerprint differs for {method}/{system}")
    source_path = Path(normalized["source_path"]).resolve()
    if not source_path.is_file() or common.sha256_file(source_path) != normalized["source_sha256"]:
        raise ValueError(f"source artifact fingerprint differs for {method}/{system}")

    if normalized["source_policy"] == "gamma_cellopt_restart":
        if normalized["source_restart"] != normalized["source_path"]:
            raise ValueError(f"Gamma source path/restart mismatch for {method}/{system}")
        if any(
            normalized[field]
            for field in (
                "structure_path",
                "structure_sha256",
                "preflight_input",
                "preflight_input_sha256",
                "preflight_output",
                "preflight_output_sha256",
                "preflight_stamp",
                "preflight_stamp_sha256",
            )
        ):
            raise ValueError(f"Gamma row contains experimental lineage for {method}/{system}")
        if normalized["structure_source"]:
            raise ValueError(f"Gamma row contains a structure lineage label for {method}/{system}")
        return normalized

    if method != "GXTB" or normalized["source_restart"]:
        raise ValueError(f"experimental_reference is GXTB-only and is not a restart: {method}/{system}")
    expected_source, expected_structure = experimental_reference_paths(system)
    if source_path != expected_source.resolve(strict=True):
        raise ValueError(f"noncanonical experimental reference input for {system}: {source_path}")
    structure_path = Path(normalized["structure_path"]).resolve()
    if structure_path != expected_structure.resolve(strict=True):
        raise ValueError(f"noncanonical experimental reference structure for {system}: {structure_path}")
    if common.sha256_file(structure_path) != normalized["structure_sha256"]:
        raise ValueError(f"experimental structure fingerprint differs for {system}")
    metadata = _manifest_metadata(system)
    if normalized["structure_source"] != str(metadata["structure_source"]):
        raise ValueError(f"experimental structure lineage label differs for {system}")
    validate_experimental_reference_source(system, source_path, structure_path)
    for path_field, hash_field in (
        ("preflight_input", "preflight_input_sha256"),
        ("preflight_output", "preflight_output_sha256"),
        ("preflight_stamp", "preflight_stamp_sha256"),
    ):
        artifact = Path(normalized[path_field]).resolve()
        if not artifact.is_file() or common.sha256_file(artifact) != normalized[hash_field]:
            raise ValueError(f"{path_field} fingerprint differs for GXTB/{system}")
    return normalized


def manifest_owned_inputs(
    output_root: Path,
    method: str | None = None,
    selected_systems: set[str] | None = None,
) -> list[Path]:
    """Select only primary inputs recorded by the k222 manifest.

    Continuation and BFGS-polish inputs live below the same directories and
    must never be interpreted as additional benchmark cases.
    """

    output_root = output_root.resolve()
    path = manifest_path(output_root)
    if not path.is_file():
        raise FileNotFoundError(path)
    methods = (method,) if method else common.PUBLISHED_METHODS
    selected: dict[tuple[str, str], Path] = {}
    with path.open(newline="") as handle:
        for raw_row in csv.DictReader(handle):
            row = _normalized_manifest_row(raw_row)
            row_method, system = row["method"], row["system"]
            if row_method not in methods:
                continue
            if selected_systems is not None and system not in selected_systems:
                continue
            key = (row_method, system)
            if key in selected:
                raise ValueError(f"duplicate k222 manifest row for {row_method}/{system}")
            row = validate_manifest_row(raw_row, output_root)
            input_path = Path(row["input"]).resolve()
            selected[key] = input_path
    return [selected[key] for key in sorted(selected)]


def manifest_record_for_input(output_root: Path, input_path: Path) -> dict[str, str]:
    input_path = input_path.resolve()
    matches: list[dict[str, str]] = []
    with manifest_path(output_root).open(newline="") as handle:
        for row in csv.DictReader(handle):
            if Path(row["input"]).resolve() == input_path:
                matches.append(validate_manifest_row(row, output_root.resolve()))
    if len(matches) != 1:
        raise ValueError(f"expected one k222 manifest row for {input_path}, found {len(matches)}")
    return matches[0]


def stamp_context(row: Mapping[str, str]) -> tuple[dict[str, object] | None, dict[str, Path] | None]:
    normalized = _normalized_manifest_row(row)
    if int(normalized["schema"]) < MANIFEST_SCHEMA:
        return None, None
    protocol = {
        "source_policy": normalized["source_policy"],
        "variant": normalized["variant"],
        "lineage": normalized["lineage"],
    }
    if normalized["source_policy"] == "experimental_reference":
        sources = {
            "reference_input": Path(normalized["source_path"]),
            "reference_structure": Path(normalized["structure_path"]),
            "preflight_input": Path(normalized["preflight_input"]),
            "preflight_output": Path(normalized["preflight_output"]),
            "preflight_stamp": Path(normalized["preflight_stamp"]),
        }
    else:
        sources = {"gamma_restart": Path(normalized["source_restart"])}
    return protocol, sources


def update_provenance(args: argparse.Namespace) -> Path:
    workflow_paths: dict[str, object] = {"k222_cellopt_root": args.output_root}
    path = manifest_path(args.output_root)
    if path.is_file():
        with path.open(newline="") as handle:
            policies = {
                _normalized_manifest_row(row)["source_policy"]
                for row in csv.DictReader(handle)
                if row.get("method") == "GXTB"
            }
        if len(policies) > 1:
            raise ValueError("GXTB k222 manifest mixes source policies")
        if policies:
            workflow_paths["k222_source_policy"] = next(iter(policies))
    provenance = ROOT / "data" / common.GXTB_PROVENANCE_NAME
    if provenance.is_file():
        payload = json.loads(provenance.read_text())
        frozen_paths = payload.get("workflow_paths", {})
        if isinstance(frozen_paths, dict):
            frozen_root = frozen_paths.get("k222_cellopt_root")
            if frozen_root and Path(str(frozen_root)).resolve() != args.output_root.resolve():
                raise ValueError(f"k222 cellopt root is already frozen as {frozen_root}")
            frozen_policy = frozen_paths.get("k222_source_policy")
            requested_policy = workflow_paths.get("k222_source_policy")
            if frozen_policy and requested_policy and frozen_policy != requested_policy:
                raise ValueError(f"k222 source policy is already frozen as {frozen_policy}")
    return common.update_gxtb_provenance(
        ROOT,
        cp2k=getattr(args, "cp2k", None),
        cp2k_source=getattr(args, "cp2k_source", None),
        save_tblite=getattr(args, "save_tblite", None),
        save_tblite_source=getattr(args, "save_tblite_source", None),
        campaign_manifest=getattr(args, "campaign_manifest", None),
        workflow_paths=workflow_paths,
    )


def systems() -> list[dict[str, object]]:
    metadata = json.loads((ROOT / "data" / "metadata.json").read_text())
    return list(metadata["systems"])


def section_bounds(lines: list[str], section: str) -> tuple[int, int]:
    target = f"&{section.upper()}"
    for start, line in enumerate(lines):
        parts = line.strip().upper().split(maxsplit=1)
        if not parts or parts[0] != target:
            continue
        depth = 0
        for end in range(start, len(lines)):
            stripped = lines[end].strip().upper()
            if stripped.startswith("&END"):
                depth -= 1
            elif stripped.startswith("&"):
                depth += 1
            if depth == 0:
                return start, end
        break
    raise ValueError(f"section {section} not found or not closed")


def _coord_elements(text: str) -> list[str]:
    lines = text.splitlines()
    start, end = section_bounds(lines, "COORD")
    elements: list[str] = []
    for line in lines[start + 1 : end]:
        stripped = line.strip()
        if not stripped or stripped.upper() == "SCALED" or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Za-z]{1,3}", fields[0]):
            elements.append(fields[0].capitalize())
    return elements


def validate_experimental_reference_source(system: str, source: Path, structure: Path) -> None:
    expected_source, expected_structure = experimental_reference_paths(system)
    source = source.resolve(strict=True)
    structure = structure.resolve(strict=True)
    if source != expected_source.resolve(strict=True) or structure != expected_structure.resolve(strict=True):
        raise ValueError(f"experimental reference source for {system} is not the frozen local pair")
    text = source.read_text()
    common.validate_method_input(text, "GXTB")
    required = {
        "RUN_TYPE ENERGY": r"^\s*RUN_TYPE\s+ENERGY\s*$",
        "analytical stress": r"^\s*STRESS_TENSOR\s+ANALYTICAL\s*$",
        "shifted k222": (
            r"^\s*SCHEME\s+MACDONALD\s+2\s+2\s+2\s+"
            r"0\.25\s+0\.25\s+0\.25\s*$"
        ),
        "SPGLIB symmetry": r"^\s*SYMMETRY\s+T\s*$",
        "reduced grid": r"^\s*FULL_GRID\s+F\s*$",
        "SPGLIB backend": r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$",
        "SPGLIB reduction": r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$",
    }
    for description, pattern in required.items():
        if len(re.findall(pattern, text, flags=re.I | re.M)) != 1:
            raise ValueError(f"frozen experimental reference input lacks unique {description}: {source}")
    if re.search(r"^\s*&MOTION\b", text, flags=re.I | re.M):
        raise ValueError(f"frozen experimental reference must be a single point: {source}")
    metadata = _manifest_metadata(system)
    if len(_coord_elements(text)) != int(metadata["n_atoms_crystal"]):
        raise ValueError(f"experimental reference atom count differs for {system}")

    # The generated input is accepted only if it is still the deterministic
    # rendering of the frozen local structure, so the two hashes form real
    # lineage rather than two unrelated labels in a manifest.
    pipeline_system = next(row for row in pipeline.SYSTEMS if str(row["id"]) == system)
    geometry = pipeline.parse_cif(structure)
    expected_text = pipeline.crystal_input(
        pipeline_system,
        geometry,
        "GXTB",
        pipeline.MESHES[2],
        "ENERGY",
    )
    if text != expected_text:
        raise ValueError(f"frozen k222 input is not reproducible from {structure}")


def experimental_reference_to_k222_input(source: Path, project: str, system: str) -> str:
    expected_source, structure = experimental_reference_paths(system)
    if source.resolve(strict=True) != expected_source.resolve(strict=True):
        raise ValueError(f"noncanonical experimental source for {system}: {source}")
    validate_experimental_reference_source(system, source, structure)
    lines = source.read_text().splitlines()
    replacements = {
        "PROJECT": (re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I), f'"{project}"'),
        "RUN_TYPE": (re.compile(r"^(\s*RUN_TYPE\s+).*$", re.I), "CELL_OPT"),
    }
    for name, (pattern, value) in replacements.items():
        matches = [index for index, line in enumerate(lines) if pattern.match(line)]
        if len(matches) != 1:
            raise ValueError(f"expected exactly one {name} in {source}")
        index = matches[0]
        match = pattern.match(lines[index])
        assert match is not None
        lines[index] = f"{match.group(1)}{value}"
    lines += [
        "",
        "&MOTION",
        "  &CELL_OPT",
        "    OPTIMIZER CG",
        "    MAX_ITER 500",
        "    EXTERNAL_PRESSURE [bar] 0.0",
        "    KEEP_ANGLES T",
        "    &CG",
        "      &LINE_SEARCH",
        "        TYPE 2PNT",
        "      &END LINE_SEARCH",
        "    &END CG",
        "  &END CELL_OPT",
        "&END MOTION",
    ]
    text = "\n".join(lines) + "\n"
    common.validate_method_input(text, "GXTB")
    return text


def restart_to_k222_input(source: Path, project: str, method: str) -> str:
    lines = source.read_text().splitlines()
    try:
        global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    except StopIteration as exc:
        raise ValueError(f"GLOBAL section missing in {source}") from exc
    lines = lines[global_start:]

    motion_start, motion_end = section_bounds(lines, "MOTION")
    del lines[motion_start : motion_end + 1]

    project_pattern = re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I)
    for index, line in enumerate(lines):
        match = project_pattern.match(line)
        if match:
            lines[index] = f'{match.group(1)}"{project}"'
            break
    else:
        raise ValueError(f"PROJECT keyword missing in {source}")

    dft_start, dft_end = section_bounds(lines, "DFT")
    if any("&KPOINTS" in line.upper() for line in lines[dft_start:dft_end]):
        raise ValueError(f"source restart already contains KPOINTS: {source}")
    indent = re.match(r"\s*", lines[dft_start]).group(0) + "  "
    kpoints = [
        f"{indent}&KPOINTS",
        f"{indent}  SCHEME MACDONALD 2 2 2 0.25 0.25 0.25",
        f"{indent}  EPS_SYMMETRY 1.0E-8",
        f"{indent}  SYMMETRY T",
        f"{indent}  FULL_GRID F",
        f"{indent}  SYMMETRY_BACKEND SPGLIB",
        f"{indent}  SYMMETRY_REDUCTION_METHOD SPGLIB",
    ]
    kpoints += [f"{indent}&END KPOINTS"]
    lines[dft_end:dft_end] = kpoints

    cell_start, cell_end = section_bounds(lines, "CELL")
    if not any("CANONICALIZE" in line.upper() for line in lines[cell_start:cell_end]):
        indent = re.match(r"\s*", lines[cell_start]).group(0) + "  "
        lines.insert(cell_start + 1, f"{indent}CANONICALIZE TRUE")

    lines += [
        "",
        "&MOTION",
        "  &CELL_OPT",
        "    OPTIMIZER CG",
        "    MAX_ITER 500",
        "    EXTERNAL_PRESSURE [bar] 0.0",
        "    KEEP_ANGLES T",
        "    &CG",
        "      &LINE_SEARCH",
        "        TYPE 2PNT",
        "      &END LINE_SEARCH",
        "    &END CG",
        "  &END CELL_OPT",
        "&END MOTION",
    ]
    text = "\n".join(lines) + "\n"
    common.validate_method_input(text, method)
    return text


def parse_overrides(values: list[str]) -> dict[tuple[str, str], Path]:
    overrides: dict[tuple[str, str], Path] = {}
    for value in values:
        key, separator, path = value.partition("=")
        if not separator or "/" not in key:
            raise ValueError(f"invalid override {value!r}; expected METHOD/system=/path/restart")
        method, system = key.split("/", 1)
        identity = (method, system)
        if identity in overrides:
            raise ValueError(f"duplicate override for {method}/{system}")
        overrides[identity] = Path(path).resolve()
    return overrides


def find_gamma_restart(source_root: Path, method: str, system: str) -> Path:
    directories = (
        source_root / "runs" / method / system / "gamma_cellopt_keep_angles",
        source_root / "runs" / "cellopt_gamma" / method / f"{system}_{method}_gamma_cellopt",
    )
    for directory in directories:
        restart = common.final_restart(directory)
        if restart is not None:
            return restart
    raise ValueError("final Gamma restart not found in: " + ", ".join(str(path) for path in directories))


def validate_gxtb_gamma_source(source_root: Path, system: str, restart: Path) -> None:
    stem = f"{system}_GXTB_gamma_cellopt"
    input_path = source_root / "inputs" / "cellopt_gamma" / "GXTB" / f"{stem}.inp"
    output = restart.parent / f"{stem}.out"
    if not output.is_file():
        raise ValueError(f"GXTB Gamma output not found for restart {restart}: {output}")
    valid, reason = common.recorded_job_stamp_matches(
        restart.parent,
        input_path,
        "GXTB",
        "x23b_cellopt_gamma",
        output,
        campaign_identity=common.load_campaign_identity(source_root),
    )
    if not valid:
        raise ValueError(f"untrusted GXTB Gamma source for {system}: {reason}")


def prepare(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / ".x23b_k222_cellopt_manifest.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _prepare(args, output_root)


def _prepare(args: argparse.Namespace, output_root: Path) -> None:
    requested_policy = getattr(args, "source_policy", None)
    if args.method == "GXTB" and requested_policy is None:
        raise ValueError(
            "GXTB k222 preparation requires an explicit --source-policy; "
            "choose gamma_cellopt_restart or experimental_reference"
        )
    source_policy = requested_policy or "gamma_cellopt_restart"
    if source_policy not in SOURCE_POLICIES:
        raise ValueError(f"unsupported k222 source policy: {source_policy}")
    if source_policy == "experimental_reference" and args.method != "GXTB":
        raise ValueError("experimental_reference is restricted to --method GXTB")
    if getattr(args, "override", []) and args.method == "GXTB":
        raise ValueError("--override is forbidden for GXTB; use an explicit verified source policy")
    if source_policy == "gamma_cellopt_restart" and getattr(args, "gamma_root", None) is None:
        raise ValueError("gamma_cellopt_restart requires --gamma-root")
    if source_policy == "experimental_reference" and getattr(args, "preflight_root", None) is None:
        raise ValueError("experimental_reference requires --preflight-root")
    if args.clean and args.method != "GXTB":
        raise ValueError("--clean is restricted to --method GXTB; published method trees are immutable")
    if (
        args.clean
        and source_policy == "gamma_cellopt_restart"
        and (output_root / args.method).exists()
    ):
        shutil.rmtree(output_root / args.method)
    overrides = parse_overrides(args.override)
    selected_methods = (args.method,) if args.method else common.PUBLISHED_METHODS
    selected_systems = set(args.system) if args.system else None
    manifest = manifest_path(output_root)
    manifest_rows: dict[tuple[str, str], dict[str, str]] = {}
    if manifest.is_file() and args.method == "GXTB":
        with manifest.open(newline="") as handle:
            existing_policies = {
                _normalized_manifest_row(row)["source_policy"]
                for row in csv.DictReader(handle)
                if row.get("method") == "GXTB"
            }
        if existing_policies and existing_policies != {source_policy}:
            raise ValueError(
                "k222 output root is already owned by a different GXTB source policy; "
                "use a new output root"
            )
    if manifest.is_file() and not args.clean:
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle):
                manifest_rows[(row["method"], row["system"])] = row
    elif manifest.is_file() and args.clean:
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["method"] != args.method or (
                    source_policy == "experimental_reference"
                    and selected_systems is not None
                    and row["system"] not in selected_systems
                ):
                    manifest_rows[(row["method"], row["system"])] = row
    prepared = 0
    for system_data in systems():
        system = str(system_data["id"])
        if selected_systems is not None and system not in selected_systems:
            continue
        for method in selected_methods:
            structure: Path | None = None
            preflight_input: Path | None = None
            preflight_output: Path | None = None
            preflight_stamp: Path | None = None
            if source_policy == "experimental_reference":
                source, structure = experimental_reference_paths(system)
                validate_experimental_reference_source(system, source, structure)
                # Import lazily because the preflight script reuses the
                # canonical-source helpers in this module.
                import x23b_experimental_k222_preflight as preflight

                preflight_record = preflight.validate_completed_case(
                    args.preflight_root.resolve(),
                    system,
                    common.load_campaign_identity(ROOT),
                )
                preflight_input = Path(str(preflight_record["input"])).resolve(strict=True)
                preflight_output = Path(str(preflight_record["output"])).resolve(strict=True)
                preflight_stamp = (
                    preflight_output.parent / common.JOB_STAMP_NAME
                ).resolve(strict=True)
            else:
                source = overrides.get((method, system))
                if source is None:
                    source = find_gamma_restart(args.gamma_root.resolve(), method, system)
                    if method == "GXTB":
                        validate_gxtb_gamma_source(args.gamma_root.resolve(), system, source)
            if not source.is_file():
                raise FileNotFoundError(source)
            variant = variant_for_policy(source_policy)
            project = f"{system}_{method}_{variant}".replace("-", "_")
            run_dir = output_root / method / system / variant
            if args.clean and source_policy == "experimental_reference" and run_dir.exists():
                shutil.rmtree(run_dir)
            if not args.clean and (run_dir / "cp2k.out").exists():
                raise ValueError(f"refusing to prepare over an existing k222 output: {run_dir / 'cp2k.out'}")
            run_dir.mkdir(parents=True, exist_ok=True)
            input_path = run_dir / f"{project}.inp"
            if source_policy == "experimental_reference":
                input_text = experimental_reference_to_k222_input(source, project, system)
            else:
                input_text = restart_to_k222_input(source, project, method)
            input_path.write_text(input_text)
            metadata = _manifest_metadata(system)
            manifest_rows[(method, system)] = {
                "schema": str(MANIFEST_SCHEMA),
                "method": method,
                "system": system,
                "variant": variant,
                "source_policy": source_policy,
                "source_kind": source_policy,
                "lineage": LINEAGE_BY_POLICY[source_policy],
                "source_path": str(source.resolve()),
                "source_sha256": common.sha256_file(source),
                "structure_path": str(structure.resolve()) if structure is not None else "",
                "structure_sha256": common.sha256_file(structure) if structure is not None else "",
                "structure_source": (
                    str(metadata["structure_source"]) if structure is not None else ""
                ),
                "preflight_input": str(preflight_input) if preflight_input is not None else "",
                "preflight_input_sha256": (
                    common.sha256_file(preflight_input) if preflight_input is not None else ""
                ),
                "preflight_output": str(preflight_output) if preflight_output is not None else "",
                "preflight_output_sha256": (
                    common.sha256_file(preflight_output) if preflight_output is not None else ""
                ),
                "preflight_stamp": str(preflight_stamp) if preflight_stamp is not None else "",
                "preflight_stamp_sha256": (
                    common.sha256_file(preflight_stamp) if preflight_stamp is not None else ""
                ),
                "source_restart": str(source.resolve()) if source_policy == "gamma_cellopt_restart" else "",
                "input": str(input_path),
                "input_sha256": common.sha256_file(input_path),
                "run_dir": str(run_dir),
            }
            prepared += 1
    rows = [manifest_rows[key] for key in sorted(manifest_rows)]
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(_normalized_manifest_row(row) for row in rows)
    print(
        f"Prepared {prepared} native-Bloch k222 cell optimizations in {output_root}; "
        f"manifest contains {len(rows)} entries"
    )


def cp2k_completed(output: Path) -> bool:
    if not output.is_file():
        return False
    text = output.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and "GEOMETRY OPTIMIZATION COMPLETED" in text


def cp2k_terminal(output: Path) -> bool:
    if not output.is_file():
        return False
    text = output.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and (
        "GEOMETRY OPTIMIZATION COMPLETED" in text
        or "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text
    )


def run_one(
    input_path: Path,
    cp2k: Path,
    method: str,
    threads: int,
    force: bool,
    prune_transients: bool,
    campaign_identity: dict[str, object] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    output = run_dir / "cp2k.out"
    common.validate_method_input(input_path.read_text(), method)
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, common.BUSY_RETURN_CODE, "BUSY"
        stamp_matches, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            method,
            "x23b_k222_cellopt",
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
        if not force and output.exists() and method == "GXTB":
            recorded_matches, _ = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                method,
                "x23b_k222_cellopt",
                output,
                campaign_identity=campaign_identity,
                accepted_status_prefixes=("converged", "max_iter"),
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
            if not stamp_matches or not recorded_matches:
                return input_path, 1, "STALE_OUTPUT"
        if not force and cp2k_terminal(output):
            action = "SKIP_CONVERGED" if cp2k_completed(output) else "SKIP_MAX_ITER"
            if action == "SKIP_CONVERGED" and prune_transients:
                common.prune_gxtb_transients(run_dir, keep_final_restart=True)
            return input_path, 0, action
        if force:
            for path in run_dir.iterdir():
                if path != input_path:
                    if path.is_dir():
                        shutil.rmtree(path)
                    else:
                        path.unlink()
        process = subprocess.run(
            [str(cp2k), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            env=common.thread_environment(threads),
            check=False,
        )
        (run_dir / "returncode.txt").write_text(f"{process.returncode}\n")
        if process.returncode != 0:
            action = "FAILED"
        elif cp2k_completed(output):
            action = "CONVERGED"
            if prune_transients:
                common.prune_gxtb_transients(run_dir, keep_final_restart=True)
        elif cp2k_terminal(output):
            action = "MAX_ITER"
        else:
            action = "INCOMPLETE"
        if method == "GXTB":
            details: dict[str, object] = {
                "returncode": process.returncode,
                "action": action,
                "output": str(output),
            }
            if output.is_file():
                details["output_sha256"] = common.sha256_file(output)
            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                method,
                "x23b_k222_cellopt",
                action.lower(),
                details=details,
                campaign_identity=campaign_identity,
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
        return input_path, process.returncode, action


def run(args: argparse.Namespace) -> None:
    if args.force and args.method != "GXTB":
        raise ValueError("--force is restricted to --method GXTB")
    if args.prune_transients and args.method != "GXTB":
        raise ValueError("--prune-transients is restricted to --method GXTB")
    all_systems = sorted(str(row["id"]) for row in systems())
    if args.system and (args.start_system or args.end_system):
        raise ValueError("--system cannot be combined with --start-system or --end-system")
    if args.system:
        selected_systems = sorted(set(args.system))
    else:
        start = all_systems.index(args.start_system) if args.start_system else 0
        end = all_systems.index(args.end_system) + 1 if args.end_system else len(all_systems)
        if start >= end:
            raise ValueError("--start-system must not follow --end-system")
        selected_systems = all_systems[start:end]
    inputs = manifest_owned_inputs(args.output_root, args.method, set(selected_systems))
    expected = len(selected_systems) * (1 if args.method else len(common.PUBLISHED_METHODS))
    if len(inputs) != expected:
        selection = f" for {args.method}" if args.method else ""
        raise ValueError(f"expected {expected} prepared inputs{selection}, found {len(inputs)}")
    campaign_identity = None
    if args.method == "GXTB":
        common.require_gxtb_build_artifacts(
            cp2k=args.cp2k,
            cp2k_source=args.cp2k_source,
            save_tblite=args.save_tblite,
            save_tblite_source=args.save_tblite_source,
            campaign_manifest=args.campaign_manifest,
        )
        update_provenance(args)
        campaign_identity = common.load_campaign_identity(ROOT)
    failed = []
    max_iter = []
    contexts = {
        path: stamp_context(manifest_record_for_input(args.output_root, path))
        for path in inputs
    }
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_one,
                path,
                args.cp2k.resolve(),
                args.method or next(part for part in path.parts if part in METHODS),
                args.threads_per_job,
                args.force,
                args.prune_transients,
                campaign_identity,
                contexts[path][0],
                contexts[path][1],
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:14s} {relative} rc={returncode}", flush=True)
            if action in {"MAX_ITER", "SKIP_MAX_ITER"}:
                max_iter.append(relative)
            elif returncode != 0 or action in {"INCOMPLETE", "BUSY", "STALE_OUTPUT"}:
                failed.append(relative)
    if args.method == "GXTB":
        update_provenance(args)
    if max_iter:
        print(f"{len(max_iter)} optimization(s) require continue-maxiter", flush=True)
    if failed:
        raise SystemExit(f"{len(failed)} CP2K jobs failed")


def latest_numbered_restart(run_dir: Path) -> tuple[int, Path] | None:
    restarts: list[tuple[int, Path]] = []
    for path in run_dir.glob("*-1_*.restart"):
        match = re.search(r"-1_(\d+)\.restart$", path.name)
        if match:
            restarts.append((int(match.group(1)), path))
    return max(restarts) if restarts else None


def optimization_records(output: Path) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    current: dict[str, float | int] = {}
    for line in output.read_text(errors="ignore").splitlines():
        step_match = CELL_OPT_PATTERNS["step"].match(line)
        if step_match:
            if "step" in current:
                records.append(current)
            current = {"step": int(step_match.group(1))}
            continue
        if "step" not in current:
            continue
        for name, pattern in CELL_OPT_PATTERNS.items():
            if name == "step":
                continue
            match = pattern.match(line)
            if match:
                current[name] = float(match.group(1))
                break
    if "step" in current:
        records.append(current)
    return records


def best_polish_restart(run_dir: Path) -> tuple[Path, dict[str, float | int]]:
    candidates: list[tuple[float, float, int, Path, dict[str, float | int]]] = []
    for output in sorted(run_dir.glob("cp2k*.out")):
        for record in optimization_records(output):
            if any(name not in record for name in CELL_OPT_LIMITS):
                continue
            step = int(record["step"])
            restarts = list(run_dir.glob(f"*-1_{step}.restart"))
            if not restarts:
                continue
            restart = max(restarts, key=lambda path: path.stat().st_mtime_ns)
            ratios = [
                abs(float(record[name])) / limit if name == "pressure" else float(record[name]) / limit
                for name, limit in CELL_OPT_LIMITS.items()
            ]
            score = max(ratios)
            norm = sum(value * value for value in ratios)
            candidates.append((score, norm, -step, restart, record))
    if not candidates:
        raise ValueError(f"no complete optimization record with a matching restart in {run_dir}")
    _, _, _, restart, record = min(candidates, key=lambda item: item[:3])
    return restart, record


def set_cell_opt_keyword(lines: list[str], keyword: str, value: str) -> None:
    start, end = section_bounds(lines, "CELL_OPT")
    pattern = re.compile(rf"^(\s*{re.escape(keyword)}\s+).*$", re.I)
    for index in range(start + 1, end):
        match = pattern.match(lines[index])
        if match:
            lines[index] = f"{match.group(1)}{value}"
            return
    indent = re.match(r"\s*", lines[start]).group(0) + "  "
    lines.insert(end, f"{indent}{keyword} {value}")


def bfgs_polish_input(source: Path, project: str, max_iter: int, trust_radius: float) -> str:
    lines = source.read_text().splitlines()
    global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    lines = lines[global_start:]

    project_pattern = re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I)
    for index, line in enumerate(lines):
        match = project_pattern.match(line)
        if match:
            lines[index] = f'{match.group(1)}"{project}"'
            break
    else:
        raise ValueError(f"PROJECT keyword missing in {source}")

    cell_start, cell_end = section_bounds(lines, "CELL_OPT")
    for subsection in ("CG", "BFGS", "LBFGS"):
        local_lines = lines[cell_start : cell_end + 1]
        try:
            local_start, local_end = section_bounds(local_lines, subsection)
        except ValueError:
            continue
        del lines[cell_start + local_start : cell_start + local_end + 1]
        cell_start, cell_end = section_bounds(lines, "CELL_OPT")

    settings = {
        "OPTIMIZER": "BFGS",
        "MAX_ITER": str(max_iter),
        "STEP_START_VAL": "0",
        "MAX_DR": "0.003",
        "RMS_DR": "0.0015",
        "MAX_FORCE": "0.00045",
        "RMS_FORCE": "0.0003",
        "PRESSURE_TOLERANCE": "[bar] 100.0",
    }
    for keyword, value in settings.items():
        set_cell_opt_keyword(lines, keyword, value)

    cell_start, cell_end = section_bounds(lines, "CELL_OPT")
    indent = re.match(r"\s*", lines[cell_start]).group(0) + "  "
    lines[cell_end:cell_end] = [
        f"{indent}&BFGS",
        f"{indent}  TRUST_RADIUS [angstrom] {trust_radius:.12g}",
        f"{indent}&END BFGS",
    ]
    return "\n".join(lines) + "\n"


def continuation_input(source: Path, project: str, additional_steps: int) -> tuple[str, int]:
    text = source.read_text()
    step_match = re.search(r"^\s*STEP_START_VAL\s+(\d+)\s*$", text, flags=re.M)
    if step_match is None:
        raise ValueError(f"STEP_START_VAL missing in {source}")
    start_step = int(step_match.group(1))
    text, project_count = re.subn(
        r'^(\s*PROJECT(?:_NAME)?\s+).+$',
        rf'\1"{project}"',
        text,
        count=1,
        flags=re.I | re.M,
    )
    text, max_iter_count = re.subn(
        r"^(\s*MAX_ITER\s+)\d+\s*$",
        rf"\g<1>{start_step + additional_steps}",
        text,
        count=1,
        flags=re.I | re.M,
    )
    if project_count != 1 or max_iter_count != 1:
        raise ValueError(f"could not update continuation input from {source}")
    return text, start_step


def archive_path(path: Path, label: str) -> Path:
    candidate = path.with_name(f"{path.stem}.{label}{path.suffix}")
    index = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.stem}.{label}.{index}{path.suffix}")
        index += 1
    return candidate


def continue_one(
    input_path: Path,
    cp2k: Path,
    additional_steps: int,
    rounds: int,
    method: str,
    threads: int,
    prune_transients: bool,
    campaign_identity: dict[str, object] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    canonical_output = run_dir / "cp2k.out"
    common.validate_method_input(input_path.read_text(), method)
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, common.BUSY_RETURN_CODE, "BUSY"
        stamp_matches, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            method,
            "x23b_k222_cellopt",
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
        if canonical_output.exists() and method == "GXTB":
            recorded_matches, _ = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                method,
                "x23b_k222_cellopt",
                canonical_output,
                campaign_identity=campaign_identity,
                accepted_status_prefixes=("converged", "max_iter"),
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
            if not stamp_matches or not recorded_matches:
                return input_path, 1, "STALE_OUTPUT"
        if cp2k_completed(canonical_output):
            if prune_transients:
                common.prune_gxtb_transients(run_dir, keep_final_restart=True)
            return input_path, 0, "SKIP"
        canonical_text = canonical_output.read_text(errors="ignore") if canonical_output.exists() else ""
        if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" not in canonical_text:
            return input_path, 0, "WAIT"

        for round_index in range(1, rounds + 1):
            latest = latest_numbered_restart(run_dir)
            if latest is None:
                return input_path, 1, "NO_RESTART"
            step, restart = latest
            project = f"{input_path.stem}_continue_{step}"
            text, start_step = continuation_input(restart, project, additional_steps)
            if start_step != step:
                raise ValueError(f"restart step mismatch in {restart}: filename={step}, input={start_step}")
            continuation = run_dir / f"{project}.inp"
            continuation.write_text(text)
            output = run_dir / f"cp2k.continue_{step}.out"
            if output.exists():
                output.unlink()
            code = subprocess.run(
                [str(cp2k), "-i", continuation.name, "-o", output.name],
                cwd=run_dir,
                env=common.thread_environment(threads),
                check=False,
            ).returncode
            if code != 0:
                return input_path, code, f"CONTINUE_{round_index}_FAILED"
            result = output.read_text(errors="ignore")
            if "GEOMETRY OPTIMIZATION COMPLETED" in result and "PROGRAM ENDED" in result:
                if canonical_output.exists():
                    canonical_output.replace(archive_path(canonical_output, f"precontinue_{step}"))
                output.replace(canonical_output)
                (run_dir / "returncode.txt").write_text("0\n")
                if method == "GXTB":
                    common.write_job_stamp(
                        run_dir,
                        input_path,
                        cp2k,
                        method,
                        "x23b_k222_cellopt",
                        "converged_after_continuation",
                        details={
                            "continuation_input": str(continuation),
                            "continuation_input_sha256": common.sha256_file(continuation),
                            "source_restart": str(restart),
                            "source_restart_sha256": common.sha256_file(restart),
                            "round": round_index,
                            "output": str(canonical_output.resolve()),
                            "output_sha256": common.sha256_file(canonical_output),
                        },
                        campaign_identity=campaign_identity,
                        protocol_identity=protocol_identity,
                        source_artifacts=source_artifacts,
                    )
                if prune_transients:
                    common.prune_gxtb_transients(run_dir, keep_final_restart=True)
                return input_path, 0, f"CONTINUE_{round_index}"
            if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" not in result:
                return input_path, 1, f"CONTINUE_{round_index}_INCOMPLETE"
        return input_path, 1, "MAX_ITER"


def continue_maxiter(args: argparse.Namespace) -> None:
    if args.prune_transients and args.method != "GXTB":
        raise ValueError("--prune-transients is restricted to --method GXTB")
    wanted = set(args.system) if args.system else None
    inputs = manifest_owned_inputs(args.output_root, args.method, wanted)
    campaign_identity = None
    if args.method == "GXTB":
        common.require_gxtb_build_artifacts(
            cp2k=args.cp2k,
            cp2k_source=args.cp2k_source,
            save_tblite=args.save_tblite,
            save_tblite_source=args.save_tblite_source,
            campaign_manifest=args.campaign_manifest,
        )
        update_provenance(args)
        campaign_identity = common.load_campaign_identity(ROOT)
    failed = []
    contexts = {
        path: stamp_context(manifest_record_for_input(args.output_root, path))
        for path in inputs
    }
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                continue_one,
                path,
                args.cp2k.resolve(),
                args.additional_steps,
                args.rounds,
                args.method or next(part for part in path.parts if part in METHODS),
                args.threads_per_job,
                args.prune_transients,
                campaign_identity,
                contexts[path][0],
                contexts[path][1],
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:24s} {relative} rc={returncode}", flush=True)
            if returncode != 0:
                failed.append(relative)
    if args.method == "GXTB":
        update_provenance(args)
    if failed:
        raise SystemExit(f"{len(failed)} continuation jobs failed")


def polish_one(
    input_path: Path,
    cp2k: Path,
    max_iter: int,
    trust_radius: float,
    force: bool,
    method: str,
    threads: int,
    prune_transients: bool,
    campaign_identity: dict[str, object] | None = None,
    protocol_identity: Mapping[str, object] | None = None,
    source_artifacts: Mapping[str, Path] | None = None,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    canonical_output = run_dir / "cp2k.out"
    common.validate_method_input(input_path.read_text(), method)
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, common.BUSY_RETURN_CODE, "BUSY"
        stamp_matches, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            method,
            "x23b_k222_cellopt",
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
        if canonical_output.exists() and method == "GXTB":
            recorded_matches, _ = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                method,
                "x23b_k222_cellopt",
                canonical_output,
                campaign_identity=campaign_identity,
                accepted_status_prefixes=("converged", "max_iter"),
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
            if not stamp_matches or not recorded_matches:
                return input_path, 1, "STALE_OUTPUT"
        if cp2k_completed(canonical_output):
            if prune_transients:
                common.prune_gxtb_transients(run_dir, keep_final_restart=True)
            return input_path, 0, "SKIP_CONVERGED"
        if not cp2k_terminal(canonical_output):
            return input_path, 0, "WAIT"

        restart, record = best_polish_restart(run_dir)
        step = int(record["step"])
        project = f"{input_path.stem}_bfgs_polish_{step}"
        polish_dir = run_dir / f"bfgs_polish_from_{step}"
        if force and polish_dir.exists():
            shutil.rmtree(polish_dir)
        polish_dir.mkdir(parents=True, exist_ok=True)
        polish_input = polish_dir / f"{project}.inp"
        polish_output = polish_dir / "cp2k.out"
        polish_input.write_text(bfgs_polish_input(restart, project, max_iter, trust_radius))
        common.validate_method_input(polish_input.read_text(), method)

        polish_stamp_matches, _ = common.job_stamp_matches(
            polish_dir,
            polish_input,
            cp2k,
            method,
            "x23b_k222_bfgs_polish",
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
        polish_recorded_matches = True
        if method == "GXTB" and polish_output.exists():
            polish_recorded_matches, _ = common.recorded_job_stamp_matches(
                polish_dir,
                polish_input,
                method,
                "x23b_k222_bfgs_polish",
                polish_output,
                campaign_identity=campaign_identity,
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
        if (
            method == "GXTB"
            and polish_output.exists()
            and not force
            and (not polish_stamp_matches or not polish_recorded_matches)
        ):
            return input_path, 1, "POLISH_STALE_OUTPUT"
        if force or not cp2k_completed(polish_output):
            code = subprocess.run(
                [str(cp2k), "-i", polish_input.name, "-o", polish_output.name],
                cwd=polish_dir,
                env=common.thread_environment(threads),
                check=False,
            ).returncode
            if method == "GXTB":
                details: dict[str, object] = {"returncode": code, "output": str(polish_output)}
                if polish_output.is_file():
                    details["output_sha256"] = common.sha256_file(polish_output)
                common.write_job_stamp(
                    polish_dir,
                    polish_input,
                    cp2k,
                    method,
                    "x23b_k222_bfgs_polish",
                    "converged" if code == 0 and cp2k_completed(polish_output) else "failed",
                    details=details,
                    campaign_identity=campaign_identity,
                    protocol_identity=protocol_identity,
                    source_artifacts=source_artifacts,
                )
        else:
            code = 0
        if code != 0 or not cp2k_completed(polish_output):
            return input_path, code or 1, "POLISH_FAILED"

        final_restarts = list(polish_dir.glob("*-1.restart"))
        if not final_restarts:
            return input_path, 1, "POLISH_NO_RESTART"
        final_restart = max(final_restarts, key=lambda path: path.stat().st_mtime_ns)
        archived_output = archive_path(canonical_output, f"prebfgs_{step}")
        canonical_output.replace(archived_output)
        shutil.copyfile(polish_output, canonical_output)
        promoted_restart = run_dir / f"{input_path.stem}_bfgs_polished-1.restart"
        shutil.copyfile(final_restart, promoted_restart)
        promoted_restart.touch()
        (run_dir / "returncode.txt").write_text("0\n")
        provenance = {
            "source_restart": str(restart),
            "source_step": step,
            "source_metrics": record,
            "optimizer": "BFGS",
            "trust_radius_angstrom": trust_radius,
            "max_iter": max_iter,
            "polish_input": str(polish_input),
            "polish_output": str(polish_output),
            "promoted_restart": str(promoted_restart),
            "archived_cg_output": str(archived_output),
        }
        (run_dir / "bfgs_polish_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        if method == "GXTB":
            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                method,
                "x23b_k222_cellopt",
                "converged_after_bfgs_polish",
                details={
                    "polish_input": str(polish_input),
                    "polish_input_sha256": common.sha256_file(polish_input),
                    "source_restart": str(restart),
                    "source_restart_sha256": common.sha256_file(restart),
                    "output": str(canonical_output.resolve()),
                    "output_sha256": common.sha256_file(canonical_output),
                    "promoted_restart_sha256": common.sha256_file(promoted_restart),
                },
                campaign_identity=campaign_identity,
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
        if prune_transients:
            common.prune_gxtb_transients(run_dir, keep_final_restart=True)
        return input_path, 0, f"POLISHED_FROM_{step}"


def polish_bfgs(args: argparse.Namespace) -> None:
    if args.force and args.method != "GXTB":
        raise ValueError("--force is restricted to --method GXTB")
    if args.prune_transients and args.method != "GXTB":
        raise ValueError("--prune-transients is restricted to --method GXTB")
    wanted = set(args.system) if args.system else None
    inputs = manifest_owned_inputs(args.output_root, args.method, wanted)
    if not inputs:
        raise ValueError("no prepared cell-optimization inputs selected")
    campaign_identity = None
    if args.method == "GXTB":
        common.require_gxtb_build_artifacts(
            cp2k=args.cp2k,
            cp2k_source=args.cp2k_source,
            save_tblite=args.save_tblite,
            save_tblite_source=args.save_tblite_source,
            campaign_manifest=args.campaign_manifest,
        )
        update_provenance(args)
        campaign_identity = common.load_campaign_identity(ROOT)
    failed = []
    contexts = {
        path: stamp_context(manifest_record_for_input(args.output_root, path))
        for path in inputs
    }
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                polish_one,
                path,
                args.cp2k.resolve(),
                args.max_iter,
                args.trust_radius,
                args.force,
                args.method or next(part for part in path.parts if part in METHODS),
                args.threads_per_job,
                args.prune_transients,
                campaign_identity,
                contexts[path][0],
                contexts[path][1],
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:24s} {relative} rc={returncode}", flush=True)
            if returncode != 0:
                failed.append(relative)
    if args.method == "GXTB":
        update_provenance(args)
    if failed:
        raise SystemExit(f"{len(failed)} BFGS polishing jobs failed")


def last_float(text: str, pattern: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.M)
    return float(matches[-1]) if matches else None


def last_int(text: str, pattern: str) -> int | None:
    matches = re.findall(pattern, text, flags=re.M)
    return int(matches[-1]) if matches else None


def format_number(value: float | int | None, digits: int = 6) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def load_gamma_rows(path: Path, methods: tuple[str, ...] = common.PUBLISHED_METHODS) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = {(row["method"], row["system"]): row for row in rows}
    expected = {(method, str(system["id"])) for system in systems() for method in methods}
    missing = expected - set(selected)
    if missing:
        raise ValueError(f"Gamma result table is missing {len(missing)} method/system rows")
    return selected


def load_molecule_rows(
    run_root: Path,
    methods: tuple[str, ...] = common.PUBLISHED_METHODS,
) -> dict[tuple[str, str], dict[str, str]]:
    rows: dict[tuple[str, str], dict[str, str]] = {}
    campaign_identity = (
        common.load_campaign_identity(run_root) if "GXTB" in methods else None
    )
    for system_data in systems():
        system = str(system_data["id"])
        for method in methods:
            stem = f"{system}_{method}_mol_geoopt"
            output = run_root / "runs" / "molecule_geoopt" / method / stem / f"{stem}.out"
            text = output.read_text(errors="ignore") if output.is_file() else ""
            energy = last_float(text, r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$")
            if "GEOMETRY OPTIMIZATION COMPLETED" not in text or "PROGRAM ENDED" not in text or energy is None:
                raise ValueError(f"completed molecule optimization not found for {method}/{system}: {output}")
            if method == "GXTB":
                input_path = run_root / "inputs" / "molecule_geoopt" / method / f"{stem}.inp"
                valid, reason = common.recorded_job_stamp_matches(
                    output.parent,
                    input_path,
                    method,
                    "x23b_molecule_geoopt",
                    output,
                    campaign_identity=campaign_identity,
                )
                if not valid:
                    raise ValueError(f"untrusted molecule result for {method}/{system}: {reason}")
            rows[(method, system)] = {"gas_energy_hartree": f"{energy:.12f}"}
    return rows


def collect(args: argparse.Namespace) -> None:
    selected_methods = (args.method,) if args.method else common.PUBLISHED_METHODS
    if "GXTB" in selected_methods and args.gamma_csv is not None:
        raise ValueError(
            "--gamma-csv is not accepted for GXTB; use the method-owned stamped "
            "--molecule-run-root"
        )
    campaign_identity = (
        common.load_campaign_identity(ROOT) if "GXTB" in selected_methods else None
    )
    output_root = args.output_root.resolve()
    primary_manifest = manifest_path(output_root)
    with primary_manifest.open(newline="") as handle:
        manifest: dict[tuple[str, str], dict[str, str]] = {}
        for raw_row in csv.DictReader(handle):
            key = (raw_row["method"], raw_row["system"])
            if key in manifest:
                raise ValueError(f"duplicate cell-optimization manifest row: {key[0]}/{key[1]}")
            if key[0] in selected_methods:
                manifest[key] = validate_manifest_row(raw_row, output_root)
    if args.gamma_csv is not None:
        gamma_rows = load_gamma_rows(args.gamma_csv.resolve(), selected_methods)
    else:
        gamma_rows = load_molecule_rows(args.molecule_run_root.resolve(), selected_methods)
    metadata = {str(row["id"]): row for row in systems()}
    rows = []
    expected_manifest = {(method, str(row["id"])) for method in selected_methods for row in systems()}
    missing_manifest = sorted(expected_manifest - set(manifest))
    if missing_manifest:
        raise ValueError(
            "cell-optimization manifest is incomplete: "
            + ", ".join(f"{method}/{system}" for method, system in missing_manifest)
        )
    for method in selected_methods:
        for system in sorted(metadata):
            source = gamma_rows[(method, system)]
            manifest_row = manifest[(method, system)]
            run_dir = Path(manifest_row["run_dir"]).resolve()
            output = run_dir / "cp2k.out"
            text = output.read_text(errors="ignore") if output.is_file() else ""
            if method == "GXTB" and output.is_file():
                input_path = Path(manifest_row["input"])
                protocol_identity, source_artifacts = stamp_context(manifest_row)
                valid, reason = common.recorded_job_stamp_matches(
                    run_dir,
                    input_path,
                    method,
                    "x23b_k222_cellopt",
                    output,
                    campaign_identity=campaign_identity,
                    protocol_identity=protocol_identity,
                    source_artifacts=source_artifacts,
                )
                if not valid:
                    raise ValueError(f"untrusted k222 result for {method}/{system}: {reason}")
            energy = last_float(text, r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$")
            volume = last_float(text, r"^\s*CELL\| Volume.*?([-+0-9.Ee]+)\s*$")
            step = last_int(text, r"^\s*OPT\| Step number\s+(\d+)\s*$")
            pressure = last_float(text, r"^\s*OPT\| (?:Internal pressure|Pressure deviation) \[?bar\]?\s+([-+0-9.Ee]+)\s*$")
            max_step = last_float(text, r"^\s*OPT\| Maximum step size\s+([-+0-9.Ee]+)\s*$")
            rms_step = last_float(text, r"^\s*OPT\| RMS step size\s+([-+0-9.Ee]+)\s*$")
            max_gradient = last_float(text, r"^\s*OPT\| Maximum gradient\s+([-+0-9.Ee]+)\s*$")
            rms_gradient = last_float(text, r"^\s*OPT\| RMS gradient\s+([-+0-9.Ee]+)\s*$")
            gas_energy = float(source["gas_energy_hartree"])
            ref_energy = float(metadata[system]["ref_energy"])
            ref_volume = float(metadata[system]["x23b_same_cell_ref_volume"])
            n_molecules = int(metadata[system]["molecules_per_cell"])
            lattice = None if energy is None else (gas_energy - energy / n_molecules) * HARTREE_TO_KJMOL
            error = None if lattice is None else lattice - ref_energy
            volume_error = None if volume is None else 100.0 * (volume - ref_volume) / ref_volume
            returncode_file = run_dir / "returncode.txt"
            returncode = int(returncode_file.read_text()) if returncode_file.is_file() else None
            row = {
                "method": method,
                "system": system,
                "variant": manifest_row["variant"],
                "mesh": "k222",
                "source_variant": (
                    "experimental_k222_preflight"
                    if manifest_row["source_policy"] == "experimental_reference"
                    else "gamma_cellopt_keep_angles"
                ),
                "source": manifest_row["source_kind"],
                "source_policy": manifest_row["source_policy"],
                "source_path": manifest_row["source_path"],
                "source_sha256": manifest_row["source_sha256"],
                "structure_path": manifest_row["structure_path"],
                "structure_sha256": manifest_row["structure_sha256"],
                "structure_source": manifest_row["structure_source"],
                "lineage": manifest_row["lineage"],
                "preflight_input": manifest_row["preflight_input"],
                "preflight_input_sha256": manifest_row["preflight_input_sha256"],
                "preflight_output": manifest_row["preflight_output"],
                "preflight_output_sha256": manifest_row["preflight_output_sha256"],
                "preflight_stamp": manifest_row["preflight_stamp"],
                "preflight_stamp_sha256": manifest_row["preflight_stamp_sha256"],
                "returncode": "" if returncode is None else returncode,
                "program_ended": "PROGRAM ENDED" in text,
                "opt_completed": "GEOMETRY OPTIMIZATION COMPLETED" in text,
                "max_iter_reached": "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text,
                "last_step": format_number(step),
                "energy_hartree": format_number(energy, 12),
                "gas_energy_hartree": format_number(gas_energy, 12),
                "lattice_energy_kJmol": format_number(lattice),
                "x23b_ref_lattice_energy_kJmol": format_number(ref_energy),
                "error_kJmol": format_number(error),
                "volume_A3": format_number(volume),
                "x23b_same_cell_ref_volume_A3": format_number(ref_volume),
                "volume_error_percent": format_number(volume_error),
                "last_pressure_bar": format_number(pressure),
                "last_max_step": format_number(max_step, 10),
                "last_rms_step": format_number(rms_step, 10),
                "last_max_gradient": format_number(max_gradient, 10),
                "last_rms_gradient": format_number(rms_gradient, 10),
                "source_restart": manifest_row["source_restart"],
                "run_dir": str(run_dir),
                "output": str(output),
            }
            rows.append(row)
    incomplete = []
    for method in selected_methods:
        selected = [row for row in rows if row["method"] == method]
        complete = sum(row["program_ended"] and row["opt_completed"] for row in selected)
        print(f"{method}: {complete}/{len(selected)} converged")
        if complete != len(systems()):
            incomplete.append(f"{method}={complete}/{len(systems())}")
    if incomplete and not args.allow_incomplete:
        raise ValueError("complete 23-system k222 coverage required: " + ", ".join(incomplete))
    if args.method == "GXTB":
        update_provenance(args)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def add_provenance_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cp2k-source", type=Path)
    parser.add_argument("--save-tblite", type=Path)
    parser.add_argument("--save-tblite-source", type=Path)
    parser.add_argument(
        "--campaign-manifest",
        type=Path,
        default=common.DEFAULT_CAMPAIGN_MANIFEST,
        help="frozen V1 build manifest (single source of truth)",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--gamma-root", type=Path)
    prepare_parser.add_argument("--preflight-root", type=Path)
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.add_argument("--override", action="append", default=[])
    prepare_parser.add_argument(
        "--source-policy",
        choices=SOURCE_POLICIES,
        help="mandatory explicit source policy for GXTB; GFN1/GFN2 retain the Gamma default",
    )
    prepare_parser.add_argument("--method", choices=METHODS)
    prepare_parser.add_argument("--system", action="append", choices=sorted(str(row["id"]) for row in systems()))
    prepare_parser.add_argument("--clean", action="store_true")
    prepare_parser.set_defaults(function=prepare)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--cp2k", type=Path, required=True)
    run_parser.add_argument("--jobs", type=int, default=8)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    run_parser.add_argument("--method", choices=METHODS)
    run_parser.add_argument("--system", action="append", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--start-system", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--end-system", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--force", action="store_true")
    run_parser.add_argument("--prune-transients", action="store_true")
    add_provenance_arguments(run_parser)
    run_parser.set_defaults(function=run)

    continue_parser = subparsers.add_parser("continue-maxiter")
    continue_parser.add_argument("--output-root", type=Path, required=True)
    continue_parser.add_argument("--cp2k", type=Path, required=True)
    continue_parser.add_argument("--jobs", type=int, default=2)
    continue_parser.add_argument("--threads-per-job", type=int, default=1)
    continue_parser.add_argument("--additional-steps", type=int, default=300)
    continue_parser.add_argument("--rounds", type=int, default=3)
    continue_parser.add_argument("--method", choices=METHODS)
    continue_parser.add_argument("--system", action="append")
    continue_parser.add_argument("--prune-transients", action="store_true")
    add_provenance_arguments(continue_parser)
    continue_parser.set_defaults(function=continue_maxiter)

    polish_parser = subparsers.add_parser("polish-bfgs")
    polish_parser.add_argument("--output-root", type=Path, required=True)
    polish_parser.add_argument("--cp2k", type=Path, required=True)
    polish_parser.add_argument("--jobs", type=int, default=2)
    polish_parser.add_argument("--threads-per-job", type=int, default=1)
    polish_parser.add_argument("--max-iter", type=int, default=300)
    polish_parser.add_argument("--trust-radius", type=float, default=0.002)
    polish_parser.add_argument("--method", choices=METHODS)
    polish_parser.add_argument("--system", action="append")
    polish_parser.add_argument("--force", action="store_true")
    polish_parser.add_argument("--prune-transients", action="store_true")
    add_provenance_arguments(polish_parser)
    polish_parser.set_defaults(function=polish_bfgs)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    energy_source = collect_parser.add_mutually_exclusive_group(required=True)
    energy_source.add_argument("--gamma-csv", type=Path)
    energy_source.add_argument("--molecule-run-root", type=Path)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.add_argument("--method", choices=METHODS)
    collect_parser.add_argument("--allow-incomplete", action="store_true")
    collect_parser.set_defaults(function=collect)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
