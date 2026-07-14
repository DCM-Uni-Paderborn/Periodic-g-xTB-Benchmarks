#!/usr/bin/env python3
"""Run native-Bloch single points on converged X23b k222 cell-opt geometries."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import math
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import x23b_kpoint_cellopt as cellopt
import x23b_common as common


FIELDS = (
    "method",
    "system",
    "target_mesh",
    "source_mesh",
    "program_ended",
    "source_energy_hartree",
    "target_energy_hartree",
    "source_lattice_energy_kJmol",
    "target_lattice_energy_kJmol",
    "delta_target_minus_source_kJmol",
    "target_error_kJmol",
    "source_restart",
    "output",
)
MANIFEST_FIELDS = (
    "method",
    "system",
    "mesh",
    "source_run_dir",
    "source_input",
    "source_input_sha256",
    "source_output_sha256",
    "source_policy",
    "source_variant",
    "source_protocol_identity",
    "source_artifacts",
    "source_artifact_hashes",
    "source_restart",
    "source_restart_sha256",
    "input",
    "input_sha256",
    "run_dir",
)


def _normalized_manifest_row(row: dict[str, str]) -> dict[str, str]:
    return {field: str(row.get(field, "") or "") for field in MANIFEST_FIELDS}


def variant(mesh: int) -> str:
    if mesh not in (3, 4):
        raise ValueError("final production mesh must be 3 or 4")
    return f"k{mesh}{mesh}{mesh}_sp_on_k222"


def gamma_centered_shift(mesh: int) -> float:
    return 0.0 if mesh % 2 else (mesh - 1) / (2.0 * mesh)


def final_restart(run_dir: Path) -> Path:
    restart = common.final_restart(run_dir)
    if restart is None:
        raise ValueError(f"final restart not found in {run_dir}")
    return restart


def restart_to_single_point(source: Path, project: str, mesh: int, method: str) -> str:
    lines = source.read_text().splitlines()
    global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    lines = lines[global_start:]

    motion_start, motion_end = cellopt.section_bounds(lines, "MOTION")
    del lines[motion_start : motion_end + 1]

    replacements = (
        (re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I), rf'\1"{project}"'),
        (re.compile(r"^(\s*RUN_TYPE\s+).*$", re.I), r"\1ENERGY"),
    )
    for pattern, replacement in replacements:
        for index, line in enumerate(lines):
            if pattern.match(line):
                lines[index] = pattern.sub(replacement, line, count=1)
                break
        else:
            raise ValueError(f"required GLOBAL keyword missing in {source}: {pattern.pattern}")

    k_start, k_end = cellopt.section_bounds(lines, "KPOINTS")
    shift = gamma_centered_shift(mesh)
    for index in range(k_start, k_end + 1):
        if lines[index].strip().upper().startswith("SCHEME"):
            indent = re.match(r"\s*", lines[index]).group(0)
            lines[index] = (
                f"{indent}SCHEME MACDONALD {mesh} {mesh} {mesh} "
                f"{shift:.12g} {shift:.12g} {shift:.12g}"
            )
            break
    else:
        raise ValueError(f"KPOINTS SCHEME missing in {source}")
    if method == "GXTB":
        block = "\n".join(lines[k_start : k_end + 1])
        required = (
            re.search(r"^\s*SYMMETRY\s+T\s*$", block, re.M | re.I),
            re.search(r"^\s*FULL_GRID\s+F\s*$", block, re.M | re.I),
            re.search(r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$", block, re.M | re.I),
            re.search(r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$", block, re.M | re.I),
        )
        if not all(required):
            raise ValueError(
                f"GXTB source restart does not contain the SPGLIB-reduced production mesh: {source}"
            )
    text = "\n".join(lines) + "\n"
    common.validate_method_input(text, method)
    return text


def manifest_path(output_root: Path) -> Path:
    return output_root / "manifest.csv"


def manifest_owned_inputs(
    output_root: Path,
    mesh: int,
    method: str | None = None,
    selected_systems: set[str] | None = None,
) -> list[Path]:
    """Return exactly the primary final-SP inputs owned by the manifest."""

    output_root = output_root.resolve()
    path = manifest_path(output_root)
    if not path.is_file():
        raise FileNotFoundError(path)
    methods = (method,) if method else common.PUBLISHED_METHODS
    selected: dict[tuple[str, str], Path] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            row_method, system = row["method"], row["system"]
            if row_method not in methods:
                continue
            if selected_systems is not None and system not in selected_systems:
                continue
            if int(row.get("mesh", mesh)) != mesh:
                raise ValueError(f"manifest mesh mismatch for {row_method}/{system}")
            if row_method == "GXTB":
                target_stamp_context(row)
            key = (row_method, system)
            if key in selected:
                raise ValueError(f"duplicate final-SP manifest row for {row_method}/{system}")
            input_path = Path(row["input"]).resolve()
            run_dir = Path(row["run_dir"]).resolve()
            if input_path.parent != run_dir:
                raise ValueError(f"manifest input/run_dir mismatch for {row_method}/{system}")
            try:
                input_path.relative_to(output_root)
            except ValueError as exc:
                raise ValueError(f"manifest input escapes output root: {input_path}") from exc
            if not input_path.is_file():
                raise FileNotFoundError(input_path)
            selected[key] = input_path
    return [selected[key] for key in sorted(selected)]


def manifest_record_for_input(output_root: Path, input_path: Path) -> dict[str, str]:
    input_path = input_path.resolve()
    matches: list[dict[str, str]] = []
    with manifest_path(output_root.resolve()).open(newline="") as handle:
        for row in csv.DictReader(handle):
            if Path(row["input"]).resolve() == input_path:
                matches.append(_normalized_manifest_row(row))
    if len(matches) != 1:
        raise ValueError(f"expected one final-SP manifest row for {input_path}, found {len(matches)}")
    return matches[0]


def target_stamp_context(
    row: dict[str, str],
) -> tuple[dict[str, object], dict[str, Path]]:
    row = _normalized_manifest_row(row)
    if not all(
        row[field]
        for field in (
            "source_policy",
            "source_variant",
            "source_protocol_identity",
            "source_artifacts",
            "source_input",
            "source_input_sha256",
            "source_output_sha256",
            "source_run_dir",
            "source_restart",
            "source_restart_sha256",
            "source_artifact_hashes",
            "input",
            "input_sha256",
        )
    ):
        raise ValueError(f"policyless GXTB final-SP manifest row: {row['system']}")
    source_protocol = json.loads(row["source_protocol_identity"])
    source_files_payload = json.loads(row["source_artifacts"])
    source_hashes = json.loads(row["source_artifact_hashes"])
    if (
        not isinstance(source_protocol, dict)
        or not isinstance(source_files_payload, dict)
        or not isinstance(source_hashes, dict)
    ):
        raise ValueError(f"invalid GXTB source lineage in final-SP manifest: {row['system']}")
    if (
        source_protocol.get("source_policy") != row["source_policy"]
        or source_protocol.get("variant") != row["source_variant"]
    ):
        raise ValueError(f"GXTB source policy/variant lineage differs: {row['system']}")
    source_files = {name: Path(path) for name, path in source_files_payload.items()}
    if set(source_files) != set(source_hashes):
        raise ValueError(f"incomplete GXTB source-artifact hashes: {row['system']}")
    for name, path in source_files.items():
        if common.sha256_file(path) != source_hashes[name]:
            raise ValueError(f"GXTB {name} fingerprint differs: {row['system']}")
    source_files.update(
        {
            "cellopt_input": Path(row["source_input"]),
            "cellopt_output": Path(row["source_run_dir"]) / "cp2k.out",
            "cellopt_restart": Path(row["source_restart"]),
        }
    )
    frozen_paths = {
        "cellopt_input": row["source_input_sha256"],
        "cellopt_output": row["source_output_sha256"],
        "cellopt_restart": row["source_restart_sha256"],
    }
    for name, expected_hash in frozen_paths.items():
        if common.sha256_file(source_files[name]) != expected_hash:
            raise ValueError(f"GXTB {name} fingerprint differs: {row['system']}")
    input_path = Path(row["input"])
    if common.sha256_file(input_path) != row["input_sha256"]:
        raise ValueError(f"GXTB target input fingerprint differs: {row['system']}")
    mesh = int(row["mesh"])
    text = input_path.read_text()
    common.validate_method_input(text, "GXTB")
    shift = gamma_centered_shift(mesh)
    expected_scheme = (
        rf"^\s*SCHEME\s+MACDONALD\s+{mesh}\s+{mesh}\s+{mesh}\s+"
        rf"{re.escape(f'{shift:.12g}')}\s+{re.escape(f'{shift:.12g}')}\s+"
        rf"{re.escape(f'{shift:.12g}')}\s*$"
    )
    if len(re.findall(expected_scheme, text, flags=re.I | re.M)) != 1:
        raise ValueError(f"GXTB target input mesh contract differs: {row['system']}")
    if len(re.findall(r"^\s*RUN_TYPE\s+ENERGY\s*$", text, flags=re.I | re.M)) != 1:
        raise ValueError(f"GXTB target input RUN_TYPE differs: {row['system']}")
    if re.search(r"^\s*&MOTION\b", text, flags=re.I | re.M):
        raise ValueError(f"GXTB target input unexpectedly contains MOTION: {row['system']}")
    protocol = {
        "source_policy": row["source_policy"],
        "source_variant": row["source_variant"],
        "source_protocol_identity": source_protocol,
        "target_mesh": f"k{int(row['mesh'])}{int(row['mesh'])}{int(row['mesh'])}",
    }
    return protocol, source_files


def update_provenance(args: argparse.Namespace) -> Path:
    mesh_id = f"k{args.mesh}{args.mesh}{args.mesh}"
    provenance = cellopt.ROOT / "data" / common.GXTB_PROVENANCE_NAME
    if provenance.is_file():
        payload = json.loads(provenance.read_text())
        workflow_paths = payload.get("workflow_paths", {})
        if isinstance(workflow_paths, dict):
            roots = workflow_paths.get("final_single_point_roots", {})
            if isinstance(roots, dict):
                target_root = args.output_root.resolve()
                for frozen_mesh, frozen_value in roots.items():
                    frozen = Path(str(frozen_value)).resolve()
                    if frozen_mesh == mesh_id and frozen != target_root:
                        raise ValueError(f"{mesh_id} final-SP root is already frozen as {frozen}")
                    if frozen_mesh != mesh_id and frozen == target_root:
                        raise ValueError(
                            f"final-SP root {target_root} is already frozen for {frozen_mesh}"
                        )
    return common.update_gxtb_provenance(
        cellopt.ROOT,
        cp2k=getattr(args, "cp2k", None),
        cp2k_source=getattr(args, "cp2k_source", None),
        save_tblite=getattr(args, "save_tblite", None),
        save_tblite_source=getattr(args, "save_tblite_source", None),
        campaign_manifest=getattr(args, "campaign_manifest", None),
        workflow_paths={"final_single_point_roots": {mesh_id: args.output_root}},
    )


def prepare(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / f".{output_root.name}.manifest.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _prepare(args, output_root)


def guard_output_root_mesh(output_root: Path, mesh: int) -> None:
    """Refuse cross-mesh reuse before prepare can delete or rewrite anything."""

    mesh_id = f"k{mesh}{mesh}{mesh}"
    provenance = cellopt.ROOT / "data" / common.GXTB_PROVENANCE_NAME
    if provenance.is_file():
        payload = json.loads(provenance.read_text())
        paths = payload.get("workflow_paths", {})
        roots = paths.get("final_single_point_roots", {}) if isinstance(paths, dict) else {}
        if isinstance(roots, dict):
            for frozen_mesh, frozen_value in roots.items():
                frozen = Path(str(frozen_value)).resolve()
                if frozen_mesh == mesh_id and frozen != output_root:
                    raise ValueError(f"{mesh_id} final-SP root is already frozen as {frozen}")
                if frozen_mesh != mesh_id and frozen == output_root:
                    raise ValueError(
                        f"final-SP root {output_root} is already frozen for {frozen_mesh}"
                    )
    path = manifest_path(output_root)
    if path.is_file():
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        observed: set[int] = set()
        for row in rows:
            try:
                observed.add(int(row.get("mesh", "3") or "3"))
            except ValueError as exc:
                raise ValueError("invalid mesh in existing final-SP manifest") from exc
        if observed and observed != {mesh}:
            raise ValueError(
                f"existing final-SP manifest belongs to meshes {sorted(observed)}, not {mesh}"
            )


def _prepare(args: argparse.Namespace, output_root: Path) -> None:
    if args.mesh not in (3, 4):
        raise ValueError("final production mesh must be 3 or 4")
    guard_output_root_mesh(output_root, args.mesh)
    if args.clean and args.method != "GXTB":
        raise ValueError("--clean is restricted to --method GXTB; published method trees are immutable")
    if args.clean and (output_root / args.method).exists():
        shutil.rmtree(output_root / args.method)
    methods = (args.method,) if args.method else common.PUBLISHED_METHODS
    selected_systems = set(args.system) if args.system else None
    campaign_identity = (
        common.load_campaign_identity(cellopt.ROOT) if "GXTB" in methods else None
    )
    cellopt_root = args.cellopt_root.resolve()
    source_manifest: dict[tuple[str, str], dict[str, str]] = {}
    with cellopt.manifest_path(cellopt_root).open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["method"], row["system"])
            if key in source_manifest:
                raise ValueError(f"duplicate k222 source manifest row for {key[0]}/{key[1]}")
            if key[0] in methods:
                source_manifest[key] = cellopt.validate_manifest_row(row, cellopt_root)
    manifest: dict[tuple[str, str], dict[str, str]] = {}
    path = manifest_path(output_root)
    if path.is_file() and not args.clean:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                row.setdefault("mesh", "3")
                manifest[(row["method"], row["system"])] = row
    elif path.is_file() and args.clean:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["method"] != args.method:
                    row.setdefault("mesh", "3")
                    manifest[(row["method"], row["system"])] = row

    prepared = 0
    for system_data in cellopt.systems():
        system = str(system_data["id"])
        if selected_systems is not None and system not in selected_systems:
            continue
        for method in methods:
            source_record = source_manifest.get((method, system))
            if source_record is None:
                raise ValueError(f"k222 source manifest row required: {method}/{system}")
            source_run_dir = Path(source_record["run_dir"]).resolve()
            source_input = Path(source_record["input"]).resolve()
            source_output = source_run_dir / "cp2k.out"
            if not cellopt.cp2k_completed(source_output):
                raise ValueError(f"converged k222 optimization required: {method}/{system}")
            if method == "GXTB":
                source_protocol, source_files = cellopt.stamp_context(source_record)
                valid, reason = common.recorded_job_stamp_matches(
                    source_run_dir,
                    source_input,
                    method,
                    "x23b_k222_cellopt",
                    source_output,
                    campaign_identity=campaign_identity,
                    protocol_identity=source_protocol,
                    source_artifacts=source_files,
                )
                if not valid:
                    raise ValueError(f"untrusted k222 source for {method}/{system}: {reason}")
            source_restart = final_restart(source_run_dir)
            target_variant = variant(args.mesh)
            project = f"{system}_{method}_{target_variant}".replace("-", "_")
            run_dir = output_root / method / system / target_variant
            if not args.clean and (run_dir / "cp2k.out").exists():
                raise ValueError(f"refusing to prepare over existing final-SP output: {run_dir / 'cp2k.out'}")
            run_dir.mkdir(parents=True, exist_ok=True)
            input_path = run_dir / f"{project}.inp"
            input_path.write_text(restart_to_single_point(source_restart, project, args.mesh, method))
            manifest[(method, system)] = {
                "method": method,
                "system": system,
                "mesh": str(args.mesh),
                "source_run_dir": str(source_run_dir),
                "source_input": str(source_input),
                "source_input_sha256": common.sha256_file(source_input),
                "source_output_sha256": common.sha256_file(source_output),
                "source_policy": source_record["source_policy"],
                "source_variant": source_record["variant"],
                "source_protocol_identity": json.dumps(
                    source_protocol if method == "GXTB" else None,
                    sort_keys=True,
                ),
                "source_artifacts": json.dumps(
                    {name: str(path.resolve()) for name, path in (source_files or {}).items()}
                    if method == "GXTB"
                    else None,
                    sort_keys=True,
                ),
                "source_artifact_hashes": json.dumps(
                    {
                        name: common.sha256_file(path)
                        for name, path in (source_files or {}).items()
                    }
                    if method == "GXTB"
                    else None,
                    sort_keys=True,
                ),
                "source_restart": str(source_restart),
                "source_restart_sha256": common.sha256_file(source_restart),
                "input": str(input_path),
                "input_sha256": common.sha256_file(input_path),
                "run_dir": str(run_dir),
            }
            prepared += 1

    rows = [manifest[key] for key in sorted(manifest)]
    output_root.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(_normalized_manifest_row(row) for row in rows)
    print(f"Prepared {prepared} k{args.mesh}{args.mesh}{args.mesh} single points; manifest has {len(rows)} entries")


def output_ok(path: Path) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and "ENERGY| Total FORCE_EVAL" in text


def run_one(
    input_path: Path,
    cp2k: Path,
    threads: int,
    mesh: int,
    method: str,
    force: bool,
    prune_transients: bool,
    campaign_identity: dict[str, object] | None = None,
    protocol_identity: dict[str, object] | None = None,
    source_artifacts: dict[str, Path] | None = None,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    output = run_dir / "cp2k.out"
    common.validate_method_input(input_path.read_text(), method)
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, common.BUSY_RETURN_CODE, "BUSY"
        stamp_phase = f"x23b_final_k{mesh}{mesh}{mesh}_on_k222"
        stamp_matches, _ = common.job_stamp_matches(
            run_dir,
            input_path,
            cp2k,
            method,
            stamp_phase,
            campaign_identity=campaign_identity,
            protocol_identity=protocol_identity,
            source_artifacts=source_artifacts,
        )
        if not force and output.exists() and method == "GXTB":
            recorded_matches, _ = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                method,
                stamp_phase,
                output,
                campaign_identity=campaign_identity,
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
            if not stamp_matches or not recorded_matches:
                return input_path, 1, "STALE_OUTPUT"
        if not force and output_ok(output):
            if prune_transients:
                common.prune_gxtb_transients(run_dir, keep_final_restart=False)
            return input_path, 0, "SKIP"
        if force:
            for path in run_dir.iterdir():
                if path != input_path:
                    path.unlink() if path.is_file() else shutil.rmtree(path)
        env = common.thread_environment(threads)
        code = subprocess.run(
            [str(cp2k), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            env=env,
            check=False,
        ).returncode
        action = "CONVERGED" if code == 0 and output_ok(output) else "FAILED"
        if method == "GXTB":
            details: dict[str, object] = {"returncode": code, "action": action, "output": str(output)}
            if output.is_file():
                details["output_sha256"] = common.sha256_file(output)
            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                method,
                stamp_phase,
                action.lower(),
                details=details,
                campaign_identity=campaign_identity,
                protocol_identity=protocol_identity,
                source_artifacts=source_artifacts,
            )
        if action == "CONVERGED" and prune_transients:
            common.prune_gxtb_transients(run_dir, keep_final_restart=False)
        return input_path, code, action


def run(args: argparse.Namespace) -> None:
    if args.force and args.method != "GXTB":
        raise ValueError("--force is restricted to --method GXTB")
    if args.prune_transients and args.method != "GXTB":
        raise ValueError("--prune-transients is restricted to --method GXTB")
    wanted = set(args.system) if args.system else None
    inputs = manifest_owned_inputs(args.output_root, args.mesh, args.method, wanted)
    system_count = len(set(args.system)) if args.system else len(cellopt.systems())
    method_count = 1 if args.method else len(common.PUBLISHED_METHODS)
    expected = system_count * method_count
    if len(inputs) != expected:
        raise ValueError(f"expected {expected} prepared single-point inputs, found {len(inputs)}")
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
        campaign_identity = common.load_campaign_identity(cellopt.ROOT)
    failed = []
    contexts = {
        path: (
            target_stamp_context(manifest_record_for_input(args.output_root, path))
            if args.method == "GXTB"
            else (None, None)
        )
        for path in inputs
    }
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_one,
                path,
                args.cp2k.resolve(),
                args.threads_per_job,
                args.mesh,
                args.method or next(part for part in path.parts if part in cellopt.METHODS),
                args.force,
                args.prune_transients,
                campaign_identity,
                contexts[path][0],
                contexts[path][1],
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, code, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:9s} {relative} rc={code}", flush=True)
            if code != 0 or action in {"FAILED", "BUSY", "STALE_OUTPUT"}:
                failed.append(relative)
    if args.method == "GXTB":
        update_provenance(args)
    if failed:
        raise SystemExit(f"{len(failed)} single-point jobs failed")


def finite(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.12f}"


def collect(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    with manifest_path(output_root).open(newline="") as handle:
        manifest = list(csv.DictReader(handle))
    selected_methods = (args.method,) if args.method else common.PUBLISHED_METHODS
    campaign_identity = (
        common.load_campaign_identity(cellopt.ROOT) if "GXTB" in selected_methods else None
    )
    manifest = [row for row in manifest if row["method"] in selected_methods]
    identities = [(row["method"], row["system"]) for row in manifest]
    if len(identities) != len(set(identities)):
        raise ValueError("duplicate method/system rows in final-kpoint manifest")
    wrong_mesh = [
        f"{row['method']}/{row['system']}=k{row.get('mesh', '?')}"
        for row in manifest
        if int(row.get("mesh", -1)) != args.mesh
    ]
    if wrong_mesh:
        raise ValueError(
            f"collect --mesh {args.mesh} differs from final-kpoint manifest: "
            + ", ".join(wrong_mesh)
        )
    expected = {
        (method, str(system["id"]))
        for method in selected_methods
        for system in cellopt.systems()
    }
    present = {(row["method"], row["system"]) for row in manifest}
    if present != expected:
        raise ValueError(
            f"complete {len(expected)}-case final-kpoint manifest required: "
            f"missing={sorted(expected - present)}, unexpected={sorted(present - expected)}"
        )
    metadata = {str(row["id"]): row for row in cellopt.systems()}
    gas = cellopt.load_molecule_rows(args.molecule_run_root.resolve(), selected_methods)
    rows = []
    for entry in manifest:
        method, system = entry["method"], entry["system"]
        target_mesh = int(entry["mesh"])
        k222_output = Path(entry["source_run_dir"]) / "cp2k.out"
        k333_output = Path(entry["run_dir"]) / "cp2k.out"
        k222_text = k222_output.read_text(errors="ignore")
        k333_text = k333_output.read_text(errors="ignore") if k333_output.is_file() else ""
        if method == "GXTB":
            if not entry.get("source_input") or not entry.get("source_protocol_identity"):
                raise ValueError(f"policyless GXTB final-SP source manifest row: {system}")
            source_input = Path(entry["source_input"])
            source_protocol = json.loads(entry["source_protocol_identity"])
            source_files_payload = json.loads(entry["source_artifacts"])
            source_files = {
                name: Path(path) for name, path in source_files_payload.items()
            }
            source_valid, source_reason = common.recorded_job_stamp_matches(
                k222_output.parent,
                source_input,
                method,
                "x23b_k222_cellopt",
                k222_output,
                campaign_identity=campaign_identity,
                protocol_identity=source_protocol,
                source_artifacts=source_files,
            )
            if not source_valid:
                raise ValueError(f"untrusted k222 source for {method}/{system}: {source_reason}")
            if k333_output.is_file():
                target_protocol, target_sources = target_stamp_context(entry)
                target_valid, target_reason = common.recorded_job_stamp_matches(
                    k333_output.parent,
                    Path(entry["input"]),
                    method,
                    f"x23b_final_k{target_mesh}{target_mesh}{target_mesh}_on_k222",
                    k333_output,
                    campaign_identity=campaign_identity,
                    protocol_identity=target_protocol,
                    source_artifacts=target_sources,
                )
                if not target_valid:
                    raise ValueError(f"untrusted final SP for {method}/{system}: {target_reason}")
        pattern = r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$"
        e222 = cellopt.last_float(k222_text, pattern)
        e333 = cellopt.last_float(k333_text, pattern)
        n_molecules = int(metadata[system]["molecules_per_cell"])
        gas_energy = float(gas[(method, system)]["gas_energy_hartree"])
        ref = float(metadata[system]["ref_energy"])
        lattice222 = None if e222 is None else (gas_energy - e222 / n_molecules) * cellopt.HARTREE_TO_KJMOL
        lattice333 = None if e333 is None else (gas_energy - e333 / n_molecules) * cellopt.HARTREE_TO_KJMOL
        rows.append(
            {
                "method": method,
                "system": system,
                "target_mesh": f"k{target_mesh}{target_mesh}{target_mesh}",
                "source_mesh": "k222_cellopt",
                "program_ended": "PROGRAM ENDED" in k333_text,
                "source_energy_hartree": finite(e222),
                "target_energy_hartree": finite(e333),
                "source_lattice_energy_kJmol": finite(lattice222),
                "target_lattice_energy_kJmol": finite(lattice333),
                "delta_target_minus_source_kJmol": finite(
                    None if lattice222 is None or lattice333 is None else lattice333 - lattice222
                ),
                "target_error_kJmol": finite(None if lattice333 is None else lattice333 - ref),
                "source_restart": entry["source_restart"],
                "output": str(k333_output),
            }
        )
    incomplete = [
        f"{row['method']}/{row['system']}"
        for row in rows
        if str(row["program_ended"]) != "True" or not row["target_energy_hartree"]
    ]
    if incomplete and not args.allow_incomplete:
        raise ValueError("complete final-kpoint coverage required: " + ", ".join(incomplete))
    if args.method == "GXTB":
        update_provenance(args)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for target_mesh in sorted({row["target_mesh"] for row in rows}):
        for method in selected_methods:
            deltas = [
                abs(float(row["delta_target_minus_source_kJmol"]))
                for row in rows
                if row["method"] == method
                and row["target_mesh"] == target_mesh
                and row["delta_target_minus_source_kJmol"]
            ]
            if deltas:
                print(
                    f"{method} {target_mesh}: {len(deltas)} complete, "
                    f"mean |target-k222|={sum(deltas) / len(deltas):.6f}, max={max(deltas):.6f} kJ/mol"
                )
            else:
                print(f"{method} {target_mesh}: 0 complete")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    choices = sorted(str(row["id"]) for row in cellopt.systems())

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--cellopt-root", type=Path, required=True)
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.add_argument("--mesh", type=int, choices=(3, 4), default=3)
    prepare_parser.add_argument("--method", choices=cellopt.METHODS)
    prepare_parser.add_argument("--system", action="append", choices=choices)
    prepare_parser.add_argument("--clean", action="store_true")
    prepare_parser.set_defaults(function=prepare)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--cp2k", type=Path, required=True)
    run_parser.add_argument("--jobs", type=int, default=4)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    run_parser.add_argument("--mesh", type=int, choices=(3, 4), default=3)
    run_parser.add_argument("--method", choices=cellopt.METHODS)
    run_parser.add_argument("--system", action="append", choices=choices)
    run_parser.add_argument("--force", action="store_true")
    run_parser.add_argument("--prune-transients", action="store_true")
    cellopt.add_provenance_arguments(run_parser)
    run_parser.set_defaults(function=run)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    collect_parser.add_argument("--mesh", type=int, choices=(3, 4), required=True)
    collect_parser.add_argument("--molecule-run-root", type=Path, required=True)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.add_argument("--method", choices=cellopt.METHODS)
    collect_parser.add_argument("--allow-incomplete", action="store_true")
    collect_parser.set_defaults(function=collect)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
