#!/usr/bin/env python3
"""Run the diagnostic LC12 g-XTB multi-start branch map for LiH and MgO.

Every scale is evaluated from a cold start.  Two additional, sequential WFN
continuation chains traverse the same grid in ascending and descending order.
The results are diagnostic and cannot be consumed as LC12 production EOS
points without a separate, hash-pinned branch classification.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

import diagnose_gxtb_wfn_hysteresis as wfn
import run_goldzak12_benchmark as base


SCHEMA_VERSION = 1
DEFAULT_PLAN = base.ROOT / "data" / "gxtb_multistart_plan.json"
CAMPAIGN_STATES = ("production_ready", "qualification_pending", "validation_in_progress")


def canonical_sha256(payload: object) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    """Replace a text artifact atomically after flushing its temporary file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextlib.contextmanager
def campaign_lock(root: Path):
    """Hold one nonblocking lock across all writes to a campaign root."""
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / ".campaign.lock"
    with lock_path.open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"multi-start campaign is already locked: {root}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def write_job_stamp_atomic(
    result: Path,
    signature: Mapping[str, object],
    *,
    completed: bool,
    return_code: int,
) -> None:
    payload = dict(signature)
    payload.update({"completed": completed, "return_code": return_code})
    atomic_write_text(
        base.job_stamp_path(result),
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )


def load_plan(path: Path) -> tuple[dict[str, object], str]:
    path = path.resolve(strict=True)
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"unsupported multi-start plan schema in {path}")
    if payload.get("policy_id") != "lc12-gxtb-multistart-v1":
        raise ValueError(f"unknown multi-start policy in {path}")
    execution = payload.get("execution")
    contract = payload.get("hamiltonian_contract")
    solids = payload.get("solids")
    if not isinstance(execution, Mapping) or execution.get("mesh") != "k444":
        raise ValueError("the LC12 multi-start map is restricted to k444")
    if execution.get("continuation_chains") != ["ascending", "descending"]:
        raise ValueError("multi-start plan must retain both ordered continuation chains")
    if execution.get("cold_starts_parallel") is not True:
        raise ValueError("multi-start plan must retain an independent cold start at every scale")
    if execution.get("production_eligible") is not False:
        raise ValueError("multi-start plan must remain diagnostic-only")
    if int(execution.get("omp_threads_per_job", 0)) != 1:
        raise ValueError("multi-start plan must pin one OpenMP thread per candidate")
    if not isinstance(contract, Mapping) or contract != {
        "accuracy": 0.05,
        "eps_scf": 1.0e-9,
        "full_grid": False,
        "kpoint_scheme": "MACDONALD 4 4 4 0.125 0.125 0.125",
        "method": "GXTB",
        "scc_mixer": "TBLITE",
        "symmetry": True,
        "symmetry_backend": "SPGLIB",
        "symmetry_reduction_method": "SPGLIB",
    }:
        raise ValueError(
            "multi-start Hamiltonian contract was changed without a policy version bump"
        )
    if not isinstance(solids, Mapping) or set(solids) != {"LiH", "MgO"}:
        raise ValueError("multi-start plan must contain exactly LiH and MgO")
    for solid, record in solids.items():
        if not isinstance(record, Mapping):
            raise ValueError(f"invalid plan record for {solid}")
        scales = tuple(float(value) for value in record.get("scales", ()))
        if len(scales) < 3 or scales != tuple(sorted(set(scales))):
            raise ValueError(f"{solid} scales must be sorted and unique")
    required = str(payload.get("required_cp2k_ancestor", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", required):
        raise ValueError("multi-start plan needs a full required_cp2k_ancestor hash")
    return payload, canonical_sha256(payload)


def require_cp2k_ancestor(cp2k_source: Path, revision: str) -> None:
    """Reject the pre-#5582 build even if somebody writes a fresh manifest for it."""
    source = cp2k_source.resolve(strict=True)
    present = subprocess.run(
        ["git", "-C", str(source), "cat-file", "-e", f"{revision}^{{commit}}"],
        check=False,
        capture_output=True,
        text=True,
    )
    if present.returncode != 0:
        raise ValueError(
            f"CP2K source does not contain required upstream #5582 commit {revision}"
        )
    ancestor = subprocess.run(
        ["git", "-C", str(source), "merge-base", "--is-ancestor", revision, "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    if ancestor.returncode != 0:
        raise ValueError(
            f"CP2K source HEAD is not descended from required upstream #5582 commit {revision}"
        )


def scale_tag(scale: float) -> str:
    return f"s{scale:.5f}".replace(".", "p")


def multistart_input(
    ref: base.Reference,
    scale: float,
    project: str,
    restart: Path | None,
) -> str:
    text = wfn.diagnostic_input(ref, "k444", scale, project, restart)
    text = text.replace(
        "        &RESTART ON\n          BACKUP_COPIES 0",
        "        &RESTART ON\n          LOG_PRINT_KEY T\n          BACKUP_COPIES 0",
        1,
    )
    validate_multistart_input(text, restart is not None)
    return text


def validate_multistart_input(text: str, continuation: bool) -> None:
    base.validate_method_input(text, "GXTB")
    required = (
        r"^\s*METHOD\s+GXTB\s*$",
        r"^\s*SCC_MIXER\s+TBLITE\s*$",
        r"^\s*ACCURACY\s+0\.05\s*$",
        r"^\s*EPS_SCF\s+1\.0E-9\s*$",
        r"^\s*SCHEME\s+MACDONALD\s+4\s+4\s+4\s+0\.125\s+0\.125\s+0\.125\s*$",
        r"^\s*SYMMETRY\s+T\s*$",
        r"^\s*FULL_GRID\s+F\s*$",
        r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$",
        r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$",
    )
    missing = [pattern for pattern in required if not re.search(pattern, text, re.I | re.M)]
    if missing:
        raise ValueError("malformed multi-start input; missing contract: " + ", ".join(missing))
    if re.search(r"^\s*FULL_GRID\s+T\s*$", text, re.I | re.M):
        raise ValueError("FULL_GRID T is forbidden in the LC12 multi-start map")
    if continuation:
        continuation_contract = (
            r"^\s*SCF_GUESS\s+RESTART\s*$",
            r"^\s*WFN_RESTART_FILE_NAME\s+\S+\s*$",
            r"^\s*LOG_PRINT_KEY\s+T\s*$",
        )
        if any(not re.search(pattern, text, re.I | re.M) for pattern in continuation_contract):
            raise ValueError("continuation input lacks explicit WFN restart evidence")
    elif re.search(r"^\s*WFN_RESTART_FILE_NAME\b", text, re.I | re.M):
        raise ValueError("cold-start input unexpectedly names a WFN restart")


def final_native_mixer_residual(text: str) -> tuple[float | None, str | None]:
    records = re.findall(
        r"^\s*(\d+)\s+GXTB-(Raw|Simple|FDIIS)\s+\S+\s+\S+\s+"
        r"([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s*$",
        text,
        flags=re.M,
    )
    if not records:
        return None, None
    _, label, residual = records[-1]
    return float(residual), label.upper()


def numerical_gates(
    output: Path,
    restart_output: Path,
    *,
    continuation: bool,
) -> tuple[dict[str, object], dict[str, object]]:
    text = output.read_text(errors="ignore") if output.exists() else ""
    obs = wfn.observables(output) if output.exists() else {}
    residual, mixer_label = final_native_mixer_residual(text)
    energy = obs.get("total_energy_extrapolated_t0_hartree")
    gates: dict[str, object] = {
        "normal_completion_without_abort": base.output_ok(output),
        "scf_converged": "SCF run converged in" in text,
        "final_native_mixer_is_fdiis": mixer_label == "FDIIS",
        "final_native_mixer_residual_at_most_eps_scf": (
            residual is not None and math.isfinite(residual) and residual <= 1.0e-9
        ),
        "finite_extrapolated_t0_energy": isinstance(energy, float) and math.isfinite(energy),
        "eight_mulliken_atoms": len(obs.get("mulliken_atoms", [])) == 8,
        "mo_occupations_printed": obs.get("mo_occupations_printed") is True,
        "restart_written": restart_output.is_file(),
        "explicit_wfn_restart_read": (
            obs.get("wfn_restart_read_explicit_log") is True if continuation else True
        ),
        "no_symmetry_or_scc_abort": not bool(
            re.search(r"(?is)\[ABORT\].{0,500}(?:g-?xtb|symmetr|SCC|FDIIS)", text)
        ),
    }
    diagnostics = {
        "final_native_mixer_label": mixer_label,
        "final_native_mixer_residual": residual,
        "observables": obs,
    }
    return gates, diagnostics


def artifact(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        return None
    return {"path": str(path.resolve()), "sha256": base.sha256(path)}


def recorded_artifact_issue(
    record: object,
    expected_path: Path,
    *,
    required: bool,
    label: str,
) -> str | None:
    if record is None:
        if required or expected_path.exists():
            return f"missing {label} artifact record for {expected_path}"
        return None
    if not isinstance(record, Mapping):
        return f"invalid {label} artifact record for {expected_path}"
    recorded_path = Path(str(record.get("path", "")))
    if recorded_path.resolve() != expected_path.resolve():
        return f"{label} artifact path collision: {recorded_path} != {expected_path}"
    if not expected_path.is_file():
        return f"missing {label} artifact {expected_path}"
    if record.get("sha256") != base.sha256(expected_path):
        return f"{label} artifact hash mismatch for {expected_path}"
    return None


def recorded_stamp_issue(
    output: Path,
    signature: Mapping[str, object],
    *,
    completed: bool,
    return_code: int,
) -> str | None:
    stamp = base.job_stamp_path(output)
    if not stamp.is_file():
        return f"missing candidate job stamp {stamp}"
    try:
        record = json.loads(stamp.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid candidate job stamp {stamp}: {exc}"
    expected = dict(signature)
    expected.update({"completed": completed, "return_code": return_code})
    if record != expected:
        return f"candidate job stamp does not exactly match {stamp}"
    return None


def existing_candidate_issue(
    recorded: Mapping[str, object],
    signature: Mapping[str, object],
    *,
    inp: Path,
    out: Path,
    restart_out: Path,
    parent_manifest: Path | None,
    parent_restart: Path | None,
) -> str | None:
    if recorded.get("job_signature") != signature:
        return "existing candidate job signature differs from the requested calculation"
    completed = recorded.get("completed") is True
    return_code = recorded.get("return_code")
    if not isinstance(return_code, int):
        return "existing candidate manifest has no integer return code"
    issues = [
        recorded_artifact_issue(
            recorded.get("input"), inp, required=True, label="input"
        ),
        recorded_artifact_issue(
            recorded.get("output"), out, required=completed, label="output"
        ),
        recorded_artifact_issue(
            recorded.get("wfn_restart"),
            restart_out,
            required=completed,
            label="WFN restart",
        ),
        recorded_stamp_issue(
            out,
            signature,
            completed=completed,
            return_code=return_code,
        ),
    ]
    parent_records = (
        ("parent_candidate_manifest", parent_manifest, "parent candidate manifest"),
        ("parent_wfn_restart", parent_restart, "parent WFN restart"),
    )
    for key, expected, label in parent_records:
        record = recorded.get(key)
        if expected is None:
            if record is not None:
                issues.append(f"cold candidate unexpectedly records a {label}")
        else:
            issues.append(
                recorded_artifact_issue(record, expected, required=True, label=label)
            )
    issue = next((value for value in issues if value), None)
    if issue:
        return issue
    if completed:
        if not base.job_stamp_matches(out, dict(signature)):
            return f"completed candidate does not have a reusable job stamp: {out}"
        gates, _ = numerical_gates(
            out,
            restart_out,
            continuation=parent_restart is not None,
        )
        if not all(value is True for value in gates.values()):
            return f"completed candidate no longer passes its numerical gates: {out}"
    return None


def archive_failed_attempt(
    inp: Path,
    out: Path,
    restart_out: Path,
    manifest_path: Path,
) -> Path:
    archive_root = manifest_path.parent / "attempt_archive"
    index = 1
    while (archive_root / f"attempt_{index:03d}").exists():
        index += 1
    destination = archive_root / f"attempt_{index:03d}"
    destination.mkdir(parents=True, exist_ok=False)
    sources = (
        inp,
        out,
        restart_out,
        base.job_stamp_path(out),
        manifest_path,
        manifest_path.parent / "mainLog.out",
    )
    copied = []
    for source in sources:
        if not source.is_file():
            continue
        target = destination / source.name
        shutil.copy2(source, target)
        copied.append(
            {
                "source": str(source.resolve()),
                "source_sha256": base.sha256(source),
                "archive": str(target.resolve()),
                "archive_sha256": base.sha256(target),
            }
        )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "archive": "lc12_gxtb_multistart_failed_attempt",
        "archived_at_utc": datetime.now(timezone.utc).isoformat(),
        "production_eligible": False,
        "artifacts": copied,
    }
    atomic_write_text(
        destination / "archive_manifest.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    for source in sources[1:]:
        if source.is_file():
            source.unlink()
    return destination


def candidate_paths(
    root: Path, solid: str, mode: str, scale: float
) -> tuple[Path, Path, Path, Path]:
    run_dir = root / solid / "k444" / mode / scale_tag(scale)
    project = f"{solid}_GXTB_multistart_{mode}_k444_{scale_tag(scale)}"
    return (
        run_dir / f"{project}.inp",
        run_dir / f"{project}.out",
        run_dir / f"{project}-RESTART.kp",
        run_dir / "candidate_manifest.json",
    )


def run_candidate(
    *,
    root: Path,
    ref: base.Reference,
    scale: float,
    mode: str,
    restart: Path | None,
    parent_manifest: Path | None,
    cp2k: Path,
    campaign_identity: dict[str, object],
    campaign_state: str,
    plan: dict[str, object],
    plan_sha256: str,
    threads: int,
    retry_failed: bool,
) -> dict[str, object]:
    inp, out, restart_out, manifest_path = candidate_paths(root, ref.solid, mode, scale)
    continuation = restart is not None or parent_manifest is not None
    if continuation != (restart is not None and parent_manifest is not None):
        raise RuntimeError("a continuation needs both its parent manifest and WFN restart")
    if continuation and (not restart.is_file() or not parent_manifest.is_file()):
        raise RuntimeError("a continuation parent artifact is missing")
    if mode == "cold" and continuation:
        raise RuntimeError("cold candidate cannot consume a continuation restart")
    if mode != "cold" and not continuation:
        raise RuntimeError("ascending/descending candidates require a continuation restart")
    parent_restart = artifact(restart) if restart is not None else None
    parent_record = artifact(parent_manifest) if parent_manifest is not None else None
    project = inp.stem
    input_text = multistart_input(ref, scale, project, restart)
    if manifest_path.is_file():
        if not inp.is_file() or inp.read_text() != input_text:
            raise RuntimeError(f"refusing changed input beside {manifest_path}")
    else:
        collisions = [
            path
            for path in (out, restart_out, base.job_stamp_path(out))
            if path.exists()
        ]
        if collisions:
            raise RuntimeError(
                "refusing orphaned candidate artifacts without a manifest: "
                + ", ".join(str(path) for path in collisions)
            )
        if inp.exists() and (not inp.is_file() or inp.read_text() != input_text):
            raise RuntimeError(f"refusing noncanonical orphaned input {inp}")
        atomic_write_text(inp, input_text)
    signature = base.job_signature(
        cp2k,
        inp,
        command_contract={
            "driver": "cp2k",
            "diagnostic": "lc12_gxtb_multistart",
            "policy_id": plan["policy_id"],
            "plan_sha256": plan_sha256,
            "solid": ref.solid,
            "mesh": "k444",
            "scale": scale,
            "mode": mode,
            "parent_restart_sha256": (
                parent_restart["sha256"] if parent_restart is not None else None
            ),
            "parent_manifest_sha256": (
                parent_record["sha256"] if parent_record is not None else None
            ),
            "omp_threads": threads,
            "production_eligible": False,
        },
        campaign_fingerprint=campaign_identity,
    )
    if manifest_path.is_file():
        try:
            recorded = json.loads(manifest_path.read_text())
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"invalid candidate manifest {manifest_path}: {exc}") from exc
        declarations = {
            "schema_version": recorded.get("schema_version"),
            "diagnostic": recorded.get("diagnostic"),
            "production_eligible": recorded.get("production_eligible"),
            "solid": recorded.get("solid"),
            "mesh": recorded.get("mesh"),
            "scale": recorded.get("scale"),
            "mode": recorded.get("mode"),
            "plan_sha256": recorded.get("plan_sha256"),
            "campaign_state_at_execution": recorded.get("campaign_state_at_execution"),
            "campaign_identity": recorded.get("campaign_identity"),
        }
        expected_declarations = {
            "schema_version": SCHEMA_VERSION,
            "diagnostic": "lc12_gxtb_multistart",
            "production_eligible": False,
            "solid": ref.solid,
            "mesh": "k444",
            "scale": scale,
            "mode": mode,
            "plan_sha256": plan_sha256,
            "campaign_state_at_execution": campaign_state,
            "campaign_identity": campaign_identity,
        }
        if declarations != expected_declarations:
            raise RuntimeError(f"candidate manifest declarations changed: {manifest_path}")
        issue = existing_candidate_issue(
            recorded,
            signature,
            inp=inp,
            out=out,
            restart_out=restart_out,
            parent_manifest=parent_manifest,
            parent_restart=restart,
        )
        if issue:
            raise RuntimeError(f"refusing stale candidate {manifest_path}: {issue}")
        if recorded.get("completed") is True or not retry_failed:
            return dict(recorded)
        archive_failed_attempt(inp, out, restart_out, manifest_path)
        atomic_write_text(inp, input_text)

    if restart_out.exists():
        restart_out.unlink()
    return_code = base.run_cp2k(cp2k, inp, out, threads)
    gates, diagnostics = numerical_gates(
        out,
        restart_out,
        continuation=restart is not None,
    )
    completed = all(value is True for value in gates.values())
    write_job_stamp_atomic(
        out,
        signature,
        completed=completed,
        return_code=return_code,
    )
    payload: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "diagnostic": "lc12_gxtb_multistart",
        "production_eligible": False,
        "completed": completed,
        "return_code": return_code,
        "solid": ref.solid,
        "mesh": "k444",
        "scale": scale,
        "lattice_a_A": ref.a_exp * scale,
        "mode": mode,
        "plan_sha256": plan_sha256,
        "policy_id": plan["policy_id"],
        "campaign_state_at_execution": campaign_state,
        "campaign_identity": campaign_identity,
        "job_signature": signature,
        "parent_candidate_manifest": parent_record,
        "parent_wfn_restart": parent_restart,
        "input": artifact(inp),
        "output": artifact(out),
        "wfn_restart": artifact(restart_out),
        "numerical_gates": gates,
        "diagnostics": diagnostics,
        "selection_status": "unclassified; never an EOS point",
    }
    atomic_write_text(
        manifest_path,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    return payload


def run_solid(
    *,
    root: Path,
    ref: base.Reference,
    scales: tuple[float, ...],
    cp2k: Path,
    campaign_identity: dict[str, object],
    campaign_state: str,
    plan: dict[str, object],
    plan_sha256: str,
    threads: int,
    workers: int,
    retry_failed: bool,
) -> dict[str, object]:
    common = {
        "root": root,
        "ref": ref,
        "cp2k": cp2k,
        "campaign_identity": campaign_identity,
        "campaign_state": campaign_state,
        "plan": plan,
        "plan_sha256": plan_sha256,
        "threads": threads,
        "retry_failed": retry_failed,
    }
    cold: dict[float, dict[str, object]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        future_to_scale = {
            pool.submit(
                run_candidate,
                scale=scale,
                mode="cold",
                restart=None,
                parent_manifest=None,
                **common,
            ): scale
            for scale in scales
        }
        for future in concurrent.futures.as_completed(future_to_scale):
            scale = future_to_scale[future]
            cold[scale] = future.result()

    chains: dict[str, list[dict[str, object]]] = {"ascending": [], "descending": []}
    for direction, ordered in (("ascending", scales), ("descending", tuple(reversed(scales)))):
        seed_scale = ordered[0]
        seed = cold[seed_scale]
        if seed.get("completed") is not True:
            chains[direction].append(
                {"scale": seed_scale, "status": "chain_not_started", "reason": "cold seed failed"}
            )
            continue
        _, _, restart, parent_manifest = candidate_paths(root, ref.solid, "cold", seed_scale)
        for scale in ordered[1:]:
            candidate = run_candidate(
                scale=scale,
                mode=direction,
                restart=restart,
                parent_manifest=parent_manifest,
                **common,
            )
            chains[direction].append(candidate)
            if candidate.get("completed") is not True:
                break
            _, _, restart, parent_manifest = candidate_paths(
                root, ref.solid, direction, scale
            )
    return {
        "solid": ref.solid,
        "requested_scales": list(scales),
        "cold_completed": sum(item.get("completed") is True for item in cold.values()),
        "cold_requested": len(scales),
        "ascending_completed": sum(item.get("completed") is True for item in chains["ascending"]),
        "ascending_requested": len(scales) - 1,
        "descending_completed": sum(item.get("completed") is True for item in chains["descending"]),
        "descending_requested": len(scales) - 1,
        "chains": chains,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--campaign-manifest", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--campaign-state", choices=CAMPAIGN_STATES, default="production_ready")
    parser.add_argument("--solid", action="append", choices=("LiH", "MgO"))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--cold-workers", type=int, default=1)
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args()
    if args.threads < 1 or args.cold_workers < 1:
        parser.error("--threads and --cold-workers must be positive")
    if args.solid and len(args.solid) != len(set(args.solid)):
        parser.error("duplicate --solid selections are not allowed")
    try:
        plan, plan_sha256 = load_plan(args.plan)
        execution = plan["execution"]
        assert isinstance(execution, Mapping)
        if args.threads != int(execution["omp_threads_per_job"]):
            raise ValueError(
                "--threads must match the versioned omp_threads_per_job contract "
                f"({execution['omp_threads_per_job']})"
            )
        require_cp2k_ancestor(args.cp2k_source, str(plan["required_cp2k_ancestor"]))
        campaign_identity, paths = base.validated_gxtb_campaign_from_manifest(
            args.campaign_manifest,
            args.cp2k_source,
            args.save_tblite_source,
            allowed_campaign_states=(args.campaign_state,),
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    root = (
        base.ROOT
        / "runs"
        / "gxtb_multistart_branches"
        / f"{campaign_identity['fingerprint_sha256'][:12]}-{plan_sha256[:12]}"
    )
    refs = {ref.solid: ref for ref in base.REFERENCES}
    selected = tuple(args.solid or ("LiH", "MgO"))
    solids = plan["solids"]
    assert isinstance(solids, Mapping)
    manifest = root / "campaign_manifest.json"
    try:
        with campaign_lock(root):
            snapshots = {
                "plan": (args.plan.resolve(strict=True), root / "plan_snapshot.json"),
                "build_manifest": (
                    args.campaign_manifest.resolve(strict=True),
                    root / "build_manifest_snapshot.json",
                ),
            }
            for source, snapshot in snapshots.values():
                source_text = source.read_text()
                if snapshot.exists() and (
                    not snapshot.is_file() or snapshot.read_text() != source_text
                ):
                    raise RuntimeError(f"frozen campaign snapshot collision: {snapshot}")
                if not snapshot.exists():
                    atomic_write_text(snapshot, source_text)

            existing_systems: dict[str, dict[str, object]] = {}
            if manifest.is_file():
                previous = json.loads(manifest.read_text())
                declarations = {
                    "schema_version": previous.get("schema_version"),
                    "diagnostic": previous.get("diagnostic"),
                    "production_eligible": previous.get("production_eligible"),
                    "campaign_state_at_execution": previous.get(
                        "campaign_state_at_execution"
                    ),
                    "campaign_identity": previous.get("campaign_identity"),
                    "plan_sha256": previous.get("plan_sha256"),
                    "required_cp2k_ancestor": previous.get("required_cp2k_ancestor"),
                }
                expected = {
                    "schema_version": SCHEMA_VERSION,
                    "diagnostic": "lc12_gxtb_multistart_campaign",
                    "production_eligible": False,
                    "campaign_state_at_execution": args.campaign_state,
                    "campaign_identity": campaign_identity,
                    "plan_sha256": plan_sha256,
                    "required_cp2k_ancestor": plan["required_cp2k_ancestor"],
                }
                if declarations != expected:
                    raise RuntimeError(f"existing campaign manifest collision: {manifest}")
                snapshot_issues = [
                    recorded_artifact_issue(
                        previous.get("plan"),
                        snapshots["plan"][1],
                        required=True,
                        label="plan snapshot",
                    ),
                    recorded_artifact_issue(
                        previous.get("build_manifest"),
                        snapshots["build_manifest"][1],
                        required=True,
                        label="build-manifest snapshot",
                    ),
                ]
                snapshot_issue = next(
                    (value for value in snapshot_issues if value), None
                )
                if snapshot_issue:
                    raise RuntimeError(snapshot_issue)
                for item in previous.get("systems", []):
                    if not isinstance(item, dict) or item.get("solid") not in solids:
                        raise RuntimeError(f"invalid prior system summary in {manifest}")
                    solid = str(item["solid"])
                    if solid in existing_systems:
                        raise RuntimeError(f"duplicate prior system summary for {solid}")
                    existing_systems[solid] = item

            for solid in selected:
                record = solids[solid]
                assert isinstance(record, Mapping)
                scales = tuple(float(value) for value in record["scales"])
                existing_systems[solid] = run_solid(
                    root=root,
                    ref=refs[solid],
                    scales=scales,
                    cp2k=paths["cp2k"],
                    campaign_identity=campaign_identity,
                    campaign_state=args.campaign_state,
                    plan=plan,
                    plan_sha256=plan_sha256,
                    threads=args.threads,
                    workers=args.cold_workers,
                    retry_failed=args.retry_failed,
                )
            summaries = [
                existing_systems[solid]
                for solid in ("LiH", "MgO")
                if solid in existing_systems
            ]
            complete = set(existing_systems) == {"LiH", "MgO"} and all(
                item["cold_completed"] == item["cold_requested"]
                and item["ascending_completed"] == item["ascending_requested"]
                and item["descending_completed"] == item["descending_requested"]
                for item in summaries
            )
            payload = {
                "schema_version": SCHEMA_VERSION,
                "diagnostic": "lc12_gxtb_multistart_campaign",
                "production_eligible": False,
                "completed": complete,
                "campaign_state_at_execution": args.campaign_state,
                "campaign_identity": campaign_identity,
                "plan": artifact(snapshots["plan"][1]),
                "plan_source": artifact(snapshots["plan"][0]),
                "plan_sha256": plan_sha256,
                "build_manifest": artifact(snapshots["build_manifest"][1]),
                "build_manifest_source": artifact(snapshots["build_manifest"][0]),
                "required_cp2k_ancestor": plan["required_cp2k_ancestor"],
                "systems": summaries,
            }
            atomic_write_text(
                manifest,
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
            )
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(f"completed={complete} manifest={manifest}", flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
