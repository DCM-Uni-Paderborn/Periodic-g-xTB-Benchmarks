#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import concurrent.futures
import fcntl
import hashlib
import json
import math
import os
import queue
import re
import signal
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import benchmark_execution as benchmark_execution  # noqa: E402


PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]
METHODS = ["GFN1", "GFN2", "GXTB"]
# Frozen production core.  Dense meshes are explicit convergence extensions
# and must never expand the default 78-job campaign implicitly.
MESHES = ["gamma", "k111", "k222", "k333", "k444", "k555"]
DENSE_EXTENSION_MESHES = [
    "k666",
    "k777",
    "k888",
    "k999",
    "k101010",
    "k111111",
    "k121212",
    "k131313",
]
SUPPORTED_MESHES = [*MESHES, *DENSE_EXTENSION_MESHES]
GXTB_PROTOCOL_ID = "dmc13-gxtb-spglib-reduced-v1"
GXTB_INPUT_DIRECTORY = "gxtb_spglib_inputs"
GXTB_RUN_DIRECTORY = "runs_gxtb_spglib"
GXTB_ANALYSIS_PREFIX = "gxtb_spglib"
DEFAULT_CAMPAIGN_MANIFEST = (
    Path(__file__).resolve().parents[1]
    / "campaigns"
    / "gxtb-pbc-v1-20260714"
    / "build_manifest.json"
)
MESH_SCHEMES = {
    "k111": "MACDONALD 1 1 1 0.0 0.0 0.0",
    "k222": "MACDONALD 2 2 2 0.25 0.25 0.25",
    "k333": "MACDONALD 3 3 3 0.0 0.0 0.0",
    "k444": "MACDONALD 4 4 4 0.375 0.375 0.375",
    "k555": "MACDONALD 5 5 5 0.0 0.0 0.0",
    "k666": (
        "MACDONALD 6 6 6 0.4166666666666667 "
        "0.4166666666666667 0.4166666666666667"
    ),
    "k777": "MACDONALD 7 7 7 0.0 0.0 0.0",
    "k888": "MACDONALD 8 8 8 0.4375 0.4375 0.4375",
    "k999": "MACDONALD 9 9 9 0.0 0.0 0.0",
    "k101010": "MACDONALD 10 10 10 0.45 0.45 0.45",
    "k111111": "MACDONALD 11 11 11 0.0 0.0 0.0",
    "k121212": (
        "MACDONALD 12 12 12 0.4583333333333333 "
        "0.4583333333333333 0.4583333333333333"
    ),
    "k131313": "MACDONALD 13 13 13 0.0 0.0 0.0",
}

BUILD_IDENTITY_FIELDS = (
    "cp2k_sha256",
    "cp2k_library_sha256",
    "tblite_static_library_sha256",
    "cp2k_source_revision",
    "tblite_source_revision",
)
VALIDATION_INDEX_SCHEMA_VERSION = 2
EXECUTION_BUILD_MANIFEST_SCHEMA_VERSION = 1
QUALIFICATION_EVIDENCE_SCHEMA_VERSION = 3
MAX_TOTAL_ENERGY_TOLERANCE_HARTREE = 1.0e-10
MAX_RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O = 1.0e-3
HARTREE_TO_KJMOL = 2625.499638
OMP_SCHEDULE = "static"
OMP_DYNAMIC = "FALSE"
OMP_WAIT_POLICY = "PASSIVE"
BLAS_THREAD_ENVIRONMENT = {
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


@dataclass(frozen=True)
class Job:
    mesh: str
    method: str
    phase: str
    input_path: Path
    run_dir: Path
    output_name: str


@dataclass(frozen=True)
class ProductionIdentity:
    campaign_id: str
    cp2k: Path
    cp2k_sha256: str
    cp2k_library: Path
    cp2k_library_sha256: str
    tblite_static_library: Path
    tblite_static_library_sha256: str
    cp2k_source_revision: str
    tblite_source_revision: str


def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def thread_configuration(
    threads_per_job: int, mpi_ranks_per_job: int = 1
) -> dict[str, object]:
    if threads_per_job < 1:
        raise ValueError("threads_per_job must be positive")
    if mpi_ranks_per_job < 1:
        raise ValueError("mpi_ranks_per_job must be positive")
    return {
        "mpi_ranks_per_job": mpi_ranks_per_job,
        "threads_per_job": threads_per_job,
        "omp_num_threads": threads_per_job,
        "omp_schedule": OMP_SCHEDULE,
        "omp_dynamic": False,
        "omp_wait_policy": OMP_WAIT_POLICY,
        "omp_proc_bind": None,
        "omp_places": None,
        "blas_threads": 1,
    }


def job_execution_configuration(
    threads_per_job: int,
    mpi_ranks_per_job: int,
    execution_contract_sha256: str | None,
) -> dict[str, object]:
    """Return the exact execution mode that a reusable stamp must match."""
    configuration = thread_configuration(threads_per_job, mpi_ranks_per_job)
    if execution_contract_sha256 is None:
        configuration.update(
            {
                "execution_mode": "direct",
                "execution_contract_sha256": None,
            }
        )
    else:
        configuration.update(
            {
                "execution_mode": "openmpi_ordered_pe_list",
                "execution_contract_sha256": execution_contract_sha256,
                "omp_proc_bind": "true",
                "omp_places": "cores",
            }
        )
    return configuration


def execution_expectation_from_args(
    args: argparse.Namespace,
) -> tuple[int, int, str | None]:
    pool = getattr(args, "execution_pool", None)
    threads = int(getattr(args, "threads_per_job", 1))
    if pool is None:
        return threads, 1, None
    return threads, int(pool.mpi_ranks_per_job), str(pool.contract_sha256)


def execution_parallelism(
    jobs: int, threads_per_job: int, mpi_ranks_per_job: int = 1
) -> dict[str, object]:
    if jobs < 1:
        raise ValueError("jobs must be positive")
    return {
        "jobs": jobs,
        "nominal_cores": jobs * threads_per_job * mpi_ranks_per_job,
        **thread_configuration(threads_per_job, mpi_ranks_per_job),
    }


def execution_build_identity(identity: ProductionIdentity) -> dict[str, str]:
    return {
        field: str(getattr(identity, field))
        for field in BUILD_IDENTITY_FIELDS
    }


def build_id(build_identity: dict[str, object]) -> str:
    canonical = {
        field: str(build_identity.get(field))
        for field in BUILD_IDENTITY_FIELDS
    }
    content = json.dumps(
        canonical, sort_keys=True, separators=(",", ":")
    ).encode()
    digest = hashlib.sha256(b"dmc13-execution-build-v1\0" + content).hexdigest()
    return f"sha256:{digest}"


def validation_record_key(record: dict[str, object]) -> tuple[str, str]:
    return str(record.get("mesh")), str(record.get("phase"))


def require_unique_selection(values: list[str], label: str) -> list[str]:
    duplicates = sorted({value for value in values if values.count(value) > 1})
    if duplicates:
        raise ValueError(f"duplicate {label} selection(s): {', '.join(duplicates)}")
    return values


def jobs(
    root: Path,
    methods: list[str],
    meshes: list[str] | None = None,
    phases: list[str] | None = None,
    gxtb_input_root: Path | None = None,
    gxtb_run_root: Path | None = None,
) -> list[Job]:
    require_unique_selection(methods, "method")
    selected_meshes = require_unique_selection(meshes or MESHES, "mesh")
    selected_phases = require_unique_selection(phases or PHASES, "phase")
    out: list[Job] = []
    for mesh in selected_meshes:
        for method in methods:
            for phase in selected_phases:
                if method == "GXTB" and gxtb_input_root is not None:
                    input_path = gxtb_input_root / mesh / f"ice_{phase}_{method}_{mesh}.inp"
                else:
                    input_path = root / "kpoint_inputs" / mesh / f"ice_{phase}_{method}_{mesh}.inp"
                if method == "GXTB" and gxtb_run_root is not None:
                    run_dir = gxtb_run_root / mesh / phase
                    output_name = f"ice_{phase}_{method}_{mesh}.out"
                elif mesh == "gamma":
                    run_dir = root / "runs" / method / phase
                    output_name = f"ice_{phase}_{method}.out"
                else:
                    run_dir = root / "runs_kpoints" / mesh / method / phase
                    output_name = f"ice_{phase}_{method}_{mesh}.out"
                out.append(Job(mesh, method, phase, input_path, run_dir, output_name))
    return out


def has_completed(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    return (
        "PROGRAM ENDED" in text
        and "ENERGY| Total FORCE_EVAL" in text
        and "SCF run converged" in text
        and "SCF run NOT converged" not in text
        and "ABORT" not in text
    )


def _gxtb_input_text_contract_errors(text: str, mesh: str) -> list[str]:
    normalised_lines = [
        line.strip().upper()
        for line in text.splitlines()
        if line.strip()
    ]
    lines = set(normalised_lines)
    errors: list[str] = []
    for required in (
        f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}".upper(),
        "METHOD XTB",
        "METHOD GXTB",
        "ACCURACY 0.1",
        "SCC_MIXER TBLITE",
        "ITERATIONS 300",
        "EPS_SCF 1.0E-9",
        "METHOD DIRECT_P_MIXING",
        "ALPHA 0.2",
        "CANONICALIZE TRUE",
    ):
        count = normalised_lines.count(required)
        if count == 0:
            errors.append(f"missing {required}")
        elif count != 1:
            errors.append(f"duplicate critical setting {required}")
    for forbidden_method in ("METHOD GFN1", "METHOD GFN2"):
        if forbidden_method in lines:
            errors.append(f"conflicting tblite method {forbidden_method}")

    if mesh == "gamma":
        if "&KPOINTS" in lines:
            errors.append("Gamma production input must use implicit Gamma without &KPOINTS")
    else:
        for required in (
            "&KPOINTS",
            f"SCHEME {MESH_SCHEMES[mesh]}",
            "SYMMETRY T",
            "FULL_GRID F",
            "SYMMETRY_BACKEND SPGLIB",
            "SYMMETRY_REDUCTION_METHOD SPGLIB",
        ):
            count = normalised_lines.count(required)
            if count == 0:
                errors.append(f"missing {required}")
            elif count != 1:
                errors.append(f"duplicate critical setting {required}")
        for forbidden in ("SYMMETRY F", "FULL_GRID T"):
            if forbidden in lines:
                errors.append(f"forbidden legacy setting {forbidden}")
    return errors


def gxtb_input_contract_errors(job: Job) -> list[str]:
    if job.method != "GXTB":
        return []
    if not job.input_path.is_file():
        return ["input file is missing"]
    return _gxtb_input_text_contract_errors(
        job.input_path.read_text(errors="ignore"), job.mesh
    )


def is_legacy_full_grid_input(path: Path) -> bool:
    if not path.is_file():
        return False
    lines = {
        line.strip().upper()
        for line in path.read_text(errors="ignore").splitlines()
        if line.strip()
    }
    return "FULL_GRID T" in lines and "SYMMETRY F" in lines


def stamp_path(job: Job) -> Path:
    return (job.run_dir / job.output_name).with_suffix(".run.json")


def frozen_input_path(job: Job) -> Path:
    return job.run_dir / job.input_path.name


def _scientifically_valid_stamp_payload(
    job: Job,
    identity: ProductionIdentity,
) -> dict[str, object] | None:
    """Return a hash-validated numerical result, independent of timing policy.

    Historical stamps predate the exact rank-binding contract.  Their numerical
    observables remain usable when every scientific artifact and build-identity
    hash revalidates, but they must not be resumed or counted as scaling data.
    ``stamp_valid`` adds that stricter execution-contract gate below.
    """
    output = job.run_dir / job.output_name
    stamp = stamp_path(job)
    frozen_input = frozen_input_path(job)
    if (
        not job.input_path.is_file()
        or not frozen_input.is_file()
        or gxtb_input_contract_errors(job)
        or not has_completed(output)
        or not stamp.is_file()
    ):
        return None
    try:
        payload = json.loads(stamp.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(payload, dict):
        return None
    source_hash = sha256(job.input_path)
    frozen_hash = sha256(frozen_input)
    stamp_input_hash = payload.get("input_sha256")
    identity_build_id = build_id(execution_build_identity(identity))
    scientifically_valid = (
        payload.get("campaign_id") == identity.campaign_id
        and payload.get("method") == job.method
        and payload.get("mesh") == job.mesh
        and payload.get("phase") == job.phase
        and stamp_input_hash == source_hash
        and stamp_input_hash == frozen_hash
        and payload.get("frozen_input", frozen_input.name) == frozen_input.name
        and payload.get("frozen_input_sha256", stamp_input_hash) == frozen_hash
        and payload.get("cp2k_sha256") == identity.cp2k_sha256
        and payload.get("cp2k_library_sha256") == identity.cp2k_library_sha256
        and payload.get("tblite_static_library_sha256")
        == identity.tblite_static_library_sha256
        and payload.get("cp2k_source_revision") == identity.cp2k_source_revision
        and payload.get("tblite_source_revision") == identity.tblite_source_revision
        and payload.get("build_id", identity_build_id) == identity_build_id
        and payload.get("output_sha256") == sha256(output)
        and (
            job.method != "GXTB"
            or (
                payload.get("gxtb_protocol_id") == GXTB_PROTOCOL_ID
                and payload.get("adopted_existing_output") is False
            )
        )
    )
    return payload if scientifically_valid else None


def scientific_stamp_valid(
    job: Job,
    identity: ProductionIdentity,
) -> bool:
    """Whether the numerical artifacts revalidate, including legacy stamps."""
    return _scientifically_valid_stamp_payload(job, identity) is not None


def stamp_valid(
    job: Job,
    identity: ProductionIdentity,
    threads_per_job: int = 1,
    mpi_ranks_per_job: int = 1,
    execution_contract_sha256: str | None = None,
) -> bool:
    """Whether a result is safe to reuse under the exact execution contract."""
    payload = _scientifically_valid_stamp_payload(job, identity)
    if payload is None:
        return False
    expected_execution = job_execution_configuration(
        threads_per_job,
        mpi_ranks_per_job,
        execution_contract_sha256,
    )
    return all(
        payload.get(field) == expected
        for field, expected in expected_execution.items()
    )


def write_stamp(
    job: Job,
    identity: ProductionIdentity,
    threads_per_job: int = 1,
    mpi_ranks_per_job: int = 1,
    execution_contract_sha256: str | None = None,
) -> None:
    output = job.run_dir / job.output_name
    frozen_input = frozen_input_path(job)
    if not frozen_input.is_file():
        raise ValueError(f"missing frozen executed input: {frozen_input}")
    source_hash = sha256(job.input_path)
    frozen_hash = sha256(frozen_input)
    if source_hash != frozen_hash:
        raise ValueError(
            f"source/frozen input mismatch for {job.mesh}/{job.phase}: "
            f"{source_hash} != {frozen_hash}"
        )
    payload = {
        "schema_version": 2,
        "campaign_id": identity.campaign_id,
        "build_id": build_id(execution_build_identity(identity)),
        "method": job.method,
        "mesh": job.mesh,
        "phase": job.phase,
        "input": job.input_path.name,
        "input_sha256": source_hash,
        "frozen_input": frozen_input.name,
        "frozen_input_sha256": frozen_hash,
        "cp2k": str(identity.cp2k),
        "cp2k_sha256": identity.cp2k_sha256,
        "cp2k_library": str(identity.cp2k_library),
        "cp2k_library_sha256": identity.cp2k_library_sha256,
        "tblite_static_library": str(identity.tblite_static_library),
        "tblite_static_library_sha256": identity.tblite_static_library_sha256,
        "cp2k_source_revision": identity.cp2k_source_revision,
        "tblite_source_revision": identity.tblite_source_revision,
        "output": output.name,
        "output_sha256": sha256(output),
        "adopted_existing_output": False,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID if job.method == "GXTB" else None,
        "input_contract_valid": not gxtb_input_contract_errors(job),
        **job_execution_configuration(
            threads_per_job,
            mpi_ranks_per_job,
            execution_contract_sha256,
        ),
    }
    atomic_write_bytes(
        stamp_path(job),
        (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode(),
    )


def archive_stale_file(path: Path) -> Path | None:
    """Move an invalid prior attempt aside before CP2K can append to it."""
    if not path.is_file():
        return None
    digest = sha256(path)[:12]
    archived = path.with_name(f"{path.name}.stale-{digest}")
    serial = 1
    while archived.exists():
        archived = path.with_name(f"{path.name}.stale-{digest}.{serial}")
        serial += 1
    path.replace(archived)
    return archived


def run_job(
    identity: ProductionIdentity,
    job: Job,
    force: bool,
    stop_event: threading.Event | None = None,
    threads_per_job: int = 1,
    execution_pool: benchmark_execution.ExecutionPool | None = None,
) -> tuple[Job, int]:
    mpi_ranks_per_job = (
        execution_pool.mpi_ranks_per_job if execution_pool is not None else 1
    )
    execution_contract_sha256 = (
        str(execution_pool.contract_sha256)
        if execution_pool is not None
        else None
    )
    thread_configuration(threads_per_job, mpi_ranks_per_job)
    if stop_event is not None and stop_event.is_set():
        return job, 130
    contract_errors = gxtb_input_contract_errors(job)
    if contract_errors:
        raise ValueError(
            f"{job.mesh}/{job.phase}: invalid g-xTB production input: "
            + "; ".join(contract_errors)
        )
    output = job.run_dir / job.output_name
    if not force and stamp_valid(
        job,
        identity,
        threads_per_job,
        mpi_ranks_per_job,
        execution_contract_sha256,
    ):
        if execution_pool is None:
            return job, 0
        if execution_pool.record_issue(output, stamp_path(job)) is None:
            return job, 0
    job.run_dir.mkdir(parents=True, exist_ok=True)
    local_input = job.run_dir / job.input_path.name
    for prior in (
        output,
        stamp_path(job),
        benchmark_execution.execution_record_path(output),
        benchmark_execution.launcher_log_path(output),
        job.run_dir / "run.log",
        local_input,
    ):
        archive_stale_file(prior)
    shutil.copy2(job.input_path, local_input)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads_per_job)
    env["OMP_SCHEDULE"] = OMP_SCHEDULE
    env["OMP_DYNAMIC"] = OMP_DYNAMIC
    env["OMP_WAIT_POLICY"] = OMP_WAIT_POLICY
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    env.update(BLAS_THREAD_ENVIRONMENT)
    observation: dict[str, object] | None = None
    if execution_pool is not None:
        if stop_event is not None and stop_event.is_set():
            return job, 130
        returncode, observation = execution_pool.run_cp2k(
            identity.cp2k, local_input, output
        )
    else:
        log_path = job.run_dir / "run.log"
        with log_path.open("w") as log:
            proc = subprocess.Popen(
                [str(identity.cp2k), "-i", local_input.name, "-o", job.output_name],
                cwd=job.run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                text=True,
                start_new_session=True,
            )
            while True:
                try:
                    returncode = proc.wait(timeout=0.2)
                    break
                except subprocess.TimeoutExpired:
                    if stop_event is None or not stop_event.is_set():
                        continue
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        returncode = proc.wait(timeout=5.0)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        returncode = proc.wait()
                    break
    completed = has_completed(output)
    if returncode == 0 and completed:
        if (
            execution_pool is not None
            and (
                observation is None
                or observation.get("runtime_affinity_gate") is not True
            )
        ):
            return job, 97
        write_stamp(
            job,
            identity,
            threads_per_job,
            mpi_ranks_per_job,
            execution_contract_sha256,
        )
        if execution_pool is not None:
            assert observation is not None
            try:
                execution_pool.write_record(output, observation, stamp_path(job))
            except Exception:
                stamp_path(job).unlink(missing_ok=True)
                raise
    elif stamp_path(job).exists():
        stamp_path(job).unlink()
    return job, returncode if returncode != 0 or completed else 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_bytes(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp-{os.getpid()}-{threading.get_ident()}"
    )
    try:
        with temporary.open("wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def acquire_runner_lock(path: Path, owner: dict[str, object]):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        handle.seek(0)
        current_owner = handle.read().strip() or "unknown owner"
        handle.close()
        raise ValueError(f"DMC13 runner is already active: {current_owner}") from error
    handle.seek(0)
    handle.truncate()
    json.dump(owner, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
    return handle


def release_runner_lock(handle) -> None:
    if handle is None or handle.closed:
        return
    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    handle.close()


def command_output(command: list[str], cwd: Path | None = None) -> str:
    process = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    return (process.stdout + process.stderr).strip()


def cp2k_rpaths(cp2k: Path) -> list[Path]:
    process = subprocess.run(
        ["otool", "-l", str(cp2k)],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(f"cannot inspect CP2K RPATHs: {process.stderr.strip()}")
    rpaths: list[Path] = []
    waiting_for_path = False
    for raw_line in process.stdout.splitlines():
        line = raw_line.strip()
        if line == "cmd LC_RPATH":
            waiting_for_path = True
            continue
        if not waiting_for_path or not line.startswith("path "):
            continue
        value = line[5:].split(" (offset ", 1)[0]
        value = value.replace("@loader_path", str(cp2k.parent))
        value = value.replace("@executable_path", str(cp2k.parent))
        rpaths.append(Path(value).resolve())
        waiting_for_path = False
    return rpaths


def resolve_cp2k_library_darwin(cp2k: Path) -> Path:
    process = subprocess.run(
        ["otool", "-L", str(cp2k)],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(f"cannot inspect CP2K dependencies: {process.stderr.strip()}")
    dependencies = [
        line.strip().split(" (", 1)[0]
        for line in process.stdout.splitlines()[1:]
        if line.strip()
    ]
    matches = [
        dependency
        for dependency in dependencies
        if Path(dependency).name.startswith("libcp2k")
        and dependency.endswith(".dylib")
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one loaded libcp2k dylib, found {matches!r}"
        )
    dependency = matches[0]
    candidates: list[Path] = []
    if dependency.startswith("@rpath/"):
        suffix = dependency.removeprefix("@rpath/")
        candidates = [rpath / suffix for rpath in cp2k_rpaths(cp2k)]
    elif dependency.startswith("@loader_path/"):
        candidates = [cp2k.parent / dependency.removeprefix("@loader_path/")]
    elif dependency.startswith("@executable_path/"):
        candidates = [cp2k.parent / dependency.removeprefix("@executable_path/")]
    else:
        candidates = [Path(dependency)]
    existing = [candidate.resolve() for candidate in candidates if candidate.is_file()]
    if len(existing) != 1:
        raise ValueError(
            "could not resolve the uniquely loaded libcp2k dylib; candidates: "
            + ", ".join(str(candidate) for candidate in candidates)
        )
    return existing[0]


def resolve_cp2k_library_linux(cp2k: Path) -> Path:
    process = subprocess.run(
        ["ldd", str(cp2k)],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(
            f"cannot inspect CP2K dependencies with ldd: {process.stderr.strip()}"
        )
    matches: list[Path] = []
    unresolved: list[str] = []
    for raw_line in process.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        left, separator, right = line.partition("=>")
        direct_path = left.strip().split(" (", 1)[0]
        dependency_name = Path(direct_path).name
        if not dependency_name.startswith("libcp2k") or ".so" not in dependency_name:
            continue
        if not separator:
            candidate = Path(direct_path)
            if candidate.is_absolute() and candidate.is_file():
                matches.append(candidate.resolve())
            else:
                unresolved.append(line)
            continue
        if right.strip().startswith("not found"):
            unresolved.append(line)
            continue
        resolved_text = right.strip().split(" (", 1)[0]
        candidate = Path(resolved_text)
        if candidate.is_file():
            matches.append(candidate.resolve())
        else:
            unresolved.append(line)
    if unresolved:
        raise ValueError(
            "could not resolve loaded libcp2k shared object: "
            + "; ".join(unresolved)
        )
    unique_matches = sorted(set(matches), key=str)
    if len(unique_matches) != 1:
        raise ValueError(
            "expected exactly one loaded libcp2k shared object, found "
            f"{[str(path) for path in unique_matches]!r}"
        )
    return unique_matches[0]


def resolve_cp2k_library(cp2k: Path) -> Path:
    if sys.platform == "darwin":
        return resolve_cp2k_library_darwin(cp2k)
    if sys.platform.startswith("linux"):
        return resolve_cp2k_library_linux(cp2k)
    raise ValueError(f"unsupported platform for CP2K library resolution: {sys.platform}")


def embedded_revision(text: str, label: str) -> str:
    match = re.search(rf"{re.escape(label)}\s*([0-9a-fA-F]{{7,40}})", text)
    if match is None:
        raise ValueError(f"cannot find {label.strip()} in build metadata")
    return match.group(1).lower()


def embedded_tblite_revision(
    text: str,
    expected_revision: str | None = None,
) -> str:
    """Read labelled metadata, allowing compiler-pooled expected revisions."""
    label = "tblite source revision:"
    token_pattern = re.compile(r"(?:[0-9a-fA-F]{7,40}|unknown)")
    lines = text.splitlines()
    candidates: list[str] = []
    label_seen = False
    for index, raw_line in enumerate(lines):
        position = raw_line.lower().find(label)
        if position < 0:
            continue
        label_seen = True
        remainder = raw_line[position + len(label) :].strip()
        if remainder and token_pattern.fullmatch(remainder):
            candidates.append(remainder.lower())
            continue
    unique = sorted(set(candidates))
    if len(unique) != 1:
        if unique:
            raise ValueError(
                "ambiguous tblite source revisions in libcp2k metadata: "
                + ", ".join(unique)
            )
        if (
            label_seen
            and expected_revision is not None
            and re.fullmatch(r"[0-9a-f]{40}", expected_revision)
            and re.search(
                rf"(?<![0-9a-fA-F]){re.escape(expected_revision)}"
                r"(?![0-9a-fA-F])",
                text,
            )
        ):
            return expected_revision
        raise ValueError("cannot find tblite source revision in libcp2k metadata")
    return unique[0]


def resolve_source_revision(source: Path, revision: str, label: str) -> str:
    source = source.resolve()
    head = command_output(["git", "rev-parse", "HEAD"], source)
    resolved = command_output(
        ["git", "rev-parse", f"{revision}^{{commit}}"],
        source,
    )
    if not re.fullmatch(r"[0-9a-f]{40}", head) or resolved != head:
        raise ValueError(
            f"{label} embedded revision {revision} does not equal source HEAD {head}"
        )
    return head


def production_identity(
    campaign_id: str,
    cp2k: Path,
    cp2k_library_expected: Path | None,
    tblite_static_library: Path,
    cp2k_source: Path,
    tblite_source: Path,
    tblite_revision_expected: str,
    require_embedded_tblite_revision: bool = False,
) -> ProductionIdentity:
    cp2k = cp2k.resolve()
    detected_library = resolve_cp2k_library(cp2k)
    if (
        cp2k_library_expected is not None
        and cp2k_library_expected.resolve() != detected_library
    ):
        raise ValueError(
            f"--cp2k-library resolves to {cp2k_library_expected.resolve()}, but "
            f"the platform dynamic loader selects {detected_library}"
        )
    tblite_static_library = tblite_static_library.resolve()
    if not tblite_static_library.is_file():
        raise ValueError(f"missing static save_tblite library: {tblite_static_library}")
    cp2k_version = command_output([str(cp2k), "--version"])
    cp2k_revision = embedded_revision(cp2k_version, "Source code revision")
    cp2k_source_revision = resolve_source_revision(
        cp2k_source,
        cp2k_revision,
        "CP2K",
    )
    library_strings = command_output(["strings", str(detected_library)])
    tblite_embedded = embedded_tblite_revision(
        library_strings, tblite_revision_expected
    )
    if tblite_embedded == "unknown" and require_embedded_tblite_revision:
        raise ValueError(
            "qualified execution build must embed the full save_tblite revision"
        )
    if tblite_embedded != "unknown":
        embedded_full = command_output(
            ["git", "rev-parse", f"{tblite_embedded}^{{commit}}"],
            tblite_source,
        )
        if embedded_full != tblite_revision_expected:
            raise ValueError(
                "libcp2k embeds save_tblite revision "
                f"{tblite_embedded}, expected {tblite_revision_expected}"
            )
    tblite_source_revision = resolve_source_revision(
        tblite_source,
        tblite_revision_expected,
        "save_tblite",
    )
    return ProductionIdentity(
        campaign_id=campaign_id,
        cp2k=cp2k,
        cp2k_sha256=sha256(cp2k),
        cp2k_library=detected_library,
        cp2k_library_sha256=sha256(detected_library),
        tblite_static_library=tblite_static_library,
        tblite_static_library_sha256=sha256(tblite_static_library),
        cp2k_source_revision=cp2k_source_revision,
        tblite_source_revision=tblite_source_revision,
    )


def read_campaign_manifest(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read campaign manifest {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"campaign manifest {path} is not a JSON object")
    return payload


def frozen_build_identity_from_manifest(
    manifest: dict[str, object],
) -> dict[str, str]:
    cp2k = manifest.get("cp2k")
    save_tblite = manifest.get("save_tblite")
    if not isinstance(cp2k, dict) or not isinstance(save_tblite, dict):
        raise ValueError("campaign manifest lacks frozen build identity")
    identity = {
        "cp2k_sha256": str(cp2k.get("binary_sha256")),
        "cp2k_library_sha256": str(cp2k.get("loaded_library_sha256")),
        "tblite_static_library_sha256": str(
            save_tblite.get("static_library_sha256")
        ),
        "cp2k_source_revision": str(cp2k.get("revision")),
        "tblite_source_revision": str(save_tblite.get("revision")),
    }
    for field in BUILD_IDENTITY_FIELDS[:3]:
        if not re.fullmatch(r"[0-9a-f]{64}", identity[field]):
            raise ValueError(f"campaign manifest has invalid frozen SHA256: {field}")
    for field in BUILD_IDENTITY_FIELDS[3:]:
        if not re.fullmatch(r"[0-9a-f]{40}", identity[field]):
            raise ValueError(f"campaign manifest has invalid frozen revision: {field}")
    return identity


def validate_campaign_identity(
    identity: ProductionIdentity,
    tblite: Path,
    manifest: dict[str, object],
) -> None:
    cp2k = manifest.get("cp2k")
    save_tblite = manifest.get("save_tblite")
    if not isinstance(cp2k, dict) or not isinstance(save_tblite, dict):
        raise ValueError("campaign manifest lacks cp2k/save_tblite identity blocks")
    checks = {
        "campaign_id": (identity.campaign_id, manifest.get("campaign_id")),
        "cp2k revision": (
            identity.cp2k_source_revision,
            cp2k.get("revision"),
        ),
        "CP2K launcher SHA256": (
            identity.cp2k_sha256,
            cp2k.get("binary_sha256"),
        ),
        "loaded libcp2k SHA256": (
            identity.cp2k_library_sha256,
            cp2k.get("loaded_library_sha256"),
        ),
        "save_tblite revision": (
            identity.tblite_source_revision,
            save_tblite.get("revision"),
        ),
        "save_tblite CLI SHA256": (
            sha256(tblite),
            save_tblite.get("cli_sha256"),
        ),
        "static libtblite SHA256": (
            identity.tblite_static_library_sha256,
            save_tblite.get("static_library_sha256"),
        ),
    }
    mismatches = [
        f"{label}: actual {actual!r}, manifest {expected!r}"
        for label, (actual, expected) in checks.items()
        if actual != expected
    ]
    if mismatches:
        raise ValueError("campaign identity mismatch: " + "; ".join(mismatches))


def _finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"qualification evidence {label} is not numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"qualification evidence {label} is not finite")
    return number


def _read_hashed_evidence_artifact(
    artifact_root: Path,
    path_value: object,
    expected_hash: object,
    label: str,
) -> tuple[Path, bytes]:
    if not isinstance(expected_hash, str) or not re.fullmatch(
        r"[0-9a-f]{64}", expected_hash
    ):
        raise ValueError(f"qualification evidence {label} has invalid SHA256")
    path = _validation_artifact_path(artifact_root, path_value, 2)
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(
            f"qualification evidence {label} cannot be read: {error}"
        ) from error
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError(f"qualification evidence {label} hash mismatch")
    return path, content


def _cp2k_input_water_count(text: str, label: str) -> int:
    in_coordinates = False
    oxygen_count = 0
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        upper = line.upper()
        if upper.startswith("&COORD"):
            in_coordinates = True
            continue
        if in_coordinates and upper.startswith("&END"):
            in_coordinates = False
            continue
        if in_coordinates and line and line.split()[0].upper() == "O":
            oxygen_count += 1
    if oxygen_count <= 0:
        raise ValueError(
            f"qualification evidence {label} input lacks explicit water oxygens"
        )
    return oxygen_count


def _validate_evidence_input(
    path: Path,
    content: bytes,
    *,
    mesh: str,
    phase: str,
    label: str,
) -> int:
    expected_project = f"ice_{phase}_GXTB_{mesh}"
    expected_name = f"{expected_project}.inp"
    if path.name != expected_name:
        raise ValueError(
            f"qualification evidence {label} input filename mismatch"
        )
    text = content.decode(errors="ignore")
    lines = {line.strip() for line in text.splitlines() if line.strip()}
    if f"PROJECT {expected_project}" not in lines:
        raise ValueError(f"qualification evidence {label} project mismatch")
    contract_errors = _gxtb_input_text_contract_errors(text, mesh)
    if contract_errors:
        raise ValueError(
            f"qualification evidence {label} violates input contract: "
            + "; ".join(contract_errors)
        )
    return _cp2k_input_water_count(text, label)


def _evidence_total_energy(
    content: bytes,
    label: str,
    *,
    expected_project: str,
    expected_source_revision: str,
    expected_tblite_source_revision: str,
    allow_unknown_tblite_revision: bool,
) -> float:
    text = content.decode(errors="ignore")
    if (
        "PROGRAM ENDED" not in text
        or "ENERGY| Total FORCE_EVAL" not in text
        or "SCF run converged" not in text
        or "SCF run NOT converged" in text
        or "ABORT" in text
    ):
        raise ValueError(f"qualification evidence {label} is not a completed output")
    expected_input = f"{expected_project}.inp"
    input_matches = re.findall(
        r"^\s*CP2K\|\s+Input file name\s+(\S+)\s*$", text, re.MULTILINE
    )
    project_matches = re.findall(
        r"^\s*GLOBAL\|\s+Project name\s+(\S+)\s*$", text, re.MULTILINE
    )
    revision_matches = re.findall(
        r"^\s*CP2K\|\s+source code revision number:\s*"
        r"([0-9a-fA-F]{7,40})\s*$",
        text,
        re.MULTILINE,
    )
    if len(input_matches) != 1 or input_matches[0] != expected_input:
        raise ValueError(f"qualification evidence {label} input header mismatch")
    if len(project_matches) != 1 or project_matches[0] != expected_project:
        raise ValueError(f"qualification evidence {label} project header mismatch")
    if "tblite_gxtb" not in text.lower():
        raise ValueError(f"qualification evidence {label} lacks tblite_gxtb header")
    if (
        len(revision_matches) != 1
        or not expected_source_revision.startswith(revision_matches[0].lower())
    ):
        raise ValueError(f"qualification evidence {label} source revision mismatch")
    tblite_revisions = re.findall(
        r"^\s*tblite source revision:\s*(\S+)\s*$",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if len(tblite_revisions) != 1:
        raise ValueError(
            f"qualification evidence {label} must contain exactly one tblite "
            "source revision"
        )
    tblite_revision = tblite_revisions[0].lower()
    if tblite_revision != expected_tblite_source_revision.lower() and not (
        allow_unknown_tblite_revision and tblite_revision == "unknown"
    ):
        raise ValueError(
            f"qualification evidence {label} tblite source revision mismatch"
        )
    energy_lines = [
        line for line in text.splitlines() if "ENERGY| Total FORCE_EVAL" in line
    ]
    if len(energy_lines) != 1:
        raise ValueError(
            f"qualification evidence {label} must contain exactly one total energy"
        )
    try:
        energy = float(energy_lines[0].split()[-1])
    except (ValueError, IndexError) as error:
        raise ValueError(
            f"qualification evidence {label} has an invalid energy"
        ) from error
    if not math.isfinite(energy):
        raise ValueError(f"qualification evidence {label} lacks a finite energy")
    return energy


def _evidence_execution_environment(
    content: bytes, label: str
) -> tuple[str, str, str, str]:
    text = content.decode(errors="ignore")
    fields = {
        "start time": re.findall(
            r"PROGRAM STARTED AT\s+(.+?)\s*$", text, re.MULTILINE
        ),
        "host": re.findall(
            r"PROGRAM STARTED ON\s+(\S+)\s*$", text, re.MULTILINE
        ),
        "compiled host": re.findall(
            r"^\s*CP2K\|\s+Program compiled on\s*(.*?)\s*$",
            text,
            re.MULTILINE,
        ),
        "platform": re.findall(
            r"^\s*CP2K\|\s+Program compiled for\s+(\S+)\s*$",
            text,
            re.MULTILINE,
        ),
    }
    for field, matches in fields.items():
        if len(matches) != 1 or (
            field != "compiled host" and not matches[0].strip()
        ):
            raise ValueError(
                f"qualification evidence {label} lacks unique execution {field}"
            )
    return (
        fields["host"][0].strip(),
        fields["compiled host"][0].strip(),
        fields["platform"][0].strip(),
        fields["start time"][0].strip(),
    )


def _validate_qualification_run_stamp(
    content: bytes,
    label: str,
    *,
    campaign_id: str,
    identity: dict[str, object],
    identity_id: str,
    mesh: str,
    phase: str,
    input_name: str,
    input_sha256: str,
    output_name: str,
    output_sha256: str,
    require_schema_v2: bool,
) -> None:
    try:
        stamp = json.loads(content)
    except json.JSONDecodeError as error:
        raise ValueError(f"qualification evidence {label} stamp is invalid: {error}") from error
    if not isinstance(stamp, dict):
        raise ValueError(f"qualification evidence {label} stamp is not an object")
    expected = {
        "campaign_id": campaign_id,
        "method": "GXTB",
        "mesh": mesh,
        "phase": phase,
        "input": input_name,
        "input_sha256": input_sha256,
        "output": output_name,
        "output_sha256": output_sha256,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "input_contract_valid": True,
        "adopted_existing_output": False,
        **{field: identity[field] for field in BUILD_IDENTITY_FIELDS},
    }
    if require_schema_v2:
        expected.update(
            {
                "schema_version": 2,
                "build_id": identity_id,
                "frozen_input": input_name,
                "frozen_input_sha256": input_sha256,
            }
        )
    else:
        stamp_identity = {
            field: str(stamp.get(field)) for field in BUILD_IDENTITY_FIELDS
        }
        if build_id(stamp_identity) != identity_id:
            raise ValueError(
                f"qualification evidence {label} derived build mismatch"
            )
        if stamp.get("build_id", identity_id) != identity_id:
            raise ValueError(f"qualification evidence {label} build mismatch")
    mismatches = [field for field, value in expected.items() if stamp.get(field) != value]
    if mismatches:
        raise ValueError(
            f"qualification evidence {label} stamp mismatch: "
            + ", ".join(mismatches)
        )


def validate_qualification_evidence(
    qualification: dict[str, object],
    artifact_root: Path,
    *,
    campaign_id: str,
    remote_identity: dict[str, object],
    reference_identity: dict[str, object],
    reference_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    remote_build_id = build_id(remote_identity)
    reference_build_id = build_id(reference_identity)
    remote_cp2k_source_revision = str(remote_identity["cp2k_source_revision"])
    reference_cp2k_source_revision = str(
        reference_identity["cp2k_source_revision"]
    )
    if qualification.get("status") != "passed":
        raise ValueError("qualification.status must be 'passed'")
    if qualification.get("evidence_schema_version") != QUALIFICATION_EVIDENCE_SCHEMA_VERSION:
        raise ValueError("qualification evidence_schema_version mismatch")
    total_tolerance = _finite_number(
        qualification.get("total_energy_tolerance_hartree"),
        "total_energy_tolerance_hartree",
    )
    relative_tolerance = _finite_number(
        qualification.get("relative_energy_tolerance_kjmol_per_h2o"),
        "relative_energy_tolerance_kjmol_per_h2o",
    )
    if not 0.0 < total_tolerance <= MAX_TOTAL_ENERGY_TOLERANCE_HARTREE:
        raise ValueError(
            "qualification total-energy tolerance is looser than 1e-10 Eh"
        )
    if not 0.0 < relative_tolerance <= MAX_RELATIVE_ENERGY_TOLERANCE_KJMOL_PER_H2O:
        raise ValueError(
            "qualification relative-energy tolerance is looser than "
            "0.001 kJ/mol/H2O"
        )
    observed_total = _finite_number(
        qualification.get("observed_max_abs_total_energy_delta_hartree"),
        "observed_max_abs_total_energy_delta_hartree",
    )
    observed_relative = _finite_number(
        qualification.get(
            "observed_max_abs_relative_energy_delta_kjmol_per_h2o"
        ),
        "observed_max_abs_relative_energy_delta_kjmol_per_h2o",
    )
    sentinels = qualification.get("same_mesh_dense_relative_sentinels")
    count = qualification.get("same_mesh_dense_relative_sentinel_count")
    if (
        not isinstance(sentinels, list)
        or not sentinels
        or isinstance(count, bool)
        or not isinstance(count, int)
        or count != len(sentinels)
    ):
        raise ValueError(
            "qualification requires at least one counted same-mesh dense sentinel"
        )
    total_deltas: list[float] = []
    relative_deltas: list[float] = []
    for index, value in enumerate(sentinels):
        label = f"sentinel[{index}]"
        if not isinstance(value, dict):
            raise ValueError(f"qualification evidence {label} is not an object")
        if value.get("kind") != "same_mesh_dense_relative_energy":
            raise ValueError(f"qualification evidence {label} kind mismatch")
        if value.get("mesh") not in DENSE_EXTENSION_MESHES:
            raise ValueError(f"qualification evidence {label} is not a dense mesh")
        if value.get("phase") not in PHASES[1:]:
            raise ValueError(f"qualification evidence {label} lacks a non-Ih phase")
        if value.get("remote_build_id") != remote_build_id:
            raise ValueError(f"qualification evidence {label} remote build mismatch")
        if value.get("reference_build_id") != reference_build_id:
            raise ValueError(
                f"qualification evidence {label} reference build mismatch"
            )
        if remote_build_id == reference_build_id:
            raise ValueError(
                f"qualification evidence {label} does not compare distinct builds"
            )
        resolved: dict[str, Path] = {}
        artifact_content: dict[str, bytes] = {}
        artifact_hashes: dict[str, str] = {}
        for artifact in (
            "phase_input",
            "ih_input",
            "remote_phase_output",
            "remote_ih_output",
            "reference_phase_output",
            "reference_ih_output",
            "remote_phase_stamp",
            "remote_ih_stamp",
            "reference_phase_stamp",
            "reference_ih_stamp",
        ):
            artifact_hash = value.get(f"{artifact}_sha256")
            artifact_path, content = _read_hashed_evidence_artifact(
                artifact_root,
                value.get(artifact),
                artifact_hash,
                f"{label}/{artifact}",
            )
            resolved[artifact] = artifact_path
            artifact_content[artifact] = content
            artifact_hashes[artifact] = str(artifact_hash)
        if (
            resolved["phase_input"] == resolved["ih_input"]
            or artifact_hashes["phase_input"] == artifact_hashes["ih_input"]
        ):
            raise ValueError(
                f"qualification evidence {label} phase/Ih inputs are not distinct"
            )
        output_artifacts = (
            "remote_phase_output",
            "remote_ih_output",
            "reference_phase_output",
            "reference_ih_output",
        )
        if len({resolved[name] for name in output_artifacts}) != len(
            output_artifacts
        ):
            raise ValueError(
                f"qualification evidence {label} output paths are not distinct"
            )
        for side in ("remote", "reference"):
            if (
                artifact_hashes[f"{side}_phase_output"]
                == artifact_hashes[f"{side}_ih_output"]
            ):
                raise ValueError(
                    f"qualification evidence {label} {side} phase/Ih outputs "
                    "are not distinct"
                )
        for system in ("phase", "ih"):
            if (
                artifact_hashes[f"remote_{system}_output"]
                == artifact_hashes[f"reference_{system}_output"]
            ):
                raise ValueError(
                    f"qualification evidence {label} copied reference {system} "
                    "output"
                )
        stamp_artifacts = (
            "remote_phase_stamp",
            "remote_ih_stamp",
            "reference_phase_stamp",
            "reference_ih_stamp",
        )
        if len({resolved[name] for name in stamp_artifacts}) != len(
            stamp_artifacts
        ):
            raise ValueError(
                f"qualification evidence {label} stamp paths are not distinct"
            )
        execution_environments = {
            (side, system): _evidence_execution_environment(
                artifact_content[f"{side}_{system}_output"],
                f"{label}/{side}_{system}_output",
            )
            for side in ("remote", "reference")
            for system in ("phase", "ih")
        }
        for side in ("remote", "reference"):
            if execution_environments[(side, "phase")][:3] != (
                execution_environments[(side, "ih")][:3]
            ):
                raise ValueError(
                    f"qualification evidence {label} {side} execution "
                    "environment mismatch"
                )
        if (
            execution_environments[("remote", "phase")][0]
            == execution_environments[("reference", "phase")][0]
            or execution_environments[("remote", "phase")][2]
            == execution_environments[("reference", "phase")][2]
        ):
            raise ValueError(
                f"qualification evidence {label} remote/reference execution "
                "environments are not distinct"
            )
        expected_remote_environment = qualification.get(
            "remote_execution_environment"
        )
        if not isinstance(expected_remote_environment, dict):
            raise ValueError(
                "qualification lacks remote_execution_environment"
            )
        expected_environment = (
            expected_remote_environment.get("program_started_on"),
            expected_remote_environment.get("program_compiled_on"),
            expected_remote_environment.get("program_compiled_for"),
        )
        if any(not isinstance(item, str) or not item for item in expected_environment):
            raise ValueError(
                "qualification remote_execution_environment is invalid"
            )
        if execution_environments[("remote", "phase")][:3] != expected_environment:
            raise ValueError(
                f"qualification evidence {label} remote execution environment "
                "does not match manifest"
            )
        counts: dict[str, int] = {}
        phase_name = str(value["phase"])
        mesh = str(value["mesh"])
        for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
            count_value = value.get(f"{system}_water_count")
            if (
                isinstance(count_value, bool)
                or not isinstance(count_value, int)
                or count_value <= 0
            ):
                raise ValueError(
                    f"qualification evidence {label}/{system}_water_count is invalid"
                )
            derived_count = _validate_evidence_input(
                resolved[f"{system}_input"],
                artifact_content[f"{system}_input"],
                mesh=mesh,
                phase=system_phase,
                label=f"{label}/{system}_input",
            )
            if count_value != derived_count:
                raise ValueError(
                    f"qualification evidence {label}/{system}_water_count "
                    "does not match input"
                )
            counts[system] = count_value
        for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
            reference_record = reference_records.get((mesh, system_phase))
            if reference_record is None:
                raise ValueError(
                    f"qualification evidence {label} lacks trusted reference "
                    f"record {mesh}/{system_phase}"
                )
            if reference_record.get("build_id") != reference_build_id:
                raise ValueError(
                    f"qualification evidence {label} reference record build mismatch"
                )
            input_hash = artifact_hashes[f"{system}_input"]
            if (
                reference_record.get("input_sha256") != input_hash
                or reference_record.get("frozen_input_sha256") != input_hash
            ):
                raise ValueError(
                    f"qualification evidence {label} {system} input is not "
                    "canonical reference bytes"
                )
            expected_project = f"ice_{system_phase}_GXTB_{mesh}"
            expected_input_name = f"{expected_project}.inp"
            expected_output_name = f"{expected_project}.out"
            remote_output = resolved[f"remote_{system}_output"]
            remote_stamp = resolved[f"remote_{system}_stamp"]
            if remote_output.name != expected_output_name:
                raise ValueError(
                    f"qualification evidence {label} remote {system} output "
                    "filename mismatch"
                )
            if resolved[f"{system}_input"] != remote_output.parent / expected_input_name:
                raise ValueError(
                    f"qualification evidence {label} remote {system} input path mismatch"
                )
            if remote_stamp != remote_output.with_suffix(".run.json"):
                raise ValueError(
                    f"qualification evidence {label} remote {system} stamp path mismatch"
                )
            reference_output_value = Path(str(reference_record.get("output")))
            reference_stamp_value = Path(str(reference_record.get("stamp")))
            reference_output_path = (
                reference_output_value
                if reference_output_value.is_absolute()
                else artifact_root / reference_output_value
            ).resolve()
            reference_stamp_path = (
                reference_stamp_value
                if reference_stamp_value.is_absolute()
                else artifact_root / reference_stamp_value
            ).resolve()
            if (
                resolved[f"reference_{system}_output"] != reference_output_path
                or artifact_hashes[f"reference_{system}_output"]
                != reference_record.get("output_sha256")
            ):
                raise ValueError(
                    f"qualification evidence {label} reference {system} output "
                    "does not match trusted record"
                )
            if (
                resolved[f"reference_{system}_stamp"] != reference_stamp_path
                or artifact_hashes[f"reference_{system}_stamp"]
                != reference_record.get("stamp_sha256")
            ):
                raise ValueError(
                    f"qualification evidence {label} reference {system} stamp "
                    "does not match trusted record"
                )
            _validate_qualification_run_stamp(
                artifact_content[f"remote_{system}_stamp"],
                f"{label}/remote_{system}_stamp",
                campaign_id=campaign_id,
                identity=remote_identity,
                identity_id=remote_build_id,
                mesh=mesh,
                phase=system_phase,
                input_name=expected_input_name,
                input_sha256=input_hash,
                output_name=expected_output_name,
                output_sha256=artifact_hashes[f"remote_{system}_output"],
                require_schema_v2=True,
            )
            _validate_qualification_run_stamp(
                artifact_content[f"reference_{system}_stamp"],
                f"{label}/reference_{system}_stamp",
                campaign_id=campaign_id,
                identity=reference_identity,
                identity_id=reference_build_id,
                mesh=mesh,
                phase=system_phase,
                input_name=expected_input_name,
                input_sha256=input_hash,
                output_name=expected_output_name,
                output_sha256=artifact_hashes[f"reference_{system}_output"],
                require_schema_v2=False,
            )
        conversion = _finite_number(
            value.get("hartree_to_kjmol"), f"{label}/hartree_to_kjmol"
        )
        if conversion != HARTREE_TO_KJMOL:
            raise ValueError(
                f"qualification evidence {label} Hartree conversion mismatch"
            )
        energies: dict[str, float] = {}
        for side, source_revision in (
            ("remote", remote_cp2k_source_revision),
            ("reference", reference_cp2k_source_revision),
        ):
            for system, system_phase in (("phase", phase_name), ("ih", "Ih")):
                output = f"{side}_{system}"
                energies[output] = _evidence_total_energy(
                    artifact_content[f"{output}_output"],
                    f"{label}/{output}_output",
                    expected_project=f"ice_{system_phase}_GXTB_{mesh}",
                    expected_source_revision=source_revision,
                    expected_tblite_source_revision=str(
                        (remote_identity if side == "remote" else reference_identity)[
                            "tblite_source_revision"
                        ]
                    ),
                    allow_unknown_tblite_revision=(side == "reference"),
                )
        phase_total_delta = abs(
            energies["remote_phase"] - energies["reference_phase"]
        )
        ih_total_delta = abs(
            energies["remote_ih"] - energies["reference_ih"]
        )
        for system, computed_total in (
            ("phase", phase_total_delta),
            ("ih", ih_total_delta),
        ):
            declared_total = _finite_number(
                value.get(f"{system}_total_energy_delta_hartree"),
                f"{label}/{system}_total_energy_delta_hartree",
            )
            if declared_total < 0.0 or not math.isclose(
                declared_total,
                computed_total,
                rel_tol=1.0e-12,
                abs_tol=1.0e-14,
            ):
                raise ValueError(
                    f"qualification evidence {label} {system} total-energy mismatch"
                )
            if declared_total > total_tolerance:
                raise ValueError(f"qualification evidence {label} exceeds tolerance")
            total_deltas.append(declared_total)
        computed_remote_relative = conversion * (
            energies["remote_phase"] / counts["phase"]
            - energies["remote_ih"] / counts["ih"]
        )
        computed_reference_relative = conversion * (
            energies["reference_phase"] / counts["phase"]
            - energies["reference_ih"] / counts["ih"]
        )
        declared_remote_relative = _finite_number(
            value.get("remote_relative_energy_kjmol_per_h2o"),
            f"{label}/remote_relative_energy_kjmol_per_h2o",
        )
        declared_reference_relative = _finite_number(
            value.get("reference_relative_energy_kjmol_per_h2o"),
            f"{label}/reference_relative_energy_kjmol_per_h2o",
        )
        for side, declared, computed in (
            ("remote", declared_remote_relative, computed_remote_relative),
            ("reference", declared_reference_relative, computed_reference_relative),
        ):
            if not math.isclose(declared, computed, rel_tol=1.0e-12, abs_tol=1.0e-9):
                raise ValueError(
                    f"qualification evidence {label} {side} relative-energy mismatch"
                )
        declared_relative = _finite_number(
            value.get("relative_energy_delta_kjmol_per_h2o"),
            f"{label}/relative_energy_delta_kjmol_per_h2o",
        )
        computed_relative = abs(
            computed_remote_relative - computed_reference_relative
        )
        if declared_relative < 0.0 or not math.isclose(
            declared_relative, computed_relative, rel_tol=1.0e-12, abs_tol=1.0e-12
        ):
            raise ValueError(f"qualification evidence {label} relative-energy mismatch")
        if declared_relative > relative_tolerance:
            raise ValueError(f"qualification evidence {label} exceeds tolerance")
        relative_deltas.append(declared_relative)
    if observed_total < 0.0 or not math.isclose(
        observed_total, max(total_deltas), rel_tol=1.0e-12, abs_tol=1.0e-14
    ):
        raise ValueError("qualification observed total-energy maximum mismatch")
    if observed_relative < 0.0 or not math.isclose(
        observed_relative,
        max(relative_deltas),
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ValueError("qualification observed relative-energy maximum mismatch")
    if observed_total > total_tolerance or observed_relative > relative_tolerance:
        raise ValueError("qualification observed maximum exceeds tolerance")


def validate_execution_manifest_payload(
    execution_manifest: dict[str, object],
    *,
    artifact_root: Path,
    campaign_id: str,
    campaign_manifest_sha256: str,
    identity_id: str,
    identity: dict[str, object],
    reference_identity: dict[str, object],
    reference_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    manifest_identity = execution_manifest.get("build_identity")
    qualification = execution_manifest.get("qualification")
    if not isinstance(manifest_identity, dict) or not isinstance(qualification, dict):
        raise ValueError(
            "execution-build manifest lacks build_identity/qualification blocks"
        )
    checks: dict[str, tuple[object, object]] = {
        "schema_version": (
            execution_manifest.get("schema_version"),
            EXECUTION_BUILD_MANIFEST_SCHEMA_VERSION,
        ),
        "campaign_id": (execution_manifest.get("campaign_id"), campaign_id),
        "campaign_manifest_sha256": (
            execution_manifest.get("campaign_manifest_sha256"),
            campaign_manifest_sha256,
        ),
        "gxtb_protocol_id": (
            execution_manifest.get("gxtb_protocol_id"),
            GXTB_PROTOCOL_ID,
        ),
        "build_id": (execution_manifest.get("build_id"), identity_id),
    }
    for field in BUILD_IDENTITY_FIELDS:
        checks[f"build_identity.{field}"] = (
            manifest_identity.get(field),
            identity.get(field),
        )
    mismatches = [
        f"{label}: actual {actual!r}, expected {expected!r}"
        for label, (actual, expected) in checks.items()
        if actual != expected
    ]
    cli_hash = manifest_identity.get("tblite_cli_sha256")
    if not isinstance(cli_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", cli_hash):
        mismatches.append("build_identity.tblite_cli_sha256 is not a SHA256")
    if mismatches:
        raise ValueError("execution-build manifest mismatch: " + "; ".join(mismatches))
    validate_qualification_evidence(
        qualification,
        artifact_root,
        campaign_id=campaign_id,
        remote_identity=identity,
        reference_identity=reference_identity,
        reference_records=reference_records,
    )


def validate_execution_build_manifest(
    identity: ProductionIdentity,
    tblite: Path,
    campaign_manifest_path: Path,
    campaign_manifest: dict[str, object],
    execution_manifest_path: Path,
    cp2k_source: Path,
    tblite_source: Path,
    artifact_root: Path,
    reference_records: dict[tuple[str, str], dict[str, object]],
) -> None:
    """Qualify an explicitly registered same-source, different-artifact build."""
    relative_manifest = _validation_relative_path(
        execution_manifest_path, artifact_root
    )
    resolved_manifest = _validation_artifact_path(
        artifact_root, relative_manifest, 2
    )
    try:
        execution_manifest = json.loads(resolved_manifest.read_bytes())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"cannot read execution-build manifest {resolved_manifest}: {error}"
        ) from error
    if not isinstance(execution_manifest, dict):
        raise ValueError("execution-build manifest is not an object")
    campaign_cp2k = campaign_manifest.get("cp2k")
    campaign_tblite = campaign_manifest.get("save_tblite")
    if not isinstance(campaign_cp2k, dict) or not isinstance(campaign_tblite, dict):
        raise ValueError("campaign manifest lacks pinned source blocks")
    actual_identity = execution_build_identity(identity)
    actual_build_id = build_id(actual_identity)
    validate_execution_manifest_payload(
        execution_manifest,
        artifact_root=artifact_root,
        campaign_id=str(campaign_manifest.get("campaign_id")),
        campaign_manifest_sha256=sha256(campaign_manifest_path),
        identity_id=actual_build_id,
        identity=actual_identity,
        reference_identity=frozen_build_identity_from_manifest(
            campaign_manifest
        ),
        reference_records=reference_records,
    )
    manifest_identity = execution_manifest["build_identity"]
    checks: dict[str, tuple[object, object]] = {
        "tblite_cli_sha256": (
            manifest_identity.get("tblite_cli_sha256"),
            sha256(tblite),
        ),
        "pinned cp2k source": (
            identity.cp2k_source_revision,
            campaign_cp2k.get("revision"),
        ),
        "pinned save_tblite source": (
            identity.tblite_source_revision,
            campaign_tblite.get("revision"),
        ),
    }
    mismatches = [
        f"{label}: actual {actual!r}, expected {expected!r}"
        for label, (actual, expected) in checks.items()
        if actual != expected
    ]
    for label, source in (("CP2K", cp2k_source), ("save_tblite", tblite_source)):
        status = command_output(["git", "status", "--porcelain"], source)
        if status:
            mismatches.append(f"{label} source worktree is not clean")
    if mismatches:
        raise ValueError("execution-build qualification mismatch: " + "; ".join(mismatches))


def git_worktree_fingerprint(source: Path, revision: str) -> dict[str, object]:
    diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD", "--"],
        cwd=source,
        capture_output=True,
        check=True,
    ).stdout
    untracked_output = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard", "-z"],
        cwd=source,
        capture_output=True,
        check=True,
    ).stdout
    untracked = sorted(
        path.decode("utf-8", errors="surrogateescape")
        for path in untracked_output.split(b"\0")
        if path
    )
    untracked_sha256 = {path: sha256(source / path) for path in untracked}

    digest = hashlib.sha256()
    digest.update(revision.encode())
    digest.update(b"\0tracked-diff\0")
    digest.update(diff)
    for path, file_digest in untracked_sha256.items():
        digest.update(b"\0untracked\0")
        digest.update(path.encode("utf-8", errors="surrogateescape"))
        digest.update(b"\0")
        digest.update(file_digest.encode())
    return {
        "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "untracked_files_sha256": untracked_sha256,
        "worktree_sha256": digest.hexdigest(),
    }


def git_metadata(source: Path | None) -> dict[str, object] | None:
    if source is None:
        return None
    source = source.resolve()
    revision = command_output(["git", "rev-parse", "HEAD"], source)
    metadata = {
        "path": str(source),
        "revision": revision,
        "branch": command_output(["git", "branch", "--show-current"], source),
        "status": command_output(["git", "status", "--short"], source),
        "remotes": command_output(["git", "remote", "-v"], source),
    }
    metadata.update(git_worktree_fingerprint(source, revision))
    return metadata


def manifest_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _validation_artifact_path(
    root: Path,
    value: object,
    schema_version: int,
) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("validation record contains an empty artifact path")
    path = Path(value)
    if schema_version == 1:
        return (path if path.is_absolute() else root / path).resolve()
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"schema-v2 validation path is not relative and safe: {value}")
    resolved_root = root.resolve()
    resolved = (resolved_root / path).resolve()
    if resolved != resolved_root and resolved_root not in resolved.parents:
        raise ValueError(f"schema-v2 validation path escapes artifact root: {value}")
    return resolved


def _read_validation_artifact_bytes(
    root: Path,
    value: object,
    expected_hash: object,
    schema_version: int,
    label: str,
) -> tuple[Path, bytes]:
    path = _validation_artifact_path(root, value, schema_version)
    if not isinstance(expected_hash, str) or not re.fullmatch(
        r"[0-9a-f]{64}", expected_hash
    ):
        raise ValueError(f"validation record {label} has invalid SHA256")
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(f"validation record {label} cannot be read: {error}") from error
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError(f"validation record {label} has invalid hash")
    return path, content


def _validation_relative_path(path: Path, root: Path) -> str:
    resolved_root = root.resolve()
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(
            f"validation artifact is outside the campaign root: {path}"
        ) from error
    if not relative.parts or ".." in relative.parts:
        raise ValueError(f"validation artifact has unsafe relative path: {path}")
    checked = _validation_artifact_path(root, str(relative), 2)
    if checked != resolved:
        raise ValueError(f"validation artifact path does not resolve canonically: {path}")
    return str(relative)


def _normalise_validation_index(
    payload: dict[str, object],
) -> tuple[dict[str, object], int]:
    schema_version = payload.get("schema_version")
    if schema_version not in (1, VALIDATION_INDEX_SCHEMA_VERSION):
        raise ValueError("validation index schema_version must be 1 or 2")
    records = payload.get("records")
    if not isinstance(records, list):
        raise ValueError("validation index lacks records")
    if schema_version == 1:
        legacy_identity = payload.get("build_identity")
        if not isinstance(legacy_identity, dict):
            raise ValueError("schema-v1 validation index lacks build_identity")
        identity = {
            field: str(legacy_identity.get(field))
            for field in BUILD_IDENTITY_FIELDS
        }
        identity_id = build_id(identity)
        normalised_records: list[dict[str, object]] = []
        for value in records:
            if not isinstance(value, dict):
                raise ValueError("validation index record is not an object")
            record = dict(value)
            record["campaign_id"] = payload.get("campaign_id")
            record["gxtb_protocol_id"] = payload.get("gxtb_protocol_id")
            record["build_id"] = identity_id
            if record.get("frozen_input") is None:
                record["frozen_input"] = str(
                    Path(str(record.get("output"))).parent
                    / Path(str(record.get("input"))).name
                )
            if record.get("frozen_input_sha256") is None:
                record["frozen_input_sha256"] = record.get("input_sha256")
            normalised_records.append(record)
        normalised = {
            "schema_version": VALIDATION_INDEX_SCHEMA_VERSION,
            "benchmark": payload.get("benchmark"),
            "method": payload.get("method"),
            "campaign_id": payload.get("campaign_id"),
            "gxtb_protocol_id": payload.get("gxtb_protocol_id"),
            "core_meshes": payload.get("core_meshes", MESHES),
            "supported_dense_extensions": payload.get(
                "supported_dense_extensions", DENSE_EXTENSION_MESHES
            ),
            "source_identity": {
                "cp2k_source_revision": identity["cp2k_source_revision"],
                "tblite_source_revision": identity["tblite_source_revision"],
            },
            "build_identities": {identity_id: identity},
            "parents": [],
            "records": normalised_records,
        }
        return normalised, 1
    identities = payload.get("build_identities")
    if not isinstance(identities, dict):
        raise ValueError("schema-v2 validation index lacks build_identities")
    return dict(payload), VALIDATION_INDEX_SCHEMA_VERSION


def read_validation_index(
    path: Path,
    root: Path,
    *,
    expected_campaign_id: str | None = None,
    expected_source_identity: dict[str, str] | None = None,
    expected_campaign_manifest_sha256: str | None = None,
    campaign_manifest_path: Path | None = None,
    expected_index_sha256: str | None = None,
) -> dict[str, object]:
    """Read v1/v2, verify every artifact/stamp, and return an in-memory v2 view."""
    try:
        index_bytes = path.read_bytes()
        index_hash = hashlib.sha256(index_bytes).hexdigest()
        if (
            expected_index_sha256 is not None
            and index_hash != expected_index_sha256
        ):
            raise ValueError(
                "validation index SHA256 pin mismatch: "
                f"actual {index_hash}, expected {expected_index_sha256}"
            )
        raw_payload = json.loads(index_bytes)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read validation index {path}: {error}") from error
    if not isinstance(raw_payload, dict):
        raise ValueError("validation index is not a JSON object")
    payload, source_schema_version = _normalise_validation_index(raw_payload)
    if payload.get("benchmark") != "DMC-ICE13" or payload.get("method") != "GXTB":
        raise ValueError("validation index benchmark/method mismatch")
    if payload.get("gxtb_protocol_id") != GXTB_PROTOCOL_ID:
        raise ValueError("validation index g-xTB protocol mismatch")
    campaign_id = str(payload.get("campaign_id"))
    if expected_campaign_id is not None and campaign_id != expected_campaign_id:
        raise ValueError("validation index campaign mismatch")

    source_identity = payload.get("source_identity")
    identities = payload.get("build_identities")
    records = payload.get("records")
    if (
        not isinstance(source_identity, dict)
        or not isinstance(identities, dict)
        or not isinstance(records, list)
    ):
        raise ValueError("validation index lacks source/build/record identity")
    source_fields = {
        "cp2k_source_revision": str(source_identity.get("cp2k_source_revision")),
        "tblite_source_revision": str(source_identity.get("tblite_source_revision")),
    }
    if expected_source_identity is not None and source_fields != expected_source_identity:
        raise ValueError("validation index source revision mismatch")
    campaign_manifest_hash = payload.get("campaign_manifest_sha256")
    frozen_identity: dict[str, str] | None = None
    trusted_campaign_manifest: Path | None = None
    if source_schema_version == 2:
        if not isinstance(campaign_manifest_hash, str) or not re.fullmatch(
            r"[0-9a-f]{64}", campaign_manifest_hash
        ):
            raise ValueError("schema-v2 campaign_manifest_sha256 is invalid")
        trusted_campaign_manifest = (
            campaign_manifest_path
            if campaign_manifest_path is not None
            else root.resolve().parent
            / "campaigns"
            / campaign_id
            / "build_manifest.json"
        )
        if not trusted_campaign_manifest.is_file():
            raise ValueError(
                "schema-v2 trusted campaign manifest does not exist: "
                f"{trusted_campaign_manifest}"
            )
        trusted_manifest_bytes = trusted_campaign_manifest.read_bytes()
        trusted_manifest_hash = hashlib.sha256(trusted_manifest_bytes).hexdigest()
        if campaign_manifest_hash != trusted_manifest_hash:
            raise ValueError("schema-v2 campaign manifest hash mismatch")
        if (
            expected_campaign_manifest_sha256 is not None
            and campaign_manifest_hash != expected_campaign_manifest_sha256
        ):
            raise ValueError("schema-v2 campaign manifest hash mismatch")
        try:
            trusted_manifest = json.loads(trusted_manifest_bytes)
        except json.JSONDecodeError as error:
            raise ValueError(
                f"cannot read campaign manifest {trusted_campaign_manifest}: {error}"
            ) from error
        if not isinstance(trusted_manifest, dict):
            raise ValueError(
                f"campaign manifest {trusted_campaign_manifest} is not a JSON object"
            )
        if trusted_manifest.get("campaign_id") != campaign_id:
            raise ValueError("schema-v2 trusted campaign manifest campaign mismatch")
        frozen_identity = frozen_build_identity_from_manifest(trusted_manifest)
        if source_fields != {
            "cp2k_source_revision": frozen_identity["cp2k_source_revision"],
            "tblite_source_revision": frozen_identity["tblite_source_revision"],
        }:
            raise ValueError("schema-v2 source revisions differ from campaign manifest")

    normalised_identities: dict[str, dict[str, object]] = {}
    frozen_identity_ids: set[str] = set()
    pending_execution_manifests: list[
        tuple[dict[str, object], str, dict[str, object]]
    ] = []
    for identity_id, value in identities.items():
        if not isinstance(identity_id, str) or not isinstance(value, dict):
            raise ValueError("validation build identity is malformed")
        identity = dict(value)
        if build_id(identity) != identity_id:
            raise ValueError(f"validation build identity digest mismatch: {identity_id}")
        for field in BUILD_IDENTITY_FIELDS[:3]:
            if not re.fullmatch(r"[0-9a-f]{64}", str(identity.get(field))):
                raise ValueError(
                    f"validation build identity has invalid SHA256: {identity_id}/{field}"
                )
        for field in BUILD_IDENTITY_FIELDS[3:]:
            if not re.fullmatch(r"[0-9a-f]{40}", str(identity.get(field))):
                raise ValueError(
                    f"validation build identity has invalid revision: {identity_id}/{field}"
                )
        for field, expected in source_fields.items():
            if identity.get(field) != expected:
                raise ValueError(
                    f"validation build identity source mismatch: {identity_id}/{field}"
                )
        manifest_path_value = identity.get("execution_build_manifest")
        manifest_hash = identity.get("execution_build_manifest_sha256")
        if (manifest_path_value is None) != (manifest_hash is None):
            raise ValueError(
                f"validation build identity has incomplete execution manifest: {identity_id}"
            )
        is_frozen_identity = (
            frozen_identity is not None
            and all(
                identity.get(field) == frozen_identity[field]
                for field in BUILD_IDENTITY_FIELDS
            )
        )
        if is_frozen_identity:
            frozen_identity_ids.add(identity_id)
        if (
            source_schema_version == 2
            and manifest_path_value is None
            and not is_frozen_identity
        ):
            raise ValueError(
                "alternate validation build identity lacks execution manifest: "
                f"{identity_id}"
            )
        if manifest_path_value is not None:
            if not isinstance(manifest_hash, str) or not re.fullmatch(
                r"[0-9a-f]{64}", manifest_hash
            ):
                raise ValueError(
                    f"validation execution manifest has invalid SHA256: {identity_id}"
                )
            resolved_manifest = _validation_artifact_path(
                root, manifest_path_value, 2
            )
            try:
                execution_manifest_bytes = resolved_manifest.read_bytes()
            except OSError as error:
                raise ValueError(
                    f"validation execution manifest cannot be read: {identity_id}: "
                    f"{error}"
                ) from error
            if hashlib.sha256(execution_manifest_bytes).hexdigest() != manifest_hash:
                raise ValueError(
                    f"validation execution manifest hash mismatch: {identity_id}"
                )
            try:
                execution_manifest = json.loads(execution_manifest_bytes)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"invalid execution manifest for {identity_id}: {error}"
                ) from error
            if not isinstance(execution_manifest, dict):
                raise ValueError("execution-build manifest is not an object")
            pending_execution_manifests.append(
                (execution_manifest, identity_id, identity)
            )
        normalised_identities[identity_id] = identity

    if (
        trusted_campaign_manifest is not None
        and sha256(trusted_campaign_manifest) != campaign_manifest_hash
    ):
        raise ValueError("schema-v2 trusted campaign manifest changed during validation")

    seen: set[tuple[str, str]] = set()
    normalised_records: list[dict[str, object]] = []
    coverage: dict[str, list[str]] = {}
    for value in records:
        if not isinstance(value, dict):
            raise ValueError("validation index record is not an object")
        record = dict(value)
        mesh, phase = validation_record_key(record)
        key = (mesh, phase)
        if mesh not in SUPPORTED_MESHES or phase not in PHASES:
            raise ValueError(f"unsupported validation record {mesh}/{phase}")
        if key in seen:
            raise ValueError(f"duplicate validation record {mesh}/{phase}")
        seen.add(key)
        if record.get("campaign_id") != campaign_id:
            raise ValueError(f"validation record campaign mismatch: {mesh}/{phase}")
        if record.get("gxtb_protocol_id") != GXTB_PROTOCOL_ID:
            raise ValueError(f"validation record protocol mismatch: {mesh}/{phase}")
        identity_id = record.get("build_id")
        identity = normalised_identities.get(str(identity_id))
        if identity is None:
            raise ValueError(f"validation record has unknown build: {mesh}/{phase}")
        if source_schema_version == 2 and (
            not isinstance(record.get("frozen_input"), str)
            or not isinstance(record.get("frozen_input_sha256"), str)
        ):
            raise ValueError(
                f"schema-v2 validation record lacks frozen input binding: {mesh}/{phase}"
            )
        input_hash = record.get("input_sha256")
        frozen_hash = record.get("frozen_input_sha256", input_hash)
        if input_hash != frozen_hash:
            raise ValueError(
                f"validation record {mesh}/{phase} source/frozen input hash mismatch"
            )
        files = {
            "input": (record.get("input"), input_hash),
            "frozen_input": (record.get("frozen_input"), frozen_hash),
            "output": (record.get("output"), record.get("output_sha256")),
            "stamp": (record.get("stamp"), record.get("stamp_sha256")),
        }
        resolved: dict[str, Path] = {}
        artifact_bytes: dict[str, bytes] = {}
        for label, (file_name, expected_hash) in files.items():
            file_path, content = _read_validation_artifact_bytes(
                root,
                file_name,
                expected_hash,
                source_schema_version,
                f"{mesh}/{phase}/{label}",
            )
            resolved[label] = file_path
            artifact_bytes[label] = content
        expected_project = f"ice_{phase}_GXTB_{mesh}"
        expected_input_name = f"{expected_project}.inp"
        expected_output_name = f"{expected_project}.out"
        if source_schema_version in (1, 2):
            canonical_input = (
                root / GXTB_INPUT_DIRECTORY / mesh / expected_input_name
            ).resolve()
            canonical_output = (
                root
                / GXTB_RUN_DIRECTORY
                / mesh
                / phase
                / expected_output_name
            ).resolve()
            if resolved["input"] != canonical_input:
                raise ValueError(
                    f"validation record {mesh}/{phase} input path is not canonical"
                )
            if resolved["output"] != canonical_output:
                raise ValueError(
                    f"validation record {mesh}/{phase} output path is not canonical"
                )
            if resolved["frozen_input"] != canonical_output.parent / expected_input_name:
                raise ValueError(
                    f"validation record {mesh}/{phase} frozen input path mismatch"
                )
            if resolved["stamp"] != canonical_output.with_suffix(".run.json"):
                raise ValueError(
                    f"validation record {mesh}/{phase} stamp path mismatch"
                )
            if resolved["input"] == resolved["frozen_input"]:
                raise ValueError(
                    f"validation record {mesh}/{phase} source/frozen paths coincide"
                )
            if artifact_bytes["input"] != artifact_bytes["frozen_input"]:
                raise ValueError(
                    f"validation record {mesh}/{phase} source/frozen bytes differ"
                )
            for input_label in ("input", "frozen_input"):
                _validate_evidence_input(
                    resolved[input_label],
                    artifact_bytes[input_label],
                    mesh=mesh,
                    phase=phase,
                    label=f"record {mesh}/{phase}/{input_label}",
                )
        validated_energy = _evidence_total_energy(
            artifact_bytes["output"],
            f"record {mesh}/{phase}/output",
            expected_project=expected_project,
            expected_source_revision=str(identity["cp2k_source_revision"]),
            expected_tblite_source_revision=str(
                identity["tblite_source_revision"]
            ),
            allow_unknown_tblite_revision=(
                source_schema_version == 1
                or str(identity_id) in frozen_identity_ids
            ),
        )
        try:
            stamp = json.loads(artifact_bytes["stamp"])
        except json.JSONDecodeError as error:
            raise ValueError(f"invalid stamp for {mesh}/{phase}: {error}") from error
        if not isinstance(stamp, dict):
            raise ValueError(f"invalid stamp for {mesh}/{phase}: not an object")
        expected_stamp_fields = {
            "campaign_id": campaign_id,
            "method": "GXTB",
            "mesh": mesh,
            "phase": phase,
            "gxtb_protocol_id": GXTB_PROTOCOL_ID,
            "input_sha256": input_hash,
            "output_sha256": record.get("output_sha256"),
            "input": expected_input_name,
            "output": expected_output_name,
            "input_contract_valid": True,
            "adopted_existing_output": False,
            **{field: identity.get(field) for field in BUILD_IDENTITY_FIELDS},
        }
        stamp_version = stamp.get("schema_version", 1)
        requires_v2_stamp = (
            source_schema_version == 2
            and str(identity_id) not in frozen_identity_ids
        )
        if requires_v2_stamp and stamp_version != 2:
            raise ValueError(
                f"stamp mismatch for {mesh}/{phase}: schema_version"
            )
        if stamp_version == 2:
            expected_stamp_fields.update(
                {
                    "schema_version": 2,
                    "build_id": identity_id,
                    "frozen_input": expected_input_name,
                    "frozen_input_sha256": frozen_hash,
                }
            )
        mismatches = [
            field
            for field, expected in expected_stamp_fields.items()
            if stamp.get(field) != expected
        ]
        if mismatches:
            raise ValueError(
                f"stamp mismatch for {mesh}/{phase}: {', '.join(mismatches)}"
            )
        if stamp.get("frozen_input", expected_input_name) != expected_input_name:
            raise ValueError(f"stamp mismatch for {mesh}/{phase}: frozen_input")
        if stamp.get("frozen_input_sha256", input_hash) != frozen_hash:
            raise ValueError(
                f"stamp mismatch for {mesh}/{phase}: frozen_input_sha256"
            )
        record["validated_energy_hartree"] = validated_energy
        normalised_records.append(record)
        coverage.setdefault(mesh, []).append(phase)

    reference_records = {
        validation_record_key(record): record
        for record in normalised_records
        if str(record.get("build_id")) in frozen_identity_ids
    }
    for execution_manifest, execution_identity_id, execution_identity in (
        pending_execution_manifests
    ):
        if frozen_identity is None:
            raise ValueError("schema-v2 execution manifest lacks frozen identity")
        validate_execution_manifest_payload(
            execution_manifest,
            artifact_root=root,
            campaign_id=campaign_id,
            campaign_manifest_sha256=str(campaign_manifest_hash),
            identity_id=execution_identity_id,
            identity=execution_identity,
            reference_identity=frozen_identity,
            reference_records=reference_records,
        )

    coverage = {
        mesh: [phase for phase in PHASES if phase in phases]
        for mesh, phases in coverage.items()
    }
    used_identity_ids = {str(record["build_id"]) for record in normalised_records}
    if (
        source_schema_version == 2
        and set(normalised_identities) != used_identity_ids
    ):
        raise ValueError("schema-v2 validation index contains unused build identities")
    if source_schema_version == 2 and payload.get("validated_phase_coverage") != coverage:
        raise ValueError("schema-v2 validation coverage is not derived from records")
    payload["build_identities"] = normalised_identities
    payload["records"] = normalised_records
    payload["validated_phase_coverage"] = coverage
    payload["source_schema_version"] = source_schema_version
    payload["source_index_sha256"] = index_hash
    if (
        trusted_campaign_manifest is not None
        and sha256(trusted_campaign_manifest) != campaign_manifest_hash
    ):
        raise ValueError("schema-v2 trusted campaign manifest changed during validation")
    return payload


def convergence_validation_index_path(args: argparse.Namespace) -> Path:
    return (
        args.root
        / "data"
        / f"dmc_ice13_{args.analysis_prefix}_validation_index.json"
    )


def convergence_validation_snapshot_path(
    args: argparse.Namespace,
    digest: str,
) -> Path:
    return (
        args.root
        / "data"
        / "validation_snapshots"
        / (
            f"dmc_ice13_{args.analysis_prefix}_validation_index."
            f"{digest}.json"
        )
    )


def validate_production_paths(
    root: Path,
    analysis_prefix: str,
    gxtb_input_root: Path,
    gxtb_run_root: Path,
) -> None:
    if analysis_prefix != GXTB_ANALYSIS_PREFIX:
        raise ValueError(
            f"production --analysis-prefix must be exactly {GXTB_ANALYSIS_PREFIX!r}"
        )
    _validation_relative_path(gxtb_input_root, root)
    _validation_relative_path(gxtb_run_root, root)
    data_root = (root / "data").resolve()
    current_index = (
        data_root / f"dmc_ice13_{analysis_prefix}_validation_index.json"
    ).resolve()
    if current_index.parent != data_root:
        raise ValueError("production validation index escapes ROOT/data")


def validated_gxtb_phase_coverage(
    args: argparse.Namespace,
    identity: ProductionIdentity,
) -> dict[str, list[str]]:
    coverage_sets: dict[str, set[str]] = {}
    base_index = getattr(args, "base_validation_index_payload", None)
    if isinstance(base_index, dict):
        for record in base_index.get("records", []):
            if isinstance(record, dict):
                mesh, phase = validation_record_key(record)
                coverage_sets.setdefault(mesh, set()).add(phase)
    for mesh in SUPPORTED_MESHES:
        valid = [
            job.phase
            for job in jobs(
                args.root,
                ["GXTB"],
                [mesh],
                gxtb_input_root=args.gxtb_input_root,
                gxtb_run_root=args.gxtb_run_root,
            )
            if scientific_stamp_valid(job, identity)
        ]
        if valid:
            coverage_sets.setdefault(mesh, set()).update(valid)
    return {
        mesh: [phase for phase in PHASES if phase in coverage_sets.get(mesh, set())]
        for mesh in SUPPORTED_MESHES
        if coverage_sets.get(mesh)
    }


def _current_validation_record(
    args: argparse.Namespace,
    job: Job,
    identity_id: str,
) -> dict[str, object]:
    output = job.run_dir / job.output_name
    stamp = stamp_path(job)
    frozen_input = frozen_input_path(job)
    return {
        "campaign_id": getattr(args, "campaign_id", None),
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "build_id": identity_id,
        "mesh": job.mesh,
        "phase": job.phase,
        "input": manifest_path(job.input_path, args.root),
        "input_sha256": sha256(job.input_path),
        "frozen_input": manifest_path(frozen_input, args.root),
        "frozen_input_sha256": sha256(frozen_input),
        "output": manifest_path(output, args.root),
        "output_sha256": sha256(output),
        "stamp": manifest_path(stamp, args.root),
        "stamp_sha256": sha256(stamp),
    }


def write_convergence_validation_index(
    args: argparse.Namespace,
    identity: ProductionIdentity,
) -> Path:
    """Freeze every currently hash-valid g-xTB phase used by analysis."""
    campaign_manifest_path = getattr(args, "campaign_manifest", None)
    if not isinstance(campaign_manifest_path, Path):
        campaign_manifest_path = (
            args.root.resolve().parent
            / "campaigns"
            / identity.campaign_id
            / "build_manifest.json"
        )
    if not campaign_manifest_path.is_file():
        raise ValueError(
            "schema-v2 trusted campaign manifest does not exist: "
            f"{campaign_manifest_path}"
        )
    campaign_manifest = read_campaign_manifest(campaign_manifest_path)
    if campaign_manifest.get("campaign_id") != identity.campaign_id:
        raise ValueError("schema-v2 trusted campaign manifest campaign mismatch")
    frozen_identity = frozen_build_identity_from_manifest(campaign_manifest)
    campaign_manifest_hash = sha256(campaign_manifest_path)
    base_path = getattr(args, "base_validation_index", None)
    expected_base_hash = getattr(args, "base_validation_index_sha256", None)
    if isinstance(base_path, Path) and isinstance(expected_base_hash, str):
        if not base_path.is_file() or sha256(base_path) != expected_base_hash:
            raise ValueError("base validation index changed after initial verification")
        args.base_validation_index_payload = read_validation_index(
            base_path,
            args.root,
            expected_campaign_id=identity.campaign_id,
            expected_source_identity={
                "cp2k_source_revision": identity.cp2k_source_revision,
                "tblite_source_revision": identity.tblite_source_revision,
            },
            expected_campaign_manifest_sha256=(
                campaign_manifest_hash
            ),
            campaign_manifest_path=campaign_manifest_path,
            expected_index_sha256=expected_base_hash,
        )
    args.campaign_id = identity.campaign_id
    current_identity = execution_build_identity(identity)
    current_identity_id = build_id(current_identity)
    identities: dict[str, dict[str, object]] = {
        current_identity_id: current_identity
    }
    base_index = getattr(args, "base_validation_index_payload", None)
    parents: list[dict[str, object]] = []
    records_by_key: dict[tuple[str, str], dict[str, object]] = {}
    if isinstance(base_index, dict):
        base_identities = base_index.get("build_identities", {})
        if not isinstance(base_identities, dict):
            raise ValueError("base validation index lacks build identities")
        for identity_id, value in base_identities.items():
            if not isinstance(value, dict):
                raise ValueError("base validation build identity is malformed")
            existing = identities.get(str(identity_id))
            if existing is not None:
                if any(
                    existing.get(field) != value.get(field)
                    for field in BUILD_IDENTITY_FIELDS
                ):
                    raise ValueError(
                        f"validation build identity collision: {identity_id}"
                    )
                identities[str(identity_id)] = {**value, **existing}
            else:
                identities[str(identity_id)] = dict(value)
        for value in base_index.get("records", []):
            if not isinstance(value, dict):
                raise ValueError("base validation record is malformed")
            record = dict(value)
            # This is an in-memory verification result, not part of the signed
            # validation-record schema and therefore must not be inherited.
            record.pop("validated_energy_hartree", None)
            key = validation_record_key(record)
            if key in records_by_key and records_by_key[key] != record:
                raise ValueError(f"validation record collision: {key[0]}/{key[1]}")
            records_by_key[key] = record
        if isinstance(base_path, Path):
            parents.append(
                {
                    "schema_version": base_index.get("source_schema_version"),
                    "sha256": sha256(base_path),
                    "record_count": len(base_index.get("records", [])),
                }
            )
    for mesh in SUPPORTED_MESHES:
        for job in jobs(
            args.root,
            ["GXTB"],
            [mesh],
            gxtb_input_root=args.gxtb_input_root,
            gxtb_run_root=args.gxtb_run_root,
        ):
            if not scientific_stamp_valid(job, identity):
                continue
            record = _current_validation_record(args, job, current_identity_id)
            key = validation_record_key(record)
            existing = records_by_key.get(key)
            if existing is not None and existing != record:
                raise ValueError(f"validation record collision: {mesh}/{job.phase}")
            records_by_key[key] = record
    records = [
        records_by_key[(mesh, phase)]
        for mesh in SUPPORTED_MESHES
        for phase in PHASES
        if (mesh, phase) in records_by_key
    ]
    coverage: dict[str, list[str]] = {}
    for record in records:
        mesh, phase = validation_record_key(record)
        coverage.setdefault(mesh, []).append(phase)
    execution_manifest = getattr(args, "execution_build_manifest", None)
    used_identity_ids = {str(record["build_id"]) for record in records}
    if isinstance(execution_manifest, Path) and current_identity_id in used_identity_ids:
        identities[current_identity_id] = {
            **identities[current_identity_id],
            "execution_build_manifest": _validation_relative_path(
                execution_manifest, args.root
            ),
            "execution_build_manifest_sha256": sha256(execution_manifest),
        }
    if (
        current_identity_id in used_identity_ids
        and not isinstance(execution_manifest, Path)
        and any(
            current_identity.get(field) != frozen_identity[field]
            for field in BUILD_IDENTITY_FIELDS
        )
    ):
        raise ValueError(
            "alternate validation build identity lacks execution manifest: "
            f"{current_identity_id}"
        )
    identities = {
        identity_id: identities[identity_id]
        for identity_id in sorted(used_identity_ids)
    }
    payload = {
        "schema_version": VALIDATION_INDEX_SCHEMA_VERSION,
        "benchmark": "DMC-ICE13",
        "method": "GXTB",
        "campaign_id": identity.campaign_id,
        "campaign_manifest_sha256": campaign_manifest_hash,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "core_meshes": MESHES,
        "supported_dense_extensions": DENSE_EXTENSION_MESHES,
        "source_identity": {
            "cp2k_source_revision": identity.cp2k_source_revision,
            "tblite_source_revision": identity.tblite_source_revision,
        },
        "build_identities": identities,
        "parents": parents,
        "validated_phase_coverage": coverage,
        "records": records,
    }
    if isinstance(base_path, Path) and isinstance(expected_base_hash, str):
        if not base_path.is_file() or sha256(base_path) != expected_base_hash:
            raise ValueError("base validation index changed before merged write")
    content = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    digest = hashlib.sha256(content).hexdigest()
    snapshot = convergence_validation_snapshot_path(args, digest)
    if snapshot.exists():
        if snapshot.read_bytes() != content:
            raise ValueError(
                f"immutable validation snapshot content mismatch: {snapshot}"
            )
    else:
        atomic_write_bytes(snapshot, content)
    atomic_write_bytes(convergence_validation_index_path(args), content)
    return snapshot


def legacy_gxtb_diagnostics(root: Path) -> list[dict[str, object]]:
    diagnostics: list[dict[str, object]] = []
    for job in jobs(root, ["GXTB"]):
        output = job.run_dir / job.output_name
        local_input = job.run_dir / job.input_path.name
        stamp = stamp_path(job)
        if not output.is_file() and not local_input.is_file() and not stamp.is_file():
            continue
        stamp_payload: dict[str, object] = {}
        if stamp.is_file():
            try:
                stamp_payload = json.loads(stamp.read_text())
            except (json.JSONDecodeError, OSError):
                stamp_payload = {}
        full_grid = is_legacy_full_grid_input(local_input)
        diagnostics.append(
            {
                "mesh": job.mesh,
                "phase": job.phase,
                "classification": (
                    "legacy_full_grid_diagnostic"
                    if full_grid
                    else "legacy_pre_protocol_diagnostic"
                ),
                "accepted_as_production": False,
                "local_input": manifest_path(local_input, root),
                "local_input_sha256": sha256(local_input) if local_input.is_file() else None,
                "legacy_full_grid": full_grid,
                "output": manifest_path(output, root),
                "output_completed": has_completed(output),
                "output_sha256": sha256(output) if output.is_file() else None,
                "old_stamp": manifest_path(stamp, root),
                "old_cp2k_sha256": stamp_payload.get("cp2k_sha256"),
            }
        )
    return diagnostics


def job_execution_evidence(
    job: Job,
    execution_pool: benchmark_execution.ExecutionPool,
    root: Path,
) -> dict[str, object]:
    """Classify one timing only after full schema-v2 artifact revalidation."""
    output = job.run_dir / job.output_name
    record_path = benchmark_execution.execution_record_path(output)
    record: dict[str, object] = {}
    if record_path.is_file():
        try:
            loaded = json.loads(record_path.read_text())
            if isinstance(loaded, dict):
                record = loaded
        except (json.JSONDecodeError, OSError):
            record = {}
    classification = benchmark_execution.execution_record_timing_classification(
        record_path,
        execution_pool.contract,
        output,
        stamp_path(job),
    )
    return {
        "execution_record": manifest_path(record_path, root),
        "execution_record_sha256": (
            sha256(record_path) if record_path.is_file() else None
        ),
        "execution_record_schema_version": record.get("schema_version"),
        "execution_return_code": record.get("return_code"),
        "timing_classification": classification,
    }


def aggregate_execution_timing(
    run_manifest: list[dict[str, object]],
    failures: list[tuple[Job, int]] | None,
) -> dict[str, object]:
    """Summarize only the explicitly included, hash-validated timing population."""
    included = [
        record
        for record in run_manifest
        if record.get("completed_and_hash_validated") is True
    ]
    counts: dict[str, int] = {}
    for record in included:
        classification = str(record.get("timing_classification", "missing"))
        counts[classification] = counts.get(classification, 0) + 1
    eligible = (
        bool(included)
        and not failures
        and all(
            record.get("timing_classification") == "production_scaling_eligible"
            and record.get("execution_record_schema_version")
            == benchmark_execution.SCHEMA_VERSION
            and record.get("execution_return_code") == 0
            for record in included
        )
    )
    return {
        "timing_population": "completed_and_hash_validated_run_manifest_jobs",
        "timing_record_count": len(included),
        "timing_classification_counts": counts,
        "timing_classification": (
            "production_scaling_eligible" if eligible else "timing_non_scaling"
        ),
        "all_included_timings_schema_v2_revalidated_and_rc0": eligible,
        "failed_jobs_present": bool(failures),
    }


def write_provenance(
    args: argparse.Namespace,
    identity: ProductionIdentity,
    methods: list[str],
    requested_meshes: list[str],
    requested_phases: list[str],
    completed: int,
    failures: list[tuple[Job, int]] | None = None,
) -> None:
    parallel_jobs = int(getattr(args, "jobs", 1))
    threads_per_job = int(getattr(args, "threads_per_job", 1))
    execution_pool = getattr(args, "execution_pool", None)
    mpi_ranks_per_job = (
        execution_pool.mpi_ranks_per_job if execution_pool is not None else 1
    )
    parallelism = execution_parallelism(
        parallel_jobs, threads_per_job, mpi_ranks_per_job
    )
    if execution_pool is not None:
        parallelism.update(
            {
                "affinity_policy": "openmpi_ordered_pe_list",
                "execution_contract": execution_pool.contract,
                "execution_contract_sha256": execution_pool.contract_sha256,
            }
        )
    phase_coverage = validated_gxtb_phase_coverage(args, identity)
    validated_gxtb_meshes = [
        mesh
        for mesh in SUPPORTED_MESHES
        if phase_coverage.get(mesh, []) == PHASES
    ]
    manifest_meshes = [
        mesh
        for mesh in SUPPORTED_MESHES
        if mesh in MESHES or (args.gxtb_input_root / mesh).is_dir()
    ]
    validation_snapshot = getattr(args, "validation_index_snapshot", None)
    verified_records: dict[tuple[str, str], dict[str, object]] = {}
    if isinstance(validation_snapshot, Path) and validation_snapshot.is_file():
        verified_index = read_validation_index(
            validation_snapshot,
            args.root,
            expected_campaign_id=identity.campaign_id,
            expected_source_identity={
                "cp2k_source_revision": identity.cp2k_source_revision,
                "tblite_source_revision": identity.tblite_source_revision,
            },
            expected_campaign_manifest_sha256=sha256(args.campaign_manifest),
            campaign_manifest_path=args.campaign_manifest,
        )
        verified_records = {
            validation_record_key(record): record
            for record in verified_index.get("records", [])
            if isinstance(record, dict)
        }
    run_manifest: list[dict[str, object]] = []
    for job in jobs(
        args.root,
        methods,
        manifest_meshes,
        gxtb_input_root=args.gxtb_input_root,
        gxtb_run_root=args.gxtb_run_root,
    ):
        output = job.run_dir / job.output_name
        verified_record = (
            verified_records.get((job.mesh, job.phase))
            if job.method == "GXTB"
            else None
        )
        completed_job = verified_record is not None
        adopted = False
        if stamp_path(job).is_file():
            try:
                adopted = bool(
                    json.loads(stamp_path(job).read_text()).get(
                        "adopted_existing_output"
                    )
                )
            except (json.JSONDecodeError, OSError):
                adopted = False
        manifest_record: dict[str, object] = {
            "method": job.method,
            "mesh": job.mesh,
            "phase": job.phase,
            "input": manifest_path(job.input_path, args.root),
            "input_sha256": (
                sha256(job.input_path) if job.input_path.is_file() else None
            ),
            "input_contract_valid": not gxtb_input_contract_errors(job),
            "gxtb_protocol_id": (
                GXTB_PROTOCOL_ID if job.method == "GXTB" else None
            ),
            "output": manifest_path(output, args.root),
            "output_sha256": (
                verified_record.get("output_sha256")
                if verified_record is not None
                else None
            ),
            "execution_build_id": (
                verified_record.get("build_id")
                if verified_record is not None
                else None
            ),
            "completed_and_hash_validated": completed_job,
            "adopted_existing_output": adopted,
        }
        if execution_pool is not None:
            manifest_record.update(job_execution_evidence(job, execution_pool, args.root))
        run_manifest.append(manifest_record)
    if execution_pool is not None:
        parallelism.update(aggregate_execution_timing(run_manifest, failures))
    payload: dict[str, object] = {
        "benchmark": "DMC-ICE13",
        "campaign": {
            "id": identity.campaign_id,
            "manifest": str(args.campaign_manifest),
            "manifest_sha256": sha256(args.campaign_manifest),
            "execution_build_manifest": (
                str(args.execution_build_manifest)
                if args.execution_build_manifest is not None
                else None
            ),
            "execution_build_manifest_sha256": (
                sha256(args.execution_build_manifest)
                if args.execution_build_manifest is not None
                else None
            ),
        },
        "methods": methods,
        "execution_parallelism": parallelism,
        "run_status": (
            "failed"
            if failures
            else "completed"
            if completed == len(PHASES) * len(MESHES)
            else "partial"
        ),
        "failed_jobs": [
            {
                "mesh": job.mesh,
                "method": job.method,
                "phase": job.phase,
                "returncode": returncode,
                "output": manifest_path(job.run_dir / job.output_name, args.root),
            }
            for job, returncode in failures or []
        ],
        "invocation": [sys.executable, *sys.argv],
        "benchmark_software": {
            "runner": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "input_and_analysis": {
                "path": str((args.root / "scripts" / "dmc_ice13_kpoint_benchmark.py").resolve()),
                "sha256": sha256(args.root / "scripts" / "dmc_ice13_kpoint_benchmark.py"),
            },
        },
        "protocol": {
            "phases": PHASES,
            "meshes": MESHES,
            "supported_dense_extensions": DENSE_EXTENSION_MESHES,
            "requested_meshes_this_invocation": requested_meshes,
            "requested_phases_this_invocation": requested_phases,
            "validated_gxtb_meshes_for_analysis": validated_gxtb_meshes,
            "validated_gxtb_phase_coverage": phase_coverage,
            "convergence_validation_index": manifest_path(
                convergence_validation_index_path(args), args.root
            ),
            "base_validation_index": (
                manifest_path(args.base_validation_index, args.root)
                if args.base_validation_index is not None
                else None
            ),
            "base_validation_index_sha256": (
                sha256(args.base_validation_index)
                if args.base_validation_index is not None
                else None
            ),
            "convergence_validation_index_sha256": (
                sha256(convergence_validation_index_path(args))
                if convergence_validation_index_path(args).is_file()
                else None
            ),
            "convergence_validation_snapshot": (
                manifest_path(validation_snapshot, args.root)
                if isinstance(validation_snapshot, Path)
                else None
            ),
            "convergence_validation_snapshot_sha256": (
                sha256(validation_snapshot)
                if isinstance(validation_snapshot, Path)
                and validation_snapshot.is_file()
                else None
            ),
            "job_count": len(methods) * len(PHASES) * len(MESHES),
            "requested_job_count": len(methods) * len(requested_phases) * len(requested_meshes),
            "completed_job_count": completed,
            "dense_extension_completed_job_count": sum(
                len(phase_coverage.get(mesh, [])) for mesh in DENSE_EXTENSION_MESHES
            ),
            "accuracy": 0.1,
            "eps_scf": 1.0e-9,
            "gxtb_protocol_id": GXTB_PROTOCOL_ID,
            "gxtb_input_root": manifest_path(args.gxtb_input_root, args.root),
            "gxtb_run_root": manifest_path(args.gxtb_run_root, args.root),
            "analysis_output_prefix": args.analysis_prefix,
            "kpoint_scheme": (
                "implicit Gamma plus explicit native-Bloch MacDonald meshes; every "
                "non-Gamma GXTB/GFN calculation uses FULL_GRID F with SPGLIB "
                "symmetry reduction; CP2K expands the GXTB density and overlap to "
                "the coupled full mesh internally before save_tblite evaluation"
            ),
            "gxtb_exchange": (
                "image-space whole-mesh exchange with inverse -ik and forward +ik "
                "Fourier transforms; energy and q-dependent shell response are "
                "Brillouin-zone contractions per primitive cell"
            ),
            "gxtb_multik_derivatives": (
                "not used in the DMC13 energy-only runs; the coupled whole-mesh "
                "reverse mode is available and validated separately by finite-difference "
                "and explicit Gamma-supercell force/stress checks"
            ),
            "gxtb_mixer": (
                "SCC_MIXER TBLITE (two damped Fock updates, then seven-vector "
                "Fock DIIS; ITERATIONS 300)"
            ),
            "density_mixing_keyword": "DIRECT_P_MIXING",
            "reference_note": (
                "legacy XI=0.16 retained for GFN1/GFN2 comparability; "
                "XI=0.15 sensitivity stored separately"
            ),
        },
        "cp2k": {
            "path": str(identity.cp2k),
            "sha256": identity.cp2k_sha256,
            "execution_build_id": build_id(execution_build_identity(identity)),
            "loaded_library": str(identity.cp2k_library),
            "loaded_library_sha256": identity.cp2k_library_sha256,
            "source_revision_validated": identity.cp2k_source_revision,
            "version": command_output([str(identity.cp2k), "--version"]),
            "source": git_metadata(args.cp2k_source),
        },
        "run_manifest": run_manifest,
        "legacy_diagnostics": legacy_gxtb_diagnostics(args.root),
    }
    if args.tblite is not None:
        tblite = args.tblite.resolve()
        payload["save_tblite"] = {
            "path": str(tblite),
            "sha256": sha256(tblite),
            "version": command_output([str(tblite), "--version"]),
            "source": git_metadata(args.tblite_source),
            "source_revision_validated": identity.tblite_source_revision,
            "static_library": str(identity.tblite_static_library),
            "static_library_sha256": identity.tblite_static_library_sha256,
        }
    provenance = (
        args.provenance
        or args.root / "data" / "build_provenance_gxtb_spglib.json"
    )
    atomic_write_bytes(
        provenance,
        (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode(),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--campaign-manifest",
        type=Path,
        default=DEFAULT_CAMPAIGN_MANIFEST,
    )
    parser.add_argument(
        "--execution-build-manifest",
        type=Path,
        help=(
            "explicit qualification manifest for a same-source build whose "
            "binary/library hashes differ from the frozen local build"
        ),
    )
    parser.add_argument(
        "--base-validation-index",
        type=Path,
        help=(
            "immutable schema-v1/v2 snapshot to verify, preserve, and merge "
            "with records produced by this invocation"
        ),
    )
    parser.add_argument(
        "--base-validation-index-sha256",
        dest="base_validation_index_expected_sha256",
        help="required exact SHA256 pin for --base-validation-index",
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument(
        "--cp2k-library",
        type=Path,
        help="optional expected libcp2k dylib; must equal the otool/RPATH result",
    )
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--tblite-static-library", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--tblite-source", type=Path, required=True)
    parser.add_argument(
        "--method",
        action="append",
        choices=["GXTB"],
        help="this additive production runner accepts GXTB only",
    )
    parser.add_argument(
        "--gxtb-input-root",
        type=Path,
        help=f"separate production input root (default: ROOT/{GXTB_INPUT_DIRECTORY})",
    )
    parser.add_argument(
        "--gxtb-run-root",
        type=Path,
        help=f"separate production run root (default: ROOT/{GXTB_RUN_DIRECTORY})",
    )
    parser.add_argument(
        "--analysis-prefix",
        default=GXTB_ANALYSIS_PREFIX,
        help="prefix for additive result and figure files",
    )
    parser.add_argument(
        "--mesh",
        action="append",
        choices=SUPPORTED_MESHES,
        help=(
            "mesh(es) to run; repeat as needed (default: frozen six-mesh core; "
            "k666 through k131313 are opt-in convergence extensions)"
        ),
    )
    parser.add_argument(
        "--phase",
        action="append",
        choices=PHASES,
        help="phase(s) to run; repeat as needed (default: every DMC13 phase)",
    )
    parser.add_argument("--provenance", type=Path)
    parser.add_argument("--jobs", type=positive_int, default=4)
    parser.add_argument(
        "--threads-per-job",
        type=positive_int,
        default=1,
        help=(
            "OpenMP threads per directly launched CP2K process (default: 1); "
            "production MPI execution requires exactly one"
        ),
    )
    parser.add_argument("--mpi-ranks-per-job", type=positive_int, default=1)
    parser.add_argument("--mpi-launcher", type=Path)
    parser.add_argument(
        "--mpi-launcher-arg",
        action="append",
        default=[],
        help=(
            "reserved compatibility option; production execution rejects every "
            "user-supplied MPI launcher argument"
        ),
    )
    parser.add_argument(
        "--pe-list",
        action="append",
        default=[],
        help=(
            "literal ordered CPU list (for example 96,97,98,99); repeat exactly "
            "once per --jobs worker"
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="retained for compatibility; completed jobs are always resumed safely",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="archive prior production files and recompute selected GXTB jobs",
    )
    args = parser.parse_args()
    if (args.base_validation_index is None) != (
        args.base_validation_index_expected_sha256 is None
    ):
        parser.error(
            "--base-validation-index and --base-validation-index-sha256 "
            "must be provided together"
        )
    if (
        args.base_validation_index_expected_sha256 is not None
        and not re.fullmatch(
            r"[0-9a-f]{64}", args.base_validation_index_expected_sha256
        )
    ):
        parser.error(
            "--base-validation-index-sha256 must be 64 lowercase hex characters"
        )
    args.root = args.root.resolve()
    args.campaign_manifest = args.campaign_manifest.resolve()
    if args.execution_build_manifest is not None:
        args.execution_build_manifest = args.execution_build_manifest.resolve()
    if args.base_validation_index is not None:
        args.base_validation_index = args.base_validation_index.resolve()
    args.cp2k = args.cp2k.resolve()
    args.tblite = args.tblite.resolve()
    args.tblite_static_library = args.tblite_static_library.resolve()
    args.cp2k_source = args.cp2k_source.resolve()
    args.tblite_source = args.tblite_source.resolve()
    args.gxtb_input_root = (
        args.gxtb_input_root or args.root / GXTB_INPUT_DIRECTORY
    ).resolve()
    args.gxtb_run_root = (
        args.gxtb_run_root or args.root / GXTB_RUN_DIRECTORY
    ).resolve()

    execution_requested = bool(
        args.mpi_launcher
        or args.mpi_launcher_arg
        or args.pe_list
        or args.mpi_ranks_per_job != 1
    )
    args.execution_pool = None
    if execution_requested:
        if args.mpi_launcher is None:
            parser.error("--mpi-launcher is required with MPI/affinity execution")
        try:
            args.execution_pool = benchmark_execution.ExecutionPool(
                concurrent_jobs=args.jobs,
                mpi_ranks_per_job=args.mpi_ranks_per_job,
                threads_per_rank=args.threads_per_job,
                mpi_launcher=args.mpi_launcher,
                mpi_launcher_args=args.mpi_launcher_arg,
                pe_lists=args.pe_list,
            )
        except (OSError, ValueError) as error:
            parser.error(str(error))

    try:
        validate_production_paths(
            args.root,
            args.analysis_prefix,
            args.gxtb_input_root,
            args.gxtb_run_root,
        )
    except ValueError as error:
        parser.error(str(error))

    script = args.root / "scripts" / "dmc_ice13_kpoint_benchmark.py"
    if not args.method:
        parser.error("select --method GXTB explicitly")
    try:
        methods = require_unique_selection(args.method, "method")
        requested_meshes = require_unique_selection(args.mesh or MESHES, "mesh")
        requested_phases = require_unique_selection(args.phase or PHASES, "phase")
    except ValueError as error:
        parser.error(str(error))
    if methods != ["GXTB"]:
        parser.error("select --method GXTB exactly once")
    if args.force and set(methods) != {"GXTB"}:
        parser.error("--force is restricted to an explicitly selected GXTB run")
    if args.tblite.parent.parent != args.tblite_static_library.parent.parent:
        parser.error(
            "--tblite and --tblite-static-library must come from the same install prefix"
        )
    try:
        lock_handle = acquire_runner_lock(
            args.gxtb_run_root / ".dmc13-runner.lock",
            {
                "pid": os.getpid(),
                "argv": [sys.executable, *sys.argv],
                "root": str(args.root),
                "gxtb_run_root": str(args.gxtb_run_root),
            },
        )
    except ValueError as error:
        parser.error(str(error))
    release_callback = lambda: release_runner_lock(lock_handle)
    atexit.register(release_callback)
    try:
        campaign_manifest = read_campaign_manifest(args.campaign_manifest)
        save_tblite_manifest = campaign_manifest.get("save_tblite")
        if not isinstance(save_tblite_manifest, dict):
            raise ValueError("campaign manifest lacks save_tblite identity block")
        identity = production_identity(
            str(campaign_manifest.get("campaign_id")),
            args.cp2k,
            args.cp2k_library,
            args.tblite_static_library,
            args.cp2k_source,
            args.tblite_source,
            str(save_tblite_manifest.get("revision")),
            require_embedded_tblite_revision=(
                args.execution_build_manifest is not None
            ),
        )
        if args.execution_build_manifest is None:
            validate_campaign_identity(identity, args.tblite, campaign_manifest)
    except ValueError as error:
        parser.error(f"build identity gate failed: {error}")
    args.base_validation_index_payload = None
    args.base_validation_index_sha256 = None
    if args.base_validation_index is not None:
        try:
            actual_base_hash = sha256(args.base_validation_index)
            if actual_base_hash != args.base_validation_index_expected_sha256:
                raise ValueError(
                    "base validation index SHA256 pin mismatch: "
                    f"actual {actual_base_hash}, expected "
                    f"{args.base_validation_index_expected_sha256}"
                )
            args.base_validation_index_sha256 = actual_base_hash
            args.base_validation_index_payload = read_validation_index(
                args.base_validation_index,
                args.root,
                expected_campaign_id=identity.campaign_id,
                expected_source_identity={
                    "cp2k_source_revision": identity.cp2k_source_revision,
                    "tblite_source_revision": identity.tblite_source_revision,
                },
                expected_campaign_manifest_sha256=sha256(args.campaign_manifest),
                campaign_manifest_path=args.campaign_manifest,
                expected_index_sha256=actual_base_hash,
            )
        except ValueError as error:
            parser.error(f"base validation index gate failed: {error}")

    if args.execution_build_manifest is not None:
        if not isinstance(args.base_validation_index_payload, dict):
            parser.error(
                "an alternate execution build requires an explicitly SHA256-pinned "
                "--base-validation-index"
            )
        frozen_reference_id = build_id(
            frozen_build_identity_from_manifest(campaign_manifest)
        )
        reference_records = {
            validation_record_key(record): record
            for record in args.base_validation_index_payload.get("records", [])
            if isinstance(record, dict)
            and record.get("build_id") == frozen_reference_id
        }
        try:
            validate_execution_build_manifest(
                identity,
                args.tblite,
                args.campaign_manifest,
                campaign_manifest,
                args.execution_build_manifest,
                args.cp2k_source,
                args.tblite_source,
                args.root,
                reference_records,
            )
        except ValueError as error:
            parser.error(f"build identity gate failed: {error}")

    failures: list[tuple[Job, int]] = []
    all_jobs = jobs(
        args.root,
        methods,
        requested_meshes,
        requested_phases,
        gxtb_input_root=args.gxtb_input_root,
        gxtb_run_root=args.gxtb_run_root,
    )
    base_keys = {
        validation_record_key(record)
        for record in (
            args.base_validation_index_payload.get("records", [])
            if isinstance(args.base_validation_index_payload, dict)
            else []
        )
        if isinstance(record, dict)
    }
    selected_base_keys = {
        (job.mesh, job.phase) for job in all_jobs if (job.mesh, job.phase) in base_keys
    }
    if args.force and selected_base_keys:
        parser.error(
            "--force cannot replace records preserved by --base-validation-index: "
            + ", ".join(f"{mesh}/{phase}" for mesh, phase in sorted(selected_base_keys))
        )
    if selected_base_keys:
        print(
            f"Preserving {len(selected_base_keys)} verified base-index job(s); "
            "only missing logical jobs will run.",
            flush=True,
        )
        all_jobs = [
            job for job in all_jobs if (job.mesh, job.phase) not in selected_base_keys
        ]
    for mesh in requested_meshes:
        missing_phases = [
            phase
            for phase in requested_phases
            if any(job.mesh == mesh and job.phase == phase for job in all_jobs)
        ]
        if not missing_phases:
            continue
        prepare_command = [
            sys.executable,
            str(script),
            "prepare",
            "--method",
            "GXTB",
            "--gxtb-input-root",
            str(args.gxtb_input_root),
            "--mesh",
            mesh,
        ]
        for phase in missing_phases:
            prepare_command += ["--phase", phase]
        subprocess.run(prepare_command, cwd=args.root, check=True)
    contract_failures = [
        (job, gxtb_input_contract_errors(job))
        for job in all_jobs
        if gxtb_input_contract_errors(job)
    ]
    if contract_failures:
        details = "\n".join(
            f"  {job.mesh}/{job.phase}: {'; '.join(errors)}"
            for job, errors in contract_failures[:20]
        )
        raise SystemExit(f"Invalid g-xTB production inputs:\n{details}")
    stop_event = threading.Event()
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs)
    futures = [
        pool.submit(
            run_job,
            identity,
            job,
            args.force,
            stop_event,
            args.threads_per_job,
            args.execution_pool,
        )
        for job in all_jobs
    ]
    try:
        done = 0
        for future in concurrent.futures.as_completed(futures):
            job, rc = future.result()
            done += 1
            if rc != 0:
                failures.append((job, rc))
            print(
                f"{done:3d}/{len(all_jobs)} {job.mesh:5s} "
                f"{job.method:4s} {job.phase:4s} rc={rc}",
                flush=True,
            )
    except KeyboardInterrupt:
        stop_event.set()
        for future in futures:
            future.cancel()
        pool.shutdown(wait=False, cancel_futures=True)
        print(
            "Interrupted: running CP2K jobs were stopped and queued jobs were cancelled.",
            file=sys.stderr,
        )
        raise SystemExit(130)
    else:
        pool.shutdown(wait=True)
    validation_index_path = write_convergence_validation_index(args, identity)
    args.validation_index_snapshot = validation_index_path
    if failures:
        for job, rc in failures[:20]:
            print(
                f"FAILED {job.mesh} {job.method} {job.phase} rc={rc} "
                f"out={job.run_dir / job.output_name}",
                file=sys.stderr,
            )
        completed = sum(
            len(validated_gxtb_phase_coverage(args, identity).get(mesh, []))
            for mesh in MESHES
        )
        write_provenance(
            args,
            identity,
            methods,
            requested_meshes,
            requested_phases,
            completed,
            failures,
        )
        raise SystemExit(f"{len(failures)} DMC13 jobs failed")

    analyse_command = [
        sys.executable,
        str(script),
        "analyse",
        "--restrict-gxtb-to-validated",
        "--validation-index",
        str(validation_index_path),
        "--validation-index-sha256",
        sha256(validation_index_path),
        "--gxtb-run-root",
        str(args.gxtb_run_root),
        "--output-prefix",
        args.analysis_prefix,
    ]
    subprocess.run(analyse_command, cwd=args.root, check=True)
    completed = sum(
        len(validated_gxtb_phase_coverage(args, identity).get(mesh, []))
        for mesh in MESHES
    )
    write_provenance(
        args,
        identity,
        methods,
        requested_meshes,
        requested_phases,
        completed,
    )
    print(
        args.root
        / "data"
        / f"dmc_ice13_{args.analysis_prefix}_kpoint_stats.csv"
    )
    atexit.unregister(release_callback)
    release_runner_lock(lock_handle)


if __name__ == "__main__":
    main()
