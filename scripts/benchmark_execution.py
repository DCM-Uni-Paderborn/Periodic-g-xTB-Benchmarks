#!/usr/bin/env python3
"""Fail-closed MPI rank binding and execution records for CP2K benchmarks.

Production MPI jobs use Open MPI's ordered PE-list mapper.  Every rank is
assigned one explicit, unique logical CPU; the assignment is then verified
from ``/proc/<pid>/environ`` and ``/proc/<pid>/status`` before a schema-v2
execution record can be finalized.  Schema-v1 taskset records remain readable
for numerical provenance, but their timings are classified as non-scaling.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import queue
import re
import signal
import shutil
import socket
import subprocess
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Mapping, Sequence


SCHEMA_VERSION = 2
LEGACY_SCHEMA_VERSION = 1
BLAS_THREAD_ENVIRONMENT = {
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
OPENMP_THREAD_ENVIRONMENT = {
    "OMP_NUM_THREADS": "1",
    "OMP_PROC_BIND": "true",
    "OMP_PLACES": "cores",
    "OMP_SCHEDULE": "static",
    "OMP_DYNAMIC": "FALSE",
    "OMP_WAIT_POLICY": "PASSIVE",
}
SANCTIONED_THREAD_ENVIRONMENT = {
    **OPENMP_THREAD_ENVIRONMENT,
    **BLAS_THREAD_ENVIRONMENT,
}
INHERITED_BINDING_ENVIRONMENT_PREFIXES = (
    "OMPI_MCA_",
    "PRTE_MCA_",
    "PMIX_MCA_",
    "HWLOC_",
    "GOMP_",
    "KMP_",
    "OMP_",
)
INHERITED_BINDING_ENVIRONMENT_KEYS = frozenset(
    {
        "LD_AUDIT",
        "LD_PRELOAD",
        "I_MPI_PIN",
        "I_MPI_PIN_DOMAIN",
        "I_MPI_PIN_ORDER",
        "I_MPI_PIN_PROCESSOR_LIST",
        "I_MPI_PIN_RESPECT_CPUSET",
        "MPICH_CPU_BINDING_POLICY",
        "MPICH_RANK_REORDER_METHOD",
        "MV2_CPU_BINDING_POLICY",
        "MV2_CPU_MAPPING",
        "MV2_ENABLE_AFFINITY",
        "SLURM_CPU_BIND",
        "SLURM_CPU_BIND_LIST",
        "SLURM_CPU_BIND_TYPE",
    }
)
THREAD_AFFINITY_EVIDENCE_SOURCE = "linux_proc_task_status"
THREAD_AFFINITY_ACCEPTED_PROCESS_STATUSES = frozenset(
    {"live", "process_disappeared_after_sample", "terminal_process"}
)


def default_cpu_reservation_lock_root() -> Path:
    return Path(f"/tmp/periodic-gxtb-cpu-reservations-{os.getuid()}")


def acquire_cpu_reservation_locks(
    cpus: Sequence[int], lock_root: Path
) -> list[IO[str]]:
    """Reserve each logical CPU across independent benchmark processes."""
    lock_root.mkdir(parents=True, exist_ok=True)
    handles: list[IO[str]] = []
    current_handle: IO[str] | None = None
    try:
        for cpu in sorted(set(cpus)):
            path = lock_root / f"cpu-{cpu}.lock"
            current_handle = path.open("a+")
            try:
                fcntl.flock(
                    current_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB
                )
            except BlockingIOError as error:
                current_handle.seek(0)
                holder = current_handle.read().strip() or "unidentified holder"
                raise ValueError(
                    f"logical CPU {cpu} is already reserved by another production "
                    f"launcher ({holder})"
                ) from error
            # Register the just-locked file before metadata I/O.  The separate
            # current_handle reference also covers a BaseException in append.
            handles.append(current_handle)
            handle = current_handle
            current_handle = None
            handle.seek(0)
            handle.truncate()
            json.dump(
                {
                    "cpu": cpu,
                    "hostname": socket.gethostname(),
                    "pid": os.getpid(),
                },
                handle,
                sort_keys=True,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        if current_handle is not None:
            current_handle.close()
        for handle in handles:
            handle.close()
        raise
    return handles


def mpi_control_environment_keys(environment: Mapping[str, str]) -> list[str]:
    """Return inherited variables that can perturb placement or preloading.

    The production command supplies its complete placement contract explicitly.
    Consequently, MPI/PMIx MCA controls, topology overrides, OpenMP runtime
    placement controls, vendor affinity controls, and loader injection are all
    removed before the narrow sanctioned thread environment is installed.
    ``LD_LIBRARY_PATH`` is deliberately not part of this set because the pinned
    CP2K/provider runtime may require it.
    """
    return sorted(
        key
        for key in environment
        if key in INHERITED_BINDING_ENVIRONMENT_KEYS
        or key.startswith(INHERITED_BINDING_ENVIRONMENT_PREFIXES)
    )


def sanitized_launch_environment(
    environment: Mapping[str, str],
) -> tuple[dict[str, str], list[str], bool]:
    """Build a hermetic rank environment while preserving library discovery."""
    sanitized = dict(environment)
    removed = mpi_control_environment_keys(sanitized)
    for key in removed:
        sanitized.pop(key, None)
    sanitized.update(SANCTIONED_THREAD_ENVIRONMENT)
    residual = [
        key
        for key in mpi_control_environment_keys(sanitized)
        if key not in SANCTIONED_THREAD_ENVIRONMENT
    ]
    if residual:
        raise RuntimeError(
            "binding/preload environment survived sanitization: "
            + ", ".join(residual)
        )
    library_path_preserved = sanitized.get("LD_LIBRARY_PATH") == environment.get(
        "LD_LIBRARY_PATH"
    ) and ("LD_LIBRARY_PATH" in sanitized) == ("LD_LIBRARY_PATH" in environment)
    if not library_path_preserved:
        raise RuntimeError("LD_LIBRARY_PATH changed during environment sanitization")
    return sanitized, removed, library_path_preserved


def binding_environment_scrub_contract() -> dict[str, object]:
    """Return the immutable environment policy embedded in execution records."""
    return {
        "removed_prefixes": list(INHERITED_BINDING_ENVIRONMENT_PREFIXES),
        "removed_exact_keys": sorted(INHERITED_BINDING_ENVIRONMENT_KEYS),
        "sanctioned_thread_environment": dict(SANCTIONED_THREAD_ENVIRONMENT),
        "preserved_exact_keys": ["LD_LIBRARY_PATH"],
    }


def _contract_requires_hardened_binding_evidence(
    contract: Mapping[str, object],
) -> bool:
    """Distinguish new schema-v2 records without invalidating older evidence."""
    return (
        contract.get("binding_environment_scrub_contract")
        == binding_environment_scrub_contract()
        and contract.get("pool_close_policy")
        == "reject_while_run_admitted_or_active"
        and contract.get("rank_affinity_observation")
        == "linux_proc_per_task_tid_starttime"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
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


def resolve_executable(value: str | Path, label: str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.parent == Path("."):
        discovered = shutil.which(str(candidate))
        if discovered is None:
            raise ValueError(f"{label} is not on PATH: {value}")
        candidate = Path(discovered)
    resolved = candidate.resolve(strict=True)
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise ValueError(f"{label} is not executable: {resolved}")
    return resolved


def parse_cpu_set(value: str) -> set[int]:
    cpus: set[int] = set()
    for raw_part in value.split(","):
        part = raw_part.strip()
        if not part:
            raise ValueError(f"invalid empty CPU-set component in {value!r}")
        if "-" in part:
            fields = part.split("-")
            if len(fields) != 2 or not all(field.isdigit() for field in fields):
                raise ValueError(f"invalid CPU-set range {part!r}")
            first, last = (int(field) for field in fields)
            if last < first:
                raise ValueError(f"descending CPU-set range {part!r}")
            selected = set(range(first, last + 1))
        else:
            if not part.isdigit():
                raise ValueError(f"invalid CPU-set index {part!r}")
            selected = {int(part)}
        overlap = cpus & selected
        if overlap:
            raise ValueError(
                f"duplicate CPU index/overlap inside {value!r}: {sorted(overlap)}"
            )
        cpus.update(selected)
    if not cpus:
        raise ValueError("CPU set must not be empty")
    return cpus


def parse_ordered_pe_list(value: str) -> tuple[int, ...]:
    """Parse a literal ordered list such as ``96,97,98,99``.

    Range notation is deliberately rejected: the rank-to-CPU mapping must be
    explicit in both the CLI and the immutable execution contract.
    """
    raw = value.split(",")
    if not raw or any(not field.strip().isdigit() for field in raw):
        raise ValueError(
            f"ordered PE list {value!r} must contain comma-separated CPU indices"
        )
    cpus = tuple(int(field.strip()) for field in raw)
    duplicates = sorted({cpu for cpu in cpus if cpus.count(cpu) > 1})
    if duplicates:
        raise ValueError(
            f"ordered PE list {value!r} contains duplicate CPUs: {duplicates}"
        )
    return cpus


def normalize_ordered_pe_list(value: str) -> str:
    return ",".join(str(cpu) for cpu in parse_ordered_pe_list(value))


def validate_mpi_launcher_args(arguments: Sequence[str]) -> None:
    """Reject every user argument that could alter the immutable MPI contract.

    Open MPI exposes many aliases, MCA control paths, appfiles, and multi-app
    separators.  A deny-list cannot prove that none of them changes mapping or
    binding, so production benchmark launchers accept no extra arguments.  The
    complete launcher policy is injected by :class:`ExecutionPool`.
    """
    if arguments:
        raise ValueError(
            "production MPI accepts no user-supplied launcher arguments; "
            "rank count, mapping, binding, and binding reports are injected "
            "by the execution contract"
        )


def validate_pe_lists(
    values: Sequence[str],
    concurrent_jobs: int,
    mpi_ranks_per_job: int,
    threads_per_rank: int,
    *,
    available_cpus: set[int] | None = None,
) -> list[tuple[int, ...]]:
    if concurrent_jobs < 1 or mpi_ranks_per_job < 1 or threads_per_rank < 1:
        raise ValueError("jobs, MPI ranks, and OpenMP threads must be positive")
    if threads_per_rank != 1:
        raise ValueError(
            "production MPI affinity requires exactly one OpenMP thread per rank"
        )
    if len(values) != concurrent_jobs:
        raise ValueError(
            "exactly one --pe-list is required per concurrent job: "
            f"got {len(values)}, expected {concurrent_jobs}"
        )
    parsed = [parse_ordered_pe_list(value) for value in values]
    for value, cpus in zip(values, parsed, strict=True):
        if len(cpus) != mpi_ranks_per_job:
            raise ValueError(
                f"ordered PE list {value!r} has {len(cpus)} CPUs; expected exactly "
                f"{mpi_ranks_per_job}, one per MPI rank"
            )
        if available_cpus is not None and not set(cpus) <= available_cpus:
            raise ValueError(
                f"ordered PE list {value!r} requests unavailable CPUs: "
                f"{sorted(set(cpus) - available_cpus)}"
            )
    for index, left in enumerate(parsed):
        for prior_index, right in enumerate(parsed[:index]):
            overlap = set(left) & set(right)
            if overlap:
                raise ValueError(
                    f"ordered PE lists {values[prior_index]!r} and {values[index]!r} "
                    f"overlap at {sorted(overlap)}"
                )
    return parsed


def require_single_pu_cores(
    pe_lists: Sequence[Sequence[int]],
    topology_root: Path = Path("/sys/devices/system/cpu"),
) -> None:
    """Fail before launch when ``--bind-to core`` can yield a non-singleton mask."""
    if not topology_root.is_dir():
        return
    for cpu in (cpu for pe_list in pe_lists for cpu in pe_list):
        siblings_path = topology_root / f"cpu{cpu}" / "topology" / "thread_siblings_list"
        try:
            siblings = parse_cpu_set(siblings_path.read_text().strip())
        except (OSError, ValueError) as error:
            raise ValueError(
                f"cannot prove singleton core topology for CPU {cpu}: {error}"
            ) from error
        if len(siblings) != 1:
            raise ValueError(
                f"CPU {cpu} belongs to SMT siblings {sorted(siblings)}; "
                "'--bind-to core' cannot prove one logical CPU per rank"
            )


MPI_RANK_ENVIRONMENT_KEYS = (
    "OMPI_COMM_WORLD_RANK",
    "PMI_RANK",
    "PMIX_RANK",
    "SLURM_PROCID",
    "MV2_COMM_WORLD_RANK",
)


def live_compute_cpu_owners(
    cpus: Sequence[int],
    proc_root: Path = Path("/proc"),
    *,
    ignore_process_identities: Mapping[int, int] | None = None,
) -> list[dict[str, object]]:
    """Find live CP2K or MPI-rank processes whose allowed masks overlap *cpus*.

    Advisory lock files prevent two cooperating benchmark launchers from using
    the same logical CPU.  Historical and third-party launchers do not take
    those locks, so Linux production runs also inspect procfs before launch.
    Launchers and daemons with host-wide masks are ignored unless their
    environment identifies them as an MPI rank; live CP2K processes are always
    considered owners.  Zombie processes own no schedulable CPU and are ignored.
    """
    selected = set(cpus)
    ignored_identities = dict(ignore_process_identities or {})
    if not selected or not proc_root.is_dir():
        return []
    owners: list[dict[str, object]] = []
    for directory in proc_root.iterdir():
        if not directory.name.isdigit():
            continue
        pid = int(directory.name)
        if pid == os.getpid():
            continue
        try:
            initial_stat_state, initial_starttime = _linux_proc_stat_identity(
                (directory / "stat").read_text(errors="replace")
            )
            status = (directory / "status").read_text(errors="replace")
        except (FileNotFoundError, ProcessLookupError):
            continue
        except (PermissionError, OSError, ValueError):
            if pid in ignored_identities:
                owners.append(
                    {
                        "pid": pid,
                        "name": "",
                        "state": "",
                        "cpus_allowed_list": "",
                        "overlap": sorted(selected),
                        "cp2k_process": False,
                        "mpi_rank_process": False,
                        "process_identity_status": "initial_identity_unreadable",
                    }
                )
            continue
        fields: dict[str, str] = {}
        for line in status.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key] = value.strip()
        state = fields.get("State", "")
        mask_text = fields.get("Cpus_allowed_list", "")
        try:
            allowed = parse_cpu_set(mask_text)
        except ValueError:
            allowed = set()
        overlap = selected & allowed
        name = fields.get("Name", "")
        is_cp2k = name.casefold().startswith("cp2k")
        is_mpi_rank = False
        try:
            environment_items = (directory / "environ").read_bytes().split(b"\0")
            environment_keys = {
                item.split(b"=", 1)[0].decode(errors="replace")
                for item in environment_items
                if item and b"=" in item
            }
            is_mpi_rank = bool(
                environment_keys.intersection(MPI_RANK_ENVIRONMENT_KEYS)
            )
        except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
            environment_keys = set()
        try:
            final_stat_state, final_starttime = _linux_proc_stat_identity(
                (directory / "stat").read_text(errors="replace")
            )
        except (FileNotFoundError, ProcessLookupError):
            continue
        except (PermissionError, OSError, ValueError):
            if overlap or pid in ignored_identities:
                owners.append(
                    {
                        "pid": pid,
                        "name": name,
                        "state": state,
                        "cpus_allowed_list": mask_text,
                        "overlap": sorted(selected),
                        "cp2k_process": is_cp2k,
                        "mpi_rank_process": is_mpi_rank,
                        "process_identity_status": "final_identity_unreadable",
                    }
                )
            continue
        if final_starttime != initial_starttime:
            owners.append(
                {
                    "pid": pid,
                    "name": name,
                    "state": state,
                    "cpus_allowed_list": mask_text,
                    "overlap": sorted(selected),
                    "cp2k_process": is_cp2k,
                    "mpi_rank_process": is_mpi_rank,
                    "process_identity_status": "pid_reused_during_scan",
                    "initial_process_starttime": initial_starttime,
                    "final_process_starttime": final_starttime,
                }
            )
            continue
        expected_starttime = ignored_identities.get(pid)
        if expected_starttime == initial_starttime:
            continue
        state = fields.get("State", "")
        if (
            state.startswith(("Z", "X"))
            or final_stat_state in {"Z", "X"}
            or initial_stat_state in {"Z", "X"}
        ):
            continue
        if not overlap:
            continue
        if not is_cp2k and not is_mpi_rank:
            continue
        owners.append(
            {
                "pid": pid,
                "name": name,
                "state": state,
                "cpus_allowed_list": mask_text,
                "overlap": sorted(overlap),
                "cp2k_process": is_cp2k,
                "mpi_rank_process": is_mpi_rank,
                "process_starttime": initial_starttime,
                "process_identity_status": "stable",
            }
        )
    return sorted(owners, key=lambda owner: int(owner["pid"]))


def require_no_live_compute_overlap(
    cpus: Sequence[int], proc_root: Path = Path("/proc")
) -> None:
    owners = live_compute_cpu_owners(cpus, proc_root)
    if not owners:
        return
    description = "; ".join(
        f"PID {owner['pid']} ({owner['name']}) mask "
        f"{owner['cpus_allowed_list']} overlaps {owner['overlap']}"
        for owner in owners
    )
    raise ValueError(
        "selected logical CPUs are already owned by live CP2K/MPI ranks: "
        + description
    )


def validate_cpu_sets(
    values: Sequence[str],
    concurrent_jobs: int,
    mpi_ranks_per_job: int,
    threads_per_rank: int,
    *,
    available_cpus: set[int] | None = None,
) -> list[set[int]]:
    if concurrent_jobs < 1 or mpi_ranks_per_job < 1 or threads_per_rank < 1:
        raise ValueError("jobs, MPI ranks, and OpenMP threads must be positive")
    if len(values) != concurrent_jobs:
        raise ValueError(
            "exactly one --cpu-set is required per concurrent job: "
            f"got {len(values)}, expected {concurrent_jobs}"
        )
    parsed = [parse_cpu_set(value) for value in values]
    required = mpi_ranks_per_job * threads_per_rank
    for value, cpus in zip(values, parsed, strict=True):
        if len(cpus) < required:
            raise ValueError(
                f"CPU set {value!r} has {len(cpus)} CPUs; execution requires at "
                f"least {required} (= ranks times threads)"
            )
        if available_cpus is not None and not cpus <= available_cpus:
            raise ValueError(
                f"CPU set {value!r} requests unavailable CPUs: "
                f"{sorted(cpus - available_cpus)}"
            )
    for index, left in enumerate(parsed):
        for prior_index, right in enumerate(parsed[:index]):
            overlap = left & right
            if overlap:
                raise ValueError(
                    f"CPU sets {values[prior_index]!r} and {values[index]!r} "
                    f"overlap at {sorted(overlap)}"
                )
    return parsed


def execution_record_path(output: Path) -> Path:
    return output.with_suffix(output.suffix + ".execution.json")


def cp2k_command(
    *,
    taskset: str | Path,
    cpu_set: str,
    mpi_launcher: str | Path,
    mpi_launcher_args: Sequence[str],
    mpi_ranks_per_job: int,
    cp2k: Path,
    inp: Path,
    out: Path,
) -> list[str]:
    """Build the exact launch command with location-independent I/O paths."""
    return [
        str(taskset),
        "-c",
        cpu_set,
        str(mpi_launcher),
        *mpi_launcher_args,
        "-np",
        str(mpi_ranks_per_job),
        str(cp2k.resolve(strict=True)),
        "-i",
        str(inp.resolve(strict=True)),
        "-o",
        str(out.resolve()),
    ]


def launcher_log_path(output: Path) -> Path:
    return output.with_suffix(output.suffix + ".launcher.log")


def execution_record_timing_classification(
    path: Path,
    expected_contract: Mapping[str, object] | None = None,
    output: Path | None = None,
    scientific_job_stamp: Path | None = None,
) -> str:
    """Classify timing reuse without modifying a historical record.

    A schema-v2 timing is eligible only after complete command, rank-mask,
    launcher-log, executable, input, output, and stamp revalidation.
    """
    try:
        record = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return "invalid_execution_record"
    if record.get("schema_version") == LEGACY_SCHEMA_VERSION:
        return "legacy_timing_non_scaling"
    if record.get("schema_version") == SCHEMA_VERSION:
        if (
            expected_contract is None
            or output is None
            or scientific_job_stamp is None
        ):
            return "timing_requires_full_revalidation"
        if (
            recorded_execution_issue(
                path, expected_contract, output, scientific_job_stamp
            )
            is None
        ):
            return "production_scaling_eligible"
    return "timing_non_scaling"


def _common_artifact_issue(
    record: Mapping[str, object],
    path: Path,
    output: Path,
    scientific_job_stamp: Path,
) -> str | None:
    cp2k_path = Path(str(record.get("cp2k", "")))
    try:
        cp2k_resolved = cp2k_path.resolve(strict=True)
    except (FileNotFoundError, OSError):
        return f"missing recorded CP2K executable in {path}"
    if not cp2k_resolved.is_file() or not os.access(cp2k_resolved, os.X_OK):
        return f"recorded CP2K executable is not executable in {path}"
    cp2k_hash = sha256(cp2k_resolved)
    if record.get("cp2k_sha256_at_launch") != cp2k_hash:
        return f"CP2K executable hash mismatch in {path}"
    input_path = Path(str(record.get("input", "")))
    try:
        input_resolved = input_path.resolve(strict=True)
    except (FileNotFoundError, OSError):
        return f"missing recorded execution input in {path}"
    input_hash = sha256(input_resolved)
    if record.get("input_sha256_at_launch") != input_hash:
        return f"execution input hash mismatch in {path}"
    if record.get("working_directory") not in (None, str(input_resolved.parent)):
        return f"execution working-directory mismatch in {path}"
    if Path(str(record.get("output", ""))).resolve() != output.resolve():
        return f"execution output path mismatch in {path}"
    if not output.is_file() or record.get("output_sha256") != sha256(output):
        return f"execution output hash mismatch in {path}"
    if (
        Path(str(record.get("scientific_job_stamp", ""))).resolve()
        != scientific_job_stamp.resolve()
    ):
        return f"scientific job-stamp path mismatch in {path}"
    if (
        not scientific_job_stamp.is_file()
        or record.get("scientific_job_stamp_sha256") != sha256(scientific_job_stamp)
    ):
        return f"scientific job-stamp hash mismatch in {path}"
    try:
        signature = json.loads(scientific_job_stamp.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid scientific job stamp for {path}: {exc}"
    if not isinstance(signature, Mapping):
        return f"invalid scientific job-stamp payload for {path}"
    stamped_cp2k_value = signature.get("executable", signature.get("cp2k"))
    stamped_cp2k_hash = signature.get(
        "executable_sha256", signature.get("cp2k_sha256")
    )
    if stamped_cp2k_value is None or stamped_cp2k_hash != cp2k_hash:
        return f"scientific job-stamp CP2K identity mismatch in {path}"
    try:
        stamped_cp2k = Path(str(stamped_cp2k_value)).resolve(strict=True)
    except (FileNotFoundError, OSError):
        return f"scientific job stamp has no reusable CP2K executable for {path}"
    if stamped_cp2k != cp2k_resolved:
        return f"recorded CP2K executable differs from scientific job stamp in {path}"
    if signature.get("input_sha256") != input_hash:
        return f"scientific job-stamp input hash mismatch in {path}"
    stamped_input_value = signature.get("input")
    if isinstance(stamped_input_value, str) and Path(stamped_input_value).is_absolute():
        try:
            stamped_input = Path(stamped_input_value).resolve(strict=True)
        except (FileNotFoundError, OSError):
            return f"scientific job stamp has no reusable input for {path}"
        if stamped_input != input_resolved:
            return f"recorded input differs from scientific job stamp in {path}"
    return None


def _legacy_recorded_execution_issue(
    record: Mapping[str, object],
    path: Path,
    expected_contract: Mapping[str, object],
    output: Path,
    scientific_job_stamp: Path,
) -> str | None:
    """Validate schema v1 without endorsing its timings for scaling."""
    if expected_contract.get("schema_version") != LEGACY_SCHEMA_VERSION:
        return f"legacy execution contract mismatch in {path}"
    if record.get("contract") != expected_contract:
        return f"execution contract mismatch in {path}"
    if record.get("contract_sha256") != canonical_sha256(dict(expected_contract)):
        return f"execution contract hash mismatch in {path}"
    cpu_sets = expected_contract.get("cpu_sets")
    assigned = str(record.get("assigned_cpu_set", ""))
    if not isinstance(cpu_sets, list) or assigned not in cpu_sets:
        return f"execution record uses an unapproved CPU set in {path}"
    launcher = str(expected_contract.get("mpi_launcher", ""))
    launcher_args = expected_contract.get("mpi_launcher_args")
    taskset = str(expected_contract.get("taskset", ""))
    ranks = expected_contract.get("mpi_ranks_per_job")
    if not isinstance(launcher_args, list) or not isinstance(ranks, int):
        return f"invalid expected execution contract for {path}"
    if record.get("runtime_affinity_gate") is not True:
        return f"runtime MPI/affinity gate failed in {path}"
    if record.get("mpiexec_internal_rebinding_detected") is not False:
        return f"MPI internal rebinding was detected in {path}"
    rank_pids = record.get("observed_cp2k_rank_pids")
    rank_masks = record.get("observed_cp2k_rank_masks")
    if not isinstance(rank_pids, list) or len(rank_pids) != ranks:
        return f"observed CP2K rank count mismatch in {path}"
    try:
        assigned_cpus = parse_cpu_set(assigned)
        normalized_rank_masks = (
            [parse_cpu_set(str(mask)) for mask in rank_masks]
            if isinstance(rank_masks, list)
            else []
        )
    except ValueError:
        return f"invalid recorded CP2K rank CPU mask in {path}"
    if len(normalized_rank_masks) != ranks or any(
        mask != assigned_cpus for mask in normalized_rank_masks
    ):
        return f"observed CP2K rank CPU-mask mismatch in {path}"
    issue = _common_artifact_issue(record, path, output, scientific_job_stamp)
    if issue is not None:
        return issue
    cp2k_path = Path(str(record.get("cp2k", ""))).resolve(strict=True)
    input_path = Path(str(record.get("input", ""))).resolve(strict=True)
    expected_command = cp2k_command(
        taskset=taskset,
        cpu_set=assigned,
        mpi_launcher=launcher,
        mpi_launcher_args=launcher_args,
        mpi_ranks_per_job=ranks,
        cp2k=cp2k_path,
        inp=input_path,
        out=output,
    )
    if record.get("command") != expected_command:
        return f"full execution command/affinity mismatch in {path}"
    return None


def _rank_process_provenance_issue(
    item: Mapping[str, object], rank: int, path: Path, cp2k: Path
) -> str | None:
    """Require one immutable, terminally resolved Linux task per rank."""
    required = {
        "pid",
        "is_cp2k_rank",
        "ompi_comm_world_rank",
        "raw_ompi_comm_world_rank",
        "rank_identity_source",
        "executable",
        "arguments",
        "process_starttime",
        "observed_process_starttimes",
        "process_starttime_changed_ever",
        "process_terminally_confirmed",
        "process_terminal_confirmation",
        "process_reappeared_after_terminal_ever",
        "executable_changed_ever",
        "cpu_mask_changed_during_sample_ever",
        "snapshot_unavailable_ever",
        "process_identity_status",
        "snapshot_consistency_status",
        "stat_state",
        "state",
        "observed_process_states",
        "rank_observation_status",
        "observed_rank_observation_statuses",
    }
    if not required.issubset(item):
        return f"missing CP2K rank process provenance in {path}"
    starttime = item.get("process_starttime")
    statuses = item.get("observed_rank_observation_statuses")
    states = item.get("observed_process_states")
    observed_starttimes = item.get("observed_process_starttimes")
    confirmation = item.get("process_terminal_confirmation")
    terminal_confirmation = confirmation == "process_disappeared" or (
        isinstance(confirmation, str)
        and confirmation in {"terminal_state_Z", "terminal_state_X"}
    )
    unavailable = item.get("rank_environment_unavailable_ever") is True
    pid = item.get("pid")
    arguments = item.get("arguments")
    executable = item.get("executable")
    cp2k_text = str(cp2k)
    argument_targets: list[str] = []
    if isinstance(arguments, list):
        for argument in arguments[:2]:
            if not isinstance(argument, str):
                continue
            candidate = Path(argument)
            if not candidate.is_absolute():
                continue
            try:
                argument_targets.append(str(candidate.resolve(strict=True)))
            except (FileNotFoundError, OSError):
                continue
    executable_name = Path(str(executable)).name
    interpreted = executable_name in {"sh", "dash", "bash"} or (
        executable_name.startswith("python")
    )
    executable_proves_cp2k = bool(
        executable == cp2k_text
        or argument_targets
        and argument_targets[0] == cp2k_text
        or interpreted
        and len(argument_targets) > 1
        and argument_targets[1] == cp2k_text
    )
    if (
        not isinstance(pid, int)
        or isinstance(pid, bool)
        or item.get("is_cp2k_rank") is not True
        or not isinstance(item.get("ompi_comm_world_rank"), int)
        or isinstance(item.get("ompi_comm_world_rank"), bool)
        or item.get("ompi_comm_world_rank") != rank
        or not isinstance(executable, str)
        or not isinstance(arguments, list)
        or any(not isinstance(argument, str) for argument in arguments)
        or not executable_proves_cp2k
        or not isinstance(starttime, int)
        or isinstance(starttime, bool)
        or not isinstance(observed_starttimes, list)
        or len(observed_starttimes) != 1
        or not isinstance(observed_starttimes[0], int)
        or isinstance(observed_starttimes[0], bool)
        or observed_starttimes[0] != starttime
        or item.get("process_starttime_changed_ever") is not False
        or item.get("process_terminally_confirmed") is not True
        or not terminal_confirmation
        or item.get("process_reappeared_after_terminal_ever") is not False
        or item.get("executable_changed_ever") is not False
        or item.get("cpu_mask_changed_during_sample_ever") is not False
        or item.get("snapshot_unavailable_ever") is not False
        or item.get("process_identity_status")
        not in {"stable", "terminal_state", "disappeared_after_sample"}
        or item.get("snapshot_consistency_status")
        not in {"consistent", "process_disappeared"}
        or not isinstance(item.get("stat_state"), str)
        or not re.fullmatch(r"[RSDZTtXIWP]", str(item.get("stat_state")))
        or not isinstance(statuses, list)
        or not statuses
        or statuses[0] != "explicit"
        or item.get("rank_observation_status") != statuses[-1]
        or not isinstance(states, list)
        or not states
        or any(
            not isinstance(state, str)
            or not re.fullmatch(r"[RSDZTtXIWP](?:\s+\([^\r\n]*\))?", state)
            for state in states
        )
        or item.get("state") not in states
    ):
        return f"invalid CP2K rank process provenance in {path}"
    if unavailable:
        if (
            item.get("raw_ompi_comm_world_rank") is not None
            or item.get("rank_identity_source")
            != "pending_terminal_environment_loss"
        ):
            return f"invalid CP2K rank process provenance in {path}"
    elif (
        statuses != ["explicit"]
        or not isinstance(item.get("raw_ompi_comm_world_rank"), int)
        or isinstance(item.get("raw_ompi_comm_world_rank"), bool)
        or item.get("raw_ompi_comm_world_rank") != rank
        or item.get("rank_identity_source") != "explicit_environment"
    ):
        return f"invalid CP2K rank process provenance in {path}"
    if (
        item.get("process_identity_status") == "disappeared_after_sample"
        and (
            confirmation != "process_disappeared"
            or item.get("snapshot_consistency_status") != "process_disappeared"
        )
    ) or (
        item.get("process_identity_status") == "terminal_state"
        and (
            confirmation != f"terminal_state_{item.get('stat_state')}"
            or item.get("snapshot_consistency_status") != "consistent"
        )
    ) or (
        item.get("process_identity_status") == "stable"
        and item.get("snapshot_consistency_status") != "consistent"
    ):
        return f"inconsistent CP2K rank terminal provenance in {path}"
    return None


def _rank_environment_evidence_issue(
    item: Mapping[str, object], assigned_mask: str, path: Path
) -> str | None:
    """Validate the narrowly allowed terminal ``/proc/environ`` loss."""
    required = {
        "rank_environment_unavailable_ever",
        "rank_environment_unavailable_sample_count",
        "rank_environment_unavailable_pending",
        "rank_environment_terminally_confirmed",
        "rank_environment_terminal_confirmation",
        "rank_environment_events",
    }
    if not required.issubset(item):
        return f"missing terminal rank-environment evidence in {path}"
    unavailable_value = item.get("rank_environment_unavailable_ever")
    unavailable_count = item.get("rank_environment_unavailable_sample_count")
    if (
        unavailable_value is not True
        and unavailable_value is not False
        or not isinstance(unavailable_count, int)
        or isinstance(unavailable_count, bool)
    ):
        return f"invalid terminal rank-environment evidence in {path}"
    unavailable = unavailable_value is True
    if not unavailable:
        if (
            item.get("rank_environment_unavailable_pending") is not False
            or item.get("rank_environment_terminally_confirmed") is not False
            or item.get("rank_environment_terminal_confirmation") is not None
            or unavailable_count != 0
            or item.get("rank_environment_events") != []
        ):
            return f"inconsistent terminal rank-environment evidence in {path}"
        return None

    pid = item.get("pid")
    sample_count = item.get("sample_count")
    starttime = item.get("process_starttime")
    starttimes = item.get("observed_process_starttimes")
    statuses = item.get("observed_rank_observation_statuses")
    events = item.get("rank_environment_events")
    confirmation = item.get("rank_environment_terminal_confirmation")
    terminal_confirmation = (
        confirmation == "process_disappeared"
        or isinstance(confirmation, str)
        and confirmation in {"terminal_state_Z", "terminal_state_X"}
    )
    if (
        not isinstance(pid, int)
        or isinstance(pid, bool)
        or not isinstance(sample_count, int)
        or isinstance(sample_count, bool)
        or not isinstance(starttime, int)
        or isinstance(starttime, bool)
        or starttimes != [starttime]
        or item.get("process_starttime_changed_ever") is not False
        or item.get("rank_environment_unavailable_pending") is not False
        or item.get("rank_environment_terminally_confirmed") is not True
        or not terminal_confirmation
        or confirmation != item.get("process_terminal_confirmation")
        or item.get("rank_identity_source")
        != "pending_terminal_environment_loss"
        or item.get("raw_ompi_comm_world_rank") is not None
        or not isinstance(statuses, list)
        or not statuses
        or statuses[0] != "explicit"
        or not isinstance(events, list)
        or not events
        or item.get("rank_environment_unavailable_sample_count") != len(events)
    ):
        return f"invalid terminal rank-environment evidence in {path}"

    expected_statuses = ["explicit"]
    previous_sample_index = 0
    allowed_statuses = {"environment_empty", "environment_unreadable"}
    for event in events:
        if not isinstance(event, dict):
            return f"invalid terminal rank-environment evidence in {path}"
        event_status = event.get("environment_status")
        event_sample_index = event.get("sample_index")
        if (
            not isinstance(event_sample_index, int)
            or isinstance(event_sample_index, bool)
            or event_sample_index <= previous_sample_index
            or event_sample_index > sample_count
            or event_sample_index < 2
            or not isinstance(event.get("pid"), int)
            or isinstance(event.get("pid"), bool)
            or event.get("pid") != pid
            or not isinstance(event.get("process_starttime"), int)
            or isinstance(event.get("process_starttime"), bool)
            or event.get("process_starttime") != starttime
            or event.get("cpus_allowed_list") != assigned_mask
            or event_status not in allowed_statuses
            or event.get("terminal_resolution") != confirmation
            or event.get("state") not in item.get("observed_process_states", [])
        ):
            return f"invalid terminal rank-environment event sequence in {path}"
        previous_sample_index = event_sample_index
        if event_status not in expected_statuses:
            expected_statuses.append(str(event_status))
    if statuses != expected_statuses:
        return f"invalid terminal rank-environment status history in {path}"
    if (
        events[-1].get("sample_index") != sample_count
        or events[-1].get("state") != item.get("state")
    ):
        return f"invalid terminal rank-environment final sample in {path}"
    return None


def _thread_affinity_evidence_issue(
    item: Mapping[str, object], assigned_mask: str, path: Path
) -> str | None:
    """Revalidate sticky per-thread Linux affinity evidence for one rank."""
    required = {
        "thread_affinity_evidence_source",
        "thread_affinity_scan_status",
        "thread_affinity_scan_issues",
        "thread_affinity_process_status",
        "live_thread_affinity",
        "thread_affinity_sample_count",
        "thread_affinity_scan_statuses",
        "thread_affinity_evidence_sources",
        "thread_affinity_scan_issues_ever",
        "thread_affinity_process_statuses",
        "observed_thread_cpu_masks",
        "observed_thread_identities",
        "current_thread_affinity_sample_exact",
        "all_thread_affinity_samples_exact",
        "thread_affinity_violation_ever",
    }
    if not required.issubset(item):
        return f"missing CP2K thread-affinity evidence in {path}"
    try:
        assigned = parse_cpu_set(assigned_mask)
    except ValueError:
        return f"invalid assigned CPU in thread-affinity evidence in {path}"
    if len(assigned) != 1:
        return f"invalid assigned CPU in thread-affinity evidence in {path}"
    expected_cpu = next(iter(assigned))
    sample_count = item.get("sample_count")
    thread_sample_count = item.get("thread_affinity_sample_count")
    statuses = item.get("thread_affinity_scan_statuses")
    sources = item.get("thread_affinity_evidence_sources")
    masks = item.get("observed_thread_cpu_masks")
    identities = item.get("observed_thread_identities")
    pid = item.get("pid")
    process_starttime = item.get("process_starttime")
    process_statuses = item.get("thread_affinity_process_statuses")
    if (
        not isinstance(sample_count, int)
        or isinstance(sample_count, bool)
        or not isinstance(thread_sample_count, int)
        or isinstance(thread_sample_count, bool)
        or thread_sample_count != sample_count
        or not isinstance(statuses, list)
        or not statuses
        or statuses != ["consistent"]
        or item.get("thread_affinity_scan_status") not in statuses
        or item.get("thread_affinity_scan_issues") != []
        or item.get("thread_affinity_scan_issues_ever") != []
        or not isinstance(process_statuses, list)
        or not process_statuses
        or process_statuses != list(dict.fromkeys(process_statuses))
        or any(
            status not in THREAD_AFFINITY_ACCEPTED_PROCESS_STATUSES
            for status in process_statuses
        )
        or item.get("thread_affinity_process_status") not in process_statuses
        or sources != [THREAD_AFFINITY_EVIDENCE_SOURCE]
        or not isinstance(masks, list)
        or masks != [assigned_mask]
        or not isinstance(identities, list)
        or not identities
        or identities != list(dict.fromkeys(identities))
        or not isinstance(pid, int)
        or isinstance(pid, bool)
        or not isinstance(process_starttime, int)
        or isinstance(process_starttime, bool)
        or f"{pid}:{process_starttime}" not in identities
        or any(
            not isinstance(identity, str)
            or re.fullmatch(r"[0-9]+:[0-9]+", identity) is None
            for identity in identities
        )
        or item.get("current_thread_affinity_sample_exact") is not True
        or item.get("all_thread_affinity_samples_exact") is not True
        or item.get("thread_affinity_violation_ever") is not False
        or not _thread_affinity_sample_matches(item, expected_cpu)
    ):
        return f"invalid CP2K thread-affinity evidence in {path}"
    return None


def recorded_execution_issue(
    path: Path,
    expected_contract: Mapping[str, object],
    output: Path,
    scientific_job_stamp: Path,
) -> str | None:
    """Revalidate a separate execution record against its immutable artifacts."""
    if not path.is_file():
        return f"missing execution record {path}"
    try:
        record = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid execution record {path}: {exc}"
    if record.get("schema_version") == LEGACY_SCHEMA_VERSION:
        return _legacy_recorded_execution_issue(
            record, path, expected_contract, output, scientific_job_stamp
        )
    if record.get("schema_version") != SCHEMA_VERSION:
        return f"execution record schema mismatch in {path}"
    if record.get("contract") != expected_contract:
        return f"execution contract mismatch in {path}"
    hardened_binding_evidence = _contract_requires_hardened_binding_evidence(
        expected_contract
    )
    hardening_declarations = {
        "binding_environment_scrub_contract",
        "pool_close_policy",
        "rank_affinity_observation",
    }
    if hardening_declarations.intersection(expected_contract) and not (
        hardened_binding_evidence
    ):
        return f"invalid hardened binding contract in {path}"
    expected_contract_sha = canonical_sha256(dict(expected_contract))
    if record.get("contract_sha256") != expected_contract_sha:
        return f"execution contract hash mismatch in {path}"
    pe_lists = expected_contract.get("ordered_pe_lists")
    assigned = str(record.get("assigned_ordered_pe_list", ""))
    if not isinstance(pe_lists, list) or assigned not in pe_lists:
        return f"execution record uses an unapproved ordered PE list in {path}"
    launcher = str(expected_contract.get("mpi_launcher", ""))
    launcher_path = Path(launcher)
    launcher_args = expected_contract.get("mpi_launcher_args")
    ranks = expected_contract.get("mpi_ranks_per_job")
    if not isinstance(launcher_args, list) or not isinstance(ranks, int):
        return f"invalid expected execution contract for {path}"
    command = record.get("command")
    input_path = Path(str(record.get("input", "")))
    cp2k_path = Path(str(record.get("cp2k", "")))
    expected_prefix = [
        launcher,
        *launcher_args,
        "--map-by",
        f"pe-list={assigned}:ordered",
        "--bind-to",
        "core",
        "--report-bindings",
        "-np",
        str(ranks),
    ]
    expected_command = [
        *expected_prefix,
        str(cp2k_path),
        "-i",
        str(input_path),
        "-o",
        str(output.resolve()),
    ]
    if not isinstance(command, list) or command != expected_command:
        return f"execution command/affinity mismatch in {path}"
    expected_launcher_sha = expected_contract.get("mpi_launcher_sha256")
    if (
        not isinstance(expected_launcher_sha, str)
        or record.get("mpi_launcher_sha256_at_launch") != expected_launcher_sha
        or not launcher_path.is_file()
        or sha256(launcher_path) != expected_launcher_sha
    ):
        return f"execution MPI-launcher hash mismatch in {path}"
    if (
        not cp2k_path.is_file()
        or record.get("cp2k_sha256_at_launch") != sha256(cp2k_path)
    ):
        return f"CP2K executable hash mismatch in {path}"
    if record.get("runtime_affinity_gate") is not True:
        return f"runtime MPI/affinity gate failed in {path}"
    if record.get("return_code") != 0:
        return f"execution return code is not zero in {path}"
    if record.get("cross_process_cpu_reservation_gate") is not True:
        return f"cross-process CPU reservation gate failed in {path}"
    if hardened_binding_evidence:
        scrub_contract = expected_contract.get(
            "binding_environment_scrub_contract"
        )
        if not isinstance(scrub_contract, Mapping):
            return f"invalid binding-environment contract in {path}"
        sanctioned_environment = scrub_contract.get(
            "sanctioned_thread_environment"
        )
        if (
            not isinstance(sanctioned_environment, dict)
            or not all(
                isinstance(key, str) and isinstance(value, str)
                for key, value in sanctioned_environment.items()
            )
            or record.get("thread_environment") != sanctioned_environment
        ):
            return f"thread environment does not reconstruct from contract in {path}"
        removed_environment = record.get(
            "removed_mpi_binding_environment_keys"
        )
        if (
            not isinstance(removed_environment, list)
            or any(not isinstance(key, str) for key in removed_environment)
            or removed_environment != sorted(set(removed_environment))
            or mpi_control_environment_keys(
                {key: "removed" for key in removed_environment}
            )
            != removed_environment
        ):
            return f"invalid removed binding-environment key evidence in {path}"
        if record.get("binding_environment_scrub_gate") is not True:
            return f"binding/preload environment scrub gate failed in {path}"
        if record.get("residual_binding_environment_keys") != []:
            return f"binding/preload environment survived scrubbing in {path}"
        if record.get("ld_library_path_preserved") is not True:
            return f"LD_LIBRARY_PATH preservation gate failed in {path}"
    if record.get("live_compute_overlap_preflight_gate") is not True:
        return f"live CP2K/MPI CPU-overlap preflight failed in {path}"
    if record.get("mpi_bind_to") != "core":
        return f"MPI core binding was not recorded in {path}"
    if record.get("timing_classification") != "production_scaling_eligible":
        return f"execution timing is not scaling-eligible in {path}"
    rank_pids = record.get("observed_cp2k_rank_pids")
    rank_ids = record.get("observed_cp2k_rank_ids")
    rank_masks = record.get("observed_cp2k_rank_masks")
    try:
        assigned_cpus = parse_ordered_pe_list(assigned)
    except ValueError:
        return f"invalid recorded CP2K rank CPU mask in {path}"
    if record.get("assigned_cpu_count") != len(assigned_cpus):
        return f"assigned CPU count mismatch in {path}"

    child_processes = record.get("observed_child_processes")
    if not isinstance(child_processes, list) or not all(
        isinstance(item, dict) for item in child_processes
    ):
        return f"invalid observed child-process evidence in {path}"
    observed_by_pid: dict[int, dict[str, object]] = {}
    for item in child_processes:
        pid = item.get("pid")
        if not isinstance(pid, int) or isinstance(pid, bool) or pid in observed_by_pid:
            return f"invalid or duplicate observed child PID in {path}"
        observed_by_pid[pid] = item

    duplicate_rank_samples = record.get("concurrent_duplicate_rank_samples")
    if not isinstance(duplicate_rank_samples, list):
        return f"invalid concurrent duplicate-rank sample evidence in {path}"
    duplicate_rank_id_set: set[int] = set()
    previous_sample_index = 0
    for sample in duplicate_rank_samples:
        if not isinstance(sample, dict):
            return f"invalid concurrent duplicate-rank sample evidence in {path}"
        sample_index = sample.get("sample_index")
        groups = sample.get("rank_pid_groups")
        if (
            not isinstance(sample_index, int)
            or isinstance(sample_index, bool)
            or sample_index <= previous_sample_index
            or not isinstance(groups, list)
            or not groups
        ):
            return f"invalid concurrent duplicate-rank sample evidence in {path}"
        previous_sample_index = sample_index
        sample_ranks: list[int] = []
        for group in groups:
            if not isinstance(group, dict):
                return f"invalid concurrent duplicate-rank sample evidence in {path}"
            rank = group.get("ompi_comm_world_rank")
            pids = group.get("pids")
            if (
                not isinstance(rank, int)
                or isinstance(rank, bool)
                or not 0 <= rank < ranks
                or not isinstance(pids, list)
                or len(pids) < 2
                or any(
                    not isinstance(pid, int) or isinstance(pid, bool) for pid in pids
                )
                or pids != sorted(set(pids))
                or any(
                    pid not in observed_by_pid
                    or observed_by_pid[pid].get("is_cp2k_rank") is not True
                    or observed_by_pid[pid].get("ompi_comm_world_rank") != rank
                    for pid in pids
                )
            ):
                return f"invalid concurrent duplicate-rank sample evidence in {path}"
            sample_ranks.append(rank)
            duplicate_rank_id_set.add(rank)
        if sample_ranks != sorted(set(sample_ranks)):
            return f"invalid concurrent duplicate-rank sample evidence in {path}"

    duplicate_rank_ids = record.get("concurrent_duplicate_rank_ids_ever")
    if (
        not isinstance(duplicate_rank_ids, list)
        or any(
            not isinstance(rank, int)
            or isinstance(rank, bool)
            or not 0 <= rank < ranks
            for rank in duplicate_rank_ids
        )
        or duplicate_rank_ids != sorted(set(duplicate_rank_ids))
        or duplicate_rank_ids != sorted(duplicate_rank_id_set)
    ):
        return f"invalid concurrent duplicate-rank evidence in {path}"
    if record.get("concurrent_duplicate_rank_processes_ever") is not bool(
        duplicate_rank_id_set
    ):
        return f"inconsistent concurrent duplicate-rank evidence in {path}"

    rank_children = [item for item in child_processes if item.get("is_cp2k_rank")]
    unranked_cp2k_process_seen = any(
        not isinstance(item.get("ompi_comm_world_rank"), int)
        or isinstance(item.get("ompi_comm_world_rank"), bool)
        for item in rank_children
    )
    if record.get("unranked_cp2k_process_seen") is not unranked_cp2k_process_seen:
        return f"inconsistent unranked CP2K process evidence in {path}"

    for item in rank_children:
        rank = item.get("ompi_comm_world_rank")
        if not isinstance(rank, int) or isinstance(rank, bool):
            continue
        if not 0 <= rank < len(assigned_cpus):
            return f"CP2K rank identity is outside the assigned PE list in {path}"
        assigned_mask = str(assigned_cpus[rank])
        provenance_issue = _rank_process_provenance_issue(
            item, rank, path, cp2k_path.resolve()
        )
        if provenance_issue is not None:
            return provenance_issue
        environment_issue = _rank_environment_evidence_issue(
            item, assigned_mask, path
        )
        if environment_issue is not None:
            return environment_issue
        if hardened_binding_evidence:
            thread_issue = _thread_affinity_evidence_issue(
                item, assigned_mask, path
            )
            if thread_issue is not None:
                return thread_issue
        observed_rank_ids = item.get("observed_rank_ids")
        if (
            not isinstance(observed_rank_ids, list)
            or len(observed_rank_ids) != 1
            or not isinstance(observed_rank_ids[0], int)
            or isinstance(observed_rank_ids[0], bool)
            or observed_rank_ids[0] != rank
            or item.get("observed_cpu_masks") != [assigned_mask]
            or item.get("cpus_allowed_list") != assigned_mask
            or not isinstance(item.get("sample_count"), int)
            or isinstance(item.get("sample_count"), bool)
            or int(item["sample_count"]) < 1
            or item.get("affinity_violation_ever") is not False
            or item.get("rank_identity_changed_ever") is not False
            or item.get("current_sample_matches_assigned_singleton") is not True
        ):
            return f"CP2K rank affinity history is incomplete or mismatched in {path}"

    recomputed_rank_evidence = _aggregate_cp2k_rank_generations(
        observed_by_pid, assigned_cpus, duplicate_rank_id_set
    )
    recomputed_rank_ids = [
        int(item["ompi_comm_world_rank"]) for item in recomputed_rank_evidence
    ]
    recomputed_rank_pids = [
        int(item["canonical_pid"]) for item in recomputed_rank_evidence
    ]
    recomputed_pid_generations = [
        item["pid_generations"] for item in recomputed_rank_evidence
    ]
    if any(len(generations) != 1 for generations in recomputed_pid_generations):
        return f"multiple CP2K PID generations are not scaling-eligible in {path}"
    recomputed_rank_masks: list[str] = []
    normalized_rank_masks: list[set[int]] = []
    for item in recomputed_rank_evidence:
        masks = item.get("observed_cpu_masks")
        if not isinstance(masks, list) or len(masks) != 1:
            recomputed_rank_masks = []
            normalized_rank_masks = []
            break
        mask_text = str(masks[0])
        try:
            normalized_rank_masks.append(parse_cpu_set(mask_text))
        except ValueError:
            recomputed_rank_masks = []
            normalized_rank_masks = []
            break
        recomputed_rank_masks.append(mask_text)

    rank_count_matches = len(recomputed_rank_evidence) == ranks
    rank_ids_exact = recomputed_rank_ids == list(range(ranks))
    masks_complete = len(normalized_rank_masks) == len(recomputed_rank_evidence)
    masks_exact = masks_complete and all(
        0 <= rank < len(assigned_cpus)
        and mask == {assigned_cpus[rank]}
        for rank, mask in zip(
            recomputed_rank_ids, normalized_rank_masks, strict=True
        )
    )
    all_rank_samples_exact = (
        rank_count_matches
        and not unranked_cp2k_process_seen
        and not duplicate_rank_id_set
        and all(
            item.get("all_samples_match_assigned_singleton") is True
            for item in recomputed_rank_evidence
        )
    )

    if record.get("observed_cp2k_rank_evidence") != recomputed_rank_evidence:
        return f"aggregated CP2K rank-generation evidence mismatch in {path}"
    if record.get("observed_cp2k_rank_pid_generations") != recomputed_pid_generations:
        return f"CP2K rank PID-generation history mismatch in {path}"
    if record.get("observed_cp2k_process_generation_count") != sum(
        len(generations) for generations in recomputed_pid_generations
    ):
        return f"CP2K process-generation count mismatch in {path}"
    if record.get("observed_cp2k_rank_count") != len(recomputed_rank_evidence):
        return f"observed CP2K logical-rank count mismatch in {path}"
    if record.get("expected_cp2k_rank_count") != ranks:
        return f"expected CP2K rank count mismatch in {path}"
    if record.get("rank_count_matches") is not rank_count_matches:
        return f"derived CP2K rank-count gate mismatch in {path}"
    if record.get("rank_ids_exactly_0_to_n_minus_1") is not rank_ids_exact:
        return f"derived CP2K rank-ID gate mismatch in {path}"
    if record.get("rank_masks_complete") is not masks_complete:
        return f"derived CP2K rank-mask completeness gate mismatch in {path}"
    if record.get("rank_masks_exactly_match_ordered_pe_list") is not masks_exact:
        return f"derived CP2K rank-mask gate mismatch in {path}"
    if (
        record.get("all_observed_rank_samples_match_ordered_pe_list")
        is not all_rank_samples_exact
    ):
        return f"derived CP2K rank-sample gate mismatch in {path}"
    if rank_pids != recomputed_rank_pids or len(recomputed_rank_pids) != ranks:
        return f"observed CP2K rank count mismatch in {path}"
    if rank_ids != recomputed_rank_ids or not rank_ids_exact:
        return f"observed CP2K MPI-rank ordering mismatch in {path}"
    if rank_masks != recomputed_rank_masks or not masks_exact:
        return f"observed CP2K rank CPU-mask mismatch in {path}"
    if duplicate_rank_id_set:
        return f"concurrently live duplicate CP2K MPI rank detected in {path}"
    if unranked_cp2k_process_seen:
        return f"unranked CP2K process generation detected in {path}"
    if not all_rank_samples_exact:
        return f"CP2K rank affinity history is incomplete or mismatched in {path}"
    log_path = Path(str(record.get("launcher_log", "")))
    if log_path.resolve() != launcher_log_path(output).resolve():
        return f"launcher-log path mismatch in {path}"
    if not log_path.is_file() or record.get("launcher_log_sha256") != sha256(log_path):
        return f"launcher-log hash mismatch in {path}"
    reported_binding_rank_ids = _reported_binding_rank_ids(
        log_path.read_text(errors="replace")
    )
    if (
        record.get("reported_binding_rank_ids") != reported_binding_rank_ids
        or reported_binding_rank_ids != list(range(ranks))
    ):
        return f"launcher binding report is incomplete in {path}"
    if record.get("binding_report_complete") is not (
        reported_binding_rank_ids == list(range(ranks))
    ):
        return f"derived launcher binding-report gate mismatch in {path}"
    if record.get("live_compute_overlap_preflight_owners") != []:
        return f"live CP2K/MPI overlap owners were recorded in {path}"
    if record.get("live_compute_overlap_runtime_gate") is not True:
        return f"runtime live CP2K/MPI overlap gate failed in {path}"
    if record.get("live_compute_overlap_runtime_samples") != []:
        return f"runtime live CP2K/MPI overlap was recorded in {path}"
    if hardened_binding_evidence:
        if record.get("local_affinity_violation_gate") is not True:
            return f"runtime local rank/thread affinity gate failed in {path}"
        if record.get("local_affinity_violation_samples") != []:
            return (
                "runtime local rank/thread affinity violation was recorded in "
                f"{path}"
            )
    return _common_artifact_issue(record, path, output, scientific_job_stamp)


def _linux_proc_stat_identity(stat_text: str) -> tuple[str, int]:
    """Return Linux task state and starttime without misparsing ``comm``."""
    closing_parenthesis = stat_text.rfind(")")
    if closing_parenthesis < 0:
        raise ValueError("malformed Linux /proc PID stat record")
    fields = stat_text[closing_parenthesis + 1 :].split()
    if len(fields) <= 19:
        raise ValueError("truncated Linux /proc PID stat record")
    return fields[0], int(fields[19])


def _linux_process_starttime(
    pid: int, proc_root: Path = Path("/proc")
) -> int | None:
    """Read one Linux task birth identity, or ``None`` if it is not provable."""
    try:
        _, starttime = _linux_proc_stat_identity(
            (proc_root / str(pid) / "stat").read_text(errors="replace")
        )
    except (
        FileNotFoundError,
        PermissionError,
        ProcessLookupError,
        OSError,
        ValueError,
    ):
        return None
    return starttime


def _status_fields(status_text: str) -> dict[str, str]:
    return {
        key: value.strip()
        for line in status_text.splitlines()
        if ":" in line
        for key, value in (line.split(":", 1),)
    }


def _linux_thread_affinity_snapshot(
    pid: int,
    process_starttime: int,
    proc_root: Path = Path("/proc"),
) -> tuple[list[dict[str, object]], str]:
    """Read every Linux task affinity with TID/starttime race detection."""
    task_root = proc_root / str(pid) / "task"
    try:
        initial_tids = sorted(
            int(entry.name) for entry in task_root.iterdir() if entry.name.isdigit()
        )
    except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
        return [], "task_directory_unreadable"
    if not initial_tids:
        return [], "task_directory_empty"

    records: list[dict[str, object]] = []
    issues: list[str] = []
    for tid in initial_tids:
        task = task_root / str(tid)
        try:
            initial_stat_state, initial_starttime = _linux_proc_stat_identity(
                (task / "stat").read_text(errors="replace")
            )
            initial_fields = _status_fields(
                (task / "status").read_text(errors="replace")
            )
        except FileNotFoundError:
            issues.append("thread_disappeared_before_sample")
            continue
        except (PermissionError, ProcessLookupError, OSError, ValueError):
            issues.append("thread_identity_unreadable_before_sample")
            continue

        initial_mask = initial_fields.get("Cpus_allowed_list", "")
        masks = [initial_mask] if initial_mask else []
        try:
            final_fields = _status_fields(
                (task / "status").read_text(errors="replace")
            )
            final_stat_state, final_starttime = _linux_proc_stat_identity(
                (task / "stat").read_text(errors="replace")
            )
        except FileNotFoundError:
            issues.append("thread_disappeared_during_sample")
            records.append(
                {
                    "tid": tid,
                    "thread_starttime": initial_starttime,
                    "observed_thread_starttimes": [initial_starttime],
                    "state": initial_fields.get("State", ""),
                    "stat_state": initial_stat_state,
                    "cpus_allowed_list": initial_mask,
                    "observed_cpu_masks": masks,
                    "live": initial_stat_state not in {"Z", "X"},
                    "identity_status": "disappeared_during_sample",
                }
            )
            continue
        except (PermissionError, ProcessLookupError, OSError, ValueError):
            issues.append("thread_identity_unreadable_after_sample")
            records.append(
                {
                    "tid": tid,
                    "thread_starttime": initial_starttime,
                    "observed_thread_starttimes": [initial_starttime],
                    "state": initial_fields.get("State", ""),
                    "stat_state": initial_stat_state,
                    "cpus_allowed_list": initial_mask,
                    "observed_cpu_masks": masks,
                    "live": initial_stat_state not in {"Z", "X"},
                    "identity_status": "identity_unreadable_after_sample",
                }
            )
            continue

        final_mask = final_fields.get("Cpus_allowed_list", "")
        if final_mask and final_mask not in masks:
            masks.append(final_mask)
        identity_status = "stable"
        if final_starttime != initial_starttime:
            identity_status = "tid_reused_during_sample"
            issues.append(identity_status)
        elif final_mask != initial_mask:
            identity_status = "cpu_mask_changed_during_sample"
            issues.append(identity_status)
        if tid == pid and initial_starttime != process_starttime:
            identity_status = "leader_starttime_mismatch"
            issues.append(identity_status)
        records.append(
            {
                "tid": tid,
                "thread_starttime": initial_starttime,
                "observed_thread_starttimes": sorted(
                    {initial_starttime, final_starttime}
                ),
                "state": final_fields.get("State", ""),
                "stat_state": final_stat_state,
                "cpus_allowed_list": final_mask,
                "observed_cpu_masks": masks,
                "live": final_stat_state not in {"Z", "X"},
                "identity_status": identity_status,
            }
        )

    try:
        final_tids = sorted(
            int(entry.name) for entry in task_root.iterdir() if entry.name.isdigit()
        )
    except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
        final_tids = []
        issues.append("task_directory_unreadable_after_sample")
    if final_tids != initial_tids:
        issues.append("thread_set_changed_during_sample")
    if pid not in initial_tids:
        issues.append("leader_thread_missing")
    status = "+".join(sorted(set(issues))) if issues else "consistent"
    return records, status


def _linux_process_snapshot(
    pid: int, cp2k: Path, proc_root: Path = Path("/proc")
) -> dict[str, object] | None:
    """Read one live process/rank affinity from procfs without external tools."""
    root = proc_root / str(pid)
    if not root.is_dir():
        return None
    try:
        stat_state, process_starttime = _linux_proc_stat_identity(
            (root / "stat").read_text(errors="replace")
        )
        initial_status = (root / "status").read_text(errors="replace")
        command = (root / "cmdline").read_bytes().split(b"\0")
        arguments = [item.decode(errors="replace") for item in command if item]
        initial_executable = str((root / "exe").resolve(strict=True))
    except (
        FileNotFoundError,
        PermissionError,
        ProcessLookupError,
        OSError,
        ValueError,
    ):
        return None
    try:
        environment_bytes = (root / "environ").read_bytes()
    except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
        environment_bytes = b""
        environment_read_status = "unreadable"
    else:
        environment_read_status = "available" if environment_bytes else "empty"
    thread_affinity, thread_affinity_scan_status = _linux_thread_affinity_snapshot(
        pid, process_starttime, proc_root
    )
    thread_affinity_scan_issues = (
        []
        if thread_affinity_scan_status == "consistent"
        else thread_affinity_scan_status.split("+")
    )
    thread_affinity_process_status = "live"
    try:
        final_status = (root / "status").read_text(errors="replace")
        final_executable = str((root / "exe").resolve(strict=True))
        final_stat_state, final_process_starttime = _linux_proc_stat_identity(
            (root / "stat").read_text(errors="replace")
        )
    except FileNotFoundError:
        process_identity_status = "disappeared_after_sample"
        snapshot_consistency_status = "process_disappeared"
        status = initial_status
        executable = initial_executable
        thread_affinity_process_status = "process_disappeared_after_sample"
    except (PermissionError, ProcessLookupError, OSError, ValueError):
        process_identity_status = "identity_unreadable_after_sample"
        snapshot_consistency_status = "final_identity_unreadable"
        status = initial_status
        executable = initial_executable
    else:
        initial_fields = {
            key: value.strip()
            for line in initial_status.splitlines()
            if ":" in line
            for key, value in (line.split(":", 1),)
        }
        final_fields = {
            key: value.strip()
            for line in final_status.splitlines()
            if ":" in line
            for key, value in (line.split(":", 1),)
        }
        status = final_status
        executable = final_executable
        if final_process_starttime != process_starttime:
            process_identity_status = "pid_reused_during_sample"
            snapshot_consistency_status = "pid_reused"
        elif final_executable != initial_executable:
            process_identity_status = "executable_changed_during_sample"
            snapshot_consistency_status = "executable_changed"
        elif final_fields.get("Cpus_allowed_list") != initial_fields.get(
            "Cpus_allowed_list"
        ):
            process_identity_status = "cpu_mask_changed_during_sample"
            snapshot_consistency_status = "cpu_mask_changed"
        else:
            stat_state = final_stat_state
            process_identity_status = (
                "terminal_state" if stat_state in {"Z", "X"} else "stable"
            )
            snapshot_consistency_status = "consistent"
            if stat_state in {"Z", "X"}:
                thread_affinity_process_status = "terminal_process"
    environment_items = environment_bytes.split(b"\0")
    environment = {
        key.decode(errors="replace"): value.decode(errors="replace")
        for item in environment_items
        if item and b"=" in item
        for key, value in (item.split(b"=", 1),)
    }
    fields: dict[str, str] = {}
    for line in status.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key] = value.strip()
    cp2k_resolved = str(cp2k.resolve())
    argument_targets: list[str] = []
    for argument in arguments[:2]:
        candidate = Path(argument)
        if not candidate.is_absolute():
            continue
        try:
            argument_targets.append(str(candidate.resolve(strict=True)))
        except (FileNotFoundError, OSError):
            continue
    executable_name = Path(executable).name
    interpreted = executable_name in {"sh", "dash", "bash"} or executable_name.startswith(
        "python"
    )
    is_cp2k_rank = (
        executable == cp2k_resolved
        or (argument_targets and argument_targets[0] == cp2k_resolved)
        or (
            interpreted
            and len(argument_targets) > 1
            and argument_targets[1] == cp2k_resolved
        )
    )
    mpi_rank: int | None = None
    rank_value = environment.get("OMPI_COMM_WORLD_RANK")
    if environment_read_status != "available":
        rank_observation_status = f"environment_{environment_read_status}"
    elif rank_value is None:
        rank_observation_status = "explicit_missing"
    else:
        try:
            mpi_rank = int(rank_value)
        except ValueError:
            rank_observation_status = "explicit_invalid"
        else:
            rank_observation_status = "explicit"
    return {
        "pid": pid,
        "ppid": int(fields.get("PPid", "0")),
        "state": fields.get("State", ""),
        "stat_state": stat_state,
        "process_starttime": process_starttime,
        "process_identity_status": process_identity_status,
        "snapshot_consistency_status": snapshot_consistency_status,
        "executable": executable,
        "arguments": arguments,
        "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
        "ompi_comm_world_rank": mpi_rank,
        "rank_observation_status": rank_observation_status,
        "is_cp2k_rank": is_cp2k_rank,
        "thread_affinity_evidence_source": THREAD_AFFINITY_EVIDENCE_SOURCE,
        "thread_affinity_scan_status": thread_affinity_scan_status,
        "thread_affinity_scan_issues": thread_affinity_scan_issues,
        "thread_affinity_process_status": thread_affinity_process_status,
        "live_thread_affinity": thread_affinity,
    }


def _linux_descendants(root_pid: int) -> set[int]:
    """Return a best-effort live descendant closure rooted at *root_pid*."""
    found = {root_pid}
    pending = [root_pid]
    while pending:
        parent = pending.pop()
        children_path = Path("/proc") / str(parent) / "task" / str(parent) / "children"
        try:
            children = [int(value) for value in children_path.read_text().split()]
        except (
            FileNotFoundError,
            PermissionError,
            ProcessLookupError,
            OSError,
            ValueError,
        ):
            children = []
        for child in children:
            if child not in found:
                found.add(child)
                pending.append(child)
    return found


def _live_process_group_members(process_group: int) -> set[int]:
    """Return non-zombie processes that can still use a reserved CPU."""
    proc_root = Path("/proc")
    if not proc_root.is_dir():
        return set()
    members: set[int] = set()
    for entry in proc_root.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            if os.getpgid(pid) != process_group:
                continue
            state = next(
                line.split(":", 1)[1].strip().split()[0]
                for line in (entry / "status").read_text().splitlines()
                if line.startswith("State:")
            )
        except (
            FileNotFoundError,
            PermissionError,
            ProcessLookupError,
            OSError,
            StopIteration,
        ):
            continue
        if state not in {"Z", "X"}:
            members.add(pid)
    return members


def _terminate_and_reap_process_group(
    process: subprocess.Popen[bytes],
    term_timeout: float = 30.0,
    tracked_rank_starttimes: Mapping[int, int] | None = None,
) -> None:
    """Drain the launcher session and every directly tracked rank task."""
    process_group = process.pid
    tracked = dict(tracked_rank_starttimes or {})

    def live_tracked_ranks() -> set[int]:
        live: set[int] = set()
        for pid, starttime in tracked.items():
            resolution = _linux_process_terminal_resolution(pid, starttime)
            if resolution is None or resolution == "identity_unreadable":
                live.add(pid)
        return live

    def signal_tracked_ranks(sig: int) -> None:
        for pid in live_tracked_ranks():
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass

    try:
        os.killpg(process_group, signal.SIGTERM)
    except ProcessLookupError:
        pass
    signal_tracked_ranks(signal.SIGTERM)
    deadline = time.monotonic() + term_timeout
    while time.monotonic() < deadline:
        process.poll()
        if process.returncode is not None and not _live_process_group_members(
            process_group
        ) and not live_tracked_ranks():
            process.wait()
            return
        time.sleep(0.05)
    try:
        os.killpg(process_group, signal.SIGKILL)
    except ProcessLookupError:
        pass
    signal_tracked_ranks(signal.SIGKILL)
    while True:
        process.poll()
        members = _live_process_group_members(process_group)
        live_ranks = live_tracked_ranks()
        if process.returncode is not None and not members and not live_ranks:
            process.wait()
            return
        if members:
            try:
                os.killpg(process_group, signal.SIGKILL)
            except ProcessLookupError:
                pass
        if live_ranks:
            signal_tracked_ranks(signal.SIGKILL)
        time.sleep(0.05)


def _tracked_rank_process_starttimes(
    observed: Mapping[int, Mapping[str, object]],
) -> dict[int, int]:
    """Return immutable Linux identities for every observed CP2K rank task."""
    identities: dict[int, int] = {}
    for pid, record in observed.items():
        starttime = record.get("process_starttime")
        if (
            record.get("is_cp2k_rank") is True
            and isinstance(pid, int)
            and not isinstance(pid, bool)
            and isinstance(starttime, int)
            and not isinstance(starttime, bool)
        ):
            identities[pid] = starttime
    return identities


def _reported_binding_rank_ids(text: str) -> list[int]:
    """Extract Open MPI ``--report-bindings`` rank identifiers."""
    matches = re.findall(
        r"\b(?:MCW\s+)?rank\s+(\d+)\s+bound\b", text, flags=re.IGNORECASE
    )
    return sorted({int(value) for value in matches})


def _ordered_cp2k_rank_processes(
    observed: Mapping[int, Mapping[str, object]], expected_ranks: int
) -> list[Mapping[str, object]]:
    """Order CP2K processes by MPI rank, never by process identifier."""
    return sorted(
        (record for record in observed.values() if record.get("is_cp2k_rank")),
        key=lambda record: (
            int(record["ompi_comm_world_rank"])
            if isinstance(record.get("ompi_comm_world_rank"), int)
            else expected_ranks,
            int(record["pid"]),
        ),
    )


def _concurrent_live_duplicate_rank_ids(
    snapshots: Sequence[Mapping[str, object]],
) -> set[int]:
    """Return MPI ranks represented by multiple live CP2K PIDs in one sample."""
    return set(_concurrent_live_rank_pid_groups(snapshots))


def _concurrent_live_rank_pid_groups(
    snapshots: Sequence[Mapping[str, object]],
) -> dict[int, list[int]]:
    """Return concurrently live CP2K PIDs grouped by duplicated MPI rank."""
    rank_pids: dict[int, set[int]] = {}
    for snapshot in snapshots:
        if snapshot.get("is_cp2k_rank") is not True:
            continue
        if str(snapshot.get("state", "")).startswith(("Z", "X")):
            continue
        rank = snapshot.get("ompi_comm_world_rank")
        pid = snapshot.get("pid")
        if isinstance(rank, int) and isinstance(pid, int):
            rank_pids.setdefault(rank, set()).add(pid)
    return {
        rank: sorted(pids)
        for rank, pids in sorted(rank_pids.items())
        if len(pids) > 1
    }


def _thread_affinity_sample_matches(
    snapshot: Mapping[str, object], expected_cpu: int
) -> bool:
    """Require one stable, singleton affinity proof for every captured task."""
    scan_issues = snapshot.get("thread_affinity_scan_issues")
    if scan_issues is None:
        scan_issues = (
            []
            if snapshot.get("thread_affinity_scan_status") == "consistent"
            else [snapshot.get("thread_affinity_scan_status")]
        )
    if (
        snapshot.get("thread_affinity_evidence_source")
        != THREAD_AFFINITY_EVIDENCE_SOURCE
        or snapshot.get("thread_affinity_scan_status") != "consistent"
        or scan_issues != []
        or snapshot.get("thread_affinity_process_status", "live")
        not in THREAD_AFFINITY_ACCEPTED_PROCESS_STATUSES
    ):
        return False
    records = snapshot.get("live_thread_affinity")
    if not isinstance(records, list) or not all(
        isinstance(record, dict) for record in records
    ):
        return False
    if not records:
        return False
    pid = snapshot.get("pid")
    process_starttime = snapshot.get("process_starttime")
    leader_proven = False
    for record in records:
        tid = record.get("tid")
        starttime = record.get("thread_starttime")
        observed_starttimes = record.get("observed_thread_starttimes")
        masks = record.get("observed_cpu_masks")
        if (
            not isinstance(tid, int)
            or isinstance(tid, bool)
            or not isinstance(starttime, int)
            or isinstance(starttime, bool)
            or not isinstance(observed_starttimes, list)
            or not observed_starttimes
            or any(
                not isinstance(value, int) or isinstance(value, bool)
                for value in observed_starttimes
            )
            or starttime not in observed_starttimes
            or not isinstance(masks, list)
            or not masks
            or not isinstance(record.get("live"), bool)
            or not isinstance(record.get("identity_status"), str)
        ):
            return False
        try:
            normalized_masks = [parse_cpu_set(str(mask)) for mask in masks]
            current_mask = parse_cpu_set(str(record.get("cpus_allowed_list", "")))
        except ValueError:
            return False
        if any(mask != {expected_cpu} for mask in normalized_masks):
            return False
        if current_mask != {expected_cpu}:
            return False
        if record.get("identity_status") != "stable":
            return False
        if (
            tid == pid
            and starttime == process_starttime
            and record.get("identity_status") == "stable"
        ):
            leader_proven = True
    return leader_proven


def _aggregate_cp2k_rank_generations(
    observed: Mapping[int, Mapping[str, object]],
    assigned_cpus: Sequence[int],
    concurrent_duplicate_rank_ids: set[int] | None = None,
) -> list[dict[str, object]]:
    """Aggregate rank tasks while making successor generations fail closed."""
    groups: dict[int, list[Mapping[str, object]]] = {}
    for record in observed.values():
        if record.get("is_cp2k_rank") is not True:
            continue
        rank = record.get("ompi_comm_world_rank")
        if isinstance(rank, int):
            groups.setdefault(rank, []).append(record)
    duplicates = concurrent_duplicate_rank_ids or set()
    aggregates: list[dict[str, object]] = []
    for rank in sorted(groups):
        generations = sorted(groups[rank], key=lambda record: int(record["pid"]))
        mask_history: set[str] = set()
        for record in generations:
            history = record.get("observed_cpu_masks")
            if isinstance(history, list):
                mask_history.update(str(mask) for mask in history)
            elif record.get("cpus_allowed_list"):
                mask_history.add(str(record["cpus_allowed_list"]))
        canonical = max(
            generations,
            key=lambda record: (int(record.get("sample_count", 0)), -int(record["pid"])),
        )
        all_samples_exact = (
            0 <= rank < len(assigned_cpus)
            and len(generations) == 1
            and mask_history == {str(assigned_cpus[rank])}
            and rank not in duplicates
            and all(
                record.get("affinity_violation_ever") is False
                and record.get("rank_identity_changed_ever") is not True
                and record.get("process_starttime_changed_ever") is not True
                and record.get("process_terminally_confirmed") is True
                and record.get("process_reappeared_after_terminal_ever") is not True
                and record.get("executable_changed_ever") is not True
                and record.get("cpu_mask_changed_during_sample_ever") is not True
                and (
                    "thread_affinity_evidence_source" not in record
                    or record.get("thread_affinity_violation_ever") is False
                    and record.get("all_thread_affinity_samples_exact") is True
                )
                and record.get("snapshot_unavailable_ever") is not True
                and record.get("rank_environment_unavailable_pending") is not True
                and (
                    record.get("rank_environment_unavailable_ever") is not True
                    or record.get("rank_environment_terminally_confirmed") is True
                )
                and record.get("current_sample_matches_assigned_singleton") is True
                for record in generations
            )
        )
        aggregates.append(
            {
                "ompi_comm_world_rank": rank,
                "canonical_pid": int(canonical["pid"]),
                "pid_generations": [int(record["pid"]) for record in generations],
                "observed_cpu_masks": sorted(mask_history),
                "all_samples_match_assigned_singleton": all_samples_exact,
                "concurrent_duplicate_pid_ever": rank in duplicates,
            }
        )
    return aggregates


def _accumulate_process_snapshot(
    previous: Mapping[str, object] | None,
    snapshot: Mapping[str, object],
    assigned_cpus: Sequence[int],
) -> dict[str, object]:
    """Accumulate affinity evidence without forgetting an earlier violation."""
    accumulated = dict(snapshot)
    sample_count = (
        int(previous.get("sample_count", 0)) + 1 if previous is not None else 1
    )
    accumulated["sample_count"] = sample_count
    if not snapshot.get("is_cp2k_rank"):
        if previous and previous.get("is_cp2k_rank") is True:
            retained = dict(previous)
            retained.update(
                {
                    "sample_count": sample_count,
                    "last_observed_executable": snapshot.get("executable"),
                    "last_observed_arguments": snapshot.get("arguments"),
                    "last_observed_process_identity_status": snapshot.get(
                        "process_identity_status"
                    ),
                    "executable_changed_ever": True,
                    "current_sample_matches_assigned_singleton": False,
                    "affinity_violation_ever": True,
                }
            )
            prior_starttime = previous.get("process_starttime")
            current_starttime = snapshot.get("process_starttime")
            if (
                isinstance(prior_starttime, int)
                and not isinstance(prior_starttime, bool)
                and isinstance(current_starttime, int)
                and not isinstance(current_starttime, bool)
                and prior_starttime != current_starttime
            ) or snapshot.get("process_identity_status") == "pid_reused_during_sample":
                retained["process_starttime_changed_ever"] = True
            return retained
        accumulated["affinity_violation_ever"] = bool(
            previous and previous.get("affinity_violation_ever")
        )
        return accumulated

    raw_rank = snapshot.get("ompi_comm_world_rank")
    observation_status = snapshot.get("rank_observation_status")
    if not isinstance(observation_status, str):
        observation_status = (
            "explicit"
            if isinstance(raw_rank, int) and not isinstance(raw_rank, bool)
            else "explicit_missing"
        )
    mask_text = str(snapshot.get("cpus_allowed_list", ""))
    try:
        mask = parse_cpu_set(mask_text)
    except ValueError:
        mask = set()
    previous_rank = previous.get("ompi_comm_world_rank") if previous else None
    previous_pending = bool(
        previous and previous.get("rank_environment_unavailable_pending")
    )
    previous_starttime = previous.get("process_starttime") if previous else None
    current_starttime = snapshot.get("process_starttime")
    starttimes_comparable = bool(
        previous
        and isinstance(previous_starttime, int)
        and not isinstance(previous_starttime, bool)
        and isinstance(current_starttime, int)
        and not isinstance(current_starttime, bool)
    )
    same_process_identity = bool(
        starttimes_comparable
        and snapshot.get("pid") == previous.get("pid")
        and current_starttime == previous_starttime
        and previous.get("process_starttime_changed_ever") is not True
        and snapshot.get("process_identity_status", "stable")
        in {"stable", "terminal_state", "disappeared_after_sample"}
    )
    unavailable_environment = observation_status in {
        "environment_empty",
        "environment_unreadable",
    }
    previous_rank_proven = bool(
        previous
        and isinstance(previous_rank, int)
        and not isinstance(previous_rank, bool)
        and previous.get("observed_rank_ids") == [previous_rank]
        and previous.get("rank_identity_changed_ever") is not True
        and previous.get("affinity_violation_ever") is not True
    )
    retain_pending_rank = bool(
        unavailable_environment
        and previous_rank_proven
        and same_process_identity
        and 0 <= int(previous_rank) < len(assigned_cpus)
        and mask == {assigned_cpus[int(previous_rank)]}
        and _thread_affinity_sample_matches(
            snapshot, assigned_cpus[int(previous_rank)]
        )
    )

    rank = previous_rank if retain_pending_rank else raw_rank
    if retain_pending_rank:
        accumulated["ompi_comm_world_rank"] = rank
        accumulated["raw_ompi_comm_world_rank"] = raw_rank
        accumulated["rank_identity_source"] = "pending_terminal_environment_loss"
    else:
        accumulated["raw_ompi_comm_world_rank"] = raw_rank
        accumulated["rank_identity_source"] = (
            "explicit_environment" if observation_status == "explicit" else "unproven"
        )
    sample_matches = (
        isinstance(rank, int)
        and not isinstance(rank, bool)
        and 0 <= rank < len(assigned_cpus)
        and mask == {assigned_cpus[rank]}
    )
    thread_sample_matches = bool(
        isinstance(rank, int)
        and not isinstance(rank, bool)
        and 0 <= rank < len(assigned_cpus)
        and _thread_affinity_sample_matches(snapshot, assigned_cpus[rank])
    )
    rank_history = list(previous.get("observed_rank_ids", [])) if previous else []
    mask_history = list(previous.get("observed_cpu_masks", [])) if previous else []
    status_history = (
        list(previous.get("observed_rank_observation_statuses", []))
        if previous
        else []
    )
    state_history = list(previous.get("observed_process_states", [])) if previous else []
    starttime_history = (
        list(previous.get("observed_process_starttimes", [])) if previous else []
    )
    thread_status_history = (
        list(previous.get("thread_affinity_scan_statuses", []))
        if previous
        else []
    )
    thread_source_history = (
        list(previous.get("thread_affinity_evidence_sources", []))
        if previous
        else []
    )
    thread_issue_history = (
        list(previous.get("thread_affinity_scan_issues_ever", []))
        if previous
        else []
    )
    thread_process_status_history = (
        list(previous.get("thread_affinity_process_statuses", []))
        if previous
        else []
    )
    thread_mask_history = (
        list(previous.get("observed_thread_cpu_masks", []))
        if previous
        else []
    )
    thread_identity_history = (
        list(previous.get("observed_thread_identities", []))
        if previous
        else []
    )
    thread_status = str(snapshot.get("thread_affinity_scan_status", "missing"))
    thread_source = str(snapshot.get("thread_affinity_evidence_source", "missing"))
    raw_thread_issues = snapshot.get("thread_affinity_scan_issues")
    if isinstance(raw_thread_issues, list) and all(
        isinstance(issue, str) for issue in raw_thread_issues
    ):
        thread_issues = list(raw_thread_issues)
    elif raw_thread_issues is None and thread_status == "consistent":
        thread_issues = []
    else:
        thread_issues = ["invalid_or_missing_thread_scan_issue_evidence"]
    thread_process_status = str(
        snapshot.get("thread_affinity_process_status", "live")
    )
    if thread_status not in thread_status_history:
        thread_status_history.append(thread_status)
    if thread_source not in thread_source_history:
        thread_source_history.append(thread_source)
    for issue in thread_issues:
        if issue not in thread_issue_history:
            thread_issue_history.append(issue)
    if thread_process_status not in thread_process_status_history:
        thread_process_status_history.append(thread_process_status)
    accumulated["thread_affinity_scan_issues"] = thread_issues
    accumulated["thread_affinity_process_status"] = thread_process_status
    thread_records = snapshot.get("live_thread_affinity")
    if isinstance(thread_records, list):
        for thread_record in thread_records:
            if not isinstance(thread_record, Mapping):
                continue
            tid = thread_record.get("tid")
            thread_starttime = thread_record.get("thread_starttime")
            if (
                isinstance(tid, int)
                and not isinstance(tid, bool)
                and isinstance(thread_starttime, int)
                and not isinstance(thread_starttime, bool)
            ):
                identity = f"{tid}:{thread_starttime}"
                if identity not in thread_identity_history:
                    thread_identity_history.append(identity)
            masks = thread_record.get("observed_cpu_masks")
            if isinstance(masks, list):
                for thread_mask in masks:
                    thread_mask_text = str(thread_mask)
                    if thread_mask_text not in thread_mask_history:
                        thread_mask_history.append(thread_mask_text)
    if observation_status == "explicit" and raw_rank not in rank_history:
        rank_history.append(raw_rank)
    if mask_text not in mask_history:
        mask_history.append(mask_text)
    if observation_status not in status_history:
        status_history.append(observation_status)
    state_text = str(snapshot.get("state", ""))
    if state_text and state_text not in state_history:
        state_history.append(state_text)
    if current_starttime is not None and current_starttime not in starttime_history:
        starttime_history.append(current_starttime)
    identity_status = snapshot.get("process_identity_status", "stable")
    process_starttime_changed = bool(
        identity_status == "pid_reused_during_sample"
        or previous
        and starttimes_comparable
        and current_starttime != previous_starttime
    )
    process_snapshot_inconsistent = identity_status not in {
        "stable",
        "terminal_state",
        "disappeared_after_sample",
    }
    explicit_rank_reappeared_after_loss = bool(
        previous_pending and observation_status == "explicit"
    )
    rank_identity_changed = bool(
        previous
        and not retain_pending_rank
        and (
            observation_status != "explicit"
            or any(
                isinstance(prior_rank, int) and prior_rank != raw_rank
                for prior_rank in previous.get("observed_rank_ids", [])
            )
        )
    )
    environment_events = (
        [dict(event) for event in previous.get("rank_environment_events", [])]
        if previous
        else []
    )
    if unavailable_environment:
        environment_events.append(
            {
                "sample_index": sample_count,
                "pid": snapshot.get("pid"),
                "process_starttime": current_starttime,
                "state": state_text,
                "cpus_allowed_list": mask_text,
                "environment_status": observation_status,
                "terminal_resolution": (
                    "pending" if retain_pending_rank else "rejected_unproven"
                ),
            }
        )
    elif explicit_rank_reappeared_after_loss:
        for event in environment_events:
            if event.get("terminal_resolution") == "pending":
                event["terminal_resolution"] = "explicit_rank_reappeared"

    terminal_at_sample = identity_status in {
        "terminal_state",
        "disappeared_after_sample",
    }
    process_reappeared_after_terminal = bool(
        previous
        and previous.get("process_terminally_confirmed") is True
        and not terminal_at_sample
    )
    if terminal_at_sample:
        process_terminal_confirmation = (
            f"terminal_state_{snapshot.get('stat_state')}"
            if identity_status == "terminal_state"
            else "process_disappeared"
        )
        process_terminally_confirmed = True
    else:
        process_terminal_confirmation = (
            previous.get("process_terminal_confirmation") if previous else None
        )
        process_terminally_confirmed = bool(
            previous and previous.get("process_terminally_confirmed")
        )
    if retain_pending_rank and terminal_at_sample:
        resolution = (
            f"terminal_state_{snapshot.get('stat_state')}"
            if identity_status == "terminal_state"
            else "process_disappeared"
        )
        for event in environment_events:
            if event.get("terminal_resolution") == "pending":
                event["terminal_resolution"] = resolution
        pending_environment_loss = False
        terminally_confirmed = True
        terminal_confirmation = resolution
    else:
        pending_environment_loss = retain_pending_rank
        terminally_confirmed = bool(
            previous and previous.get("rank_environment_terminally_confirmed")
        )
        terminal_confirmation = (
            previous.get("rank_environment_terminal_confirmation")
            if previous
            else None
        )
    accumulated.update(
        {
            "observed_rank_ids": rank_history,
            "observed_cpu_masks": mask_history,
            "observed_rank_observation_statuses": status_history,
            "observed_process_states": state_history,
            "observed_process_starttimes": starttime_history,
            "current_sample_matches_assigned_singleton": sample_matches,
            "thread_affinity_sample_count": int(
                previous.get("thread_affinity_sample_count", 0) if previous else 0
            )
            + 1,
            "thread_affinity_scan_statuses": thread_status_history,
            "thread_affinity_evidence_sources": thread_source_history,
            "thread_affinity_scan_issues_ever": thread_issue_history,
            "thread_affinity_process_statuses": thread_process_status_history,
            "observed_thread_cpu_masks": thread_mask_history,
            "observed_thread_identities": thread_identity_history,
            "current_thread_affinity_sample_exact": thread_sample_matches,
            "all_thread_affinity_samples_exact": bool(
                (previous is None or previous.get("all_thread_affinity_samples_exact"))
                and thread_sample_matches
            ),
            "thread_affinity_violation_ever": bool(
                (previous and previous.get("thread_affinity_violation_ever"))
                or not thread_sample_matches
            ),
            "process_starttime_changed_ever": bool(
                (previous and previous.get("process_starttime_changed_ever"))
                or process_starttime_changed
            ),
            "process_terminally_confirmed": process_terminally_confirmed,
            "process_terminal_confirmation": process_terminal_confirmation,
            "process_reappeared_after_terminal_ever": bool(
                (previous and previous.get("process_reappeared_after_terminal_ever"))
                or process_reappeared_after_terminal
            ),
            "executable_changed_ever": bool(
                (previous and previous.get("executable_changed_ever"))
                or identity_status == "executable_changed_during_sample"
            ),
            "cpu_mask_changed_during_sample_ever": bool(
                (previous and previous.get("cpu_mask_changed_during_sample_ever"))
                or identity_status == "cpu_mask_changed_during_sample"
            ),
            "snapshot_unavailable_ever": bool(
                previous and previous.get("snapshot_unavailable_ever")
            ),
            "rank_identity_changed_ever": bool(
                (previous and previous.get("rank_identity_changed_ever"))
                or rank_identity_changed
                or explicit_rank_reappeared_after_loss
            ),
            "rank_environment_unavailable_ever": bool(
                (previous and previous.get("rank_environment_unavailable_ever"))
                or unavailable_environment
            ),
            "rank_environment_unavailable_sample_count": int(
                previous.get("rank_environment_unavailable_sample_count", 0)
                if previous
                else 0
            )
            + int(unavailable_environment),
            "rank_environment_unavailable_pending": pending_environment_loss,
            "rank_environment_terminally_confirmed": terminally_confirmed,
            "rank_environment_terminal_confirmation": terminal_confirmation,
            "rank_environment_events": environment_events,
            "affinity_violation_ever": bool(
                (previous and previous.get("affinity_violation_ever"))
                or not sample_matches
                or not thread_sample_matches
                or rank_identity_changed
                or process_starttime_changed
                or process_snapshot_inconsistent
                or process_reappeared_after_terminal
                or explicit_rank_reappeared_after_loss
            ),
        }
    )
    return accumulated


def _linux_process_terminal_resolution(
    pid: int, expected_starttime: int, proc_root: Path = Path("/proc")
) -> str | None:
    """Resolve one pending rank-environment loss without trusting PID alone."""
    root = proc_root / str(pid)
    if not root.is_dir():
        return "process_disappeared"
    try:
        state, starttime = _linux_proc_stat_identity(
            (root / "stat").read_text(errors="replace")
        )
    except FileNotFoundError:
        return "process_disappeared"
    except (PermissionError, ProcessLookupError, OSError, ValueError):
        return "identity_unreadable"
    if starttime != expected_starttime:
        return "pid_reused"
    if state in {"Z", "X"}:
        return f"terminal_state_{state}"
    return None


def _resolve_pending_rank_environment(
    record: dict[str, object], resolution: str
) -> None:
    """Resolve or reject one terminal-tail environment-loss observation."""
    terminal = resolution == "process_disappeared" or resolution.startswith(
        "terminal_state_"
    )
    events = [dict(event) for event in record.get("rank_environment_events", [])]
    for event in events:
        if event.get("terminal_resolution") == "pending":
            event["terminal_resolution"] = resolution
    record["rank_environment_events"] = events
    record["rank_environment_unavailable_pending"] = False
    record["rank_environment_terminally_confirmed"] = terminal
    record["rank_environment_terminal_confirmation"] = resolution
    _resolve_rank_process_lifetime(record, resolution)
    if not terminal:
        record["rank_identity_changed_ever"] = True


def _resolve_rank_process_lifetime(
    record: dict[str, object], resolution: str
) -> None:
    """Persist the terminal proof for every observed CP2K rank process."""
    terminal = resolution == "process_disappeared" or resolution.startswith(
        "terminal_state_"
    )
    record["process_terminally_confirmed"] = terminal
    record["process_terminal_confirmation"] = resolution
    if not terminal:
        record["affinity_violation_ever"] = True
        record["current_sample_matches_assigned_singleton"] = False
        if resolution == "pid_reused":
            record["process_starttime_changed_ever"] = True


def _observed_rank_process_is_still_live(
    pid: int,
    record: dict[str, object],
    proc_root: Path = Path("/proc"),
) -> bool:
    """Track one proven rank task directly, even outside descendant scans."""
    starttime = record.get("process_starttime")
    if isinstance(starttime, int) and not isinstance(starttime, bool):
        resolution = _linux_process_terminal_resolution(pid, starttime, proc_root)
    else:
        resolution = "identity_unreadable"
    if resolution is None:
        record["snapshot_unavailable_ever"] = True
        record["current_sample_matches_assigned_singleton"] = False
        record["affinity_violation_ever"] = True
        return True
    if record.get("rank_environment_unavailable_pending"):
        _resolve_pending_rank_environment(record, resolution)
    else:
        _resolve_rank_process_lifetime(record, resolution)
    return False


class ExecutionPool:
    """Allocate disjoint ordered PE lists to concurrent CP2K MPI launchers."""

    def __init__(
        self,
        *,
        concurrent_jobs: int,
        mpi_ranks_per_job: int,
        threads_per_rank: int,
        mpi_launcher: str | Path,
        mpi_launcher_args: Sequence[str],
        pe_lists: Sequence[str],
        check_current_affinity: bool = True,
        cpu_reservation_lock_root: Path | None = None,
    ) -> None:
        validate_mpi_launcher_args(mpi_launcher_args)
        available: set[int] | None = None
        if check_current_affinity and hasattr(os, "sched_getaffinity"):
            available = set(os.sched_getaffinity(0))
        parsed_pe_lists = validate_pe_lists(
            pe_lists,
            concurrent_jobs,
            mpi_ranks_per_job,
            threads_per_rank,
            available_cpus=available,
        )
        if check_current_affinity:
            require_single_pu_cores(parsed_pe_lists)
        resolved_launcher = resolve_executable(mpi_launcher, "MPI launcher")
        if "taskset" in resolved_launcher.name.lower():
            raise ValueError("outer taskset launchers are forbidden for production MPI")
        normalized_launcher_args = tuple(mpi_launcher_args)
        normalized_pe_lists = tuple(
            ",".join(str(cpu) for cpu in cpus) for cpus in parsed_pe_lists
        )
        reservation_lock_root = (
            cpu_reservation_lock_root or default_cpu_reservation_lock_root()
        ).resolve()
        selected_cpus = [cpu for pe_list in parsed_pe_lists for cpu in pe_list]
        reservation_handles: list[IO[str]] = []
        try:
            reservation_handles = acquire_cpu_reservation_locks(
                selected_cpus, reservation_lock_root
            )
            require_no_live_compute_overlap(selected_cpus)
            available_queue: queue.Queue[str] = queue.Queue()
            for value in normalized_pe_lists:
                available_queue.put(value)
            contract: dict[str, object] = {
                "schema_version": SCHEMA_VERSION,
                "mode": "openmpi_ordered_pe_list",
                "concurrent_jobs": concurrent_jobs,
                "mpi_ranks_per_job": mpi_ranks_per_job,
                "threads_per_rank": threads_per_rank,
                "mpi_launcher": str(resolved_launcher),
                "mpi_launcher_sha256": sha256(resolved_launcher),
                "mpi_launcher_args": list(normalized_launcher_args),
                "ordered_pe_lists": list(normalized_pe_lists),
                "ordered_pe_lists_disjoint": True,
                "cross_process_cpu_reservation": "flock_per_logical_cpu",
                "live_compute_overlap_preflight": (
                    "linux_procfs_live_cp2k_or_mpi_rank_allowed_cpu_masks"
                ),
                "live_compute_overlap_runtime_monitor": (
                    "every_affinity_sample_excluding_verified_own_pid_starttimes"
                ),
                "cpu_reservation_lock_root": str(reservation_lock_root),
                "exact_cpus_per_job": mpi_ranks_per_job,
                "exact_cpus_per_rank": 1,
                "single_pu_core_preflight": True,
                "mpi_map_by": "pe-list=<ordered-list>:ordered",
                "mpi_bind_to": "core",
                "mpi_report_bindings": True,
                "outer_taskset": False,
                "openmp_environment": dict(OPENMP_THREAD_ENVIRONMENT),
                "blas_environment": dict(BLAS_THREAD_ENVIRONMENT),
                "removed_mpi_binding_environment_key_policy": (
                    "all inherited placement, topology, threading, and preload "
                    "override variables"
                ),
                "binding_environment_scrub_contract": (
                    binding_environment_scrub_contract()
                ),
                "pool_close_policy": "reject_while_run_admitted_or_active",
                "rank_affinity_observation": "linux_proc_per_task_tid_starttime",
            }
            contract_sha256 = canonical_sha256(contract)

            # Publish owned handles only after every potentially failing
            # initialization step.  Until then the local list is the single
            # cleanup authority for the whole post-acquisition region.
            self.concurrent_jobs = concurrent_jobs
            self.mpi_ranks_per_job = mpi_ranks_per_job
            self.threads_per_rank = threads_per_rank
            self.mpi_launcher = resolved_launcher
            self.mpi_launcher_args = normalized_launcher_args
            self.pe_lists = normalized_pe_lists
            self.cpu_reservation_lock_root = reservation_lock_root
            self._available = available_queue
            self._active: set[str] = set()
            self._lifecycle_lock = threading.Lock()
            self._admitted_runs = 0
            self.contract = contract
            self.contract_sha256 = contract_sha256
            self._reservation_handles = reservation_handles
            self._closed = False
        except BaseException:
            for handle in reservation_handles:
                handle.close()
            self._reservation_handles = []
            self._closed = True
            raise

    def close(self) -> None:
        lifecycle_lock = getattr(self, "_lifecycle_lock", None)
        if lifecycle_lock is None:
            return
        with lifecycle_lock:
            if getattr(self, "_closed", True) and not getattr(
                self, "_reservation_handles", []
            ):
                return
            if self._admitted_runs or self._active:
                raise RuntimeError(
                    "cannot release CPU reservations while CP2K runs are "
                    "admitted or active"
                )
            # Make admission fail closed before the first individual lock is
            # released.  Close every handle even if one close reports an error;
            # failed handles remain owned so a later close can retry them.
            self._closed = True
            handles = self._reservation_handles
            self._reservation_handles = []
            first_error: BaseException | None = None
            for handle in handles:
                try:
                    handle.close()
                except BaseException as error:
                    self._reservation_handles.append(handle)
                    if first_error is None:
                        first_error = error
            if first_error is not None:
                raise RuntimeError(
                    "failed to release every CPU reservation; the execution "
                    "pool remains closed"
                ) from first_error

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            # Destructors cannot safely report lifecycle misuse.  In particular,
            # never release reservations from here while a run is admitted.
            pass

    def run_cp2k(
        self,
        cp2k: Path,
        inp: Path,
        out: Path,
    ) -> tuple[int, dict[str, object]]:
        admission_owned = False
        pe_list: str | None = None
        proc: subprocess.Popen[bytes] | None = None
        observed: dict[int, dict[str, object]] = {}
        with self._lifecycle_lock:
            if self._closed:
                raise RuntimeError(
                    "execution pool CPU reservations were already released"
                )
            self._admitted_runs += 1
            admission_owned = True
        try:
            cp2k = cp2k.resolve(strict=True)
            if not cp2k.is_file() or not os.access(cp2k, os.X_OK):
                raise ValueError(f"CP2K executable is not executable: {cp2k}")
            pe_list = self._available.get()
            with self._lifecycle_lock:
                if pe_list in self._active:
                    raise RuntimeError(
                        f"ordered PE list was allocated twice: {pe_list}"
                    )
                self._active.add(pe_list)
            assigned_cpus = parse_ordered_pe_list(pe_list)
            require_no_live_compute_overlap(assigned_cpus)
            live_overlap_preflight_gate = True
            cp2k_resolved = resolve_executable(cp2k, "CP2K executable")
            cp2k_sha256_at_launch = sha256(cp2k_resolved)
            input_resolved = inp.resolve(strict=True)
            input_sha256_at_launch = sha256(input_resolved)
            out.parent.mkdir(parents=True, exist_ok=True)
            launcher_sha256_at_launch = sha256(self.mpi_launcher)
            if launcher_sha256_at_launch != self.contract["mpi_launcher_sha256"]:
                raise RuntimeError("MPI launcher changed after execution-pool creation")
            log_path = launcher_log_path(out)
            for stale in (out, execution_record_path(out), log_path):
                if stale.exists():
                    stale.unlink()
            main_log = inp.parent / "mainLog.out"
            if main_log.exists():
                main_log.unlink()
            (
                env,
                removed_binding_environment,
                library_path_preserved,
            ) = sanitized_launch_environment(os.environ)
            residual_binding_environment = [
                key
                for key in mpi_control_environment_keys(env)
                if key not in SANCTIONED_THREAD_ENVIRONMENT
            ]
            binding_environment_scrub_gate = bool(
                not residual_binding_environment and library_path_preserved
            )
            command = [
                str(self.mpi_launcher),
                *self.mpi_launcher_args,
                "--map-by",
                f"pe-list={pe_list}:ordered",
                "--bind-to",
                "core",
                "--report-bindings",
                "-np",
                str(self.mpi_ranks_per_job),
                str(cp2k_resolved),
                "-i",
                str(input_resolved),
                "-o",
                str(out.resolve()),
            ]
            started = datetime.now(timezone.utc).isoformat()
            concurrent_duplicate_rank_ids_ever: set[int] = set()
            concurrent_duplicate_rank_samples: list[dict[str, object]] = []
            live_overlap_runtime_samples: list[dict[str, object]] = []
            local_affinity_violation_samples: list[dict[str, object]] = []
            affinity_sample_index = 0
            with log_path.open("wb") as launcher_log:
                proc = subprocess.Popen(
                    command,
                    cwd=inp.parent,
                    stdout=launcher_log,
                    stderr=subprocess.STDOUT,
                    env=env,
                    start_new_session=True,
                )
                launcher_starttime = _linux_process_starttime(proc.pid)
                if Path("/proc").is_dir() and launcher_starttime is None:
                    raise RuntimeError(
                        "could not prove the MPI launcher PID/starttime identity"
                    )
                while True:
                    affinity_sample_index += 1
                    snapshots_this_sample: list[dict[str, object]] = []
                    descendant_pids = {proc.pid}
                    if Path("/proc").is_dir():
                        descendant_pids = _linux_descendants(proc.pid)
                        tracked_rank_pids = {
                            pid
                            for pid, record in observed.items()
                            if record.get("is_cp2k_rank") is True
                            and record.get("process_terminally_confirmed") is not True
                        }
                        candidate_pids = (
                            descendant_pids | tracked_rank_pids
                        )
                        for pid in candidate_pids:
                            snapshot = _linux_process_snapshot(pid, cp2k_resolved)
                            if snapshot is None:
                                tracked = observed.get(pid)
                                if tracked and tracked.get("is_cp2k_rank") is True:
                                    if _observed_rank_process_is_still_live(
                                        pid, tracked
                                    ):
                                        snapshots_this_sample.append(tracked)
                                continue
                            observed[pid] = _accumulate_process_snapshot(
                                observed.get(pid), snapshot, assigned_cpus
                            )
                            snapshots_this_sample.append(observed[pid])
                    verified_own_identities: dict[int, int] = {}
                    if launcher_starttime is not None:
                        verified_own_identities[proc.pid] = launcher_starttime
                    for pid, record in observed.items():
                        starttime = record.get("process_starttime")
                        if (
                            record.get("is_cp2k_rank") is True
                            and record.get("process_terminally_confirmed") is not True
                            and isinstance(starttime, int)
                            and not isinstance(starttime, bool)
                            and _linux_process_terminal_resolution(pid, starttime)
                            is None
                        ):
                            verified_own_identities[pid] = starttime
                    external_owners = live_compute_cpu_owners(
                        assigned_cpus,
                        ignore_process_identities=verified_own_identities,
                    )
                    if external_owners:
                        live_overlap_runtime_samples.append(
                            {
                                "sample_index": affinity_sample_index,
                                "owners": external_owners,
                            }
                        )
                    duplicate_groups = _concurrent_live_rank_pid_groups(
                        snapshots_this_sample
                    )
                    concurrent_duplicate_rank_ids_ever.update(duplicate_groups)
                    if duplicate_groups:
                        concurrent_duplicate_rank_samples.append(
                            {
                                "sample_index": affinity_sample_index,
                                "rank_pid_groups": [
                                    {
                                        "ompi_comm_world_rank": rank,
                                        "pids": pids,
                                    }
                                    for rank, pids in duplicate_groups.items()
                                ],
                            }
                        )
                    local_affinity_violations = [
                        {
                            "pid": record.get("pid"),
                            "rank": record.get("ompi_comm_world_rank"),
                            "leader_mask": record.get("cpus_allowed_list"),
                            "thread_scan_status": record.get(
                                "thread_affinity_scan_status"
                            ),
                            "thread_masks": record.get(
                                "observed_thread_cpu_masks", []
                            ),
                        }
                        for record in snapshots_this_sample
                        if record.get("is_cp2k_rank") is True
                        and (
                            record.get("current_sample_matches_assigned_singleton")
                            is not True
                            or record.get("current_thread_affinity_sample_exact")
                            is not True
                            or record.get("affinity_violation_ever") is True
                        )
                    ]
                    if local_affinity_violations:
                        local_affinity_violation_samples.append(
                            {
                                "sample_index": affinity_sample_index,
                                "ranks": local_affinity_violations,
                            }
                        )
                    if external_owners or local_affinity_violations:
                        _terminate_and_reap_process_group(
                            proc, tracked_rank_starttimes=(
                                _tracked_rank_process_starttimes(observed)
                            )
                        )
                        return_code = proc.returncode if proc.returncode is not None else 97
                        break
                    try:
                        return_code = proc.wait(timeout=0.05)
                        break
                    except subprocess.TimeoutExpired:
                        continue
            rank_live_after_launcher = False
            for pid, record in observed.items():
                if (
                    record.get("is_cp2k_rank") is not True
                    or record.get("process_terminally_confirmed") is True
                ):
                    continue
                starttime = record.get("process_starttime")
                if isinstance(starttime, int) and not isinstance(starttime, bool):
                    resolution = _linux_process_terminal_resolution(pid, starttime)
                else:
                    resolution = "identity_unreadable"
                final_resolution = resolution or "launcher_ended_while_process_live"
                rank_live_after_launcher = rank_live_after_launcher or resolution is None
                if record.get("rank_environment_unavailable_pending"):
                    _resolve_pending_rank_environment(record, final_resolution)
                else:
                    _resolve_rank_process_lifetime(record, final_resolution)
            if rank_live_after_launcher or _live_process_group_members(proc.pid):
                _terminate_and_reap_process_group(
                    proc,
                    tracked_rank_starttimes=(
                        _tracked_rank_process_starttimes(observed)
                    ),
                )
                if return_code == 0:
                    return_code = 97
            finished = datetime.now(timezone.utc).isoformat()
            if main_log.exists() and (
                not out.exists()
                or "PROGRAM ENDED" not in out.read_text(errors="ignore")
            ):
                shutil.copyfile(main_log, out)
            rank_generations = _aggregate_cp2k_rank_generations(
                observed,
                assigned_cpus,
                concurrent_duplicate_rank_ids_ever,
            )
            rank_ids = [
                int(record["ompi_comm_world_rank"]) for record in rank_generations
            ]
            rank_mask_texts: list[str] = []
            rank_masks: list[set[int]] = []
            for record in rank_generations:
                mask_history = record.get("observed_cpu_masks")
                if not isinstance(mask_history, list) or len(mask_history) != 1:
                    rank_mask_texts = []
                    rank_masks = []
                    break
                mask_text = str(mask_history[0])
                try:
                    rank_masks.append(parse_cpu_set(mask_text))
                except ValueError:
                    rank_mask_texts = []
                    rank_masks = []
                    break
                rank_mask_texts.append(mask_text)
            unranked_cp2k_process_seen = any(
                record.get("is_cp2k_rank") is True
                and not isinstance(record.get("ompi_comm_world_rank"), int)
                for record in observed.values()
            )
            rank_count_matches = len(rank_generations) == self.mpi_ranks_per_job
            rank_ids_exact = rank_ids == list(range(self.mpi_ranks_per_job))
            masks_complete = len(rank_masks) == len(rank_generations)
            masks_exact = masks_complete and all(
                0 <= rank < len(assigned_cpus)
                and mask == {assigned_cpus[rank]}
                for rank, mask in zip(rank_ids, rank_masks, strict=True)
            )
            all_rank_samples_exact = (
                rank_count_matches
                and not unranked_cp2k_process_seen
                and not concurrent_duplicate_rank_ids_ever
                and all(
                    record.get("all_samples_match_assigned_singleton") is True
                    for record in rank_generations
                )
            )
            launcher_text = log_path.read_text(errors="replace")
            reported_binding_rank_ids = _reported_binding_rank_ids(launcher_text)
            binding_report_complete = reported_binding_rank_ids == list(
                range(self.mpi_ranks_per_job)
            )
            cross_process_reservation_gate = (
                not self._closed
                and len(self._reservation_handles)
                == sum(len(parse_ordered_pe_list(value)) for value in self.pe_lists)
                and all(not handle.closed for handle in self._reservation_handles)
            )
            runtime_affinity_gate = (
                rank_count_matches
                and rank_ids_exact
                and masks_exact
                and all_rank_samples_exact
                and binding_report_complete
                and cross_process_reservation_gate
                and binding_environment_scrub_gate
                and live_overlap_preflight_gate
                and not live_overlap_runtime_samples
                and not local_affinity_violation_samples
            )
            observation: dict[str, object] = {
                "schema_version": SCHEMA_VERSION,
                "contract": self.contract,
                "contract_sha256": self.contract_sha256,
                "assigned_ordered_pe_list": pe_list,
                "assigned_cpu_count": len(assigned_cpus),
                "command": command,
                "working_directory": str(inp.parent.resolve()),
                "mpi_launcher_sha256_at_launch": launcher_sha256_at_launch,
                "cp2k": str(cp2k_resolved),
                "cp2k_sha256_at_launch": cp2k_sha256_at_launch,
                "input": str(input_resolved),
                "input_sha256_at_launch": input_sha256_at_launch,
                "output": str(out.resolve()),
                "return_code": return_code,
                "started_at_utc": started,
                "finished_at_utc": finished,
                "hostname": socket.gethostname(),
                "pid": os.getpid(),
                "thread_environment": {
                    **OPENMP_THREAD_ENVIRONMENT,
                    **BLAS_THREAD_ENVIRONMENT,
                },
                "removed_mpi_binding_environment_keys": removed_binding_environment,
                "binding_environment_scrub_gate": binding_environment_scrub_gate,
                "residual_binding_environment_keys": residual_binding_environment,
                "ld_library_path_preserved": library_path_preserved,
                "launcher_log": str(log_path.resolve()),
                "launcher_log_sha256": sha256(log_path),
                "reported_binding_rank_ids": reported_binding_rank_ids,
                "binding_report_complete": binding_report_complete,
                "cross_process_cpu_reservation_gate": (
                    cross_process_reservation_gate
                ),
                "live_compute_overlap_preflight_gate": live_overlap_preflight_gate,
                "live_compute_overlap_preflight_owners": [],
                "live_compute_overlap_runtime_gate": not live_overlap_runtime_samples,
                "live_compute_overlap_runtime_samples": live_overlap_runtime_samples,
                "local_affinity_violation_gate": not local_affinity_violation_samples,
                "local_affinity_violation_samples": local_affinity_violation_samples,
                "observed_child_processes": sorted(
                    observed.values(), key=lambda record: int(record["pid"])
                ),
                "observed_cp2k_rank_pids": [
                    int(record["canonical_pid"]) for record in rank_generations
                ],
                "observed_cp2k_rank_pid_generations": [
                    record["pid_generations"] for record in rank_generations
                ],
                "observed_cp2k_rank_evidence": rank_generations,
                "observed_cp2k_rank_ids": rank_ids,
                "observed_cp2k_rank_masks": rank_mask_texts,
                "observed_cp2k_rank_count": len(rank_generations),
                "observed_cp2k_process_generation_count": sum(
                    len(record["pid_generations"]) for record in rank_generations
                ),
                "concurrent_duplicate_rank_ids_ever": sorted(
                    concurrent_duplicate_rank_ids_ever
                ),
                "concurrent_duplicate_rank_samples": (
                    concurrent_duplicate_rank_samples
                ),
                "concurrent_duplicate_rank_processes_ever": bool(
                    concurrent_duplicate_rank_ids_ever
                ),
                "unranked_cp2k_process_seen": unranked_cp2k_process_seen,
                "expected_cp2k_rank_count": self.mpi_ranks_per_job,
                "rank_count_matches": rank_count_matches,
                "rank_ids_exactly_0_to_n_minus_1": rank_ids_exact,
                "rank_masks_complete": masks_complete,
                "rank_masks_exactly_match_ordered_pe_list": masks_exact,
                "all_observed_rank_samples_match_ordered_pe_list": (
                    all_rank_samples_exact
                ),
                "mpi_bind_to": "core",
                "runtime_affinity_gate": runtime_affinity_gate,
                "timing_classification": (
                    "production_scaling_eligible"
                    if return_code == 0 and runtime_affinity_gate
                    else "timing_non_scaling"
                ),
            }
            return return_code, observation
        except BaseException:
            if proc is not None:
                _terminate_and_reap_process_group(
                    proc,
                    tracked_rank_starttimes=(
                        _tracked_rank_process_starttimes(observed)
                    ),
                )
            raise
        finally:
            if admission_owned:
                # The admission count keeps close() from releasing the process-
                # wide reservation locks until cleanup is complete.  Reclaim
                # Queue ownership and Active membership under the same lock so
                # another waiter cannot observe a half-transition, including
                # when BaseException interrupts Queue -> Active.
                with self._lifecycle_lock:
                    try:
                        if pe_list is not None:
                            self._active.discard(pe_list)
                            # Once put_nowait returns, Queue ownership itself is
                            # authoritative; no second local ownership flag can
                            # introduce a post-transfer update window.
                            self._available.put_nowait(pe_list)
                    finally:
                        self._admitted_runs -= 1
                        admission_owned = False

    def write_record(
        self,
        output: Path,
        observation: Mapping[str, object],
        scientific_job_stamp: Path,
    ) -> dict[str, object]:
        if observation.get("contract") != self.contract:
            raise ValueError("execution observation contract differs from this pool")
        if observation.get("runtime_affinity_gate") is not True:
            raise ValueError(
                "execution observation did not prove the requested child-rank CPU masks"
            )
        if Path(str(observation.get("output", ""))).resolve() != output.resolve():
            raise ValueError("execution observation is bound to another output")
        if not output.is_file():
            raise ValueError(f"cannot finalize execution record without output {output}")
        if not scientific_job_stamp.is_file():
            raise ValueError(
                f"cannot finalize execution record without scientific stamp {scientific_job_stamp}"
            )
        payload = dict(observation)
        payload.update(
            {
                "output_sha256": sha256(output),
                "scientific_job_stamp": str(scientific_job_stamp.resolve()),
                "scientific_job_stamp_sha256": sha256(scientific_job_stamp),
            }
        )
        path = execution_record_path(output)
        atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return {"path": str(path.resolve()), "sha256": sha256(path)}

    def record_issue(self, output: Path, scientific_job_stamp: Path) -> str | None:
        path = execution_record_path(output)
        return recorded_execution_issue(
            path,
            self.contract,
            output,
            scientific_job_stamp,
        )
