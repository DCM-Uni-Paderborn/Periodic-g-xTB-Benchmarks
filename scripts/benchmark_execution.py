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
import shutil
import socket
import subprocess
import tempfile
import threading
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
    """Return every inherited Open MPI/PRRTE MCA control variable.

    Scrubbing only names containing ``hwloc`` or ``rmaps`` is insufficient:
    MCA parameter-file selectors can inject those settings indirectly.  The
    production command supplies its complete placement contract explicitly, so
    no inherited OMPI/PRRTE MCA variable is permitted.
    """
    return sorted(
        key
        for key in environment
        if key.startswith(("OMPI_MCA_", "PRTE_MCA_"))
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
    ignore_pids: set[int] | None = None,
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
    ignored = set(ignore_pids or ()) | {os.getpid()}
    if not selected or not proc_root.is_dir():
        return []
    owners: list[dict[str, object]] = []
    for directory in proc_root.iterdir():
        if not directory.name.isdigit():
            continue
        pid = int(directory.name)
        if pid in ignored:
            continue
        try:
            status = (directory / "status").read_text(errors="replace")
        except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
            continue
        fields: dict[str, str] = {}
        for line in status.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key] = value.strip()
        state = fields.get("State", "")
        if state.startswith(("Z", "X")):
            continue
        mask_text = fields.get("Cpus_allowed_list", "")
        try:
            allowed = parse_cpu_set(mask_text)
        except ValueError:
            continue
        overlap = selected & allowed
        if not overlap:
            continue
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
        if (
            item.get("observed_rank_ids") != [rank]
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
    return _common_artifact_issue(record, path, output, scientific_job_stamp)


def _linux_process_snapshot(pid: int, cp2k: Path) -> dict[str, object] | None:
    """Read one live process/rank affinity from procfs without external tools."""
    root = Path("/proc") / str(pid)
    if not root.is_dir():
        return None
    try:
        status = (root / "status").read_text(errors="replace")
        command = (root / "cmdline").read_bytes().split(b"\0")
        arguments = [item.decode(errors="replace") for item in command if item]
        environment_items = (root / "environ").read_bytes().split(b"\0")
        environment = {
            key.decode(errors="replace"): value.decode(errors="replace")
            for item in environment_items
            if item and b"=" in item
            for key, value in (item.split(b"=", 1),)
        }
        executable = str((root / "exe").resolve(strict=True))
    except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
        return None
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
    if rank_value is not None:
        try:
            mpi_rank = int(rank_value)
        except ValueError:
            mpi_rank = None
    return {
        "pid": pid,
        "ppid": int(fields.get("PPid", "0")),
        "state": fields.get("State", ""),
        "executable": executable,
        "arguments": arguments,
        "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
        "ompi_comm_world_rank": mpi_rank,
        "is_cp2k_rank": is_cp2k_rank,
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
        except (FileNotFoundError, PermissionError, ProcessLookupError, OSError, ValueError):
            children = []
        for child in children:
            if child not in found:
                found.add(child)
                pending.append(child)
    return found


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


def _aggregate_cp2k_rank_generations(
    observed: Mapping[int, Mapping[str, object]],
    assigned_cpus: Sequence[int],
    concurrent_duplicate_rank_ids: set[int] | None = None,
) -> list[dict[str, object]]:
    """Aggregate sequential same-rank PID generations without weakening gates."""
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
            and mask_history == {str(assigned_cpus[rank])}
            and rank not in duplicates
            and all(
                record.get("affinity_violation_ever") is False
                and record.get("rank_identity_changed_ever") is not True
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
    accumulated["sample_count"] = (
        int(previous.get("sample_count", 0)) + 1 if previous is not None else 1
    )
    if not snapshot.get("is_cp2k_rank"):
        accumulated["affinity_violation_ever"] = bool(
            previous and previous.get("affinity_violation_ever")
        )
        return accumulated

    rank = snapshot.get("ompi_comm_world_rank")
    mask_text = str(snapshot.get("cpus_allowed_list", ""))
    try:
        mask = parse_cpu_set(mask_text)
    except ValueError:
        mask = set()
    sample_matches = (
        isinstance(rank, int)
        and 0 <= rank < len(assigned_cpus)
        and mask == {assigned_cpus[rank]}
    )
    rank_history = list(previous.get("observed_rank_ids", [])) if previous else []
    mask_history = list(previous.get("observed_cpu_masks", [])) if previous else []
    if rank not in rank_history:
        rank_history.append(rank)
    if mask_text not in mask_history:
        mask_history.append(mask_text)
    rank_identity_changed = bool(
        previous
        and any(
            isinstance(prior_rank, int) and prior_rank != rank
            for prior_rank in previous.get("observed_rank_ids", [])
        )
    )
    accumulated.update(
        {
            "observed_rank_ids": rank_history,
            "observed_cpu_masks": mask_history,
            "current_sample_matches_assigned_singleton": sample_matches,
            "rank_identity_changed_ever": bool(
                (previous and previous.get("rank_identity_changed_ever"))
                or rank_identity_changed
            ),
            "affinity_violation_ever": bool(
                (previous and previous.get("affinity_violation_ever"))
                or not sample_matches
                or rank_identity_changed
            ),
        }
    )
    return accumulated


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
                    "all inherited OMPI_MCA_*/PRTE_MCA_* control variables"
                ),
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
            self._active_lock = threading.Lock()
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
        if getattr(self, "_closed", True):
            return
        for handle in self._reservation_handles:
            handle.close()
        self._reservation_handles = []
        self._closed = True

    def __del__(self) -> None:
        self.close()

    def run_cp2k(
        self,
        cp2k: Path,
        inp: Path,
        out: Path,
    ) -> tuple[int, dict[str, object]]:
        if self._closed:
            raise RuntimeError("execution pool CPU reservations were already released")
        cp2k = cp2k.resolve(strict=True)
        if not cp2k.is_file() or not os.access(cp2k, os.X_OK):
            raise ValueError(f"CP2K executable is not executable: {cp2k}")
        pe_list = self._available.get()
        with self._active_lock:
            if pe_list in self._active:
                self._available.put(pe_list)
                raise RuntimeError(f"ordered PE list was allocated twice: {pe_list}")
            self._active.add(pe_list)
        try:
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
            env = os.environ.copy()
            env.update(BLAS_THREAD_ENVIRONMENT)
            env.update(OPENMP_THREAD_ENVIRONMENT)
            removed_binding_environment = mpi_control_environment_keys(env)
            for key in removed_binding_environment:
                env.pop(key, None)
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
            observed: dict[int, dict[str, object]] = {}
            concurrent_duplicate_rank_ids_ever: set[int] = set()
            concurrent_duplicate_rank_samples: list[dict[str, object]] = []
            affinity_sample_index = 0
            with log_path.open("wb") as launcher_log:
                proc = subprocess.Popen(
                    command,
                    cwd=inp.parent,
                    stdout=launcher_log,
                    stderr=subprocess.STDOUT,
                    env=env,
                )
                while True:
                    affinity_sample_index += 1
                    snapshots_this_sample: list[dict[str, object]] = []
                    if Path("/proc").is_dir():
                        for pid in _linux_descendants(proc.pid):
                            snapshot = _linux_process_snapshot(pid, cp2k_resolved)
                            if snapshot is None:
                                continue
                            snapshots_this_sample.append(snapshot)
                            observed[pid] = _accumulate_process_snapshot(
                                observed.get(pid), snapshot, assigned_cpus
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
                    try:
                        return_code = proc.wait(timeout=0.05)
                        break
                    except subprocess.TimeoutExpired:
                        continue
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
                and live_overlap_preflight_gate
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
                "launcher_log": str(log_path.resolve()),
                "launcher_log_sha256": sha256(log_path),
                "reported_binding_rank_ids": reported_binding_rank_ids,
                "binding_report_complete": binding_report_complete,
                "cross_process_cpu_reservation_gate": (
                    cross_process_reservation_gate
                ),
                "live_compute_overlap_preflight_gate": live_overlap_preflight_gate,
                "live_compute_overlap_preflight_owners": [],
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
        finally:
            with self._active_lock:
                self._active.remove(pe_list)
            self._available.put(pe_list)

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
