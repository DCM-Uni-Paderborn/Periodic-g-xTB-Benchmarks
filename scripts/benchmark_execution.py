#!/usr/bin/env python3
"""Fail-closed MPI/affinity execution records for benchmark CP2K jobs.

The scientific job stamp remains owned by the benchmark driver.  This
module writes a separate, additive record that binds the exact launcher,
``--bind-to none`` policy, taskset CPU allocation, input, output, and completed
scientific stamp without changing the scientific stamp schema or matcher.
"""

from __future__ import annotations

import hashlib
import json
import os
import queue
import shutil
import socket
import subprocess
import tempfile
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence


SCHEMA_VERSION = 1
BLAS_THREAD_ENVIRONMENT = {
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
MPI_NO_REBIND_ENVIRONMENT = {
    "OMPI_MCA_hwloc_base_binding_policy": "none",
    "PRTE_MCA_hwloc_base_binding_policy": "none",
}


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


def require_bind_to_none(arguments: Sequence[str]) -> None:
    """Reject absent, conflicting, or internally rebinding MPI policies."""
    bindings: list[str] = []
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument == "--bind-to":
            if index + 1 >= len(arguments):
                raise ValueError("MPI launcher --bind-to has no value")
            bindings.append(arguments[index + 1].lower())
            index += 2
            continue
        if argument.startswith("--bind-to="):
            bindings.append(argument.split("=", 1)[1].lower())
        index += 1
    if bindings != ["none"]:
        raise ValueError(
            "MPI launcher arguments must contain exactly one '--bind-to none' "
            f"and no internal rebinding policy; observed {bindings or '<none>'}"
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
    if record.get("schema_version") != SCHEMA_VERSION:
        return f"execution record schema mismatch in {path}"
    if record.get("contract") != expected_contract:
        return f"execution contract mismatch in {path}"
    expected_contract_sha = canonical_sha256(dict(expected_contract))
    if record.get("contract_sha256") != expected_contract_sha:
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
    cp2k_path = Path(str(record.get("cp2k", "")))
    try:
        cp2k_resolved = cp2k_path.resolve(strict=True)
    except (FileNotFoundError, OSError):
        return f"missing recorded CP2K executable in {path}"
    if not cp2k_resolved.is_file() or not os.access(cp2k_resolved, os.X_OK):
        return f"recorded CP2K executable is not executable in {path}"
    cp2k_sha256 = sha256(cp2k_resolved)
    if record.get("cp2k_sha256_at_launch") != cp2k_sha256:
        return f"CP2K executable hash mismatch in {path}"
    input_path = Path(str(record.get("input", "")))
    try:
        input_resolved = input_path.resolve(strict=True)
    except (FileNotFoundError, OSError):
        return f"missing recorded execution input in {path}"
    if record.get("input_sha256_at_launch") != sha256(input_resolved):
        return f"execution input hash mismatch in {path}"
    if Path(str(record.get("output", ""))).resolve() != output.resolve():
        return f"execution output path mismatch in {path}"
    if not output.is_file() or record.get("output_sha256") != sha256(output):
        return f"execution output hash mismatch in {path}"
    if Path(str(record.get("scientific_job_stamp", ""))).resolve() != scientific_job_stamp.resolve():
        return f"scientific job-stamp path mismatch in {path}"
    if (
        not scientific_job_stamp.is_file()
        or record.get("scientific_job_stamp_sha256") != sha256(scientific_job_stamp)
    ):
        return f"scientific job-stamp hash mismatch in {path}"
    try:
        scientific_signature = json.loads(scientific_job_stamp.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return f"invalid scientific job stamp for {path}: {exc}"
    if not isinstance(scientific_signature, Mapping):
        return f"invalid scientific job-stamp payload for {path}"
    try:
        stamped_cp2k = Path(str(scientific_signature.get("executable", ""))).resolve(
            strict=True
        )
    except (FileNotFoundError, OSError):
        return f"scientific job stamp has no reusable CP2K executable for {path}"
    if stamped_cp2k != cp2k_resolved:
        return f"recorded CP2K executable differs from scientific job stamp in {path}"
    if scientific_signature.get("executable_sha256") != cp2k_sha256:
        return f"scientific job-stamp CP2K hash mismatch in {path}"
    try:
        stamped_input = Path(str(scientific_signature.get("input", ""))).resolve(
            strict=True
        )
    except (FileNotFoundError, OSError):
        return f"scientific job stamp has no reusable input for {path}"
    if stamped_input != input_resolved:
        return f"recorded input differs from scientific job stamp in {path}"
    if scientific_signature.get("input_sha256") != sha256(input_resolved):
        return f"scientific job-stamp input hash mismatch in {path}"
    expected_command = cp2k_command(
        taskset=taskset,
        cpu_set=assigned,
        mpi_launcher=launcher,
        mpi_launcher_args=launcher_args,
        mpi_ranks_per_job=ranks,
        cp2k=cp2k_resolved,
        inp=input_resolved,
        out=output,
    )
    if record.get("command") != expected_command:
        return f"full execution command/affinity mismatch in {path}"
    return None


def _linux_process_snapshot(pid: int, cp2k: Path) -> dict[str, object] | None:
    """Read one live process/rank affinity from procfs without external tools."""
    root = Path("/proc") / str(pid)
    if not root.is_dir():
        return None
    try:
        status = (root / "status").read_text(errors="replace")
        command = (root / "cmdline").read_bytes().split(b"\0")
        arguments = [item.decode(errors="replace") for item in command if item]
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
    return {
        "pid": pid,
        "ppid": int(fields.get("PPid", "0")),
        "executable": executable,
        "arguments": arguments,
        "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
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


class ExecutionPool:
    """Allocate disjoint taskset masks to concurrent CP2K MPI launchers."""

    def __init__(
        self,
        *,
        concurrent_jobs: int,
        mpi_ranks_per_job: int,
        threads_per_rank: int,
        mpi_launcher: str | Path,
        mpi_launcher_args: Sequence[str],
        cpu_sets: Sequence[str],
        taskset: str | Path,
        check_current_affinity: bool = True,
    ) -> None:
        require_bind_to_none(mpi_launcher_args)
        available: set[int] | None = None
        if check_current_affinity and hasattr(os, "sched_getaffinity"):
            available = set(os.sched_getaffinity(0))
        validate_cpu_sets(
            cpu_sets,
            concurrent_jobs,
            mpi_ranks_per_job,
            threads_per_rank,
            available_cpus=available,
        )
        self.concurrent_jobs = concurrent_jobs
        self.mpi_ranks_per_job = mpi_ranks_per_job
        self.threads_per_rank = threads_per_rank
        self.mpi_launcher = resolve_executable(mpi_launcher, "MPI launcher")
        self.mpi_launcher_args = tuple(mpi_launcher_args)
        self.taskset = resolve_executable(taskset, "taskset launcher")
        self.cpu_sets = tuple(cpu_sets)
        self._available: queue.Queue[str] = queue.Queue()
        for value in self.cpu_sets:
            self._available.put(value)
        self._active: set[str] = set()
        self._active_lock = threading.Lock()
        self.contract: dict[str, object] = {
            "schema_version": SCHEMA_VERSION,
            "mode": "taskset_mpi",
            "concurrent_jobs": concurrent_jobs,
            "mpi_ranks_per_job": mpi_ranks_per_job,
            "threads_per_rank": threads_per_rank,
            "mpi_launcher": str(self.mpi_launcher),
            "mpi_launcher_sha256": sha256(self.mpi_launcher),
            "mpi_launcher_args": list(self.mpi_launcher_args),
            "mpi_bind_to": "none",
            "taskset": str(self.taskset),
            "taskset_sha256": sha256(self.taskset),
            "cpu_sets": list(self.cpu_sets),
            "cpu_sets_disjoint": True,
            "minimum_cpu_set_size": mpi_ranks_per_job * threads_per_rank,
            "mpi_no_rebind_environment": dict(MPI_NO_REBIND_ENVIRONMENT),
        }
        self.contract_sha256 = canonical_sha256(self.contract)

    def run_cp2k(
        self,
        cp2k: Path,
        inp: Path,
        out: Path,
        ) -> tuple[int, dict[str, object]]:
        cpu_set = self._available.get()
        with self._active_lock:
            if cpu_set in self._active:
                self._available.put(cpu_set)
                raise RuntimeError(f"CPU set was allocated twice: {cpu_set}")
            self._active.add(cpu_set)
        try:
            cp2k_resolved = resolve_executable(cp2k, "CP2K executable")
            cp2k_sha256_at_launch = sha256(cp2k_resolved)
            input_sha256_at_launch = sha256(inp.resolve(strict=True))
            out.parent.mkdir(parents=True, exist_ok=True)
            for stale in (out, execution_record_path(out)):
                if stale.exists():
                    stale.unlink()
            main_log = inp.parent / "mainLog.out"
            if main_log.exists():
                main_log.unlink()
            env = os.environ.copy()
            env.update(BLAS_THREAD_ENVIRONMENT)
            env.update(MPI_NO_REBIND_ENVIRONMENT)
            env["OMP_NUM_THREADS"] = str(self.threads_per_rank)
            env["OMP_PROC_BIND"] = "false"
            env["OMP_WAIT_POLICY"] = "PASSIVE"
            command = cp2k_command(
                taskset=self.taskset,
                cpu_set=cpu_set,
                mpi_launcher=self.mpi_launcher,
                mpi_launcher_args=self.mpi_launcher_args,
                mpi_ranks_per_job=self.mpi_ranks_per_job,
                cp2k=cp2k_resolved,
                inp=inp,
                out=out,
            )
            started = datetime.now(timezone.utc).isoformat()
            proc = subprocess.Popen(
                command,
                cwd=inp.parent,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
                env=env,
            )
            observed: dict[int, dict[str, object]] = {}
            while True:
                if Path("/proc").is_dir():
                    for pid in _linux_descendants(proc.pid):
                        snapshot = _linux_process_snapshot(pid, cp2k_resolved)
                        if snapshot is None:
                            continue
                        previous = observed.get(pid)
                        snapshot["sample_count"] = (
                            int(previous.get("sample_count", 0)) + 1
                            if previous is not None
                            else 1
                        )
                        observed[pid] = snapshot
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
            assigned_cpus = parse_cpu_set(cpu_set)
            rank_processes = sorted(
                (record for record in observed.values() if record["is_cp2k_rank"]),
                key=lambda record: int(record["pid"]),
            )
            rank_masks = [
                parse_cpu_set(str(record["cpus_allowed_list"]))
                for record in rank_processes
                if record.get("cpus_allowed_list")
            ]
            rank_count_matches = len(rank_processes) == self.mpi_ranks_per_job
            masks_complete = len(rank_masks) == len(rank_processes)
            masks_exact = masks_complete and all(mask == assigned_cpus for mask in rank_masks)
            runtime_affinity_gate = rank_count_matches and masks_exact
            observation: dict[str, object] = {
                "schema_version": SCHEMA_VERSION,
                "contract": self.contract,
                "contract_sha256": self.contract_sha256,
                "assigned_cpu_set": cpu_set,
                "assigned_cpu_count": len(parse_cpu_set(cpu_set)),
                "command": command,
                "working_directory": str(inp.parent.resolve()),
                "cp2k": str(cp2k_resolved),
                "cp2k_sha256_at_launch": cp2k_sha256_at_launch,
                "input": str(inp.resolve()),
                "input_sha256_at_launch": input_sha256_at_launch,
                "output": str(out.resolve()),
                "return_code": return_code,
                "started_at_utc": started,
                "finished_at_utc": finished,
                "hostname": socket.gethostname(),
                "pid": os.getpid(),
                "thread_environment": {
                    "OMP_NUM_THREADS": env["OMP_NUM_THREADS"],
                    "OMP_PROC_BIND": env["OMP_PROC_BIND"],
                    "OMP_WAIT_POLICY": env["OMP_WAIT_POLICY"],
                    **BLAS_THREAD_ENVIRONMENT,
                    **MPI_NO_REBIND_ENVIRONMENT,
                },
                "observed_child_processes": sorted(
                    observed.values(), key=lambda record: int(record["pid"])
                ),
                "observed_cp2k_rank_pids": [
                    int(record["pid"]) for record in rank_processes
                ],
                "observed_cp2k_rank_masks": [
                    str(record["cpus_allowed_list"]) for record in rank_processes
                ],
                "observed_cp2k_rank_count": len(rank_processes),
                "expected_cp2k_rank_count": self.mpi_ranks_per_job,
                "rank_count_matches": rank_count_matches,
                "rank_masks_complete": masks_complete,
                "rank_masks_exactly_match_taskset": masks_exact,
                "mpiexec_internal_rebinding_detected": not masks_exact,
                "runtime_affinity_gate": runtime_affinity_gate,
            }
            return return_code, observation
        finally:
            with self._active_lock:
                self._active.remove(cpu_set)
            self._available.put(cpu_set)

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
