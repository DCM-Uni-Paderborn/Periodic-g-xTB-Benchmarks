#!/usr/bin/env python3
"""Run the fixed mixer symmetry-star storage qualification matrix."""

from __future__ import annotations

import concurrent.futures
import fcntl
import hashlib
import json
import os
import re
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import IO, Sequence


ROOT = Path(__file__).resolve().parent
MATRIX = json.loads((ROOT / "test_matrix.json").read_text())
INPUT_ROOT = ROOT / "test_inputs"
CP2K = Path(os.environ["CP2K_EXE"]).resolve()
CP2K_LIB = Path(os.environ["CP2K_LIB"]).resolve()
MPIEXEC = Path(os.environ["MPIEXEC_EXE"]).resolve()
CPU_SLOTS = int(os.environ.get("CPU_SLOTS", "8"))
CPUS_PER_SLOT = int(os.environ.get("CPUS_PER_SLOT", "4"))
RUN_ROOT = ROOT / os.environ.get("RUN_ROOT", "runs_v2_exact_binding")
SLOT_CPUS: tuple[tuple[int, ...], ...] = ()
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OMP_PROC_BIND": "true",
    "OMP_PLACES": "cores",
    "OMP_DYNAMIC": "FALSE",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "MKL_DYNAMIC": "FALSE",
    "BLIS_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
}
MPI_RANK_ENVIRONMENT_KEYS = {
    "OMPI_COMM_WORLD_RANK",
    "PMI_RANK",
    "PMIX_RANK",
    "SLURM_PROCID",
    "MV2_COMM_WORLD_RANK",
}


def parse_ordered_pe_list(value: str) -> tuple[int, ...]:
    fields = value.split(",")
    if not fields or any(not field.strip().isdigit() for field in fields):
        raise ValueError("ORDERED_PE_RESERVATION must be a literal comma-separated list")
    cpus = tuple(int(field.strip()) for field in fields)
    if len(set(cpus)) != len(cpus):
        raise ValueError("ORDERED_PE_RESERVATION contains duplicate logical CPUs")
    return cpus


def parse_linux_cpu_list(value: str) -> set[int]:
    cpus: set[int] = set()
    for component in value.split(","):
        component = component.strip()
        if component.isdigit():
            selected = {int(component)}
        elif "-" in component:
            bounds = component.split("-", 1)
            if len(bounds) != 2 or not all(bound.isdigit() for bound in bounds):
                raise ValueError(f"invalid Linux CPU list {value!r}")
            first, last = (int(bound) for bound in bounds)
            if last < first:
                raise ValueError(f"descending Linux CPU range {component!r}")
            selected = set(range(first, last + 1))
        else:
            raise ValueError(f"invalid Linux CPU list {value!r}")
        if cpus & selected:
            raise ValueError(f"overlapping Linux CPU list {value!r}")
        cpus.update(selected)
    return cpus


def require_single_pu_cores(
    cpus: Sequence[int], topology_root: Path = Path("/sys/devices/system/cpu")
) -> None:
    if not topology_root.is_dir():
        return
    for cpu in cpus:
        path = topology_root / f"cpu{cpu}" / "topology" / "thread_siblings_list"
        try:
            siblings = parse_linux_cpu_list(path.read_text().strip())
        except (OSError, ValueError) as error:
            raise RuntimeError(
                f"cannot prove singleton core topology for CPU {cpu}: {error}"
            ) from error
        if len(siblings) != 1:
            raise RuntimeError(
                f"CPU {cpu} has SMT siblings {sorted(siblings)}; "
                "--bind-to core would not be singleton"
            )


def mpi_control_environment_keys(environment: dict[str, str]) -> list[str]:
    """Remove direct and indirect inherited Open MPI/PRRTE MCA controls."""
    return sorted(
        key
        for key in environment
        if key.startswith(("OMPI_MCA_", "PRTE_MCA_"))
    )


def acquire_cpu_locks(
    cpus: Sequence[int], lock_root: Path | None = None
) -> list[IO[str]]:
    lock_root = lock_root or Path(
        f"/tmp/periodic-gxtb-cpu-reservations-{os.getuid()}"
    )
    lock_root.mkdir(parents=True, exist_ok=True)
    handles: list[IO[str]] = []
    current_handle: IO[str] | None = None
    try:
        for cpu in sorted(set(cpus)):
            current_handle = (lock_root / f"cpu-{cpu}.lock").open("a+")
            try:
                fcntl.flock(
                    current_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB
                )
            except BlockingIOError as error:
                current_handle.seek(0)
                holder = current_handle.read().strip() or "unidentified holder"
                raise RuntimeError(
                    f"logical CPU {cpu} is already reserved ({holder})"
                ) from error
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
                    "source": "run_test_matrix.py",
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


def live_compute_cpu_owners(
    cpus: Sequence[int],
    proc_root: Path = Path("/proc"),
    *,
    ignore_process_identities: dict[int, int] | None = None,
) -> list[dict]:
    selected = set(cpus)
    ignored_identities = dict(ignore_process_identities or {})
    if not selected or not proc_root.is_dir():
        return []
    owners = []
    for entry in proc_root.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid == os.getpid():
            continue
        try:
            initial_stat_state, initial_starttime = linux_proc_stat_identity(
                (entry / "stat").read_text(errors="replace")
            )
            status = (entry / "status").read_text(errors="replace")
        except (FileNotFoundError, ProcessLookupError):
            continue
        except (PermissionError, OSError, ValueError):
            if pid in ignored_identities:
                owners.append({
                    "pid": pid,
                    "name": "",
                    "state": "",
                    "cpus_allowed_list": "",
                    "overlap": sorted(selected),
                    "cp2k_process": False,
                    "mpi_rank_process": False,
                    "process_identity_status": "initial_identity_unreadable",
                })
            continue
        fields = {
            key: value.strip()
            for line in status.splitlines()
            if ":" in line
            for key, value in (line.split(":", 1),)
        }
        state = fields.get("State", "")
        try:
            overlap = selected & parse_linux_cpu_list(
                fields.get("Cpus_allowed_list", "")
            )
        except ValueError:
            overlap = set()
        name = fields.get("Name", "")
        is_cp2k = name.casefold().startswith("cp2k")
        try:
            environment_keys = {
                item.split(b"=", 1)[0].decode(errors="replace")
                for item in (entry / "environ").read_bytes().split(b"\0")
                if item and b"=" in item
            }
        except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
            environment_keys = set()
        is_mpi_rank = bool(environment_keys & MPI_RANK_ENVIRONMENT_KEYS)
        try:
            final_stat_state, final_starttime = linux_proc_stat_identity(
                (entry / "stat").read_text(errors="replace")
            )
        except (FileNotFoundError, ProcessLookupError):
            continue
        except (PermissionError, OSError, ValueError):
            if overlap or pid in ignored_identities:
                owners.append({
                    "pid": pid,
                    "name": name,
                    "state": state,
                    "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
                    "overlap": sorted(selected),
                    "cp2k_process": is_cp2k,
                    "mpi_rank_process": is_mpi_rank,
                    "process_identity_status": "final_identity_unreadable",
                })
            continue
        if final_starttime != initial_starttime:
            owners.append({
                "pid": pid,
                "name": name,
                "state": state,
                "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
                "overlap": sorted(selected),
                "cp2k_process": is_cp2k,
                "mpi_rank_process": is_mpi_rank,
                "process_identity_status": "pid_reused_during_scan",
                "initial_process_starttime": initial_starttime,
                "final_process_starttime": final_starttime,
            })
            continue
        if ignored_identities.get(pid) == initial_starttime:
            continue
        if (
            state.startswith(("Z", "X"))
            or initial_stat_state in {"Z", "X"}
            or final_stat_state in {"Z", "X"}
            or not overlap
            or not (is_cp2k or is_mpi_rank)
        ):
            continue
        owners.append({
            "pid": pid,
            "name": name,
            "state": state,
            "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
            "overlap": sorted(overlap),
            "cp2k_process": is_cp2k,
            "mpi_rank_process": is_mpi_rank,
            "process_starttime": initial_starttime,
            "process_identity_status": "stable",
        })
    return sorted(owners, key=lambda owner: owner["pid"])


def require_no_live_compute_overlap(
    cpus: Sequence[int], proc_root: Path = Path("/proc")
) -> None:
    owners = live_compute_cpu_owners(cpus, proc_root)
    if owners:
        description = "; ".join(
            f"PID {owner['pid']} ({owner['name']}) mask "
            f"{owner['cpus_allowed_list']} overlaps {owner['overlap']}"
            for owner in owners
        )
        raise RuntimeError(
            "selected CPUs are already owned by live CP2K/MPI ranks: "
            + description
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jobs() -> list[tuple[dict, int, str]]:
    result = []
    for case in MATRIX["cases"]:
        for ranks in case["ranks"]:
            for variant in ("DENSE", "STREAMED", "QUALIFY"):
                result.append((case, int(ranks), variant))
    return result


def descendants(root_pid: int) -> set[int]:
    """Return the current Linux /proc descendant set of one launcher."""
    found: set[int] = set()
    pending = [root_pid]
    while pending:
        parent = pending.pop()
        children_path = Path(f"/proc/{parent}/task/{parent}/children")
        try:
            children = [int(value) for value in children_path.read_text().split()]
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        for child in children:
            if child not in found:
                found.add(child)
                pending.append(child)
    return found


def live_process_group_members(process_group: int) -> set[int]:
    """Return non-zombie Linux processes still able to use reserved CPUs."""
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
        except (FileNotFoundError, PermissionError, ProcessLookupError, StopIteration):
            continue
        if state not in {"Z", "X"}:
            members.add(pid)
    return members


def _signal_process_group(process_group: int, sig: signal.Signals) -> None:
    try:
        os.killpg(process_group, sig)
    except ProcessLookupError:
        pass


def terminate_and_reap_process_group(
    process: subprocess.Popen,
    term_timeout: float = 30.0,
    tracked_rank_starttimes: dict[int, int] | None = None,
) -> None:
    """TERM, then KILL, and retain the caller's locks until no rank is live.

    The post-KILL wait intentionally has no timeout.  Releasing CPU reservation
    locks while an uninterruptible rank remains alive would be less safe than
    retaining the launcher process and its locks until the kernel can reap it.
    """
    process_group = process.pid
    tracked = dict(tracked_rank_starttimes or {})

    def live_tracked_ranks() -> set[int]:
        live: set[int] = set()
        for pid, starttime in tracked.items():
            resolution = process_terminal_resolution(pid, starttime)
            if resolution is None or resolution == "identity_unreadable":
                live.add(pid)
        return live

    def signal_tracked_ranks(sig: signal.Signals) -> None:
        for pid in live_tracked_ranks():
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass

    _signal_process_group(process_group, signal.SIGTERM)
    signal_tracked_ranks(signal.SIGTERM)
    deadline = time.monotonic() + term_timeout
    while time.monotonic() < deadline:
        process.poll()
        if process.returncode is not None and not live_process_group_members(
            process_group
        ) and not live_tracked_ranks():
            process.wait()
            return
        time.sleep(0.05)

    _signal_process_group(process_group, signal.SIGKILL)
    signal_tracked_ranks(signal.SIGKILL)
    while True:
        process.poll()
        live_members = live_process_group_members(process_group)
        live_ranks = live_tracked_ranks()
        if process.returncode is not None and not live_members and not live_ranks:
            process.wait()
            return
        if live_members:
            _signal_process_group(process_group, signal.SIGKILL)
        if live_ranks:
            signal_tracked_ranks(signal.SIGKILL)
        time.sleep(0.05)


def linux_proc_stat_identity(stat_text: str) -> tuple[str, int]:
    closing_parenthesis = stat_text.rfind(")")
    if closing_parenthesis < 0:
        raise ValueError("malformed Linux /proc PID stat record")
    fields = stat_text[closing_parenthesis + 1 :].split()
    if len(fields) <= 19:
        raise ValueError("truncated Linux /proc PID stat record")
    return fields[0], int(fields[19])


def linux_process_starttime(
    pid: int, proc_root: Path = Path("/proc")
) -> int | None:
    try:
        _, starttime = linux_proc_stat_identity(
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


def rank_snapshot(pid: int, proc_root: Path = Path("/proc")) -> dict | None:
    root = proc_root / str(pid)
    try:
        stat_state, process_starttime = linux_proc_stat_identity(
            (root / "stat").read_text(errors="replace")
        )
        initial_executable = (root / "exe").resolve(strict=True)
        initial_status = (root / "status").read_text(errors="replace")
        arguments = [
            item.decode(errors="replace")
            for item in (root / "cmdline").read_bytes().split(b"\0")
            if item
        ]
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
    try:
        final_status = (root / "status").read_text(errors="replace")
        final_executable = (root / "exe").resolve(strict=True)
        final_stat_state, final_process_starttime = linux_proc_stat_identity(
            (root / "stat").read_text(errors="replace")
        )
    except FileNotFoundError:
        process_identity_status = "disappeared_after_sample"
        snapshot_consistency_status = "process_disappeared"
        status = initial_status
        executable = initial_executable
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
    environment = {
        key.decode(errors="replace"): value.decode(errors="replace")
        for item in environment_bytes.split(b"\0")
        if item and b"=" in item
        for key, value in (item.split(b"=", 1),)
    }
    allowed_text = next(
        (
            line.split(":", 1)[1].strip()
            for line in status.splitlines()
            if line.startswith("Cpus_allowed_list:")
        ),
        "",
    )
    state = next(
        (
            line.split(":", 1)[1].strip()
            for line in status.splitlines()
            if line.startswith("State:")
        ),
        "",
    )
    rank = None
    rank_value = environment.get("OMPI_COMM_WORLD_RANK")
    if environment_read_status != "available":
        rank_observation_status = f"environment_{environment_read_status}"
    elif rank_value is None:
        rank_observation_status = "explicit_missing"
    else:
        try:
            rank = int(rank_value)
        except ValueError:
            rank_observation_status = "explicit_invalid"
        else:
            rank_observation_status = "explicit"
    return {
        "pid": pid,
        "rank": rank,
        "raw_rank": rank,
        "state": state,
        "stat_state": stat_state,
        "process_starttime": process_starttime,
        "process_identity_status": process_identity_status,
        "snapshot_consistency_status": snapshot_consistency_status,
        "rank_observation_status": rank_observation_status,
        "cpus_allowed_list": allowed_text,
        "executable": str(executable),
        "arguments": arguments,
        "is_cp2k_rank": executable == CP2K,
    }


def accumulate_rank_snapshot(
    previous: dict | None, snapshot: dict, expected_cpus: Sequence[int]
) -> dict:
    accumulated = dict(snapshot)
    sample_count = int(previous.get("sample_count", 0)) + 1 if previous else 1
    if snapshot.get("is_cp2k_rank") is False:
        if previous and previous.get("is_cp2k_rank") is True:
            retained = dict(previous)
            retained.update({
                "sample_count": sample_count,
                "last_observed_executable": snapshot.get("executable"),
                "last_observed_process_identity_status": snapshot.get(
                    "process_identity_status"
                ),
                "executable_changed_ever": True,
                "current_sample_matches_assigned_singleton": False,
                "affinity_violation_ever": True,
            })
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
        return accumulated
    raw_rank = snapshot.get("rank")
    observation_status = snapshot.get("rank_observation_status")
    if not isinstance(observation_status, str):
        observation_status = (
            "explicit"
            if isinstance(raw_rank, int) and not isinstance(raw_rank, bool)
            else "explicit_missing"
        )
    mask_text = str(snapshot.get("cpus_allowed_list", ""))
    try:
        mask = parse_linux_cpu_list(mask_text)
    except ValueError:
        mask = set()
    previous_rank = previous.get("rank") if previous else None
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
    identity_status = snapshot.get("process_identity_status", "stable")
    same_process_identity = bool(
        starttimes_comparable
        and snapshot.get("pid") == previous.get("pid")
        and current_starttime == previous_starttime
        and previous.get("process_starttime_changed_ever") is not True
        and identity_status
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
        and 0 <= int(previous_rank) < len(expected_cpus)
        and mask == {expected_cpus[int(previous_rank)]}
    )
    rank = previous_rank if retain_pending_rank else raw_rank
    accumulated["rank"] = rank
    accumulated["raw_rank"] = raw_rank
    accumulated["rank_identity_source"] = (
        "pending_terminal_environment_loss"
        if retain_pending_rank
        else "explicit_environment"
        if observation_status == "explicit"
        else "unproven"
    )
    sample_matches = (
        isinstance(rank, int)
        and not isinstance(rank, bool)
        and 0 <= rank < len(expected_cpus)
        and mask == {expected_cpus[rank]}
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
        and any(
            observation_status != "explicit"
            or isinstance(prior_rank, int) and prior_rank != raw_rank
            for prior_rank in previous.get("observed_rank_ids", [])
        )
    )
    environment_events = (
        [dict(event) for event in previous.get("rank_environment_events", [])]
        if previous
        else []
    )
    if unavailable_environment:
        environment_events.append({
            "sample_index": sample_count,
            "pid": snapshot.get("pid"),
            "process_starttime": current_starttime,
            "state": state_text,
            "cpus_allowed_list": mask_text,
            "environment_status": observation_status,
            "terminal_resolution": (
                "pending" if retain_pending_rank else "rejected_unproven"
            ),
        })
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
    accumulated.update({
        "sample_count": sample_count,
        "observed_rank_ids": rank_history,
        "observed_cpu_masks": mask_history,
        "observed_rank_observation_statuses": status_history,
        "observed_process_states": state_history,
        "observed_process_starttimes": starttime_history,
        "current_sample_matches_assigned_singleton": sample_matches,
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
        ) + int(unavailable_environment),
        "rank_environment_unavailable_pending": pending_environment_loss,
        "rank_environment_terminally_confirmed": terminally_confirmed,
        "rank_environment_terminal_confirmation": terminal_confirmation,
        "rank_environment_events": environment_events,
        "affinity_violation_ever": bool(
            (previous and previous.get("affinity_violation_ever"))
            or not sample_matches
            or rank_identity_changed
            or process_starttime_changed
            or process_snapshot_inconsistent
            or process_reappeared_after_terminal
            or explicit_rank_reappeared_after_loss
        ),
    })
    return accumulated


def process_terminal_resolution(
    pid: int, expected_starttime: int, proc_root: Path = Path("/proc")
) -> str | None:
    root = proc_root / str(pid)
    if not root.is_dir():
        return "process_disappeared"
    try:
        state, starttime = linux_proc_stat_identity(
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


def resolve_pending_rank_environment(record: dict, resolution: str) -> None:
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
    resolve_rank_process_lifetime(record, resolution)
    if not terminal:
        record["rank_identity_changed_ever"] = True


def resolve_rank_process_lifetime(record: dict, resolution: str) -> None:
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


def observed_rank_process_is_still_live(
    pid: int, record: dict, proc_root: Path = Path("/proc")
) -> bool:
    starttime = record.get("process_starttime")
    if isinstance(starttime, int) and not isinstance(starttime, bool):
        resolution = process_terminal_resolution(pid, starttime, proc_root)
    else:
        resolution = "identity_unreadable"
    if resolution is None:
        record["snapshot_unavailable_ever"] = True
        record["current_sample_matches_assigned_singleton"] = False
        record["affinity_violation_ever"] = True
        return True
    if record.get("rank_environment_unavailable_pending"):
        resolve_pending_rank_environment(record, resolution)
    else:
        resolve_rank_process_lifetime(record, resolution)
    return False


def concurrent_rank_pid_groups(current: Sequence[dict]) -> dict[int, list[int]]:
    rank_pids: dict[int, set[int]] = {}
    for snapshot in current:
        if str(snapshot.get("state", "")).startswith(("Z", "X")):
            continue
        rank = snapshot.get("rank")
        if isinstance(rank, int) and not isinstance(rank, bool):
            rank_pids.setdefault(rank, set()).add(int(snapshot["pid"]))
    return {
        rank: sorted(pids)
        for rank, pids in sorted(rank_pids.items())
        if len(pids) > 1
    }


def runtime_external_owners(
    launcher_pid: int,
    launcher_starttime: int,
    observed: dict[int, dict],
    expected_cpus: Sequence[int],
) -> list[dict]:
    own_identities = {launcher_pid: launcher_starttime}
    for pid, item in observed.items():
        starttime = item.get("process_starttime")
        if (
            item.get("is_cp2k_rank") is True
            and item.get("process_terminally_confirmed") is not True
            and isinstance(starttime, int)
            and not isinstance(starttime, bool)
            and process_terminal_resolution(pid, starttime) is None
        ):
            own_identities[pid] = starttime
    return live_compute_cpu_owners(
        expected_cpus, ignore_process_identities=own_identities
    )


def tracked_rank_starttimes(observed: dict[int, dict]) -> dict[int, int]:
    identities: dict[int, int] = {}
    for pid, item in observed.items():
        starttime = item.get("process_starttime")
        if (
            item.get("is_cp2k_rank") is True
            and isinstance(starttime, int)
            and not isinstance(starttime, bool)
        ):
            identities[pid] = starttime
    return identities


def sample_rank_affinity(
    launcher_pid: int, observed: dict[int, dict], expected_cpus: Sequence[int]
) -> dict[int, list[int]]:
    current = []
    tracked_rank_pids = {
        pid
        for pid, item in observed.items()
        if item.get("is_cp2k_rank") is True
        and item.get("process_terminally_confirmed") is not True
    }
    for pid in sorted(descendants(launcher_pid) | tracked_rank_pids):
        snapshot = rank_snapshot(pid)
        if snapshot is None:
            tracked = observed.get(pid)
            if tracked and tracked.get("is_cp2k_rank") is True:
                if observed_rank_process_is_still_live(pid, tracked):
                    current.append(tracked)
            continue
        if snapshot.get("is_cp2k_rank") is not True and pid not in observed:
            continue
        observed[pid] = accumulate_rank_snapshot(
            observed.get(pid), snapshot, expected_cpus
        )
        current.append(observed[pid])
    return concurrent_rank_pid_groups(current)


def ordered_rank_proof(
    observed: dict[int, dict],
    ranks: int,
    expected_cpus: Sequence[int],
    concurrent_duplicate_rank_ids: set[int] | None = None,
    *,
    final: bool = True,
) -> list[dict]:
    groups: dict[int, list[dict]] = {}
    for item in observed.values():
        rank = item.get("rank")
        if isinstance(rank, int):
            groups.setdefault(rank, []).append(item)
    if sorted(groups) != list(range(ranks)):
        raise RuntimeError(f"could not prove exactly ranks 0..{ranks - 1}: {observed}")
    duplicates = concurrent_duplicate_rank_ids or set()
    ordered = []
    for rank in range(ranks):
        generations = sorted(groups[rank], key=lambda item: int(item["pid"]))
        mask_history = sorted(
            {
                str(mask)
                for item in generations
                for mask in item.get("observed_cpu_masks", [])
            }
        )
        exact = (
            len(generations) == 1
            and mask_history == [str(expected_cpus[rank])]
            and rank not in duplicates
            and all(
                item.get("affinity_violation_ever") is False
                and item.get("rank_identity_changed_ever") is not True
                and item.get("process_starttime_changed_ever") is not True
                and (
                    not final
                    or item.get("process_terminally_confirmed") is True
                )
                and item.get("process_reappeared_after_terminal_ever") is not True
                and item.get("executable_changed_ever") is not True
                and item.get("cpu_mask_changed_during_sample_ever") is not True
                and item.get("snapshot_unavailable_ever") is not True
                and (
                    not final
                    or item.get("rank_environment_unavailable_pending") is not True
                )
                and (
                    not final
                    or item.get("rank_environment_unavailable_ever") is not True
                    or item.get("rank_environment_terminally_confirmed") is True
                )
                and item.get("current_sample_matches_assigned_singleton") is True
                for item in generations
            )
        )
        canonical = max(
            generations,
            key=lambda item: (int(item.get("sample_count", 0)), -int(item["pid"])),
        )
        ordered.append(
            {
                "rank": rank,
                "pid": int(canonical["pid"]),
                "pid_generations": [int(item["pid"]) for item in generations],
                "cpus_allowed_list": (
                    mask_history[0] if len(mask_history) == 1 else ""
                ),
                "observed_cpu_masks": mask_history,
                "current_sample_matches_assigned_singleton": exact,
                "affinity_violation_ever": not exact,
                "concurrent_duplicate_pid_ever": rank in duplicates,
            }
        )
    if any(item["affinity_violation_ever"] for item in ordered):
        raise RuntimeError(f"rank affinity did not remain singleton-exact: {ordered}")
    return ordered


def prove_initial_rank_affinity(
    process: subprocess.Popen,
    launcher_starttime: int,
    observed: dict[int, dict],
    ranks: int,
    expected_cpus: Sequence[int],
    runtime_overlap_samples: list[dict],
) -> tuple[dict[int, dict], set[int], list[dict], int]:
    duplicate_rank_ids: set[int] = set()
    duplicate_rank_samples: list[dict] = []
    sample_index = 0
    deadline = time.monotonic() + 30.0
    while time.monotonic() < deadline and process.poll() is None:
        sample_index += 1
        duplicate_groups = sample_rank_affinity(
            process.pid, observed, expected_cpus
        )
        external_owners = runtime_external_owners(
            process.pid, launcher_starttime, observed, expected_cpus
        )
        if external_owners:
            runtime_overlap_samples.append({
                "sample_index": sample_index,
                "owners": external_owners,
            })
            raise RuntimeError(
                f"external live CP2K/MPI CPU overlap: {external_owners}"
            )
        duplicate_rank_ids.update(duplicate_groups)
        if duplicate_groups:
            duplicate_rank_samples.append({
                "sample_index": sample_index,
                "rank_pid_groups": [
                    {"rank": rank, "pids": pids}
                    for rank, pids in duplicate_groups.items()
                ],
            })
        if any(item.get("affinity_violation_ever") for item in observed.values()):
            raise RuntimeError(f"rank affinity mismatch before initial proof: {observed}")
        try:
            ordered_rank_proof(
                observed,
                ranks,
                expected_cpus,
                duplicate_rank_ids,
                final=False,
            )
            return (
                observed,
                duplicate_rank_ids,
                duplicate_rank_samples,
                sample_index,
            )
        except RuntimeError:
            time.sleep(0.01)
    raise RuntimeError(f"could not prove initial singleton affinity for {ranks} ranks")


def reported_binding_rank_ids(text: str) -> list[int]:
    return sorted(
        {
            int(value)
            for value in re.findall(
                r"\b(?:MCW\s+)?rank\s+(\d+)\s+bound\b",
                text,
                flags=re.IGNORECASE,
            )
        }
    )


def _run_one_inner(
    job: tuple[dict, int, str], slot: int, lifecycle: dict
) -> str:
    case, ranks, variant = job
    run_id = f"{case['name']}_p{ranks}_{variant.lower()}"
    run_dir = RUN_ROOT / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"refusing to reuse nonempty run directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = (INPUT_ROOT / case["input"]).resolve()
    if not input_path.is_file():
        raise RuntimeError(f"missing input: {input_path}")
    if ranks > CPUS_PER_SLOT:
        raise RuntimeError(f"rank count exceeds fixed CPU slot: {run_id}")
    if not SLOT_CPUS or slot >= len(SLOT_CPUS):
        raise RuntimeError("ordered PE slots were not initialized")
    pe_cpus = SLOT_CPUS[slot][:ranks]
    pe_list = ",".join(str(cpu) for cpu in pe_cpus)
    require_no_live_compute_overlap(pe_cpus)

    env = os.environ.copy()
    env.update(THREAD_ENV)
    removed_binding_environment = mpi_control_environment_keys(env)
    for key in removed_binding_environment:
        env.pop(key, None)
    env["CP2K_GXTB_SYMMETRY_STAR_CONTRACTION"] = variant
    env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_ROOT"
    env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = "3"
    env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
    command = [
        str(MPIEXEC), "--map-by", f"pe-list={pe_list}:ordered",
        "--bind-to", "core", "--report-bindings", "-np", str(ranks),
        str(CP2K), "-i", str(input_path),
    ]
    metadata = {
        "schema_version": 2,
        "run_id": run_id,
        "case": case["name"],
        "features": case["features"],
        "expected_nfull": case["nfull"],
        "ranks": ranks,
        "variant": variant,
        "ordered_pe_list": pe_list,
        "exact_cpus_per_rank": 1,
        "mpi_map_by": "pe-list=<ordered-list>:ordered",
        "mpi_bind_to": "core",
        "mpi_report_bindings": True,
        "outer_taskset": False,
        "cross_process_cpu_reservation_gate": True,
        "live_compute_overlap_preflight_gate": True,
        "live_compute_overlap_preflight_owners": [],
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "working_directory": str(run_dir.resolve()),
        "cp2k": str(CP2K),
        "cp2k_sha256": sha256(CP2K),
        "cp2k_lib": str(CP2K_LIB),
        "cp2k_lib_sha256": sha256(CP2K_LIB),
        "mpi_launcher": str(MPIEXEC),
        "mpi_launcher_sha256": sha256(MPIEXEC),
        "command": command,
        "removed_mpi_binding_environment_keys": removed_binding_environment,
        "environment": {key: env[key] for key in sorted(set(THREAD_ENV) | {
            "CP2K_GXTB_SYMMETRY_STAR_CONTRACTION",
            "CP2K_GXTB_EXCHANGE_STREAM_MODE",
            "CP2K_GXTB_EXCHANGE_GRADIENT_MODE",
            "CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE",
            "CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION",
        })},
        "started_unix": time.time(),
        "timing_classification": "timing_pending_full_revalidation",
    }
    lifecycle["run_dir"] = run_dir
    lifecycle["runtime_overlap_samples"] = []
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    output_path = run_dir / "cp2k.out"
    launcher_log_path = run_dir / "launcher.log"
    started = time.perf_counter()
    with output_path.open("wb") as stdout, launcher_log_path.open("wb") as launcher_log:
        process = subprocess.Popen(
            command,
            cwd=run_dir,
            env=env,
            stdout=stdout,
            stderr=launcher_log,
            start_new_session=True,
        )
        lifecycle["process"] = process
        launcher_starttime = linux_process_starttime(process.pid)
        if Path("/proc").is_dir() and launcher_starttime is None:
            raise RuntimeError(
                "could not prove the MPI launcher PID/starttime identity"
            )
        lifecycle["launcher_starttime"] = launcher_starttime
        runtime_overlap_samples = lifecycle["runtime_overlap_samples"]
        observed: dict[int, dict] = {}
        lifecycle["observed"] = observed
        (
            observed,
            duplicate_rank_ids,
            duplicate_rank_samples,
            affinity_sample_index,
        ) = prove_initial_rank_affinity(
            process,
            int(launcher_starttime or 0),
            observed,
            ranks,
            pe_cpus,
            runtime_overlap_samples,
        )
        while process.poll() is None:
            affinity_sample_index += 1
            duplicate_groups = sample_rank_affinity(
                process.pid, observed, pe_cpus
            )
            external_owners = runtime_external_owners(
                process.pid,
                int(launcher_starttime or 0),
                observed,
                pe_cpus,
            )
            if external_owners:
                runtime_overlap_samples.append({
                    "sample_index": affinity_sample_index,
                    "owners": external_owners,
                })
                raise RuntimeError(
                    f"external live CP2K/MPI CPU overlap: {external_owners}"
                )
            duplicate_rank_ids.update(duplicate_groups)
            if duplicate_groups:
                duplicate_rank_samples.append({
                    "sample_index": affinity_sample_index,
                    "rank_pid_groups": [
                        {"rank": rank, "pids": pids}
                        for rank, pids in duplicate_groups.items()
                    ],
                })
            ordered_rank_proof(
                observed,
                ranks,
                pe_cpus,
                duplicate_rank_ids,
                final=False,
            )
            time.sleep(0.01)
        returncode = process.wait()
    wall = time.perf_counter() - started
    for pid, item in observed.items():
        if item.get("process_terminally_confirmed") is True:
            continue
        starttime = item.get("process_starttime")
        if isinstance(starttime, int) and not isinstance(starttime, bool):
            resolution = process_terminal_resolution(pid, starttime)
        else:
            resolution = "identity_unreadable"
        resolution = resolution or "launcher_ended_while_process_live"
        if item.get("rank_environment_unavailable_pending"):
            resolve_pending_rank_environment(item, resolution)
        else:
            resolve_rank_process_lifetime(item, resolution)
    live_after_launcher = {
        pid
        for pid, item in observed.items()
        if item.get("process_terminally_confirmed") is not True
    }
    live_group_members = live_process_group_members(process.pid)
    if live_after_launcher or live_group_members:
        terminate_and_reap_process_group(
            process,
            tracked_rank_starttimes=tracked_rank_starttimes(observed),
        )
        returncode = 97
        for pid, item in observed.items():
            if item.get("process_terminally_confirmed") is True:
                continue
            starttime = item.get("process_starttime")
            resolution = (
                process_terminal_resolution(pid, starttime)
                if isinstance(starttime, int) and not isinstance(starttime, bool)
                else "identity_unreadable"
            )
            resolution = resolution or "process_still_live_after_group_cleanup"
            if item.get("rank_environment_unavailable_pending"):
                resolve_pending_rank_environment(item, resolution)
            else:
                resolve_rank_process_lifetime(item, resolution)
    affinity = ordered_rank_proof(
        observed, ranks, pe_cpus, duplicate_rank_ids, final=True
    )
    launcher_text = launcher_log_path.read_text(errors="replace")
    binding_ranks = reported_binding_rank_ids(launcher_text)
    runtime_affinity_gate = (
        binding_ranks == list(range(ranks))
        and len(affinity) == ranks
        and all(item.get("affinity_violation_ever") is False for item in affinity)
        and not lifecycle["runtime_overlap_samples"]
    )
    if returncode == 0 and not runtime_affinity_gate:
        returncode = 97
    metadata.update({
        "finished_unix": time.time(),
        "wall_seconds": wall,
        "returncode": returncode,
        "affinity_proof": affinity,
        "observed_cp2k_rank_pid_generations": [
            item["pid_generations"] for item in affinity
        ],
        "observed_cp2k_process_generation_count": sum(
            len(item["pid_generations"]) for item in affinity
        ),
        "observed_child_processes": sorted(
            observed.values(), key=lambda item: int(item["pid"])
        ),
        "concurrent_duplicate_rank_ids_ever": sorted(duplicate_rank_ids),
        "concurrent_duplicate_rank_samples": duplicate_rank_samples,
        "concurrent_duplicate_rank_processes_ever": bool(duplicate_rank_ids),
        "all_observed_rank_samples_match_ordered_pe_list": all(
            item.get("affinity_violation_ever") is False for item in affinity
        ),
        "reported_binding_rank_ids": binding_ranks,
        "launcher_log": launcher_log_path.name,
        "launcher_log_sha256": sha256(launcher_log_path),
        "runtime_affinity_gate": runtime_affinity_gate,
        "live_compute_overlap_runtime_gate": not lifecycle[
            "runtime_overlap_samples"
        ],
        "live_compute_overlap_runtime_samples": lifecycle[
            "runtime_overlap_samples"
        ],
        "timing_classification": (
            "production_scaling_eligible"
            if returncode == 0 and runtime_affinity_gate
            else "timing_non_scaling"
        ),
        "output_sha256": sha256(output_path),
    })
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "returncode.txt").write_text(f"{returncode}\n")
    text = output_path.read_text(errors="replace")
    if returncode != 0:
        raise RuntimeError(f"{run_id}: return code {returncode}")
    if text.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"{run_id}: CP2K did not terminate exactly once")
    if variant == "STREAMED" and "GXTB-MIXER-STAR-STREAMED" not in text:
        raise RuntimeError(f"{run_id}: missing streamed selector marker")
    if variant == "QUALIFY" and "GXTB-QUALIFICATION_ONLY MIXER-STAR" not in text:
        raise RuntimeError(f"{run_id}: missing qualification marker")
    return f"COMPLETE\t{run_id}\t{wall:.3f}s\tPE={pe_list}"


def run_one(job: tuple[dict, int, str], slot: int) -> str:
    """Run one isolated case and drain every spawned task on any failure."""
    lifecycle: dict = {}
    try:
        return _run_one_inner(job, slot, lifecycle)
    except BaseException as error:
        process = lifecycle.get("process")
        if process is not None:
            terminate_and_reap_process_group(
                process,
                tracked_rank_starttimes=tracked_rank_starttimes(
                    lifecycle.get("observed", {})
                ),
            )
        run_dir = lifecycle.get("run_dir")
        if isinstance(run_dir, Path):
            record_path = run_dir / "run.json"
            try:
                record = json.loads(record_path.read_text())
            except (FileNotFoundError, OSError, json.JSONDecodeError):
                record = {}
            samples = lifecycle.get("runtime_overlap_samples", [])
            record.update({
                "returncode": 97,
                "live_compute_overlap_runtime_gate": not samples,
                "live_compute_overlap_runtime_samples": samples,
                "timing_classification": "timing_non_scaling",
                "monitor_failure": f"{type(error).__name__}: {error}",
            })
            record_path.write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
        raise


def main() -> int:
    global SLOT_CPUS
    if not CP2K.is_file():
        raise RuntimeError(f"missing CP2K executable: {CP2K}")
    if not CP2K_LIB.is_file():
        raise RuntimeError(f"missing CP2K shared library: {CP2K_LIB}")
    if not MPIEXEC.is_file():
        raise RuntimeError(f"missing MPI launcher: {MPIEXEC}")
    reservation_text = os.environ.get("ORDERED_PE_RESERVATION")
    if reservation_text is None:
        raise RuntimeError(
            "ORDERED_PE_RESERVATION is required as an explicit literal CPU list"
        )
    reservation = parse_ordered_pe_list(reservation_text)
    required = CPU_SLOTS * CPUS_PER_SLOT
    if len(reservation) != required:
        raise RuntimeError(
            f"ORDERED_PE_RESERVATION has {len(reservation)} CPUs; expected {required}"
        )
    if hasattr(os, "sched_getaffinity"):
        unavailable = set(reservation) - set(os.sched_getaffinity(0))
        if unavailable:
            raise RuntimeError(f"reservation contains unavailable CPUs: {sorted(unavailable)}")
    require_single_pu_cores(reservation)
    SLOT_CPUS = tuple(
        tuple(reservation[index:index + CPUS_PER_SLOT])
        for index in range(0, required, CPUS_PER_SLOT)
    )
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    writer_lock = (RUN_ROOT / ".writer.lock").open("a+")
    try:
        fcntl.flock(writer_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        writer_lock.close()
        raise RuntimeError(f"another writer owns {RUN_ROOT}") from error
    cpu_locks = acquire_cpu_locks(reservation)
    try:
        require_no_live_compute_overlap(reservation)
        all_jobs = jobs()
        if len(all_jobs) != 48:
            raise RuntimeError(
                f"fixed rerun matrix changed unexpectedly: {len(all_jobs)} != 48"
            )
        slots = list(range(CPU_SLOTS))
        with concurrent.futures.ThreadPoolExecutor(max_workers=CPU_SLOTS) as pool:
            running: dict[concurrent.futures.Future[str], int] = {}
            next_job = 0
            while next_job < len(all_jobs) or running:
                while slots and next_job < len(all_jobs):
                    slot = slots.pop(0)
                    future = pool.submit(run_one, all_jobs[next_job], slot)
                    running[future] = slot
                    next_job += 1
                done, _ = concurrent.futures.wait(
                    running, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in done:
                    slot = running.pop(future)
                    print(future.result(), flush=True)
                    slots.append(slot)
                    slots.sort()
    finally:
        for handle in cpu_locks:
            handle.close()
        writer_lock.close()
    print(f"ALL_COMPLETE\t{len(all_jobs)}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
