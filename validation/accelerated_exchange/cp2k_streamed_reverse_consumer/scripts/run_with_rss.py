#!/usr/bin/env python3
"""Run a command, sample RSS, and optionally prove MPI rank affinity."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import subprocess
import time
from pathlib import Path


def children(pid: int) -> list[int]:
    try:
        raw = Path(f"/proc/{pid}/task/{pid}/children").read_text().strip()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return []
    return [int(item) for item in raw.split()] if raw else []


def process_tree(root: int) -> set[int]:
    found: set[int] = set()
    pending = [root]
    while pending:
        pid = pending.pop()
        if pid in found:
            continue
        found.add(pid)
        pending.extend(children(pid))
    return found


def rss_kib(pid: int) -> int:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
        pass
    return 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ordered_pe_list(value: str) -> tuple[int, ...]:
    fields = value.split(",")
    if not fields or any(not field.strip().isdigit() for field in fields):
        raise ValueError("ordered PE list must be a literal comma-separated CPU list")
    cpus = tuple(int(field.strip()) for field in fields)
    if len(set(cpus)) != len(cpus):
        raise ValueError("ordered PE list contains duplicate CPUs")
    return cpus


def parse_linux_cpu_list(value: str) -> set[int]:
    cpus: set[int] = set()
    for field in value.split(","):
        if "-" in field:
            bounds = field.split("-", 1)
            if len(bounds) != 2 or not all(bound.isdigit() for bound in bounds):
                raise ValueError("invalid Linux CPU list")
            first, last = (int(bound) for bound in bounds)
            if last < first:
                raise ValueError("invalid descending Linux CPU range")
            cpus.update(range(first, last + 1))
        elif field.isdigit():
            cpus.add(int(field))
        else:
            raise ValueError("invalid Linux CPU list")
    return cpus


def rank_snapshot(pid: int, cp2k: Path) -> dict[str, object] | None:
    root = Path("/proc") / str(pid)
    try:
        if (root / "exe").resolve(strict=True) != cp2k:
            return None
        status = (root / "status").read_text(errors="replace")
        environment = {
            key.decode(errors="replace"): value.decode(errors="replace")
            for item in (root / "environ").read_bytes().split(b"\0")
            if item and b"=" in item
            for key, value in (item.split(b"=", 1),)
        }
    except (FileNotFoundError, PermissionError, ProcessLookupError, OSError):
        return None
    fields = {
        key: value.strip()
        for line in status.splitlines()
        if ":" in line
        for key, value in (line.split(":", 1),)
    }
    try:
        rank = int(environment["OMPI_COMM_WORLD_RANK"])
    except (KeyError, ValueError):
        rank = None
    return {
        "pid": pid,
        "rank": rank,
        "cpus_allowed_list": fields.get("Cpus_allowed_list", ""),
    }


def parse_singleton_mask(value: str) -> int | None:
    if value.isdigit():
        return int(value)
    return None


def accumulate_rank_snapshot(
    previous: dict[str, object] | None,
    snapshot: dict[str, object],
    expected_cpus: tuple[int, ...],
) -> dict[str, object]:
    """Keep a sticky failure bit across every observed rank-affinity sample."""
    accumulated = dict(snapshot)
    rank = snapshot.get("rank")
    mask_text = str(snapshot.get("cpus_allowed_list", ""))
    singleton = parse_singleton_mask(mask_text)
    sample_matches = (
        isinstance(rank, int)
        and 0 <= rank < len(expected_cpus)
        and singleton == expected_cpus[rank]
    )
    rank_history = list(previous.get("observed_rank_ids", [])) if previous else []
    mask_history = list(previous.get("observed_cpu_masks", [])) if previous else []
    if rank not in rank_history:
        rank_history.append(rank)
    if mask_text not in mask_history:
        mask_history.append(mask_text)
    accumulated.update(
        {
            "sample_count": (
                int(previous.get("sample_count", 0)) + 1 if previous else 1
            ),
            "observed_rank_ids": rank_history,
            "observed_cpu_masks": mask_history,
            "current_sample_matches_assigned_singleton": sample_matches,
            "affinity_violation_ever": bool(
                (previous and previous.get("affinity_violation_ever"))
                or not sample_matches
            ),
        }
    )
    return accumulated


def reported_binding_ranks(text: str) -> list[int]:
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


def mpi_control_environment_keys(environment: dict[str, str]) -> list[str]:
    """Remove direct and indirect inherited Open MPI/PRRTE MCA controls."""
    return sorted(
        key
        for key in environment
        if key.startswith(("OMPI_MCA_", "PRTE_MCA_"))
    )


def require_single_pu_cores(
    cpus: tuple[int, ...],
    topology_root: Path = Path("/sys/devices/system/cpu"),
) -> None:
    if not topology_root.is_dir():
        return
    for cpu in cpus:
        path = topology_root / f"cpu{cpu}" / "topology" / "thread_siblings_list"
        try:
            siblings = parse_linux_cpu_list(path.read_text().strip())
        except (OSError, ValueError) as error:
            raise ValueError(
                f"cannot prove singleton core topology for CPU {cpu}: {error}"
            ) from error
        if len(siblings) != 1:
            raise ValueError(
                f"CPU {cpu} belongs to SMT siblings {siblings}; "
                "--bind-to core is not singleton"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpi-ranks", type=int)
    parser.add_argument("--ordered-pe-list")
    parser.add_argument("--cp2k", type=Path)
    parser.add_argument("--launcher-log", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("COMMAND is required")
    affinity_values = (
        args.mpi_ranks,
        args.ordered_pe_list,
        args.cp2k,
        args.launcher_log,
    )
    affinity_requested = any(value is not None for value in affinity_values)
    if affinity_requested and any(value is None for value in affinity_values):
        parser.error("all MPI affinity proof options must be supplied together")
    expected_cpus: tuple[int, ...] = ()
    cp2k: Path | None = None
    if affinity_requested:
        assert args.mpi_ranks is not None
        assert args.ordered_pe_list is not None
        assert args.cp2k is not None
        if args.mpi_ranks < 1:
            parser.error("--mpi-ranks must be positive")
        try:
            expected_cpus = parse_ordered_pe_list(args.ordered_pe_list)
        except ValueError as error:
            parser.error(str(error))
        if len(expected_cpus) != args.mpi_ranks:
            parser.error("ordered PE-list length must equal --mpi-ranks")
        if hasattr(os, "sched_getaffinity"):
            unavailable = set(expected_cpus) - set(os.sched_getaffinity(0))
            if unavailable:
                parser.error(f"ordered PE list contains unavailable CPUs: {sorted(unavailable)}")
        try:
            require_single_pu_cores(expected_cpus)
        except ValueError as error:
            parser.error(str(error))
        cp2k = args.cp2k.resolve(strict=True)

    result_path = args.result
    command = args.command
    launch_environment = os.environ.copy()
    removed_binding_environment: list[str] = []
    if affinity_requested:
        removed_binding_environment = mpi_control_environment_keys(launch_environment)
        for key in removed_binding_environment:
            launch_environment.pop(key, None)
    started_wall = time.time()
    started_mono = time.monotonic()
    proc = subprocess.Popen(command, env=launch_environment)

    peak_tree_rss_kib = 0
    peak_single_process_rss_kib = 0
    peak_single_process_pid = proc.pid
    max_process_count = 0
    samples = 0
    observed_ranks: dict[int, dict[str, object]] = {}
    while proc.poll() is None:
        pids = process_tree(proc.pid)
        values = [(pid, rss_kib(pid)) for pid in pids]
        tree_rss = sum(value for _, value in values)
        if tree_rss > peak_tree_rss_kib:
            peak_tree_rss_kib = tree_rss
        if values:
            pid, value = max(values, key=lambda item: item[1])
            if value > peak_single_process_rss_kib:
                peak_single_process_rss_kib = value
                peak_single_process_pid = pid
        max_process_count = max(max_process_count, len(pids))
        if cp2k is not None:
            for pid in pids:
                snapshot = rank_snapshot(pid, cp2k)
                if snapshot is not None:
                    observed_ranks[pid] = accumulate_rank_snapshot(
                        observed_ranks.get(pid), snapshot, expected_cpus
                    )
        samples += 1
        time.sleep(0.02)

    returncode = proc.wait()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    affinity: dict[str, object] = {
        "requested": affinity_requested,
        "runtime_affinity_gate": None,
        "timing_classification": "not_an_mpi_affinity_run",
    }
    if affinity_requested:
        assert args.mpi_ranks is not None
        assert args.ordered_pe_list is not None
        assert args.launcher_log is not None
        ordered = sorted(
            observed_ranks.values(),
            key=lambda item: (
                int(item["rank"])
                if isinstance(item.get("rank"), int)
                else args.mpi_ranks,
                int(item["pid"]),
            ),
        )
        rank_ids = [item.get("rank") for item in ordered]
        masks = [str(item.get("cpus_allowed_list", "")) for item in ordered]
        singleton_cpus = [parse_singleton_mask(mask) for mask in masks]
        all_rank_samples_exact = len(ordered) == args.mpi_ranks and all(
            item.get("affinity_violation_ever") is False
            and item.get("current_sample_matches_assigned_singleton") is True
            for item in ordered
        )
        launcher_log = args.launcher_log.resolve()
        launcher_text = (
            launcher_log.read_text(errors="replace")
            if launcher_log.is_file()
            else ""
        )
        report_ranks = reported_binding_ranks(launcher_text)
        gate = (
            len(ordered) == args.mpi_ranks
            and rank_ids == list(range(args.mpi_ranks))
            and singleton_cpus == list(expected_cpus)
            and all_rank_samples_exact
            and report_ranks == list(range(args.mpi_ranks))
        )
        affinity = {
            "requested": True,
            "mpi_ranks": args.mpi_ranks,
            "ordered_pe_list": args.ordered_pe_list,
            "observed_rank_pids": [int(item["pid"]) for item in ordered],
            "observed_rank_ids": rank_ids,
            "observed_rank_cpu_masks": masks,
            "observed_rank_histories": ordered,
            "all_observed_rank_samples_match_ordered_pe_list": (
                all_rank_samples_exact
            ),
            "reported_binding_rank_ids": report_ranks,
            "removed_mpi_binding_environment_keys": removed_binding_environment,
            "launcher_log": str(launcher_log),
            "launcher_log_sha256": (
                sha256(launcher_log) if launcher_log.is_file() else None
            ),
            "runtime_affinity_gate": gate,
            "timing_classification": (
                "production_scaling_eligible"
                if gate and returncode == 0
                else "timing_non_scaling"
            ),
        }
        if returncode == 0 and not gate:
            returncode = 97
    result = {
        "affinity": affinity,
        "command": command,
        "elapsed_seconds": time.monotonic() - started_mono,
        "finished_unix_seconds": time.time(),
        "hostname": os.uname().nodename,
        "max_process_count": max_process_count,
        "peak_single_process_pid": peak_single_process_pid,
        "peak_single_process_rss_kib": peak_single_process_rss_kib,
        "peak_tree_rss_kib": peak_tree_rss_kib,
        "resource_children_maxrss_kib": usage.ru_maxrss,
        "returncode": returncode,
        "samples": samples,
        "started_unix_seconds": started_wall,
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
