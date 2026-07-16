#!/usr/bin/env python3
"""Run a command and sample the resident set of its complete Linux process tree."""

from __future__ import annotations

import json
import os
import resource
import subprocess
import sys
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


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit(f"usage: {sys.argv[0]} RESULT.json COMMAND [ARG ...]")

    result_path = Path(sys.argv[1])
    command = sys.argv[2:]
    started_wall = time.time()
    started_mono = time.monotonic()
    proc = subprocess.Popen(command)

    peak_tree_rss_kib = 0
    peak_single_process_rss_kib = 0
    peak_single_process_pid = proc.pid
    max_process_count = 0
    samples = 0
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
        samples += 1
        time.sleep(0.02)

    returncode = proc.wait()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {
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
