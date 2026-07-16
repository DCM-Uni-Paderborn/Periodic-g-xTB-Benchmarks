#!/usr/bin/env python3
"""Run one command and record wall time plus child rusage on Linux."""

from __future__ import annotations

import argparse
import resource
import subprocess
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("missing command after --")

    started = time.perf_counter()
    completed = subprocess.run(command, check=False)
    wall_seconds = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    args.metrics.write_text(
        f"wall_seconds={wall_seconds:.9f}\n"
        f"child_user_seconds={usage.ru_utime:.9f}\n"
        f"child_system_seconds={usage.ru_stime:.9f}\n"
        f"child_maxrss_kb={usage.ru_maxrss}\n"
        f"returncode={completed.returncode}\n"
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
