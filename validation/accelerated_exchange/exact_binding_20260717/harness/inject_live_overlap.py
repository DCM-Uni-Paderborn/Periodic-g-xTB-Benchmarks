#!/usr/bin/env python3
"""Inject one live CP2K-named process and prove the overlap preflight blocks."""

from __future__ import annotations

import ctypes
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("CP2K_EXE", "/usr/bin/true")
os.environ.setdefault("CP2K_LIB", "/usr/bin/true")
os.environ.setdefault("MPIEXEC_EXE", "/usr/bin/true")

import run_test_matrix as runner


def child(cpu: int) -> int:
    os.sched_setaffinity(0, {cpu})
    libc = ctypes.CDLL(None)
    if libc.prctl(15, ctypes.c_char_p(b"cp2k.inject"), 0, 0, 0) != 0:
        raise OSError("PR_SET_NAME failed")
    print(json.dumps({"pid": os.getpid(), "cpu": cpu}), flush=True)
    time.sleep(60.0)
    return 0


def parent(cpu: int) -> int:
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--child", str(cpu)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        assert process.stdout is not None
        ready = json.loads(process.stdout.readline())
        try:
            runner.require_no_live_compute_overlap((cpu,))
        except RuntimeError as error:
            print(
                json.dumps(
                    {
                        "status": "PASS",
                        "injected_pid": ready["pid"],
                        "selected_cpu": cpu,
                        "preflight_error": str(error),
                    },
                    sort_keys=True,
                )
            )
            return 0
        print("FAIL: live overlap was accepted", file=sys.stderr)
        return 97
    finally:
        process.terminate()
        process.wait(timeout=5.0)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--child":
        raise SystemExit(child(int(sys.argv[2])))
    selected = int(sys.argv[1]) if len(sys.argv) == 2 else min(os.sched_getaffinity(0))
    raise SystemExit(parent(selected))
