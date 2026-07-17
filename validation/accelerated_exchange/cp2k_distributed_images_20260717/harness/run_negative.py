#!/usr/bin/env python3
"""Run fail-closed qualification faults for distributed-image importers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import signal
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = json.loads((ROOT / "negative_cases.json").read_text())["cases"]
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "MKL_DYNAMIC": "FALSE",
    "BLIS_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_case(args: argparse.Namespace, case: dict) -> None:
    run_dir = ROOT / "negative_runs" / case["name"]
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"refusing to reuse nonempty negative run: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    source = (ROOT / "inputs" / case["input"]).resolve()
    if not source.is_file():
        raise RuntimeError(f"missing negative input: {source}")
    env = os.environ.copy()
    env.update(THREAD_ENV)
    env.update({
        "CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE": "3",
        "CP2K_GXTB_EXCHANGE_STREAM_MODE": "KGROUP_PARTIAL_DISTRIBUTED_IMAGES",
        "CP2K_GXTB_EXCHANGE_GRADIENT_MODE": "QUALIFY",
        "CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION": "1",
        "CP2K_GXTB_PARTIAL_QUALIFY_INJECT": case["injection"],
    })
    command = [args.mpiexec, *args.mpiexec_arg, "-np", str(case["ranks"])]
    if args.rank_prefix:
        command.extend(shlex.split(args.rank_prefix))
    command.extend([str(args.cp2k.resolve()), "-i", str(source)])
    output_path = run_dir / "cp2k.out"
    stderr_path = run_dir / "cp2k.err"
    started = time.time()
    timed_out = False
    with output_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(
            command, env=env, stdout=stdout, stderr=stderr, start_new_session=True
        )
        try:
            returncode = process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                returncode = process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                returncode = process.wait(timeout=10)
    combined = output_path.read_text(errors="replace") + "\n" + stderr_path.read_text(errors="replace")
    metadata = {
        "schema": 1,
        "case": case["name"],
        "ranks": case["ranks"],
        "injection": case["injection"],
        "expected_diagnostic": case["diagnostic"],
        "input": source.name,
        "input_sha256": sha256(source),
        "cp2k": str(args.cp2k.resolve()),
        "cp2k_sha256": sha256(args.cp2k.resolve()),
        "command": command,
        "returncode": returncode,
        "timed_out": timed_out,
        "wall_seconds": time.time() - started,
        "output_sha256": sha256(output_path),
        "stderr_sha256": sha256(stderr_path),
        "diagnostic_count": combined.count(case["diagnostic"]),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if timed_out:
        raise RuntimeError(f"{case['name']}: timed out")
    if returncode == 0:
        raise RuntimeError(f"{case['name']}: fault unexpectedly returned zero")
    if case["diagnostic"] not in combined:
        raise RuntimeError(f"{case['name']}: expected diagnostic missing")
    if "PROGRAM ENDED" in combined:
        raise RuntimeError(f"{case['name']}: failed run printed PROGRAM ENDED")
    print(f"PASS_NEGATIVE\t{case['name']}\trc={returncode}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", required=True, type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--mpiexec-arg", action="append", default=[])
    parser.add_argument("--rank-prefix", default="")
    parser.add_argument("--timeout", type=float, default=180.0)
    args = parser.parse_args()
    if not args.cp2k.is_file():
        parser.error(f"CP2K executable not found: {args.cp2k}")
    for case in CASES:
        run_case(args, case)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
