#!/usr/bin/env python3
"""Run a fixed DENSE/KGROUP_PARTIAL_DISTRIBUTED_IMAGES CP2K oracle matrix.

The driver is intentionally fail-closed: a missing input, unknown case/rank,
pre-existing nonempty run directory, failed subprocess, or incomplete CP2K
termination stops the campaign.  Numerical acceptance is delegated to
``verify_matrix.py`` so raw output is never silently accepted by this runner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MATRIX = json.loads((ROOT / "matrix.json").read_text())
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


def cases() -> list[dict]:
    expanded = []
    for case in MATRIX["cases"]:
        source = ROOT / "inputs" / case["input"]
        if not source.is_file():
            raise RuntimeError(f"missing frozen input: {source}")
        for ranks in case["ranks"]:
            expanded.append({**case, "ranks": int(ranks)})
    return expanded


def complete(text: str) -> bool:
    # CP2K's normal footer contains both "PROGRAM ENDED AT" and the provenance
    # line "PROGRAM STOPPED IN <directory>"; the latter is not an error.
    return text.count("PROGRAM ENDED") == 1


def run_one(args: argparse.Namespace, case: dict, variant: str) -> None:
    ranks = case["ranks"]
    run_id = f"{case['name']}_p{ranks}_{variant.lower()}"
    run_dir = ROOT / "runs" / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"refusing to reuse nonempty run directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)

    input_path = (ROOT / "inputs" / case["input"]).resolve()
    output_path = run_dir / "cp2k.out"
    stderr_path = run_dir / "cp2k.err"

    env = os.environ.copy()
    env.update(THREAD_ENV)
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = str(args.batch_size)
    if variant == "DENSE":
        env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "LEGACY"
        env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "DENSE"
        env.pop("CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION", None)
    elif variant == "PARTIAL_DISTRIBUTED_IMAGES":
        env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_DISTRIBUTED_IMAGES"
        env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
        env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
    else:
        raise RuntimeError(f"unknown variant: {variant}")

    command = [args.mpiexec, *args.mpiexec_arg, "-np", str(ranks)]
    if args.rank_prefix:
        command.extend(shlex.split(args.rank_prefix))
    command.extend([str(args.cp2k.resolve()), "-i", str(input_path)])

    metadata = {
        "schema": 1,
        "run_id": run_id,
        "case": case["name"],
        "features": case["features"],
        "expected_nfull": case["nfull"],
        "ranks": ranks,
        "variant": variant,
        "batch_size": args.batch_size,
        "input": input_path.name,
        "input_sha256": sha256(input_path),
        "cp2k": str(args.cp2k.resolve()),
        "cp2k_sha256": sha256(args.cp2k.resolve()),
        "command": command,
        "environment": {key: env[key] for key in sorted(set(THREAD_ENV) | {
            "CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE",
            "CP2K_GXTB_EXCHANGE_STREAM_MODE",
            "CP2K_GXTB_EXCHANGE_GRADIENT_MODE",
            "CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION",
        }) if key in env},
        "started_unix": time.time(),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    with output_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(command, env=env, stdout=stdout, stderr=stderr, check=False)
    wall = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)

    metadata.update({
        "finished_unix": time.time(),
        "wall_seconds": wall,
        "returncode": completed.returncode,
        "child_user_seconds_delta": after.ru_utime - before.ru_utime,
        "child_system_seconds_delta": after.ru_stime - before.ru_stime,
        "child_maxrss_kb_campaign_highwater": after.ru_maxrss,
        "output_sha256": sha256(output_path),
        "stderr_sha256": sha256(stderr_path),
    })
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (run_dir / "returncode.txt").write_text(f"{completed.returncode}\n")

    text = output_path.read_text(errors="replace")
    if completed.returncode != 0:
        raise RuntimeError(f"{run_id}: subprocess return code {completed.returncode}")
    if not complete(text):
        raise RuntimeError(f"{run_id}: CP2K did not terminate exactly once and cleanly")
    print(f"COMPLETE\t{run_id}\t{wall:.3f}s", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", required=True, type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--mpiexec-arg", action="append", default=[])
    parser.add_argument("--rank-prefix", default="")
    parser.add_argument("--batch-size", type=int, default=3)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument(
        "--variant",
        choices=("DENSE", "PARTIAL_DISTRIBUTED_IMAGES", "BOTH"),
        default="BOTH",
    )
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    if not args.cp2k.is_file():
        parser.error(f"CP2K executable not found: {args.cp2k}")

    selected = cases()
    if args.case:
        requested = set(args.case)
        selected = [case for case in selected if f"{case['name']}_p{case['ranks']}" in requested]
        found = {f"{case['name']}_p{case['ranks']}" for case in selected}
        if found != requested:
            raise RuntimeError(f"unknown case(s): {sorted(requested - found)}")

    variants = (
        ("DENSE", "PARTIAL_DISTRIBUTED_IMAGES")
        if args.variant == "BOTH"
        else (args.variant,)
    )
    for case in selected:
        for variant in variants:
            run_one(args, case, variant)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
