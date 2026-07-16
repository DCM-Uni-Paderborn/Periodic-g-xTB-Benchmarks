#!/usr/bin/env python3
"""Run the fixed mixer symmetry-star storage qualification matrix."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MATRIX = json.loads((ROOT / "test_matrix.json").read_text())
CP2K = Path(os.environ["CP2K_EXE"]).resolve()
CP2K_LIB = Path(os.environ["CP2K_LIB"]).resolve()
MPIEXEC = Path(os.environ["MPIEXEC_EXE"]).resolve()
CPU_FIRST = int(os.environ.get("CPU_FIRST", "192"))
CPU_SLOTS = int(os.environ.get("CPU_SLOTS", "8"))
CPUS_PER_SLOT = int(os.environ.get("CPUS_PER_SLOT", "4"))
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "MKL_DYNAMIC": "FALSE",
    "BLIS_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
}


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


def prove_rank_affinity(launcher_pid: int, ranks: int, allowed: set[int]) -> list[dict]:
    """Fail unless all live CP2K ranks inherit the exact requested CPU set."""
    deadline = time.monotonic() + 30.0
    while time.monotonic() < deadline:
        snapshots = []
        for pid in sorted(descendants(launcher_pid)):
            try:
                executable = Path(f"/proc/{pid}/exe").resolve()
                if executable != CP2K:
                    continue
                status = Path(f"/proc/{pid}/status").read_text()
                allowed_text = next(
                    line.split(":", 1)[1].strip()
                    for line in status.splitlines()
                    if line.startswith("Cpus_allowed_list:")
                )
                actual = set(os.sched_getaffinity(pid))
                stat_text = Path(f"/proc/{pid}/stat").read_text()
                tail = stat_text[stat_text.rfind(")") + 2:].split()
                processor = int(tail[36])  # Linux /proc PID stat field 39.
            except (FileNotFoundError, PermissionError, ProcessLookupError):
                continue
            snapshots.append({
                "pid": pid,
                "cpus_allowed_list": allowed_text,
                "cpus_allowed": sorted(actual),
                "processor": processor,
            })
        if len(snapshots) == ranks:
            if all(set(item["cpus_allowed"]) == allowed and item["processor"] in allowed
                   for item in snapshots):
                return snapshots
            raise RuntimeError(f"rank affinity escaped requested CPU set: {snapshots}")
        if not Path(f"/proc/{launcher_pid}").exists():
            break
        time.sleep(0.01)
    raise RuntimeError(f"could not prove live affinity for {ranks} CP2K ranks")


def run_one(job: tuple[dict, int, str], slot: int) -> str:
    case, ranks, variant = job
    run_id = f"{case['name']}_p{ranks}_{variant.lower()}"
    run_dir = ROOT / "runs" / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(f"refusing to reuse nonempty run directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = (ROOT / "inputs" / case["input"]).resolve()
    if not input_path.is_file():
        raise RuntimeError(f"missing input: {input_path}")
    first_cpu = CPU_FIRST + slot * CPUS_PER_SLOT
    last_cpu = first_cpu + CPUS_PER_SLOT - 1
    if ranks > CPUS_PER_SLOT:
        raise RuntimeError(f"rank count exceeds fixed CPU slot: {run_id}")

    env = os.environ.copy()
    env.update(THREAD_ENV)
    env["CP2K_GXTB_SYMMETRY_STAR_CONTRACTION"] = variant
    env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_ROOT"
    env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = "3"
    env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
    command = [
        "taskset", "-c", f"{first_cpu}-{last_cpu}", str(MPIEXEC),
        "--bind-to", "none", "-np", str(ranks),
        str(CP2K), "-i", str(input_path),
    ]
    metadata = {
        "schema": 1,
        "run_id": run_id,
        "case": case["name"],
        "features": case["features"],
        "expected_nfull": case["nfull"],
        "ranks": ranks,
        "variant": variant,
        "cpu_set": f"{first_cpu}-{last_cpu}",
        "input": input_path.name,
        "input_sha256": sha256(input_path),
        "cp2k": str(CP2K),
        "cp2k_sha256": sha256(CP2K),
        "cp2k_lib": str(CP2K_LIB),
        "cp2k_lib_sha256": sha256(CP2K_LIB),
        "command": command,
        "environment": {key: env[key] for key in sorted(set(THREAD_ENV) | {
            "CP2K_GXTB_SYMMETRY_STAR_CONTRACTION",
            "CP2K_GXTB_EXCHANGE_STREAM_MODE",
            "CP2K_GXTB_EXCHANGE_GRADIENT_MODE",
            "CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE",
            "CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION",
        })},
        "started_unix": time.time(),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    output_path = run_dir / "cp2k.out"
    stderr_path = run_dir / "cp2k.err"
    started = time.perf_counter()
    with output_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(
            command, env=env, stdout=stdout, stderr=stderr, start_new_session=True
        )
        try:
            affinity = prove_rank_affinity(process.pid, ranks, set(range(first_cpu, last_cpu + 1)))
        except Exception:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            process.wait(timeout=30)
            raise
        returncode = process.wait()
    wall = time.perf_counter() - started
    metadata.update({
        "finished_unix": time.time(),
        "wall_seconds": wall,
        "returncode": returncode,
        "affinity_proof": affinity,
        "output_sha256": sha256(output_path),
        "stderr_sha256": sha256(stderr_path),
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
    return f"COMPLETE\t{run_id}\t{wall:.3f}s\tCPU={first_cpu}-{last_cpu}"


def main() -> int:
    if not CP2K.is_file():
        raise RuntimeError(f"missing CP2K executable: {CP2K}")
    if not CP2K_LIB.is_file():
        raise RuntimeError(f"missing CP2K shared library: {CP2K_LIB}")
    if not MPIEXEC.is_file():
        raise RuntimeError(f"missing MPI launcher: {MPIEXEC}")
    all_jobs = jobs()
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
    print(f"ALL_COMPLETE\t{len(all_jobs)}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
