#!/usr/bin/env python3
"""Run the fail-closed CP2K distributed-image qualification campaign."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[1]
MATRIX_PATH = CAMPAIGN / "harness" / "campaign_matrix.json"
RUNS = CAMPAIGN / "formal_runs"
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "MKL_DYNAMIC": "FALSE",
    "BLIS_NUM_THREADS": "1",
    "GOTO_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
LANES = [list(range(first, first + 4)) for first in range(120, 152, 4)]
PRINT_LOCK = threading.Lock()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def status_fields(pid: int) -> dict[str, str]:
    result: dict[str, str] = {}
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                result[key] = value.strip()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        pass
    return result


def session_processes(session_id: int) -> list[dict]:
    rows = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            stat = (entry / "stat").read_text()
            close = stat.rfind(")")
            fields = stat[close + 2 :].split()
            if len(fields) < 37 or int(fields[3]) != session_id:
                continue
            processor = int(fields[36])
            command = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode(
                errors="replace"
            ).strip()
            executable = os.readlink(entry / "exe")
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
            continue
        status = status_fields(pid)
        rows.append(
            {
                "pid": pid,
                "processor": processor,
                "allowed": status.get("Cpus_allowed_list", ""),
                "rss_kb": int(status.get("VmRSS", "0 kB").split()[0]),
                "hwm_kb": int(status.get("VmHWM", "0 kB").split()[0]),
                "command": command,
                "is_cp2k": Path(executable).name == "cp2k.psmp",
            }
        )
    return rows


def run_process(
    command: list[str], env: dict[str, str], run_dir: Path, expected_cores: list[int]
) -> tuple[int, dict]:
    output_path = run_dir / "cp2k.out"
    stderr_path = run_dir / "cp2k.err"
    affinity_path = run_dir / "proc_affinity.tsv"
    observed_masks: dict[int, str] = {}
    peak_rank_rss = 0
    peak_rank_hwm = 0
    affinity_violation = ""
    started = time.perf_counter()
    with output_path.open("wb") as stdout, stderr_path.open("wb") as stderr, affinity_path.open(
        "w"
    ) as affinity:
        affinity.write("elapsed_s\tpid\tprocessor\tallowed\tis_cp2k\trss_kb\thwm_kb\tcommand\n")
        process = subprocess.Popen(
            command,
            env=env,
            stdout=stdout,
            stderr=stderr,
            cwd=run_dir,
            start_new_session=True,
        )
        while process.poll() is None:
            elapsed = time.perf_counter() - started
            rows = session_processes(process.pid)
            rank_rss = 0
            live_masks = []
            for row in rows:
                affinity.write(
                    f"{elapsed:.6f}\t{row['pid']}\t{row['processor']}\t{row['allowed']}\t"
                    f"{int(row['is_cp2k'])}\t{row['rss_kb']}\t{row['hwm_kb']}\t"
                    f"{row['command']}\n"
                )
                if row["is_cp2k"]:
                    observed_masks[row["pid"]] = row["allowed"]
                    live_masks.append(row["allowed"])
                    rank_rss += row["rss_kb"]
                    peak_rank_hwm = max(peak_rank_hwm, row["hwm_kb"])
            expected = {str(core) for core in expected_cores}
            if any("," in mask or "-" in mask or mask not in expected for mask in live_masks):
                affinity_violation = f"invalid live rank masks: {live_masks}, expected={sorted(expected)}"
            elif len(live_masks) != len(set(live_masks)):
                affinity_violation = f"duplicate live rank masks: {live_masks}"
            elif len(live_masks) > len(expected_cores):
                affinity_violation = f"too many live CP2K ranks: {len(live_masks)}"
            if affinity_violation:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            peak_rank_rss = max(peak_rank_rss, rank_rss)
            affinity.flush()
            time.sleep(0.02)
        returncode = process.wait()
    wall = time.perf_counter() - started

    expected = {str(core) for core in expected_cores}
    actual = set(observed_masks.values())
    if affinity_violation:
        raise RuntimeError(f"{run_dir.name}: {affinity_violation}")
    if actual and actual != expected:
        raise RuntimeError(
            f"rank affinity mismatch in {run_dir.name}: observed={sorted(actual)} "
            f"expected={sorted(expected)}"
        )
    if any("," in mask or "-" in mask for mask in actual):
        raise RuntimeError(f"non-singleton CP2K rank mask in {run_dir.name}: {sorted(actual)}")

    proof_rows = []
    for rank, core in enumerate(expected_cores):
        proof_path = run_dir / f"preexec_rank_{rank}.tsv"
        if not proof_path.is_file():
            raise RuntimeError(f"missing pre-exec affinity proof: {proof_path}")
        lines = proof_path.read_text().splitlines()
        if len(lines) != 2:
            raise RuntimeError(f"malformed pre-exec affinity proof: {proof_path}")
        fields = lines[1].split("\t")
        if len(fields) != 5 or int(fields[0]) != rank or fields[3] != str(core) or fields[4] != str(core):
            raise RuntimeError(f"wrong pre-exec affinity proof: {proof_path}: {fields}")
        proof_rows.append(lines[1])
    proof_path = run_dir / "preexec_affinity.tsv"
    proof_path.write_text("rank\tpid\tprocessor\tallowed\texpected\n" + "\n".join(proof_rows) + "\n")

    report = stderr_path.read_text(errors="replace")
    for rank in range(len(expected_cores)):
        if not re.search(rf"\brank\s+{rank}\b.*\bbound to\b", report, re.IGNORECASE):
            raise RuntimeError(f"missing Open MPI binding report for rank {rank}: {run_dir.name}")
    return returncode, {
        "wall_seconds": wall,
        "peak_sampled_cp2k_rank_rss_kb": peak_rank_rss,
        "peak_observed_cp2k_rank_vmhwm_kb": peak_rank_hwm,
        "observed_rank_pids": sorted(observed_masks),
        "observed_singleton_masks": sorted(actual, key=int),
        "preexec_singleton_masks": [str(core) for core in expected_cores],
        "preexec_affinity_sha256": sha256(proof_path),
        "affinity_sha256": sha256(affinity_path),
        "bindings_sha256": sha256(stderr_path),
    }


def make_environment(batch: int, variant: str, injection: str = "") -> dict[str, str]:
    env = os.environ.copy()
    env.update(THREAD_ENV)
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = str(batch)
    if variant == "DENSE":
        env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "LEGACY"
        env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "DENSE"
        env.pop("CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION", None)
        env.pop("CP2K_GXTB_PARTIAL_QUALIFY_INJECT", None)
    else:
        env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_DISTRIBUTED_IMAGES"
        env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
        env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
        if injection:
            env["CP2K_GXTB_PARTIAL_QUALIFY_INJECT"] = injection
        else:
            env.pop("CP2K_GXTB_PARTIAL_QUALIFY_INJECT", None)
    return env


def mpi_command(mpirun: Path, cp2k: Path, ranks: int, cores: list[int], input_path: Path) -> list[str]:
    if len(cores) != ranks:
        raise RuntimeError("core/rank cardinality mismatch")
    pe_list = ",".join(str(core) for core in cores)
    wrapper = CAMPAIGN / "scripts" / "gxtb_rank_affinity_proof.sh"
    if not wrapper.is_file():
        raise RuntimeError(f"missing rank affinity wrapper: {wrapper}")
    return [
        str(mpirun),
        "--map-by",
        f"pe-list={pe_list}:ordered",
        "--bind-to",
        "core",
        "--report-bindings",
        "-np",
        str(ranks),
        str(wrapper),
        str(cp2k),
        "-i",
        str(input_path),
    ]


def run_variant(args: argparse.Namespace, job: dict, variant: str, lane: list[int]) -> None:
    ranks = int(job["ranks"])
    batch = int(job["batch"])
    run_id = f"{job['case']}_p{ranks}_b{batch}_{variant.lower()}"
    run_dir = RUNS / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        if args.resume:
            meta_path = run_dir / "run.json"
            output_path = run_dir / "cp2k.out"
            if meta_path.is_file() and output_path.is_file():
                meta = json.loads(meta_path.read_text())
                if meta.get("returncode") == 0 and output_path.read_text(errors="replace").count(
                    "PROGRAM ENDED"
                ) == 1:
                    with PRINT_LOCK:
                        print(f"SKIP_COMPLETE\t{run_id}", flush=True)
                    return
        raise RuntimeError(f"refusing to reuse nonempty run directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = (CAMPAIGN / "harness" / "inputs" / job["input"]).resolve()
    if not input_path.is_file():
        raise RuntimeError(f"missing frozen input: {input_path}")
    cores = lane[:ranks]
    command = mpi_command(args.mpirun, args.cp2k, ranks, cores, input_path)
    env = make_environment(batch, variant)
    env["GXTB_AFFINITY_DIR"] = str(run_dir)
    env["GXTB_EXPECTED_CORES"] = ",".join(str(core) for core in cores)
    env["GXTB_EXPECTED_RANKS"] = str(ranks)
    metadata = {
        "schema": 1,
        "kind": "positive",
        "run_id": run_id,
        "case": job["case"],
        "features": job["features"],
        "expected_nfull": int(job["nfull"]),
        "ranks": ranks,
        "batch_size": batch,
        "variant": variant,
        "cores": cores,
        "input": input_path.name,
        "input_sha256": sha256(input_path),
        "cp2k": str(args.cp2k),
        "cp2k_sha256": sha256(args.cp2k),
        "command": command,
        "environment": {key: env[key] for key in sorted(set(THREAD_ENV) | {
            "CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE",
            "CP2K_GXTB_EXCHANGE_STREAM_MODE",
            "CP2K_GXTB_EXCHANGE_GRADIENT_MODE",
            "CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION",
            "GXTB_AFFINITY_DIR",
            "GXTB_EXPECTED_CORES",
            "GXTB_EXPECTED_RANKS",
        }) if key in env},
        "started_unix": time.time(),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    returncode, runtime = run_process(command, env, run_dir, cores)
    metadata.update(runtime)
    metadata.update(
        {
            "finished_unix": time.time(),
            "returncode": returncode,
            "output_sha256": sha256(run_dir / "cp2k.out"),
            "stderr_sha256": sha256(run_dir / "cp2k.err"),
        }
    )
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    output = (run_dir / "cp2k.out").read_text(errors="replace")
    if returncode != 0 or output.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"{run_id}: CP2K did not terminate exactly once and cleanly")
    with PRINT_LOCK:
        print(f"COMPLETE\t{run_id}\t{runtime['wall_seconds']:.3f}s\tcores={cores}", flush=True)


def run_pair(args: argparse.Namespace, job: dict, lane: list[int]) -> None:
    run_variant(args, job, "DENSE", lane)
    run_variant(args, job, "DISTRIBUTED_IMAGES", lane)


def run_fault(args: argparse.Namespace, fault: dict, lane: list[int]) -> None:
    ranks = int(fault["ranks"])
    cores = lane[:ranks]
    run_dir = RUNS / f"fault_{fault['name']}"
    if run_dir.exists() and any(run_dir.iterdir()):
        if args.resume:
            meta_path = run_dir / "run.json"
            output_path = run_dir / "cp2k.out"
            stderr_path = run_dir / "cp2k.err"
            if meta_path.is_file() and output_path.is_file() and stderr_path.is_file():
                meta = json.loads(meta_path.read_text())
                combined = output_path.read_text(errors="replace") + "\n" + stderr_path.read_text(
                    errors="replace"
                )
                if (
                    meta.get("returncode") != 0
                    and fault["diagnostic"] in combined
                    and "PROGRAM ENDED" not in combined
                ):
                    with PRINT_LOCK:
                        print(f"SKIP_COMPLETE_FAULT\t{fault['name']}", flush=True)
                    return
        raise RuntimeError(f"refusing to reuse nonempty fault directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = (CAMPAIGN / "harness" / "inputs" / fault["input"]).resolve()
    env = make_environment(3, "DISTRIBUTED_IMAGES", fault["injection"])
    env["GXTB_AFFINITY_DIR"] = str(run_dir)
    env["GXTB_EXPECTED_CORES"] = ",".join(str(core) for core in cores)
    env["GXTB_EXPECTED_RANKS"] = str(ranks)
    command = mpi_command(args.mpirun, args.cp2k, ranks, cores, input_path)
    metadata = {
        "schema": 1,
        "kind": "fault",
        "case": fault["name"],
        "ranks": ranks,
        "cores": cores,
        "injection": fault["injection"],
        "expected_diagnostic": fault["diagnostic"],
        "input": input_path.name,
        "input_sha256": sha256(input_path),
        "cp2k_sha256": sha256(args.cp2k),
        "command": command,
        "started_unix": time.time(),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    returncode, runtime = run_process(command, env, run_dir, cores)
    combined = (run_dir / "cp2k.out").read_text(errors="replace") + "\n" + (
        run_dir / "cp2k.err"
    ).read_text(errors="replace")
    metadata.update(runtime)
    metadata.update(
        {
            "finished_unix": time.time(),
            "returncode": returncode,
            "diagnostic_count": combined.count(fault["diagnostic"]),
            "output_sha256": sha256(run_dir / "cp2k.out"),
            "stderr_sha256": sha256(run_dir / "cp2k.err"),
        }
    )
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if returncode == 0 or fault["diagnostic"] not in combined or "PROGRAM ENDED" in combined:
        raise RuntimeError(f"fault gate failed: {fault['name']}")
    with PRINT_LOCK:
        print(f"PASS_FAULT\t{fault['name']}\trc={returncode}\tcores={cores}", flush=True)


def distribute(items: list[dict]) -> list[list[dict]]:
    queues = [[] for _ in LANES]
    loads = [0 for _ in LANES]
    for item in sorted(items, key=lambda row: (int(row["ranks"]), int(row.get("nfull", 1))), reverse=True):
        lane = min(range(len(LANES)), key=lambda index: loads[index])
        queues[lane].append(item)
        loads[lane] += int(item["ranks"]) * int(item.get("nfull", 1))
    return queues


def run_lane(args: argparse.Namespace, lane_index: int, jobs: list[dict]) -> None:
    lane = LANES[lane_index]
    for job in jobs:
        run_pair(args, job, lane)


def run_fault_lane(args: argparse.Namespace, lane_index: int, faults: list[dict]) -> None:
    lane = LANES[lane_index]
    for fault in faults:
        run_fault(args, fault, lane)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", required=True, type=Path)
    parser.add_argument("--mpirun", required=True, type=Path)
    parser.add_argument("--only-job", action="append", default=[])
    parser.add_argument("--skip-faults", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.cp2k = args.cp2k.resolve()
    args.mpirun = args.mpirun.resolve()
    if not args.cp2k.is_file() or not args.mpirun.is_file():
        parser.error("CP2K executable or mpirun is missing")
    matrix = json.loads(MATRIX_PATH.read_text())
    RUNS.mkdir(parents=True, exist_ok=True)
    if any(RUNS.iterdir()) and not args.resume:
        raise RuntimeError("campaign run directory is not empty")

    api_gate = CAMPAIGN / "scripts" / "check_partial_api.py"
    source = CAMPAIGN.parents[0] / "src" / "src" / "tblite_interface.F"
    subprocess.run([sys.executable, str(api_gate), str(source)], check=True)

    jobs = matrix["jobs"]
    if args.only_job:
        requested = set(args.only_job)
        jobs = [
            job
            for job in jobs
            if f"{job['case']}_p{job['ranks']}_b{job['batch']}" in requested
        ]
        found = {f"{job['case']}_p{job['ranks']}_b{job['batch']}" for job in jobs}
        if found != requested:
            raise RuntimeError(f"unknown requested job(s): {sorted(requested - found)}")
    queues = distribute(jobs)
    with ThreadPoolExecutor(max_workers=len(LANES)) as pool:
        futures = [
            pool.submit(run_lane, args, lane_index, jobs)
            for lane_index, jobs in enumerate(queues)
            if jobs
        ]
        for future in as_completed(futures):
            future.result()

    if not args.skip_faults:
        fault_queues = distribute(matrix["faults"])
        with ThreadPoolExecutor(max_workers=len(LANES)) as pool:
            futures = [
                pool.submit(run_fault_lane, args, lane_index, faults)
                for lane_index, faults in enumerate(fault_queues)
                if faults
            ]
            for future in as_completed(futures):
                future.result()

    marker = "SMOKE_FINISHED.txt" if args.only_job or args.skip_faults else "RUN_FINISHED.txt"
    (CAMPAIGN / marker).write_text(time.strftime("%Y-%m-%dT%H:%M:%S%z\n"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
