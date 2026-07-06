#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]
METHODS = ["GFN1", "GFN2"]
MESHES = ["gamma", "k111", "k222", "k333", "k444", "k555"]


@dataclass(frozen=True)
class Job:
    mesh: str
    method: str
    phase: str
    input_path: Path
    run_dir: Path
    output_name: str


def jobs(root: Path) -> list[Job]:
    out: list[Job] = []
    for mesh in MESHES:
        for method in METHODS:
            for phase in PHASES:
                input_path = root / "kpoint_inputs" / mesh / f"ice_{phase}_{method}_{mesh}.inp"
                if mesh == "gamma":
                    run_dir = root / "runs" / method / phase
                    output_name = f"ice_{phase}_{method}.out"
                else:
                    run_dir = root / "runs_kpoints" / mesh / method / phase
                    output_name = f"ice_{phase}_{method}_{mesh}.out"
                out.append(Job(mesh, method, phase, input_path, run_dir, output_name))
    return out


def has_completed(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and "ENERGY| Total FORCE_EVAL" in text


def run_job(cp2k: Path, job: Job, resume: bool) -> tuple[Job, int]:
    output = job.run_dir / job.output_name
    if resume and has_completed(output):
        return job, 0
    if job.run_dir.exists() and not resume:
        shutil.rmtree(job.run_dir)
    job.run_dir.mkdir(parents=True, exist_ok=True)
    local_input = job.run_dir / job.input_path.name
    shutil.copy2(job.input_path, local_input)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    log_path = job.run_dir / "run.log"
    with log_path.open("w") as log:
        proc = subprocess.run(
            [str(cp2k), "-i", local_input.name, "-o", job.output_name],
            cwd=job.run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            text=True,
        )
    return job, proc.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    script = args.root / "scripts" / "dmc_ice13_kpoint_benchmark.py"
    subprocess.run([sys.executable, str(script), "prepare"], cwd=args.root, check=True)

    failures: list[tuple[Job, int]] = []
    all_jobs = jobs(args.root)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_job, args.cp2k, job, args.resume) for job in all_jobs]
        done = 0
        for future in concurrent.futures.as_completed(futures):
            job, rc = future.result()
            done += 1
            if rc != 0:
                failures.append((job, rc))
            print(f"{done:3d}/{len(all_jobs)} {job.mesh:5s} {job.method:4s} {job.phase:4s} rc={rc}", flush=True)
    if failures:
        for job, rc in failures[:20]:
            print(f"FAILED {job.mesh} {job.method} {job.phase} rc={rc} out={job.run_dir / job.output_name}", file=sys.stderr)
        raise SystemExit(f"{len(failures)} DMC13 jobs failed")

    subprocess.run([sys.executable, str(script), "analyse"], cwd=args.root, check=True)
    print(args.root / "data" / "dmc_ice13_kpoint_stats.csv")


if __name__ == "__main__":
    main()
