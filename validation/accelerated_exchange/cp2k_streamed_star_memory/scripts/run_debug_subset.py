#!/usr/bin/env python3
"""Run a small bounds/FPE/leak-check selector qualification subset."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
import re
import signal
import subprocess
import time
from pathlib import Path

from run_test_matrix import prove_rank_affinity


ROOT = Path(__file__).resolve().parent
CP2K = Path(os.environ["CP2K_EXE"]).resolve()
CP2K_LIB = Path(os.environ["CP2K_LIB"]).resolve()
MPIEXEC = Path(os.environ["MPIEXEC_EXE"]).resolve()
INPUT_ROOT = Path(os.environ["CP2K_INPUT_ROOT"]).resolve()
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
QUALIFY_RE = re.compile(
    rf"GXTB-QUALIFICATION_ONLY MIXER-STAR iter=(\d+)"
    rf"\s+denseCov=\s*({FLOAT})\s+streamCov=\s*({FLOAT})"
    rf"\s+streamRoundtrip=\s*({FLOAT})\s+covDelta=\s*({FLOAT})"
)
CASES = [
    ("ch4_k290_p2", "ch4_k290_energy_force.inp", 2, "192-195"),
    ("o2_uks_tr_p4", "o2_uks_k311_tr_energy_force.inp", 4, "196-199"),
    ("ar2_1d_tr_p4", "ar2_1d_energy_force_print.inp", 4, "200-203"),
    ("ar4_2d_tr_p4", "ar4_2d_energy_force_print.inp", 4, "204-207"),
]
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


def run_case(case: tuple[str, str, int, str]) -> dict:
    name, input_name, ranks, cpus = case
    run_dir = ROOT / "runs" / name
    run_dir.mkdir(parents=True, exist_ok=False)
    input_path = INPUT_ROOT / input_name
    env = os.environ.copy()
    env.update(THREAD_ENV)
    env["CP2K_GXTB_SYMMETRY_STAR_CONTRACTION"] = "QUALIFY"
    env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_ROOT"
    env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = "3"
    env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
    # Open MPI/PMIx intentionally retain process-global allocations at exit;
    # disable that third-party LSan end-of-process report while retaining the
    # Debug build's bounds, FPE, undefined-behavior, and runtime checks.
    env["LSAN_OPTIONS"] = "detect_leaks=0"
    command = [
        "taskset", "-c", cpus, str(MPIEXEC), "--bind-to", "none",
        "-np", str(ranks), str(CP2K), "-i", str(input_path),
    ]
    started = time.time()
    with (run_dir / "cp2k.out").open("wb") as stdout, (run_dir / "cp2k.err").open("wb") as stderr:
        process = subprocess.Popen(
            command, env=env, stdout=stdout, stderr=stderr, start_new_session=True
        )
        first_cpu, last_cpu = map(int, cpus.split("-"))
        try:
            affinity = prove_rank_affinity(
                process.pid, ranks, set(range(first_cpu, last_cpu + 1))
            )
            returncode = process.wait(timeout=300)
        except Exception:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            process.wait(timeout=30)
            raise
    output = (run_dir / "cp2k.out").read_text(errors="replace")
    stderr = (run_dir / "cp2k.err").read_text(errors="replace")
    if returncode != 0 or output.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"{name}: debug execution failed ({returncode})")
    matches = QUALIFY_RE.findall(output)
    if not matches:
        raise RuntimeError(f"{name}: no mixer-star qualification marker")
    dense_cov = max(float(match[1]) for match in matches)
    stream_cov = max(float(match[2]) for match in matches)
    roundtrip = max(float(match[3]) for match in matches)
    covariance_delta = max(float(match[4]) for match in matches)
    if max(dense_cov, stream_cov, roundtrip, covariance_delta) > 1.0e-10:
        raise RuntimeError(f"{name}: internal selector gate failed")
    bad_stderr = [
        token for token in (
            "AddressSanitizer", "LeakSanitizer", "runtime error:",
            "Fortran runtime error", "SIGFPE", "SIGSEGV",
        ) if token in stderr
    ]
    if bad_stderr:
        raise RuntimeError(f"{name}: debug diagnostic(s): {bad_stderr}")
    result = {
        "schema": 1,
        "name": name,
        "ranks": ranks,
        "cpu_set": cpus,
        "returncode": returncode,
        "affinity_proof": affinity,
        "wall_seconds": time.time() - started,
        "dense_covariance": dense_cov,
        "stream_covariance": stream_cov,
        "stream_roundtrip": roundtrip,
        "covariance_delta": covariance_delta,
        "cp2k": str(CP2K),
        "cp2k_sha256": sha256(CP2K),
        "cp2k_lib": str(CP2K_LIB),
        "cp2k_lib_sha256": sha256(CP2K_LIB),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "output_sha256": sha256(run_dir / "cp2k.out"),
        "stderr_sha256": sha256(run_dir / "cp2k.err"),
        "command": command,
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(CASES)) as pool:
        results = list(pool.map(run_case, CASES))
    (ROOT / "summary.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {len(results)}/{len(CASES)} debug selector qualifications")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
