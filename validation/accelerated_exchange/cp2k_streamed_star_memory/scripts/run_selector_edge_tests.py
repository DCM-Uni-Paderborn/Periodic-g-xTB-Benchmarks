#!/usr/bin/env python3
"""Exercise fail-closed and default behavior of the symmetry-star selector."""

from __future__ import annotations

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
INPUT = Path(os.environ["CP2K_INPUT"]).resolve()
REFERENCE = Path(os.environ["DENSE_REFERENCE"]).resolve()
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
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


def observables(text: str) -> tuple[float, list[float], list[float]]:
    energies = [float(value) for value in re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )]
    forces = [tuple(map(float, row)) for row in re.findall(
        rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
        text,
        re.MULTILINE,
    )]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    if not energies or not forces or not stress_blocks:
        raise RuntimeError("missing energy, force, or stress observable")
    stress_rows = re.findall(
        rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        stress_blocks[-1],
        re.MULTILINE,
    )
    return (
        energies[-1],
        [value for row in forces for value in row],
        [float(value) for row in stress_rows for value in row],
    )


def run_case(name: str, selector: str | None, expected_message: str | None) -> dict:
    run_dir = ROOT / "runs" / name
    run_dir.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(THREAD_ENV)
    env["CP2K_GXTB_EXCHANGE_STREAM_MODE"] = "KGROUP_PARTIAL_ROOT"
    env["CP2K_GXTB_EXCHANGE_GRADIENT_MODE"] = "QUALIFY"
    env["CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE"] = "3"
    env["CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION"] = "1"
    env.pop("CP2K_GXTB_SYMMETRY_STAR_CONTRACTION", None)
    if selector is not None:
        env["CP2K_GXTB_SYMMETRY_STAR_CONTRACTION"] = selector
    command = [
        "taskset", "-c", "192-195", str(MPIEXEC), "--bind-to", "none",
        "-np", "1", str(CP2K), "-i", str(INPUT),
    ]
    started = time.time()
    process = subprocess.Popen(
        command,
        env=env,
        stdout=(run_dir / "cp2k.out").open("wb"),
        stderr=(run_dir / "cp2k.err").open("wb"),
        start_new_session=True,
    )
    try:
        affinity = prove_rank_affinity(process.pid, 1, set(range(192, 196)))
        returncode = process.wait(timeout=120)
    except Exception:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        process.wait(timeout=30)
        raise
    output = (run_dir / "cp2k.out").read_text(errors="replace")
    stderr = (run_dir / "cp2k.err").read_text(errors="replace")
    combined = output + "\n" + stderr
    if expected_message is None:
        if returncode != 0 or output.count("PROGRAM ENDED") != 1:
            raise RuntimeError(f"{name}: default selector did not complete cleanly")
        if "GXTB-MIXER-STAR-STREAMED" in output or "MIXER-STAR iter=" in output:
            raise RuntimeError(f"{name}: absent selector did not remain DENSE")
        ref = observables(REFERENCE.read_text(errors="replace"))
        got = observables(output)
        deltas = [
            abs(got[0] - ref[0]),
            max(abs(a - b) for a, b in zip(got[1], ref[1])),
            max(abs(a - b) for a, b in zip(got[2], ref[2])),
        ]
        if any(delta != 0.0 for delta in deltas):
            raise RuntimeError(f"{name}: default differs from explicit DENSE: {deltas}")
    else:
        if returncode == 0:
            raise RuntimeError(f"{name}: invalid selector was accepted")
        if expected_message not in combined:
            raise RuntimeError(f"{name}: missing fail-closed diagnostic")
        if "PROGRAM ENDED" in output:
            raise RuntimeError(f"{name}: invalid selector reached clean termination")
        deltas = None
    metadata = {
        "schema": 1,
        "name": name,
        "selector": selector,
        "expected_message": expected_message,
        "returncode": returncode,
        "affinity_proof": affinity,
        "wall_seconds": time.time() - started,
        "cp2k": str(CP2K),
        "cp2k_sha256": sha256(CP2K),
        "cp2k_lib": str(CP2K_LIB),
        "cp2k_lib_sha256": sha256(CP2K_LIB),
        "input": str(INPUT),
        "input_sha256": sha256(INPUT),
        "output_sha256": sha256(run_dir / "cp2k.out"),
        "stderr_sha256": sha256(run_dir / "cp2k.err"),
        "observable_deltas": deltas,
        "command": command,
    }
    (run_dir / "result.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def main() -> int:
    cases = [
        ("default_dense", None, None),
        ("unknown_value", "BOGUS", "Unknown CP2K_GXTB_SYMMETRY_STAR_CONTRACTION value"),
        (
            "overlong_value",
            "X" * 33,
            "Invalid or overlong CP2K_GXTB_SYMMETRY_STAR_CONTRACTION value",
        ),
    ]
    results = [run_case(*case) for case in cases]
    (ROOT / "summary.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print("PASS: default DENSE and 2/2 fail-closed selector cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
