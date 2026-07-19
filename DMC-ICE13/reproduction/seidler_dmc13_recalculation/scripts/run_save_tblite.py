#!/usr/bin/env python3
"""Run one frozen DMC-ICE13 BvK cell with a selected save_tblite CLI."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

from build_bvk_from_poscar import build


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("phase", choices=PHASES)
    parser.add_argument("mesh", type=int)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--require-binary-sha256")
    parser.add_argument("--accuracy", type=float, default=0.1)
    args = parser.parse_args()
    if args.mesh < 1 or args.accuracy <= 0.0:
        parser.error("mesh and accuracy must be positive")
    executable = args.executable.resolve()
    if not executable.is_file():
        parser.error(f"missing executable: {executable}")
    binary_hash = sha256(executable)
    if args.require_binary_sha256 and binary_hash != args.require_binary_sha256.lower():
        parser.error("executable SHA-256 differs from --require-binary-sha256")

    mesh_id = f"k{args.mesh}{args.mesh}{args.mesh}"
    structure = args.output_root / "structures" / mesh_id / args.phase / "POSCAR"
    build(PACKAGE / "structures/primitive" / args.phase / "POSCAR", structure, args.mesh)
    run = args.output_root / "runs" / mesh_id / args.phase
    run.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable), "run", "--method", "gxtb", "--acc", str(args.accuracy),
        "--iterations", "300", "--no-restart", "--json", "result.json", str(structure.resolve()),
    ]
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    completed = subprocess.run(
        command, cwd=run, env=environment, text=True, capture_output=True, check=False
    )
    (run / "process.out").write_text(completed.stdout, encoding="utf-8")
    (run / "process.err").write_text(completed.stderr, encoding="utf-8")
    (run / "exit_status").write_text(f"{completed.returncode}\n", encoding="utf-8")
    (run / "binary.sha256").write_text(f"{binary_hash}  {executable}\n", encoding="utf-8")
    (run / "input.sha256").write_text(f"{sha256(structure)}  {structure.resolve()}\n", encoding="utf-8")
    (run / "command.json").write_text(json.dumps(command, indent=2) + "\n", encoding="utf-8")
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
