#!/usr/bin/env python3
"""Run the complete 2^3/3^3 DMC-ICE13 CLI matrix for one executable."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--meshes", type=int, nargs="+", default=(2, 3))
    parser.add_argument("--accuracy", type=float, default=0.1)
    parser.add_argument("--require-binary-sha256")
    args = parser.parse_args()

    runner = HERE / "run_save_tblite.py"
    for mesh in args.meshes:
        for phase in PHASES:
            command = [
                sys.executable,
                str(runner),
                str(args.executable),
                phase,
                str(mesh),
                str(args.output_root),
                "--accuracy",
                str(args.accuracy),
            ]
            if args.require_binary_sha256:
                command.extend(("--require-binary-sha256", args.require_binary_sha256))
            print(f"running mesh={mesh} phase={phase}", flush=True)
            subprocess.run(command, check=True)
    print("status=PASS", flush=True)


if __name__ == "__main__":
    main()
