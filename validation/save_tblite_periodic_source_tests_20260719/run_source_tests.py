#!/usr/bin/env python3
"""Run and archive the periodic save_tblite source-level qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


CURRENT_REVISION = "15915c9435644eb257178ca8f8bf7220c38b1a84"
PBC_REVISION = "c932120d25"


def run(command: list[str], cwd: Path) -> tuple[int, str]:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "GFORTRAN_UNBUFFERED_ALL": "1",
        }
    )
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return completed.returncode, completed.stdout


def git(root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", *arguments], cwd=root, text=True
    ).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_differences(output: str) -> list[float]:
    marker = output.find("Difference:")
    if marker < 0:
        raise RuntimeError("CeCl3 output contains no Difference block")
    values: list[float] = []
    for line in output[marker:].splitlines()[1:5]:
        values.extend(
            float(value)
            for value in re.findall(r"[-+]?\d+\.\d+E[-+]\d+", line)
        )
    if len(values) != 12:
        raise RuntimeError(f"expected 12 CeCl3 differences, found {len(values)}")
    return values


def extract_residual(output: str, label: str) -> float:
    match = re.search(rf"{re.escape(label)}\s+([-+0-9.Ee]+)", output)
    if not match:
        raise RuntimeError(f"missing residual {label}")
    return float(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--pbc-root", type=Path, required=True)
    arguments = parser.parse_args()

    archive = Path(__file__).resolve().parent
    current_binary = arguments.current_root / "build-tests/test/unit/tblite-tester"
    pbc_binary = arguments.pbc_root / "build-tests-packages/test/unit/tblite-tester"
    cases = {
        "current_h0_diamond": (
            arguments.current_root,
            [str(current_binary), "hamiltonian", "anisotropy-gxtb-diamond-pbc"],
        ),
        "current_h0_supercell": (
            arguments.current_root,
            [str(current_binary), "hamiltonian", "anisotropy-gxtb-supercell-pbc"],
        ),
        "current_h0_gradient": (
            arguments.current_root,
            [str(current_binary), "hamiltonian", "anisotropy-gradient-gxtb-diamond-pbc"],
        ),
        "current_wignerseitz": (
            arguments.current_root,
            [str(current_binary), "wignerseitz"],
        ),
        "current_exchange": (
            arguments.current_root,
            [str(current_binary), "exchange"],
        ),
        "current_acp": (
            arguments.current_root,
            [str(current_binary), "acp"],
        ),
        "current_coulomb_charge": (
            arguments.current_root,
            [str(current_binary), "coulomb-charge"],
        ),
        "current_coulomb_multipole": (
            arguments.current_root,
            [str(current_binary), "coulomb-multipole"],
        ),
        "current_dispersion": (
            arguments.current_root,
            [str(current_binary), "dispersion"],
        ),
        "current_repulsion": (
            arguments.current_root,
            [str(current_binary), "repulsion"],
        ),
        "current_gxtb": (
            arguments.current_root,
            [str(current_binary), "gxtb"],
        ),
        "current_hamiltonian": (
            arguments.current_root,
            [str(current_binary), "hamiltonian"],
        ),
        "current_cecl3": (
            arguments.current_root,
            [str(current_binary), "hamiltonian", "hamiltonian-gradient-gxtb-cecl3"],
        ),
        "pbc_cecl3": (
            arguments.pbc_root,
            [str(pbc_binary), "hamiltonian", "hamiltonian-gradient-gxtb-cecl3"],
        ),
    }

    results: dict[str, object] = {}
    for name, (root, command) in cases.items():
        returncode, output = run(command, root)
        log_name = f"{name}.log"
        (archive / log_name).write_text(output)
        results[name] = {
            "command": command,
            "log": log_name,
            "returncode": returncode,
            "passed_count": output.count("[PASSED]"),
            "failed_count": output.count("[FAILED]"),
        }

    current_cecl3 = (archive / "current_cecl3.log").read_text()
    pbc_cecl3 = (archive / "pbc_cecl3.log").read_text()
    exchange = (archive / "current_exchange.log").read_text()
    current_cecl3_differences = extract_differences(current_cecl3)
    pbc_cecl3_differences = extract_differences(pbc_cecl3)
    summary = {
        "schema": "save_tblite-periodic-source-tests-v1",
        "current": {
            "revision": git(arguments.current_root, "rev-parse", "HEAD"),
            "binary": str(current_binary),
            "binary_sha256": sha256(current_binary),
            "git_status": git(arguments.current_root, "status", "--short"),
        },
        "pbc_baseline": {
            "revision": git(arguments.pbc_root, "rev-parse", "HEAD"),
            "binary": str(pbc_binary),
            "binary_sha256": sha256(pbc_binary),
            "build_only_repair_files": git(
                arguments.pbc_root, "diff", "--name-only"
            ).splitlines(),
        },
        "compiler": subprocess.check_output(
            ["gfortran", "--version"], text=True
        ).splitlines()[0],
        "tests": results,
        "exchange_residuals_hartree": {
            "whole_mesh_vs_explicit_bvk": extract_residual(
                exchange, "BVK_EXCHANGE_SUPERCELL_MAX_RESID"
            ),
            "mixed_radix_fft_vs_dense": extract_residual(
                exchange, "MIXED_RADIX_FFT_EXCHANGE_MAX_RESID"
            ),
        },
        "cecl3_nonperiodic_finite_difference": {
            "current_differences_hartree_per_bohr": current_cecl3_differences,
            "pbc_differences_hartree_per_bohr": pbc_cecl3_differences,
            "current_max_abs_hartree_per_bohr": max(
                abs(value) for value in current_cecl3_differences
            ),
            "pbc_max_abs_hartree_per_bohr": max(
                abs(value) for value in pbc_cecl3_differences
            ),
            "test_tolerance_hartree_per_bohr": 10.0e5 * 2.220446049250313e-16,
            "identical_difference_components": current_cecl3_differences
            == pbc_cecl3_differences,
        },
    }
    (archive / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
