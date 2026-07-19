#!/usr/bin/env python3
"""Qualification tests for the native/full/Gamma/direct energy oracles."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
SYMMETRY = TOOLS / "compare_native_symmetry_cli.py"
GAMMA = TOOLS / "compare_gamma_supercell_oracle.py"
BINARY_DIGEST = "a" * 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_cp2k(run_dir: Path, energy: float, input_path: Path) -> Path:
    run_dir.mkdir(parents=True)
    output = run_dir / "cp2k.out"
    output.write_text(
        f" ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.): {energy:.15f}\n"
        " PROGRAM ENDED AT 2026-07-19 00:00:00\n",
        encoding="utf-8",
    )
    (run_dir / "exit_status").write_text("0\n", encoding="utf-8")
    (run_dir / "binary.sha256").write_text(
        f"{BINARY_DIGEST}  cp2k.psmp\n", encoding="utf-8"
    )
    (run_dir / "input.sha256").write_text(
        f"{digest(input_path)}  {input_path}\n", encoding="utf-8"
    )
    return output


def write_cli(run_dir: Path, energy: float) -> Path:
    run_dir.mkdir(parents=True)
    result = run_dir / "result.json"
    result.write_text(json.dumps({"energy": energy}) + "\n", encoding="utf-8")
    (run_dir / "exit_status").write_text("0\n", encoding="utf-8")
    return result


class OracleQualificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.reduced_input = self.root / "reduced.inp"
        self.full_input = self.root / "full.inp"
        self.gamma_input = self.root / "gamma.inp"
        self.cli_input = self.root / "POSCAR"
        self.reduced_input.write_text("reduced\n", encoding="utf-8")
        self.full_input.write_text("full\n", encoding="utf-8")
        self.gamma_input.write_text("gamma\n", encoding="utf-8")
        self.cli_input.write_text("structure\n", encoding="utf-8")
        self.reduced = write_cp2k(self.root / "reduced", -1.0, self.reduced_input)
        self.full = write_cp2k(self.root / "full", -1.0, self.full_input)
        self.gamma = write_cp2k(self.root / "gamma", -8.0, self.gamma_input)
        self.cli = write_cli(self.root / "cli", -8.0)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_command(self, arguments: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, *arguments],
            text=True,
            capture_output=True,
            check=False,
        )

    def symmetry_arguments(self) -> list[str]:
        return [
            str(SYMMETRY),
            str(self.reduced),
            str(self.full),
            str(self.cli),
            "--replicas",
            "8",
            "--require-binary-sha256",
            BINARY_DIGEST,
            "--reduced-input",
            str(self.reduced_input),
            "--full-input",
            str(self.full_input),
            "--cli-input",
            str(self.cli_input),
            "--require-reduced-input-sha256",
            digest(self.reduced_input),
            "--require-full-input-sha256",
            digest(self.full_input),
            "--require-cli-input-sha256",
            digest(self.cli_input),
        ]

    def test_symmetry_oracle_accepts_fully_qualified_parity(self) -> None:
        completed = self.run_command(self.symmetry_arguments())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(
            payload["provenance"]["reduced_cp2k"]["binary_sha256"],
            BINARY_DIGEST,
        )

    def test_symmetry_oracle_rejects_unsuccessful_cli(self) -> None:
        (self.cli.parent / "exit_status").write_text("1\n", encoding="utf-8")
        completed = self.run_command(self.symmetry_arguments())
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("nonzero CLI exit status", completed.stderr)

    def test_symmetry_oracle_rejects_changed_input(self) -> None:
        arguments = self.symmetry_arguments()
        self.reduced_input.write_text("changed\n", encoding="utf-8")
        completed = self.run_command(arguments)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("wrong CP2K input", completed.stderr)

    def test_gamma_oracle_accepts_fully_qualified_parity(self) -> None:
        completed = self.run_command(
            [
                str(GAMMA),
                str(self.reduced),
                str(self.gamma),
                str(self.cli),
                "--replicas",
                "8",
                "--require-binary-sha256",
                BINARY_DIGEST,
                "--native-input",
                str(self.reduced_input),
                "--gamma-input",
                str(self.gamma_input),
                "--cli-input",
                str(self.cli_input),
                "--require-native-input-sha256",
                digest(self.reduced_input),
                "--require-gamma-input-sha256",
                digest(self.gamma_input),
                "--require-cli-input-sha256",
                digest(self.cli_input),
            ]
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(
            payload["provenance"]["gamma_supercell_cp2k"]["input_sha256"],
            digest(self.gamma_input),
        )


if __name__ == "__main__":
    unittest.main()
