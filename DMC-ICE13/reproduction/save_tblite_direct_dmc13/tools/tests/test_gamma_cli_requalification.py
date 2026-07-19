#!/usr/bin/env python3
"""Qualification tests for a fresh direct Gamma CLI repetition."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_gamma_cli_requalification.py"
BINARY = "a" * 64
SOURCE_REVISION = "b" * 40
PHASES = ("Ih", "VII")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class GammaCliRequalificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = self.root / "archive"
        self.current = self.root / "current"
        self.affinity = self.root / "affinity_preexec.txt"
        self.controller_status = self.root / "exit_status"
        self.source_identity = self.root / "source_identity.txt"
        self.affinity.write_text(
            "pid=1 expected_cpu=74 allowed=74\nCpus_allowed_list:\t74\n",
            encoding="utf-8",
        )
        self.controller_status.write_text("0\n", encoding="utf-8")
        self.source_identity.write_text(
            f"commit={SOURCE_REVISION}\nbranch=cp2k-integration\n",
            encoding="utf-8",
        )
        (self.current / "run_metadata.txt").parent.mkdir(parents=True)
        (self.current / "run_metadata.txt").write_text(
            "executable=/synthetic/tblite\n"
            f"executable_sha256={BINARY}\n"
            "meshes=1\n"
            "phases=Ih VII\n"
            "accuracy=0.1\n"
            "iterations=300\n",
            encoding="utf-8",
        )
        for index, phase in enumerate(PHASES):
            structure = self.archive / "structures" / "k111" / phase / "POSCAR"
            structure.parent.mkdir(parents=True)
            structure.write_text(
                f"{phase}\n1.0\n1 0 0\n0 1 0\n0 0 1\nO H\n1 2\nCartesian\n"
                "0 0 0\n0.8 0 0\n0 0.8 0\n",
                encoding="utf-8",
            )
            archived_energy = -10.0 + 0.1 * index
            archived_json = (
                self.archive
                / "results"
                / "current_save_tblite_cli"
                / "k111"
                / phase
                / "result.json"
            )
            archived_json.parent.mkdir(parents=True)
            archived_json.write_text(
                json.dumps({"energy": archived_energy}) + "\n", encoding="utf-8"
            )
            native_output = (
                self.archive
                / "results"
                / "current_cp2k_native"
                / "k111"
                / phase
                / "cp2k.out"
            )
            native_output.parent.mkdir(parents=True)
            native_output.write_text(
                " ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.): "
                f"{archived_energy + 1.0e-8:.15f}\n"
                " PROGRAM ENDED AT synthetic\n",
                encoding="utf-8",
            )
            current_dir = self.current / "k111" / phase
            current_dir.mkdir(parents=True)
            (current_dir / "result.json").write_text(
                json.dumps({"energy": archived_energy}) + "\n", encoding="utf-8"
            )
            (current_dir / "process.out").write_text(
                "total energy -1.0\nJSON dump of results written\n",
                encoding="utf-8",
            )
            (current_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (current_dir / "binary.sha256").write_text(
                f"{BINARY}  /synthetic/tblite\n", encoding="utf-8"
            )
            (current_dir / "input.sha256").write_text(
                f"{digest(structure)}  {structure}\n", encoding="utf-8"
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_verifier(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERIFY),
                str(self.current),
                "--archive-root",
                str(self.archive),
                "--require-binary-sha256",
                BINARY,
                "--source-identity",
                str(self.source_identity),
                "--require-source-revision",
                SOURCE_REVISION,
                "--affinity-proof",
                str(self.affinity),
                "--controller-exit-status",
                str(self.controller_status),
                "--phases",
                ",".join(PHASES),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_complete_repetition_passes(self) -> None:
        completed = self.run_verifier()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["coverage"], 2)
        self.assertEqual(payload["provenance"]["affinity"]["allowed"], "74")

    def test_changed_structure_is_rejected(self) -> None:
        structure = self.archive / "structures" / "k111" / "VII" / "POSCAR"
        structure.write_text(structure.read_text() + "\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("input hash mismatch", completed.stderr)

    def test_nonzero_phase_exit_is_rejected(self) -> None:
        (self.current / "k111" / "VII" / "exit_status").write_text(
            "9\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("nonzero direct CLI exit status", completed.stderr)

    def test_non_singleton_affinity_is_rejected(self) -> None:
        self.affinity.write_text(
            "pid=1 expected_cpu=74 allowed=74-75\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("invalid singleton affinity proof", completed.stderr)

    def test_wrong_source_revision_is_rejected(self) -> None:
        self.source_identity.write_text(
            f"commit={'c' * 40}\nbranch=cp2k-integration\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("source revision mismatch", completed.stderr)


if __name__ == "__main__":
    unittest.main()
