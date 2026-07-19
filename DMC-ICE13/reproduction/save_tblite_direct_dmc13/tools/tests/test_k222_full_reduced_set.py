#!/usr/bin/env python3
"""Tests for the all-phase full/reduced native 2x2x2 gate."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_k222_full_reduced_set.py"
PHASES = (
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
)
BINARY = "a" * 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class K222FullReducedSetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.reduced_runs = self.root / "reduced-runs"
        self.full_runs = self.root / "full-runs"
        self.reduced_inputs = self.root / "reduced-inputs"
        self.full_inputs = self.root / "full-inputs"
        for index, phase in enumerate(PHASES):
            reduced_input = self.reduced_inputs / phase / "input.inp"
            full_input = self.full_inputs / phase / "input.inp"
            reduced_input.parent.mkdir(parents=True)
            full_input.parent.mkdir(parents=True)
            common = (
                "&GLOBAL\n"
                f"  PROJECT phase_{phase}\n"
                "&END GLOBAL\n"
                "&KPOINTS\n"
                "  SCHEME MACDONALD 2 2 2 0.25 0.25 0.25\n"
            )
            reduced_input.write_text(
                common + "  SYMMETRY T\n  FULL_GRID F\n&END KPOINTS\n",
                encoding="utf-8",
            )
            full_input.write_text(
                common + "  SYMMETRY F\n  FULL_GRID T\n&END KPOINTS\n",
                encoding="utf-8",
            )
            energy = -100.0 - index * 0.01
            for run_root, input_path, cpu in (
                (self.reduced_runs, reduced_input, 80 + index),
                (self.full_runs, full_input, 100 + index),
            ):
                run = run_root / phase
                run.mkdir(parents=True)
                (run / "cp2k.out").write_text(
                    " Non-self consistent dispersion energy: 0.01250000000000\n"
                    " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] "
                    f"{energy:.15f}\nPROGRAM ENDED AT now\n",
                    encoding="utf-8",
                )
                (run / "exit_status").write_text("0\n", encoding="utf-8")
                (run / "binary.sha256").write_text(
                    f"{BINARY}  /build/cp2k.psmp\n", encoding="utf-8"
                )
                (run / "input.sha256").write_text(
                    f"{digest(input_path)}  {input_path}\n", encoding="utf-8"
                )
                (run / "affinity_preexec.txt").write_text(
                    f"pid=42 expected_cpu={cpu} allowed={cpu}\n", encoding="utf-8"
                )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_verifier(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERIFY),
                str(self.reduced_runs),
                str(self.full_runs),
                str(self.reduced_inputs),
                str(self.full_inputs),
                "--expected-binary",
                BINARY,
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_complete_set_passes(self) -> None:
        completed = self.run_verifier()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["phase_count"], 13)

    def test_energy_mismatch_is_rejected(self) -> None:
        output = self.full_runs / "VII" / "cp2k.out"
        output.write_text(
            output.read_text(encoding="utf-8").replace(
                "PROGRAM ENDED AT now",
                " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] "
                "-100.050000001000000\nPROGRAM ENDED AT now",
            ),
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("full/reduced energy mismatch", completed.stderr)

    def test_dispersion_mismatch_is_rejected(self) -> None:
        output = self.full_runs / "VIII" / "cp2k.out"
        output.write_text(
            output.read_text(encoding="utf-8").replace(
                "Non-self consistent dispersion energy: 0.01250000000000",
                "Non-self consistent dispersion energy: 0.01250000001000",
            ),
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("full/reduced dispersion mismatch", completed.stderr)

    def test_non_symmetry_input_change_is_rejected(self) -> None:
        input_path = self.full_inputs / "XIV" / "input.inp"
        input_path.write_text(
            input_path.read_text(encoding="utf-8") + "# unrelated change\n",
            encoding="utf-8",
        )
        manifest = self.full_runs / "XIV" / "input.sha256"
        manifest.write_text(
            f"{digest(input_path)}  {input_path}\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("non-symmetry input difference", completed.stderr)

    def test_non_singleton_affinity_is_rejected(self) -> None:
        affinity = self.full_runs / "II" / "affinity_preexec.txt"
        affinity.write_text(
            "pid=42 expected_cpu=101 allowed=101-102\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("non-singleton or wrong affinity", completed.stderr)

    def test_wrong_mesh_flags_are_rejected(self) -> None:
        input_path = self.full_inputs / "III" / "input.inp"
        input_path.write_text(
            input_path.read_text(encoding="utf-8").replace(
                "SYMMETRY F", "SYMMETRY T"
            ),
            encoding="utf-8",
        )
        manifest = self.full_runs / "III" / "input.sha256"
        manifest.write_text(
            f"{digest(input_path)}  {input_path}\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("wrong symmetry/full-grid flags", completed.stderr)

    def test_wrong_binary_is_rejected(self) -> None:
        manifest = self.reduced_runs / "XVII" / "binary.sha256"
        manifest.write_text(f"{'b' * 64}  /build/cp2k.psmp\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("binary mismatch", completed.stderr)


if __name__ == "__main__":
    unittest.main()
