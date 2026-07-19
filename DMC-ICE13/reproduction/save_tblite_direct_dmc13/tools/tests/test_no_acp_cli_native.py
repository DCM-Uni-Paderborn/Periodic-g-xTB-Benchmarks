#!/usr/bin/env python3
"""Qualification tests for the No-ACP direct-CLI/native-CP2K gate."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_no_acp_cli_native.py"
DIRECT_BINARY = "a" * 64
NATIVE_BINARY = "b" * 64
SOURCE_REVISION = "c" * 40
NO_ACP_PARAMETER_TEXT = (
    "[hamiltonian.xtb]\n"
    "[element.H.acp]\nacp_levels = [-0.1, -0.2, -0.3, -0.4]\n"
    "[element.O.acp]\nacp_levels = [-0.5, -0.6, -0.7, -0.8]\n"
)
FULL_PARAMETER_TEXT = NO_ACP_PARAMETER_TEXT.replace(
    "[element.H.acp]", "[acp]\n[element.H.acp]"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class NoAcpCliNativeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.direct = self.root / "direct"
        self.native = self.root / "native"
        self.structures = self.root / "structures"
        self.parameter = self.root / "gxtb_no_acp.toml"
        self.full_parameter = self.root / "gxtb_full.toml"
        self.controller_status = self.root / "controller_exit_status"
        self.source_identity = self.root / "source_identity.txt"
        self.parameter.write_text(NO_ACP_PARAMETER_TEXT, encoding="utf-8")
        self.full_parameter.write_text(FULL_PARAMETER_TEXT, encoding="utf-8")
        self.controller_status.write_text("0\n", encoding="utf-8")
        self.source_identity.write_text(
            f"commit={SOURCE_REVISION}\nbranch=cp2k-integration\n",
            encoding="utf-8",
        )

        primitive = {"Ih": -10.0, "XVII": -9.9}
        for index, phase in enumerate(("Ih", "XVII")):
            structure = self.structures / phase / "POSCAR"
            structure.parent.mkdir(parents=True)
            structure.write_text(
                f"{phase}\n1.0\n2 0 0\n0 2 0\n0 0 2\nO H\n8 16\nCartesian\n"
                + "\n".join(f"{i / 10:.1f} 0 0" for i in range(24))
                + "\n",
                encoding="utf-8",
            )
            direct_dir = self.direct / "k222" / phase
            direct_dir.mkdir(parents=True)
            direct_energy = 8.0 * primitive[phase]
            (direct_dir / "result.json").write_text(
                json.dumps({"energy": direct_energy}) + "\n", encoding="utf-8"
            )
            (direct_dir / "process.out").write_text(
                "total energy -1.0\nJSON dump of results written\n",
                encoding="utf-8",
            )
            (direct_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (direct_dir / "binary.sha256").write_text(
                f"{DIRECT_BINARY}  /synthetic/tblite\n", encoding="utf-8"
            )
            (direct_dir / "input.sha256").write_text(
                f"{digest(structure)}  {structure}\n", encoding="utf-8"
            )
            (direct_dir / "parameter.sha256").write_text(
                f"{digest(self.parameter)}  {self.parameter}\n", encoding="utf-8"
            )
            (direct_dir / "affinity_preexec.txt").write_text(
                f"pid={index + 1} expected_cpu=73 allowed=73\n"
                "Cpus_allowed_list:\t73\n",
                encoding="utf-8",
            )

            native_dir = self.native / phase
            native_dir.mkdir(parents=True)
            (native_dir / "cp2k.out").write_text(
                " ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.): "
                f"{primitive[phase] + 1.0e-8:.15f}\n"
                " PROGRAM ENDED AT synthetic\n",
                encoding="utf-8",
            )
            (native_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (native_dir / "binary.sha256").write_text(
                f"{NATIVE_BINARY}  /synthetic/cp2k.psmp\n", encoding="utf-8"
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_verifier(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERIFY),
                str(self.direct),
                "--native-root",
                str(self.native),
                "--structure-root",
                str(self.structures),
                "--parameter-file",
                str(self.parameter),
                "--full-parameter-file",
                str(self.full_parameter),
                "--controller-exit-status",
                str(self.controller_status),
                "--source-identity",
                str(self.source_identity),
                "--require-source-revision",
                SOURCE_REVISION,
                "--require-binary-sha256",
                DIRECT_BINARY,
                "--require-native-binary-sha256",
                NATIVE_BINARY,
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_complete_pair_passes(self) -> None:
        completed = self.run_verifier()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertLess(
            payload["statistics"]["max_abs_native_minus_direct_Ha"], 2.0e-7
        )
        self.assertEqual(payload["rows"][0]["affinity"]["allowed"], "73")

    def test_changed_structure_is_rejected(self) -> None:
        structure = self.structures / "XVII" / "POSCAR"
        structure.write_text(structure.read_text() + "\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("input hash mismatch", completed.stderr)

    def test_changed_parameter_is_rejected(self) -> None:
        self.parameter.write_text(
            NO_ACP_PARAMETER_TEXT + "# changed\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("parameter hash mismatch", completed.stderr)

    def test_global_acp_is_rejected_even_with_matching_hashes(self) -> None:
        self.parameter.write_text(
            FULL_PARAMETER_TEXT,
            encoding="utf-8",
        )
        parameter_hash = digest(self.parameter)
        for phase in ("Ih", "XVII"):
            manifest = self.direct / "k222" / phase / "parameter.sha256"
            manifest.write_text(
                f"{parameter_hash}  {self.parameter}\n", encoding="utf-8"
            )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("still activates the global ACP table", completed.stderr)

    def test_unrelated_parameter_change_is_rejected_with_matching_hashes(self) -> None:
        self.parameter.write_text(
            NO_ACP_PARAMETER_TEXT.replace("-0.8", "-0.9"), encoding="utf-8"
        )
        parameter_hash = digest(self.parameter)
        for phase in ("Ih", "XVII"):
            manifest = self.direct / "k222" / phase / "parameter.sha256"
            manifest.write_text(
                f"{parameter_hash}  {self.parameter}\n", encoding="utf-8"
            )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("changes content beyond", completed.stderr)

    def test_nonzero_controller_exit_is_rejected(self) -> None:
        self.controller_status.write_text("9\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("nonzero direct CLI controller", completed.stderr)

    def test_non_singleton_affinity_is_rejected(self) -> None:
        affinity = self.direct / "k222" / "XVII" / "affinity_preexec.txt"
        affinity.write_text(
            "pid=2 expected_cpu=73 allowed=73-74\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("invalid singleton affinity proof", completed.stderr)

    def test_wrong_native_binary_is_rejected(self) -> None:
        manifest = self.native / "Ih" / "binary.sha256"
        manifest.write_text(f"{'d' * 64}  /synthetic/cp2k.psmp\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("native CP2K binary hash mismatch", completed.stderr)

    def test_energy_mismatch_is_rejected(self) -> None:
        output = self.native / "XVII" / "cp2k.out"
        output.write_text(
            " ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.): -9.899000000000000\n"
            " PROGRAM ENDED AT synthetic\n",
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("No-ACP native/direct mismatch", completed.stderr)


if __name__ == "__main__":
    unittest.main()
