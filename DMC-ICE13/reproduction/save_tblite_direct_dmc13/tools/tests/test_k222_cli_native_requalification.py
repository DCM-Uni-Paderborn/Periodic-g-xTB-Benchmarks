#!/usr/bin/env python3
"""Qualification tests for the fresh all-phase 2x2x2 parity gate."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_k222_cli_native_requalification.py"
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
DIRECT_BINARY = "a" * 64
NATIVE_BINARY = "b" * 64
SOURCE_REVISION = "c" * 40
NATIVE_CP2K_REVISION = "d" * 40
NATIVE_PROVIDER_ARCHIVE = "e" * 64
NATIVE_CMAKE_CACHE = "f" * 64
NATIVE_BUILD_NINJA = "1" * 64
CPU = 72


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class K222CliNativeRequalificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = self.root / "archive"
        self.direct = self.root / "direct"
        self.native_runs = self.root / "native_runs"
        self.native_inputs = self.root / "native_inputs"
        self.controller = self.root / "controller_exit_status"
        self.source = self.root / "source_identity.txt"
        self.controller.write_text("0\n", encoding="utf-8")
        self.source.write_text(
            f"repository=https://example.invalid/provider.git\n"
            f"branch=cp2k-integration\ncommit={SOURCE_REVISION}\n"
            f"executable_sha256={DIRECT_BINARY}\n"
            f"native_provider_commit={SOURCE_REVISION}\n"
            f"native_provider_archive_sha256={NATIVE_PROVIDER_ARCHIVE}\n"
            f"native_cp2k_commit={NATIVE_CP2K_REVISION}\n"
            f"native_cp2k_binary_sha256={NATIVE_BINARY}\n"
            f"native_cmake_provider=SAVE\n"
            f"native_cmake_provider_revision={SOURCE_REVISION}\n"
            f"native_cmake_cache_sha256={NATIVE_CMAKE_CACHE}\n"
            f"native_build_ninja_sha256={NATIVE_BUILD_NINJA}\n",
            encoding="utf-8",
        )
        for index, phase in enumerate(PHASES):
            structure = self.archive / "structures" / "k222" / phase / "POSCAR"
            structure.parent.mkdir(parents=True)
            structure.write_text(
                "test\n1.0\n8 0 0\n0 8 0\n0 0 8\nH O\n16 8\nCartesian\n",
                encoding="utf-8",
            )

            direct_dir = self.direct / "k222" / phase
            direct_dir.mkdir(parents=True)
            primitive = -100.0 - 0.01 * index
            (direct_dir / "result.json").write_text(
                json.dumps({"energy": 8.0 * primitive}) + "\n",
                encoding="utf-8",
            )
            (direct_dir / "process.out").write_text(
                "total energy\nJSON dump of results written\n", encoding="utf-8"
            )
            (direct_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (direct_dir / "binary.sha256").write_text(
                f"{DIRECT_BINARY}  /build/tblite\n", encoding="utf-8"
            )
            (direct_dir / "input.sha256").write_text(
                f"{digest(structure)}  {structure}\n", encoding="utf-8"
            )
            (direct_dir / "affinity_preexec.txt").write_text(
                f"pid=42 expected_cpu={CPU} allowed={CPU}\n",
                encoding="utf-8",
            )

            native_input = self.native_inputs / phase / "input.inp"
            native_input.parent.mkdir(parents=True)
            native_input.write_text(
                "&KPOINTS\n"
                "  SCHEME MACDONALD 2 2 2 0.25 0.25 0.25\n"
                "  SYMMETRY T\n"
                "  FULL_GRID F\n"
                "&END KPOINTS\n",
                encoding="utf-8",
            )
            native_dir = self.native_runs / phase
            native_dir.mkdir(parents=True)
            delta = (index - 6) * 1.0e-10
            (native_dir / "cp2k.out").write_text(
                " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] "
                f"{primitive + delta:.15f}\nPROGRAM ENDED AT now\n",
                encoding="utf-8",
            )
            (native_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (native_dir / "binary.sha256").write_text(
                f"{NATIVE_BINARY}  /build/cp2k.psmp\n", encoding="utf-8"
            )
            (native_dir / "input.sha256").write_text(
                f"{digest(native_input)}  {native_input}\n", encoding="utf-8"
            )
            (native_dir / "affinity_preexec.txt").write_text(
                f"pid=84 expected_cpu={80 + index} allowed={80 + index}\n",
                encoding="utf-8",
            )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_verifier(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERIFY),
                str(self.archive),
                str(self.direct),
                str(self.native_runs),
                str(self.native_inputs),
                "--direct-controller-status",
                str(self.controller),
                "--source-identity",
                str(self.source),
                "--expected-source-revision",
                SOURCE_REVISION,
                "--expected-direct-binary",
                DIRECT_BINARY,
                "--expected-native-binary",
                NATIVE_BINARY,
                "--expected-native-provider-archive",
                NATIVE_PROVIDER_ARCHIVE,
                "--expected-native-cp2k-revision",
                NATIVE_CP2K_REVISION,
                "--expected-native-cmake-cache",
                NATIVE_CMAKE_CACHE,
                "--expected-native-build-ninja",
                NATIVE_BUILD_NINJA,
                "--expected-direct-cpu",
                str(CPU),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_complete_fresh_pair_passes(self) -> None:
        completed = self.run_verifier()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["phase_count"], 13)
        self.assertLess(
            payload["statistics"]["max_abs_native_minus_direct_Ha"], 1.0e-8
        )

    def test_changed_direct_binary_is_rejected(self) -> None:
        manifest = self.direct / "k222" / "VII" / "binary.sha256"
        manifest.write_text(f"{'d' * 64}  /build/tblite\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("direct binary mismatch", completed.stderr)

    def test_noncanonical_native_mesh_is_rejected(self) -> None:
        native_input = self.native_inputs / "XVII" / "input.inp"
        native_input.write_text(
            native_input.read_text(encoding="utf-8").replace("0.25", "0.0"),
            encoding="utf-8",
        )
        native_manifest = self.native_runs / "XVII" / "input.sha256"
        native_manifest.write_text(
            f"{digest(native_input)}  {native_input}\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("noncanonical 2x2x2 MacDonald mesh", completed.stderr)

    def test_non_singleton_native_affinity_is_rejected(self) -> None:
        affinity = self.native_runs / "II" / "affinity_preexec.txt"
        affinity.write_text(
            "pid=84 expected_cpu=81 allowed=80-81\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("non-singleton or wrong affinity", completed.stderr)

    def test_source_identity_binary_mismatch_is_rejected(self) -> None:
        self.source.write_text(
            self.source.read_text(encoding="utf-8").replace(
                DIRECT_BINARY, "e" * 64
            ),
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("direct provider executable mismatch", completed.stderr)

    def test_native_provider_revision_mismatch_is_rejected(self) -> None:
        self.source.write_text(
            self.source.read_text(encoding="utf-8").replace(
                f"native_provider_commit={SOURCE_REVISION}",
                f"native_provider_commit={'2' * 40}",
            ),
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("native provider revision mismatch", completed.stderr)

    def test_native_link_plan_mismatch_is_rejected(self) -> None:
        self.source.write_text(
            self.source.read_text(encoding="utf-8").replace(
                f"native_build_ninja_sha256={NATIVE_BUILD_NINJA}",
                f"native_build_ninja_sha256={'3' * 64}",
            ),
            encoding="utf-8",
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("native link plan mismatch", completed.stderr)


if __name__ == "__main__":
    unittest.main()
