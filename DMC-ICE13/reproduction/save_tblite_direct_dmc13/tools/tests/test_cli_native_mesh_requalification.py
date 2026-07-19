#!/usr/bin/env python3
"""Tests for the mesh-independent direct/native qualification gate."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_cli_native_mesh_requalification.py"
PHASES = ("Ih", "VII")
DIRECT_BINARY = "a" * 64
NATIVE_BINARY = "b" * 64
SOURCE_REVISION = "c" * 40
NATIVE_CP2K_REVISION = "d" * 40
NATIVE_PROVIDER_ARCHIVE = "e" * 64
NATIVE_CMAKE_CACHE = "f" * 64
NATIVE_BUILD_NINJA = "1" * 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class CliNativeMeshRequalificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = self.root / "archive"
        self.direct = self.root / "direct"
        self.native_runs = self.root / "native_runs"
        self.native_inputs = self.root / "native_inputs"
        self.direct_status = self.root / "direct_status"
        self.native_status = self.root / "native_status"
        self.source = self.root / "source_identity.txt"
        self.direct_status.write_text("0\n", encoding="utf-8")
        self.native_status.write_text("0\n", encoding="utf-8")
        self.source.write_text(
            f"commit={SOURCE_REVISION}\n"
            f"executable_sha256={DIRECT_BINARY}\n"
            f"native_provider_commit={SOURCE_REVISION}\n"
            f"native_provider_archive_sha256={NATIVE_PROVIDER_ARCHIVE}\n"
            f"native_cp2k_commit={NATIVE_CP2K_REVISION}\n"
            f"native_cp2k_binary_sha256={NATIVE_BINARY}\n"
            "native_cmake_provider=SAVE\n"
            f"native_cmake_provider_revision={SOURCE_REVISION}\n"
            f"native_cmake_cache_sha256={NATIVE_CMAKE_CACHE}\n"
            f"native_build_ninja_sha256={NATIVE_BUILD_NINJA}\n",
            encoding="utf-8",
        )
        primitive_positions = (
            ("H", (0.10, 0.20, 0.30)),
            ("H", (0.20, 0.30, 0.40)),
            ("O", (0.15, 0.25, 0.35)),
        )
        for index, phase in enumerate(PHASES):
            structure = self.archive / "structures" / "k333" / phase / "POSCAR"
            structure.parent.mkdir(parents=True)
            positions: list[str] = []
            for element in ("H", "O"):
                for primitive_element, fractional in primitive_positions:
                    if primitive_element != element:
                        continue
                    for iz in range(3):
                        for iy in range(3):
                            for ix in range(3):
                                positions.append(
                                    " ".join(
                                        f"{3.0 * (fractional[axis] + shift):.12f}"
                                        for axis, shift in enumerate((ix, iy, iz))
                                    )
                                )
            structure.write_text(
                "test\n1.0\n9 0 0\n0 9 0\n0 0 9\nH O\n54 27\nCartesian\n"
                + "\n".join(positions)
                + "\n",
                encoding="utf-8",
            )

            primitive_energy = -100.0 - 0.01 * index
            direct_dir = self.direct / "k333" / phase
            direct_dir.mkdir(parents=True)
            (direct_dir / "result.json").write_text(
                json.dumps({"energy": 27.0 * primitive_energy}) + "\n",
                encoding="utf-8",
            )
            (direct_dir / "process.out").write_text(
                "dispersion energy              2.7000000000000E-01 Eh\n"
                "total energy\nJSON dump of results written\n",
                encoding="utf-8",
            )
            (direct_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (direct_dir / "binary.sha256").write_text(
                f"{DIRECT_BINARY}  /build/tblite\n", encoding="utf-8"
            )
            (direct_dir / "input.sha256").write_text(
                f"{digest(structure)}  {structure}\n", encoding="utf-8"
            )
            direct_cpu = 72 + index
            (direct_dir / "affinity_preexec.txt").write_text(
                f"pid=42 expected_cpu={direct_cpu} allowed={direct_cpu}\n",
                encoding="utf-8",
            )

            native_input = self.native_inputs / phase / "input.inp"
            native_input.parent.mkdir(parents=True)
            native_input.write_text(
                "&KPOINTS\n"
                "  SCHEME MACDONALD 3 3 3 0 0 0\n"
                "  SYMMETRY T\n"
                "  FULL_GRID F\n"
                "&END KPOINTS\n"
                "&SUBSYS\n"
                "  &CELL\n"
                "    PERIODIC XYZ\n"
                "    A 3 0 0\n"
                "    B 0 3 0\n"
                "    C 0 0 3\n"
                "  &END CELL\n"
                "  &COORD\n"
                "    SCALED\n"
                "    H 0.10 0.20 0.30\n"
                "    H 0.20 0.30 0.40\n"
                "    O 0.15 0.25 0.35\n"
                "  &END COORD\n"
                "&END SUBSYS\n",
                encoding="utf-8",
            )
            native_dir = self.native_runs / phase
            native_dir.mkdir(parents=True)
            delta = (index + 1) * 1.0e-10
            (native_dir / "cp2k.out").write_text(
                " Non-self consistent dispersion energy: "
                "0.01000000000000\n"
                " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] "
                f"{primitive_energy + delta:.15f}\nPROGRAM ENDED AT now\n",
                encoding="utf-8",
            )
            (native_dir / "exit_status").write_text("0\n", encoding="utf-8")
            (native_dir / "binary.sha256").write_text(
                f"{NATIVE_BINARY}  /build/cp2k.psmp\n", encoding="utf-8"
            )
            (native_dir / "input.sha256").write_text(
                f"{digest(native_input)}  {native_input}\n", encoding="utf-8"
            )
            native_cpu = 80 + index
            (native_dir / "affinity_preexec.txt").write_text(
                f"pid=84 expected_cpu={native_cpu} allowed={native_cpu}\n",
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
                "--mesh-size",
                "3",
                "--phases",
                ",".join(PHASES),
                "--direct-controller-status",
                str(self.direct_status),
                "--native-controller-status",
                str(self.native_status),
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
                "--expected-direct-cpus",
                "72,73",
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_three_mesh_pair_passes(self) -> None:
        completed = self.run_verifier()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["mesh"], "3x3x3")
        self.assertEqual(payload["replicas"], 27)
        self.assertEqual(payload["phase_count"], 2)

    def test_even_mesh_shift_formula(self) -> None:
        module_root = str(TOOLS)
        sys.path.insert(0, module_root)
        try:
            from verify_cli_native_mesh_requalification import macdonald_shift

            self.assertEqual(macdonald_shift(3), 0.0)
            self.assertEqual(macdonald_shift(4), 0.375)
            self.assertEqual(macdonald_shift(8), 0.4375)
        finally:
            sys.path.remove(module_root)

    def test_noncanonical_shift_is_rejected(self) -> None:
        path = self.native_inputs / "VII" / "input.inp"
        path.write_text(
            path.read_text(encoding="utf-8").replace("3 3 3 0 0 0", "3 3 3 0.1 0 0"),
            encoding="utf-8",
        )
        (self.native_runs / "VII" / "input.sha256").write_text(
            f"{digest(path)}  {path}\n", encoding="utf-8"
        )
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("noncanonical 3x3x3 MacDonald mesh", completed.stderr)

    def test_unlisted_direct_cpu_is_rejected(self) -> None:
        affinity = self.direct / "k333" / "VII" / "affinity_preexec.txt"
        affinity.write_text("pid=42 expected_cpu=74 allowed=74\n", encoding="utf-8")
        completed = self.run_verifier()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("unexpected direct CPU", completed.stderr)


if __name__ == "__main__":
    unittest.main()
