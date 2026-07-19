#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
SCRIPT = TOOLS / "verify_cli_native_sentinels.py"
DIRECT_BINARY = "a" * 64
NATIVE_BINARY = "b" * 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_run_metadata(run: Path, binary: str, input_path: Path, cpu: int) -> None:
    run.mkdir(parents=True)
    (run / "exit_status").write_text("0\n")
    (run / "binary.sha256").write_text(f"{binary}  executable\n")
    (run / "input.sha256").write_text(f"{digest(input_path)}  {input_path}\n")
    (run / "affinity_preexec.txt").write_text(
        f"pid=1 expected_cpu={cpu} allowed={cpu}\n"
        f"Cpus_allowed_list:\t{cpu}\n"
    )


class CliNativeSentinelTest(unittest.TestCase):
    def test_passes_exact_hash_and_relative_energy_gates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mesh = 2
            for index, phase in enumerate(("Ih", "II")):
                direct_input = root / "direct-input" / phase / "POSCAR"
                native_input = root / "native-input" / phase / "input.inp"
                direct_input.parent.mkdir(parents=True)
                native_input.parent.mkdir(parents=True)
                direct_input.write_text(
                    "test\n1.0\n1 0 0\n0 1 0\n0 0 1\nH O\n16 8\nDirect\n"
                )
                native_input.write_text(
                    "&SUBSYS\n&COORD\nSCALED\nH 0 0 0\nH 0 0 1\nO 0 0 0.5\n&END COORD\n&END SUBSYS\n"
                )
                direct_run = root / "direct" / phase
                native_run = root / "native" / phase
                write_run_metadata(direct_run, DIRECT_BINARY, direct_input, 10 + index)
                write_run_metadata(native_run, NATIVE_BINARY, native_input, 20 + index)
                native_energy = -10.0 + index * 0.01
                direct_energy = (native_energy - (index + 1) * 1.0e-9) * mesh**3
                (direct_run / "tblite.json").write_text(
                    json.dumps({"energy": direct_energy}) + "\n"
                )
                (native_run / "cp2k.out").write_text(
                    f" ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] {native_energy:.16f}\n"
                    " PROGRAM ENDED AT test\n"
                )

            output = root / "result.json"
            command = [
                "python3",
                str(SCRIPT),
                "--mesh",
                str(mesh),
                "--phase",
                "Ih",
                "--phase",
                "II",
                "--direct-root",
                str(root / "direct"),
                "--native-root",
                str(root / "native"),
                "--direct-input-root",
                str(root / "direct-input"),
                "--native-input-root",
                str(root / "native-input"),
                "--direct-binary-sha256",
                DIRECT_BINARY,
                "--native-binary-sha256",
                NATIVE_BINARY,
                "--output",
                str(output),
            ]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(output.read_text())
            self.assertEqual(payload["status"], "pass")
            self.assertEqual(payload["phases"], ["Ih", "II"])
            self.assertLess(payload["summary"]["max_abs_native_minus_direct_Ha"], 3e-9)

    def test_rejects_binary_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            direct_input = root / "direct-input/Ih/POSCAR"
            native_input = root / "native-input/Ih/input.inp"
            direct_input.parent.mkdir(parents=True)
            native_input.parent.mkdir(parents=True)
            direct_input.write_text("test\n1\n1 0 0\n0 1 0\n0 0 1\nH O\n2 1\nDirect\n")
            native_input.write_text("&COORD\nH 0 0 0\nH 0 0 1\nO 0 0 0.5\n&END COORD\n")
            write_run_metadata(root / "direct/Ih", "c" * 64, direct_input, 10)
            write_run_metadata(root / "native/Ih", NATIVE_BINARY, native_input, 20)
            (root / "direct/Ih/tblite.json").write_text('{"energy": -1.0}\n')
            (root / "native/Ih/cp2k.out").write_text(
                " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] -1.0\n PROGRAM ENDED AT test\n"
            )
            result = subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--mesh",
                    "1",
                    "--phase",
                    "Ih",
                    "--direct-root",
                    str(root / "direct"),
                    "--native-root",
                    str(root / "native"),
                    "--direct-input-root",
                    str(root / "direct-input"),
                    "--native-input-root",
                    str(root / "native-input"),
                    "--direct-binary-sha256",
                    DIRECT_BINARY,
                    "--native-binary-sha256",
                    NATIVE_BINARY,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("direct binary mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
