#!/usr/bin/env python3
"""End-to-end qualification tests for the direct save_tblite driver."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


SOURCE_DRIVER = Path(__file__).resolve().parents[2] / "run_save_tblite.sh"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class DirectCliDriverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.driver = self.root / "run_save_tblite.sh"
        shutil.copy2(SOURCE_DRIVER, self.driver)
        self.structure = self.root / "structures" / "k111" / "Ih" / "POSCAR"
        self.structure.parent.mkdir(parents=True)
        self.structure.write_text("synthetic structure\n", encoding="utf-8")
        self.executable = self.root / "tblite"
        self.executable.write_text(
            """#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then
  printf 'qualified test executable\\n'
  exit 0
fi
printf 'run\\n' >> "${FAKE_COUNT:?}"
if [[ ${FAIL_RUN:-0} == 1 ]]; then
  printf 'intentional failure\\n' >&2
  exit 7
fi
json=''
while [[ $# -gt 0 ]]; do
  if [[ $1 == --json ]]; then
    json=$2
    shift 2
  else
    shift
  fi
done
printf 'total energy      -1.000000000000\\n'
if [[ ${NO_JSON:-0} != 1 ]]; then
  printf '{"energy": -1.0}\\n' > "$json"
  printf 'JSON dump of results written\\n'
fi
""",
            encoding="utf-8",
        )
        self.executable.chmod(0o755)
        self.result_root = self.root / "results"
        self.count = self.root / "count.txt"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_driver(self, **overrides: str) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment.update(
            {
                "TBLITE_EXE": str(self.executable),
                "MESHES": "1",
                "PHASES": "Ih",
                "RESULT_ROOT": str(self.result_root),
                "FAKE_COUNT": str(self.count),
            }
        )
        environment.update(overrides)
        return subprocess.run(
            ["bash", str(self.driver)],
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )

    def output(self) -> Path:
        return self.result_root / "k111" / "Ih"

    def test_success_records_complete_provenance_and_safe_resume(self) -> None:
        completed = self.run_driver()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        output = self.output()
        self.assertEqual((output / "exit_status").read_text().strip(), "0")
        self.assertEqual(
            (output / "binary.sha256").read_text().split()[0],
            digest(self.executable),
        )
        self.assertEqual(
            (output / "input.sha256").read_text().split()[0],
            digest(self.structure),
        )
        self.assertIn("total energy", (output / "process.out").read_text())
        self.assertIn(
            "JSON dump of results written", (output / "process.out").read_text()
        )
        for line in (output / "SHA256SUMS").read_text().splitlines():
            expected, path = line.split(maxsplit=1)
            self.assertEqual(digest(Path(path)), expected)
        self.assertEqual(self.count.read_text().splitlines(), ["run"])

        resumed = self.run_driver(SKIP_EXISTING="1")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(self.count.read_text().splitlines(), ["run"])

        self.structure.write_text("changed structure\n", encoding="utf-8")
        repeated = self.run_driver(SKIP_EXISTING="1")
        self.assertEqual(repeated.returncode, 0, repeated.stderr)
        self.assertEqual(self.count.read_text().splitlines(), ["run", "run"])

    def test_failure_status_is_preserved(self) -> None:
        completed = self.run_driver(FAIL_RUN="1")
        self.assertEqual(completed.returncode, 7)
        self.assertEqual((self.output() / "exit_status").read_text().strip(), "7")
        self.assertFalse((self.output() / "result.json").exists())

    def test_missing_json_is_a_failed_result(self) -> None:
        completed = self.run_driver(NO_JSON="1")
        self.assertEqual(completed.returncode, 66)
        self.assertEqual((self.output() / "exit_status").read_text().strip(), "66")


if __name__ == "__main__":
    unittest.main()
