#!/usr/bin/env python3
"""Dry-run the adaptive Terok controller with a synthetic pinned launcher."""

from __future__ import annotations

import hashlib
import io
import json
import os
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
CONTROLLER = HERE / "run_strict_adaptive_completion.sh"
PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ALL = PHASES + ("Ih",)


def executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip(), encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


class AdaptiveControllerTest(unittest.TestCase):
    def test_priority_order_and_phase_local_pruning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "tools").mkdir()
            (root / "status").mkdir()
            events = root / "events.log"
            binary = root / "cp2k.psmp"
            binary.write_bytes(b"synthetic-qualified-cp2k\n")
            digest = hashlib.sha256(binary.read_bytes()).hexdigest()
            (root / "tools/dmc_ice13_relative_energies.csv").write_text(
                "phase,DMC_relative_kJmol\n", encoding="utf-8"
            )
            for phase in ALL:
                directory = root / "inputs/k444-reduced" / phase
                directory.mkdir(parents=True)
                (directory / "input.inp").write_text(
                    "SCHEME MACDONALD 4 4 4 0.375 0.375 0.375\n",
                    encoding="utf-8",
                )

            fake_bin = root / "fake-bin"
            fake_bin.mkdir()
            executable(
                fake_bin / "date",
                r"""
                #!/usr/bin/env bash
                printf '2026-07-19T12:00:00+02:00\n'
                """,
            )
            executable(
                root / "tools/build_native_mesh_input.py",
                r"""
                #!/usr/bin/env python3
                import json, pathlib, sys
                source, target, mesh = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
                provenance = pathlib.Path(sys.argv[sys.argv.index('--provenance') + 1])
                shift = 0.0 if mesh % 2 else (mesh - 1) / (2 * mesh)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(f'SCHEME MACDONALD {mesh} {mesh} {mesh} {shift} {shift} {shift}\n')
                provenance.write_text(json.dumps({'source': str(source), 'mesh': mesh}) + '\n')
                """,
            )
            executable(
                root / "launch_pinned_cp2k.sh",
                r"""
                #!/usr/bin/env bash
                set -euo pipefail
                job=$1 cpu=$2 binary=$3 input=$4 result=$5
                printf 'launch %s cpu=%s\n' "$job" "$cpu" >>"$EVENTS"
                mkdir -p "$result"
                printf ' ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.0\n PROGRAM ENDED AT synthetic\n' >"$result/cp2k.out"
                printf '0\n' >"$result/exit_status"
                sha256sum "$binary" >"$result/binary.sha256"
                sha256sum "$input" >"$result/input.sha256"
                printf 'pid=1 expected_cpu=%s allowed=%s\n' "$cpu" "$cpu" >"$result/affinity_preexec.txt"
                """,
            )
            executable(
                root / "tools/dmc_phase_convergence.py",
                r"""
                #!/usr/bin/env python3
                import pathlib, sys
                root, previous, current, phase = pathlib.Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
                threshold = sys.argv[sys.argv.index('--threshold') + 1]
                if threshold != '0.10':
                    raise SystemExit(2)
                endpoints = {'II': 5, 'III': 6, 'IV': 7}
                endpoint = endpoints.get(phase, 8)
                print(f'{phase}\t{previous}\t{current}\tthreshold={threshold}')
                raise SystemExit(0 if current >= endpoint else 1)
                """,
            )
            executable(
                root / "tools/select_adaptive_endpoints.py",
                r"""
                #!/usr/bin/env python3
                import json, pathlib, sys
                if sys.argv[sys.argv.index('--threshold') + 1] != '0.10':
                    raise SystemExit(2)
                output_json = pathlib.Path(sys.argv[sys.argv.index('--output-json') + 1])
                output_csv = pathlib.Path(sys.argv[sys.argv.index('--output-csv') + 1])
                output_json.write_text(json.dumps({'complete': True, 'status': 'PASS'}) + '\n')
                output_csv.write_text('phase,status\n')
                print(json.dumps({'complete': True, 'status': 'PASS'}))
                """,
            )
            executable(
                root / "run_gamma_supercell_oracle.sh",
                r"""
                #!/usr/bin/env bash
                set -euo pipefail
                for phase in VII Ih; do
                  test "$(tr -d '\n' <"$DMC_ROOT/runs/k888-reduced/$phase/exit_status")" = 0
                  test "$(awk 'NR==1{print $1}' "$DMC_ROOT/runs/k888-reduced/$phase/binary.sha256")" = "$REQUIRED_BINARY_SHA256"
                done
                printf 'oracle\n' >>"$EVENTS"
                """,
            )

            environment = os.environ.copy()
            environment.update(
                {
                    "AVAILABLE_GIB_OVERRIDE": "999",
                    "CONVERGENCE_THRESHOLD": "0.10",
                    "CP2K_BINARY": str(binary),
                    "DMC_ROOT": str(root),
                    "EVENTS": str(events),
                    "GAMMA_ORACLE_CONTROLLER": str(root / "run_gamma_supercell_oracle.sh"),
                    "MAXIMUM_MESH": "8",
                    "MINIMUM_AVAILABLE_GIB": "0",
                    "PATH": f"{fake_bin}:{environment['PATH']}",
                    "PINNED_LAUNCHER": str(root / "launch_pinned_cp2k.sh"),
                    "REQUIRED_BINARY_SHA256": digest,
                }
            )
            completed = subprocess.run(
                [str(CONTROLLER)],
                text=True,
                capture_output=True,
                check=False,
                env=environment,
                timeout=30,
            )
            self.assertEqual(
                completed.returncode,
                0,
                msg=f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
            )
            lines = events.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[:3], [
                "launch dmc-k888-strict-VII cpu=91",
                "launch dmc-k888-strict-Ih cpu=85",
                "oracle",
            ])
            launched = {line.split()[1] for line in lines if line.startswith("launch ")}
            self.assertIn("dmc-k555-strict-II", launched)
            self.assertNotIn("dmc-k666-strict-II", launched)
            self.assertIn("dmc-k666-strict-III", launched)
            self.assertNotIn("dmc-k777-strict-III", launched)
            self.assertIn("dmc-k777-strict-IV", launched)
            self.assertNotIn("dmc-k888-strict-IV", launched)
            self.assertEqual(
                (root / "status/strict-adaptive-completion.status")
                .read_text(encoding="utf-8")
                .splitlines()[0],
                f"required_binary_sha256={digest}",
            )
            self.assertIn(
                "status=PASS",
                (root / "status/strict-adaptive-completion.status").read_text(),
            )


if __name__ == "__main__":
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(AdaptiveControllerTest)
    )
    print(
        json.dumps(
            {
                "errors": len(result.errors),
                "failures": len(result.failures),
                "status": "PASS" if result.wasSuccessful() else "FAIL",
                "tests_run": result.testsRun,
            },
            sort_keys=True,
        )
    )
    if not result.wasSuccessful():
        print(stream.getvalue(), file=sys.stderr)
        raise SystemExit(1)
