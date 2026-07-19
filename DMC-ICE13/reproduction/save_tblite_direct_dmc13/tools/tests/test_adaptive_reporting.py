#!/usr/bin/env python3
"""End-to-end tests for qualified adaptive DMC-ICE13 reporting."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
PHASES = ("II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ALL_STRUCTURES = PHASES + ("Ih",)
CURRENT_DIGEST = "a" * 64
OLDER_DIGEST = "b" * 64
HARTREE_TO_KJMOL = 2625.4996394799


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class AdaptiveReportingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.reference = self.root / "reference.csv"
        with self.reference.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("phase", "DMC_relative_kJmol")
            )
            writer.writeheader()
            for index, phase in enumerate(PHASES, start=1):
                writer.writerow(
                    {"phase": phase, "DMC_relative_kJmol": f"{index / 10:.12f}"}
                )
        for mesh in (1, 2):
            for phase in ALL_STRUCTURES:
                self.write_structure(mesh, phase, CURRENT_DIGEST)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_structure(self, mesh: int, phase: str, digest: str) -> None:
        mesh_name = f"k{mesh}{mesh}{mesh}-reduced"
        input_dir = self.root / "inputs" / mesh_name / phase
        run_dir = self.root / "runs" / mesh_name / phase
        input_dir.mkdir(parents=True, exist_ok=True)
        run_dir.mkdir(parents=True, exist_ok=True)
        input_path = input_dir / "input.inp"
        shift = "0.0" if mesh % 2 else repr((mesh - 1) / (2 * mesh))
        input_path.write_text(
            "&FORCE_EVAL\n"
            "  &DFT\n"
            "    &KPOINTS\n"
            f"      SCHEME MACDONALD {mesh} {mesh} {mesh} {shift} {shift} {shift}\n"
            "    &END KPOINTS\n"
            "  &END DFT\n"
            "  &SUBSYS\n"
            "    &COORD\n"
            "      O 0.0 0.0 0.0\n"
            "      H 0.8 0.0 0.0\n"
            "      H 0.0 0.8 0.0\n"
            "    &END COORD\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n",
            encoding="utf-8",
        )
        relative = 0.0 if phase == "Ih" else (PHASES.index(phase) + 1) / 10
        energy = -10.0 + relative / HARTREE_TO_KJMOL
        (run_dir / "cp2k.out").write_text(
            f" ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy:.15f}\n"
            " PROGRAM ENDED AT synthetic-test\n",
            encoding="utf-8",
        )
        (run_dir / "exit_status").write_text("0\n", encoding="utf-8")
        (run_dir / "binary.sha256").write_text(
            f"{digest}  /synthetic/cp2k.psmp\n", encoding="utf-8"
        )
        (run_dir / "input.sha256").write_text(
            f"{sha256(input_path)}  {input_path}\n", encoding="utf-8"
        )

    def run_tool(
        self, name: str, *arguments: object, expected_returncode: int = 0
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, str(TOOLS / name), *(str(value) for value in arguments)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            expected_returncode,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    @staticmethod
    def mesh_vector(output: str) -> dict[str, int]:
        rows: dict[str, int] = {}
        for line in output.splitlines():
            fields = line.split("\t")
            if fields[0] in PHASES:
                rows[fields[0]] = int(fields[1])
        return rows

    def set_mesh_digest(self, mesh: int, digest: str) -> None:
        mesh_name = f"k{mesh}{mesh}{mesh}-reduced"
        for phase in ALL_STRUCTURES:
            (self.root / "runs" / mesh_name / phase / "binary.sha256").write_text(
                f"{digest}  /synthetic/cp2k.psmp\n", encoding="utf-8"
            )

    def test_selector_and_independent_verifier_agree(self) -> None:
        endpoints = self.root / "endpoints.json"
        self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--threshold",
            "0.05",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            "--output-json",
            endpoints,
        )
        verification = self.root / "independent-verification.json"
        self.run_tool(
            "verify_adaptive_dmc13.py",
            self.root,
            endpoints,
            self.reference,
            "--meshes",
            "1,2",
            "--threshold",
            "0.05",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            "--output-json",
            verification,
        )
        selected = json.loads(endpoints.read_text(encoding="utf-8"))
        verified = json.loads(verification.read_text(encoding="utf-8"))
        self.assertTrue(selected["complete"])
        self.assertEqual(verified["status"], "PASS")
        self.assertEqual({row["endpoint_mesh"] for row in verified["rows"]}, {2})
        self.assertAlmostEqual(
            selected["statistics"]["mae_kj_mol_per_water"],
            verified["statistics"]["mae_kj_mol_per_water"],
            places=14,
        )

    def test_default_threshold_is_user_selected_point_one(self) -> None:
        mesh_name = "k222-reduced"
        for phase in PHASES:
            run_dir = self.root / "runs" / mesh_name / phase
            relative = (PHASES.index(phase) + 1) / 10 + 0.075
            energy = -10.0 + relative / HARTREE_TO_KJMOL
            (run_dir / "cp2k.out").write_text(
                f" ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy:.15f}\n"
                " PROGRAM ENDED AT synthetic-test\n",
                encoding="utf-8",
            )

        selected = self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        payload = json.loads(selected.stdout)
        self.assertTrue(payload["complete"])
        self.assertEqual(payload["threshold_kj_mol_per_water"], 0.10)

        phase_default = self.run_tool(
            "dmc_phase_convergence.py",
            self.root,
            1,
            2,
            "II",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        self.assertIn("threshold_kj_mol=0.100000000000", phase_default.stdout)
        self.assertIn("status=converged", phase_default.stdout)

        phase_strict = self.run_tool(
            "dmc_phase_convergence.py",
            self.root,
            1,
            2,
            "II",
            "--threshold",
            "0.05",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            expected_returncode=1,
        )
        self.assertIn("status=unresolved", phase_strict.stdout)

        strict = self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--threshold",
            "0.05",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            expected_returncode=1,
        )
        strict_payload = json.loads(strict.stdout)
        self.assertEqual(strict_payload["unresolved_phase_count"], len(PHASES))

    def test_mixed_report_uses_highest_qualified_pair(self) -> None:
        self.set_mesh_digest(2, OLDER_DIGEST)
        current = self.run_tool(
            "dmc_mixed_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "2,1",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        older = self.run_tool(
            "dmc_mixed_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "2,1",
            "--require-binary-sha256",
            OLDER_DIGEST,
        )
        self.assertEqual(set(self.mesh_vector(current.stdout).values()), {1})
        self.assertEqual(set(self.mesh_vector(older.stdout).values()), {2})

    def test_mixed_report_skips_bad_input_hash(self) -> None:
        bad_hash = self.root / "runs" / "k222-reduced" / "II" / "input.sha256"
        bad_hash.write_text(f"{'0' * 64}  /synthetic/input.inp\n", encoding="utf-8")
        result = self.run_tool(
            "dmc_mixed_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "2,1",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        meshes = self.mesh_vector(result.stdout)
        self.assertEqual(meshes["II"], 1)
        self.assertEqual(meshes["III"], 2)

        endpoints = self.root / "bad-hash-endpoints.json"
        selection = self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            "--output-json",
            endpoints,
            expected_returncode=2,
        )
        self.assertIn("input hash mismatch", selection.stdout)

    def test_wrong_macdonald_shift_is_not_qualified(self) -> None:
        input_path = self.root / "inputs" / "k222-reduced" / "II" / "input.inp"
        input_path.write_text(
            input_path.read_text(encoding="utf-8").replace(
                "0.25 0.25 0.25", "0.375 0.375 0.375"
            ),
            encoding="utf-8",
        )
        run_dir = self.root / "runs" / "k222-reduced" / "II"
        (run_dir / "input.sha256").write_text(
            f"{sha256(input_path)}  {input_path}\n", encoding="utf-8"
        )
        result = self.run_tool(
            "dmc_mixed_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "2,1",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        meshes = self.mesh_vector(result.stdout)
        self.assertEqual(meshes["II"], 1)
        self.assertEqual(meshes["III"], 2)

        endpoints = self.root / "wrong-shift-endpoints.json"
        selection = self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            "--output-json",
            endpoints,
            expected_returncode=2,
        )
        self.assertIn("noncanonical Gamma-centred BvK shift", selection.stdout)

    def test_monitor_reports_ready_and_not_ready_states(self) -> None:
        paper = self.root / "paper.json"
        relative_values = {"Ih": 0.0}
        relative_values.update(
            {phase: (index + 1) / 10 for index, phase in enumerate(PHASES)}
        )
        paper.write_text(
            json.dumps(
                {
                    "results": {
                        f"k{mesh}{mesh}{mesh}": {
                            "GXTB": {"relative_kjmol": relative_values}
                        }
                        for mesh in (1, 2)
                    }
                }
            ),
            encoding="utf-8",
        )
        environment = os.environ.copy()
        environment.update({"ONCE": "1", "MESHES": "2,1"})
        status = self.root / "monitor-ready"
        result = subprocess.run(
            [
                str(TOOLS / "monitor_qualified_mixed_mae.sh"),
                str(self.root),
                str(self.reference),
                str(paper),
                CURRENT_DIGEST,
                str(status),
                "1",
            ],
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("status=READY", (status / "dmc_mixed_qualified.status").read_text())
        history = (status / "dmc_mixed_qualified_history.tsv").read_text().splitlines()
        self.assertEqual(len(history), 2)
        self.assertEqual(len(history[0].split("\t")), 7)
        self.assertEqual(len(history[1].split("\t")), 7)

        for mesh in (1, 2):
            (self.root / "runs" / f"k{mesh}{mesh}{mesh}-reduced" / "II" / "cp2k.out").unlink()
        not_ready = self.root / "monitor-not-ready"
        result = subprocess.run(
            [
                str(TOOLS / "monitor_qualified_mixed_mae.sh"),
                str(self.root),
                str(self.reference),
                str(paper),
                CURRENT_DIGEST,
                str(not_ready),
                "1",
            ],
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )
        self.assertEqual(result.returncode, 1, msg=result.stderr)
        readiness = (not_ready / "dmc_mixed_qualified.status").read_text()
        self.assertIn("status=NOT_READY", readiness)
        self.assertIn("no complete qualified same-mesh result for II", readiness)

    def test_fixed_mesh_mae_requires_complete_same_build_series(self) -> None:
        result = self.run_tool(
            "dmc_fixed_mesh_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
        )
        rows = result.stdout.splitlines()
        self.assertEqual(
            rows[0],
            "mesh\tMAE_kJ_mol_H2O\tRMSE_kJ_mol_H2O\tMaxAE_kJ_mol_H2O",
        )
        self.assertEqual(len(rows), 3)
        self.assertTrue(rows[1].startswith("Gamma\t"))
        self.assertTrue(rows[2].startswith("2x2x2\t"))

        (self.root / "runs/k222-reduced/II/binary.sha256").write_text(
            f"{OLDER_DIGEST}  /synthetic/cp2k.psmp\n", encoding="utf-8"
        )
        failed = self.run_tool(
            "dmc_fixed_mesh_mae.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            expected_returncode=1,
        )
        self.assertIn("binary mismatch: mesh=2 phase=II", failed.stderr)

    def test_verifier_rejects_nonfirst_passing_pair(self) -> None:
        endpoints = self.root / "endpoints.json"
        self.run_tool(
            "select_adaptive_endpoints.py",
            self.root,
            self.reference,
            "--meshes",
            "1,2",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            "--output-json",
            endpoints,
        )
        payload = json.loads(endpoints.read_text(encoding="utf-8"))
        for row in payload["rows"]:
            row["previous_mesh"] = 2
            row["endpoint_mesh"] = 3
        invalid = self.root / "not-first-pair.json"
        invalid.write_text(json.dumps(payload), encoding="utf-8")
        result = self.run_tool(
            "verify_adaptive_dmc13.py",
            self.root,
            invalid,
            self.reference,
            "--meshes",
            "1,2,3",
            "--require-binary-sha256",
            CURRENT_DIGEST,
            expected_returncode=2,
        )
        self.assertIn("does not record its first passing adjacent pair", result.stderr)


if __name__ == "__main__":
    unittest.main()
