from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "X23b" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import x23b_common as common  # noqa: E402
import x23b_k222_fd_gate as fd_gate  # noqa: E402


def fake_executable(directory: Path) -> Path:
    path = directory / "cp2k"
    path.write_text("#!/bin/sh\nexit 99\n")
    path.chmod(path.stat().st_mode | 0o111)
    return path


def fake_campaign(cp2k: Path) -> dict[str, object]:
    return common.make_campaign_identity(
        campaign_id="fd-test-campaign",
        cp2k_executable_sha256=common.sha256_file(cp2k),
        cp2k_loaded_library_sha256="2" * 64,
        cp2k_cmake_cache_sha256="3" * 64,
        cp2k_embedded_source_revision="0123456789",
        cp2k_source_revision="0123456789" + "0" * 30,
        save_tblite_executable_sha256="4" * 64,
        save_tblite_source_revision="5" * 40,
        save_tblite_library_sha256="6" * 64,
        save_tblite_cmake_cache_sha256="7" * 64,
        dependency_lock_sha256="8" * 64,
    )


def kpoint_report() -> list[str]:
    lines = [
        f"       Number of Special K-points: {8:45d}",
        "       K-point Mesh:                                               2     2     2",
        "                   Wavevector Basis                   Special Points    Rotation",
    ]
    for index in range(8):
        x, y, z = ((index >> shift) & 1 for shift in (2, 1, 0))
        lines.append(
            f" {index + 1:10d} {0.5 * x:9.5f} {0.5 * y:9.5f} {0.5 * z:9.5f}"
            f" {index + 1:14d} {index + 1:11d} {1:11d}"
        )
    return lines


def energy_output(energy: float) -> str:
    return "\n".join(
        [
            *kpoint_report(),
            f" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] {energy:.16f}",
            " PROGRAM ENDED AT 2026-07-14 12:00:00",
        ]
    ) + "\n"


def baseline_output(
    energy: float,
    elements: list[str],
    forces: list[list[float]],
    stress: list[list[float]],
) -> str:
    lines = [
        *kpoint_report(),
        f" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] {energy:.16f}",
        " FORCES| Atomic forces [hartree/bohr]",
        " FORCES|   Atom     x                   y                   z                   |f|",
    ]
    for atom, vector in enumerate(forces, start=1):
        norm = math.sqrt(sum(value * value for value in vector))
        lines.append(
            f" FORCES| {atom:6d} {vector[0]:19.12E} {vector[1]:19.12E} "
            f"{vector[2]:19.12E} {norm:19.12E}"
        )
    lines += [
        " FORCES| Sum     0.0 0.0 0.0",
        " STRESS| Analytical stress tensor [GPa]",
        " STRESS|                         x                   y                   z",
    ]
    for axis, row in zip("xyz", stress):
        lines.append(
            f" STRESS|      {axis} {row[0]:19.11E} {row[1]:19.11E} {row[2]:19.11E}"
        )
    lines.append(" PROGRAM ENDED AT 2026-07-14 12:00:00")
    return "\n".join(lines) + "\n"


class X23bK222FiniteDifferenceGateTests(unittest.TestCase):
    def _prepared(self, base: Path) -> tuple[Path, Path, dict[str, object]]:
        cp2k = fake_executable(base)
        campaign = fake_campaign(cp2k)
        output_root = base / "fd"
        fd_gate.prepare(output_root, campaign)
        return output_root, cp2k, campaign

    def _write_exact_measurements(
        self,
        output_root: Path,
        cp2k: Path,
        campaign: dict[str, object],
    ) -> None:
        cases = fd_gate.load_manifest(output_root, campaign)
        stress = [[1.0, 0.2, 0.3], [0.2, 2.0, 0.4], [0.3, 0.4, 3.0]]
        reference_energy = -10.0
        for case in cases.values():
            forces = [
                [atom * 1.0e-4, -atom * 0.7e-4, atom * 0.3e-4]
                for atom in range(1, int(case["atom_count"]) + 1)
            ]
            directions = {
                str(row["id"]): row["vector_cartesian"]
                for row in case["coordinate_directions"]
            }
            for job in case["jobs"]:
                output = Path(str(job["output"]))
                if job["job_type"] == "baseline_energy_force":
                    text = baseline_output(reference_energy, case["elements"], forces, stress)
                elif job["job_type"] == "coordinate_energy":
                    vector = directions[str(job["direction_id"])]
                    derivative = -sum(
                        force_component * direction_component
                        for force, direction_row in zip(forces, vector)
                        for force_component, direction_component in zip(force, direction_row)
                    )
                    energy = (
                        reference_energy
                        + int(job["sign"]) * float(job["step"]) * derivative
                    )
                    text = energy_output(energy)
                else:
                    generator = job["generator"]
                    conjugation = sum(
                        stress[i][j] * generator[i][j]
                        for i in range(3)
                        for j in range(3)
                    )
                    derivative = (
                        -float(case["reference_volume_A3"])
                        * conjugation
                        * fd_gate.GPA_ANGSTROM3_TO_HARTREE
                    )
                    energy = (
                        reference_energy
                        + int(job["sign"]) * float(job["step"]) * derivative
                    )
                    text = energy_output(energy)
                output.write_text(text)
                fd_gate.parse_job_output(case, job)
                common.write_job_stamp(
                    Path(str(job["run_dir"])),
                    Path(str(job["input"])),
                    cp2k,
                    "GXTB",
                    fd_gate.PHASE,
                    "converged_measured_not_approved",
                    details={
                        "returncode": 0,
                        "output": str(output),
                        "output_sha256": common.sha256_file(output),
                    },
                    campaign_identity=campaign,
                    protocol_identity=fd_gate.protocol_identity(output_root, job),
                    source_artifacts=fd_gate.source_artifacts(output_root, case),
                )

    def test_four_system_pilot_is_small_exact_and_includes_triclinic_case(self) -> None:
        self.assertEqual(
            fd_gate.DEFAULT_SYSTEMS,
            ("ammonia", "14-cyclohexanedione", "acetic_acid", "ethylcarbamate"),
        )
        source, _ = fd_gate.cellopt.experimental_reference_paths("ethylcarbamate")
        cell = fd_gate._cell(source.read_text())
        # The expanded P1 reference has a fully oblique (triclinic) cell.
        self.assertNotAlmostEqual(cell[0][1], cell[1][0])
        self.assertNotEqual(cell[1][0], 0.0)
        self.assertNotEqual(cell[2][0], 0.0)
        self.assertNotEqual(cell[2][1], 0.0)

    def test_prepare_freezes_36_jobs_and_normalized_translation_free_directions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, _, campaign = self._prepared(Path(tmp))
            cases = fd_gate.load_manifest(output_root, campaign)
            self.assertEqual(sum(len(case["jobs"]) for case in cases.values()), 36)
            for case in cases.values():
                for direction in case["coordinate_directions"]:
                    vector = direction["vector_cartesian"]
                    self.assertAlmostEqual(
                        sum(value * value for row in vector for value in row), 1.0, places=13
                    )
                    for axis in range(3):
                        self.assertAlmostEqual(sum(row[axis] for row in vector), 0.0, places=13)
            # Exact prepare is resumable; a protocol change cannot reuse the root.
            fd_gate.prepare(output_root, campaign)
            with self.assertRaisesRegex(ValueError, "already freezes"):
                fd_gate.prepare(output_root, campaign, coordinate_step_bohr=2.0e-3)

    def test_coordinate_pair_is_a_cartesian_central_displacement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, _, campaign = self._prepared(Path(tmp))
            case = fd_gate.load_manifest(output_root, campaign)["ethylcarbamate"]
            direction = case["coordinate_directions"][0]["vector_cartesian"]
            jobs = fd_gate._jobs_by_id(case)
            plus = Path(str(jobs["coord_collective_cartesian_1_plus"]["input"])).read_text()
            minus = Path(str(jobs["coord_collective_cartesian_1_minus"]["input"])).read_text()
            _, plus_scaled = fd_gate._scaled_coordinates(plus)
            _, minus_scaled = fd_gate._scaled_coordinates(minus)
            cell = fd_gate._cell(plus)
            step = fd_gate.DEFAULT_COORDINATE_STEP_BOHR
            for atom in range(len(direction)):
                delta_fractional = [
                    (plus_scaled[atom][axis] - minus_scaled[atom][axis]) / 2.0
                    for axis in range(3)
                ]
                delta_cart = fd_gate._row_times_matrix(delta_fractional, cell)
                for observed, expected in zip(delta_cart, direction[atom]):
                    self.assertAlmostEqual(
                        observed / (step * fd_gate.BOHR_TO_ANGSTROM), expected, places=10
                    )

    def test_manifest_rejects_changed_generated_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, _, campaign = self._prepared(Path(tmp))
            case = fd_gate.load_manifest(output_root, campaign)["ammonia"]
            input_path = Path(str(case["jobs"][0]["input"]))
            input_path.write_text(input_path.read_text() + "# changed\n")
            with self.assertRaisesRegex(ValueError, "FD input changed"):
                fd_gate.load_manifest(output_root, campaign)

    def test_manifest_rejects_self_consistent_protocol_rewrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, _, campaign = self._prepared(Path(tmp))
            path = fd_gate.manifest_path(output_root)
            payload = json.loads(path.read_text())
            payload["protocol"]["stress_derivative_identity"] = "wrong sign"
            payload_without_digest = dict(payload)
            payload_without_digest.pop("payload_sha256")
            payload["payload_sha256"] = fd_gate._fingerprint(payload_without_digest)
            path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "not the deterministic rendering"):
                fd_gate.load_manifest(output_root, campaign)

    def test_existing_unstamped_output_is_stale_and_never_executed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, cp2k, campaign = self._prepared(Path(tmp))
            case = fd_gate.load_manifest(output_root, campaign)["ammonia"]
            job = case["jobs"][0]
            Path(str(job["output"])).write_text("partial foreign output\n")
            system, job_id, code, action = fd_gate.run_one(
                output_root, case, job, cp2k, 1, campaign
            )
            self.assertEqual((system, job_id), ("ammonia", "baseline"))
            self.assertNotEqual(code, 0)
            self.assertEqual(action, "STALE_OUTPUT")

    def test_missing_output_from_completed_stamp_is_stale_not_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root, cp2k, campaign = self._prepared(Path(tmp))
            case = fd_gate.load_manifest(output_root, campaign)["ammonia"]
            job = case["jobs"][0]
            output = Path(str(job["output"]))
            output.write_text(
                baseline_output(
                    -10.0,
                    case["elements"],
                    [[0.0, 0.0, 0.0] for _ in range(int(case["atom_count"]))],
                    [[0.0, 0.0, 0.0] for _ in range(3)],
                )
            )
            common.write_job_stamp(
                Path(str(job["run_dir"])),
                Path(str(job["input"])),
                cp2k,
                "GXTB",
                fd_gate.PHASE,
                "converged_measured_not_approved",
                details={
                    "returncode": 0,
                    "output": str(output),
                    "output_sha256": common.sha256_file(output),
                },
                campaign_identity=campaign,
                protocol_identity=fd_gate.protocol_identity(output_root, job),
                source_artifacts=fd_gate.source_artifacts(output_root, case),
            )
            output.unlink()
            _, _, code, action = fd_gate.run_one(output_root, case, job, cp2k, 1, campaign)
            self.assertNotEqual(code, 0)
            self.assertEqual(action, "STALE_STAMP")

    def test_collect_reports_force_and_stress_conjugations_without_approval(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root, cp2k, campaign = self._prepared(base)
            self._write_exact_measurements(output_root, cp2k, campaign)
            csv_path = base / "measured.csv"
            json_path = base / "measured.json"
            report = fd_gate.collect(output_root, csv_path, json_path, campaign)
            self.assertEqual(report["row_count"], 16)
            self.assertFalse(report["approved"])
            self.assertEqual(report["scientific_status"], "measured_not_approved")
            coordinate = next(
                row
                for row in report["rows"]
                if row["measurement_type"] == "coordinate_directional_derivative"
            )
            self.assertAlmostEqual(
                coordinate["finite_difference_energy_derivative_hartree_per_parameter"],
                coordinate["force_projection_minus_F_dot_d_hartree_per_bohr"],
                places=10,
            )
            isotropic = next(
                row
                for row in report["rows"]
                if row["system"] == "ammonia"
                and row["direction_id"] == "isotropic_linear"
            )
            self.assertAlmostEqual(isotropic["stress_conjugation_sigma_colon_G_GPa"], 6.0)
            self.assertAlmostEqual(isotropic["finite_difference_stress_conjugation_GPa"], 6.0, places=8)
            self.assertIn("-V*(sigma:G)", report["formulae"]["strain"])

            report_sha = common.sha256_file(json_path)
            approval_path = base / "approval.json"
            approval = fd_gate.approve_report(
                json_path,
                approval_path,
                reviewer="unit-test reviewer",
                coordinate_abs_tolerance_hartree_per_bohr=1.0e-8,
                stress_abs_tolerance_gpa=1.0e-6,
            )
            self.assertEqual(approval["decision"], "approved")
            self.assertEqual(common.sha256_file(json_path), report_sha)
            self.assertFalse(json.loads(json_path.read_text())["approved"])
            self.assertTrue(approval_path.is_file())


if __name__ == "__main__":
    unittest.main()
