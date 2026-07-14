from __future__ import annotations

import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "X23b" / "scripts"
sys.path.insert(0, str(SCRIPTS))


def load_script():
    path = SCRIPTS / "x23b_adaptive_kpoint_convergence.py"
    spec = importlib.util.spec_from_file_location(
        "x23b_adaptive_kpoint_convergence_test", path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


kconv = load_script()


def fake_observation(
    energy: float,
    *,
    lengths: tuple[float, float, float] = (10.0, 11.0, 12.0),
    angles: tuple[float, float, float] = (80.0, 90.0, 100.0),
    volume: float = 1000.0,
    digest: str = "a",
) -> dict[str, object]:
    return {
        "energy_hartree": energy,
        "lengths_A": list(lengths),
        "angles_deg": list(angles),
        "volume_A3": volume,
        "output_sha256": digest * 64,
    }


class X23bAdaptiveProtocolTests(unittest.TestCase):
    def test_exact_scope_thresholds_and_native_input_contract(self) -> None:
        self.assertEqual(len(kconv.systems()), 23)
        self.assertEqual(kconv.METHODS, ("GFN1", "GFN2", "GXTB"))
        self.assertEqual(kconv.MESH_IDS[0], "k111")
        self.assertEqual(kconv.MESH_IDS[-1], "k888")
        self.assertEqual(kconv.ENERGY_TOLERANCE_KJMOL, 0.05)
        self.assertEqual(kconv.LENGTH_TOLERANCE_PERCENT, 0.05)
        self.assertEqual(kconv.VOLUME_TOLERANCE_PERCENT, 0.10)
        self.assertEqual(kconv.ANGLE_TOLERANCE_DEG, 0.05)
        self.assertEqual(kconv.gamma_centered_shift(3), 0.0)
        self.assertEqual(kconv.gamma_centered_shift(4), 0.375)

        for series in kconv.SERIES:
            for method in kconv.METHODS:
                text = kconv.render_input(series, method, "ammonia", 4)
                self.assertIn("SCHEME MACDONALD 4 4 4 0.375 0.375 0.375", text)
                self.assertIn("SYMMETRY_BACKEND SPGLIB", text)
                self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", text)
                self.assertIn("FULL_GRID F", text)
                self.assertNotIn("MULTIPLE_UNIT_CELL", text)
                self.assertNotIn("SUPERCELL", text)
                if series == "fixed_experimental_sp":
                    self.assertIn("RUN_TYPE ENERGY", text)
                    self.assertNotIn("&MOTION", text)
                else:
                    self.assertIn("RUN_TYPE CELL_OPT", text)
                    self.assertIn("KEEP_ANGLES F", text)
                    self.assertIn("&KPOINTS", text)

    def test_single_step_energy_gate_selects_the_denser_mesh(self) -> None:
        molecules = int(kconv.system_metadata("ammonia")["molecules_per_cell"])
        exact_cell_delta = (
            kconv.ENERGY_TOLERANCE_KJMOL
            * molecules
            / kconv.HARTREE_TO_KJMOL
        )
        observations = {
            1: fake_observation(-100.0, digest="a"),
            2: fake_observation(-100.0 + exact_cell_delta, digest="b"),
        }
        result = kconv.assess_series(
            "fixed_experimental_sp", "GFN2", "ammonia", observations
        )
        self.assertEqual(result["status"], "converged")
        self.assertEqual(result["selected_mesh"], "k222")
        self.assertIs(result["selected_observation"], observations[2])
        row = result["rows"][0]
        self.assertTrue(row["energy_pass"])
        self.assertTrue(row["selected_step"])
        self.assertAlmostEqual(
            row["energy_abs_delta_kJmol"],
            kconv.ENERGY_TOLERANCE_KJMOL,
            places=9,
        )

    def test_cell_gate_starts_at_k333_and_records_each_criterion(self) -> None:
        observations = {
            1: fake_observation(-10.0, digest="a"),
            2: fake_observation(-10.1, digest="b"),
            3: fake_observation(
                -10.2,
                lengths=(10.004, 11.004, 12.004),
                angles=(80.04, 90.04, 100.04),
                volume=1000.9,
                digest="c",
            ),
        }
        result = kconv.assess_series(
            "independent_cellopt", "GXTB", "ammonia", observations
        )
        self.assertEqual(result["selected_mesh"], "k333")
        first, selected = result["rows"]
        self.assertFalse(first["eligible_for_stopping"])
        self.assertTrue(selected["length_pass"])
        self.assertTrue(selected["volume_pass"])
        self.assertTrue(selected["angle_pass"])
        self.assertTrue(selected["all_required_criteria_pass"])
        self.assertTrue(selected["selected_step"])

        failed = copy.deepcopy(observations)
        failed[3] = fake_observation(
            -10.2,
            lengths=(10.006, 11.0, 12.0),
            angles=(80.04, 90.04, 100.04),
            volume=1000.9,
            digest="d",
        )
        result = kconv.assess_series(
            "independent_cellopt", "GXTB", "ammonia", failed
        )
        row = result["rows"][-1]
        self.assertFalse(row["length_pass"])
        self.assertTrue(row["volume_pass"])
        self.assertTrue(row["angle_pass"])
        self.assertFalse(row["all_required_criteria_pass"])
        self.assertEqual(result["required_meshes"], [4])

    def test_maximum_mesh_and_post_convergence_extras_are_fail_closed(self) -> None:
        observations = {
            mesh: fake_observation(-float(mesh), digest=chr(96 + mesh))
            for mesh in kconv.MESH_NUMBERS
        }
        result = kconv.assess_series(
            "fixed_experimental_sp", "GFN1", "ammonia", observations
        )
        self.assertEqual(result["status"], "maximum_mesh_unconverged")
        passing = {
            1: fake_observation(-10.0, digest="a"),
            2: fake_observation(-10.0, digest="b"),
            3: fake_observation(-10.0, digest="c"),
        }
        with self.assertRaisesRegex(ValueError, "after first converged step"):
            kconv.assess_series(
                "fixed_experimental_sp", "GFN1", "ammonia", passing
            )

    def test_cell_parser_returns_lengths_angles_and_volume(self) -> None:
        text = """
 CELL| Vector a [angstrom]:  4.0 0.0 0.0 |a| = 4.0
 CELL| Vector b [angstrom]:  1.0 5.0 0.0 |b| = 5.1
 CELL| Vector c [angstrom]:  0.5 0.7 6.0 |c| = 6.1
"""
        cell = kconv.parse_cell(text)
        self.assertAlmostEqual(cell["lengths_A"][0], 4.0)
        self.assertAlmostEqual(cell["volume_A3"], 120.0)
        self.assertEqual(len(cell["angles_deg"]), 3)


class X23bAdaptiveManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name) / "workflow"
        self.root.mkdir()
        build_manifest = Path(self.temporary.name) / "build_manifest.json"
        build_manifest.write_text('{"campaign_state":"production_ready"}\n')
        identity = kconv.common.make_campaign_identity(
            campaign_id="test-post-5582",
            cp2k_executable_sha256="1" * 64,
            cp2k_loaded_library_sha256="2" * 64,
            cp2k_cmake_cache_sha256="3" * 64,
            cp2k_embedded_source_revision="4" * 12,
            cp2k_source_revision="4" * 40,
            save_tblite_executable_sha256="5" * 64,
            save_tblite_source_revision="6" * 40,
            save_tblite_library_sha256="7" * 64,
            save_tblite_cmake_cache_sha256="8" * 64,
            dependency_lock_sha256="9" * 64,
        )
        self.binding = {
            "campaign_identity": identity,
            "campaign_manifest": {
                "path": str(build_manifest),
                "file_sha256": kconv.sha256(build_manifest),
                "campaign_id": "test-post-5582",
                "campaign_state": "production_ready",
                "authority": "test",
            },
            "cp2k": {
                "path": "/test/cp2k.psmp",
                "sha256": "1" * 64,
                "source_revision": "4" * 40,
                "loaded_library_sha256": "2" * 64,
            },
            "save_tblite": {
                "path": "/test/tblite",
                "sha256": "5" * 64,
                "source_revision": "6" * 40,
                "static_library_sha256": "7" * 64,
            },
        }
        for series in kconv.SERIES:
            for method in kconv.METHODS:
                for system in (str(row["id"]) for row in kconv.systems()):
                    for mesh in kconv.MESH_NUMBERS:
                        mesh_id = f"k{mesh}{mesh}{mesh}"
                        run_dir = self.root / kconv.job_relative_dir(
                            series, method, system, mesh_id
                        )
                        run_dir.mkdir(parents=True, exist_ok=True)
                        (run_dir / "input.inp").write_text(
                            kconv.render_input(series, method, system, mesh)
                        )
        self.payload = kconv.workflow_payload(self.root, self.binding)
        self.manifest = self.root / "workflow_manifest.json"
        self.write_payload(self.payload)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_payload(self, payload: dict[str, object]) -> str:
        self.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return kconv.sha256(self.manifest)

    def test_exact_1104_job_matrix_loads_and_reordering_is_rejected(self) -> None:
        digest = self.write_payload(self.payload)
        loaded = kconv.load_workflow(self.manifest, digest)
        self.assertEqual(len(loaded["jobs"]), 2 * 3 * 23 * 8)

        reordered = copy.deepcopy(self.payload)
        reordered["jobs"][0], reordered["jobs"][1] = (
            reordered["jobs"][1],
            reordered["jobs"][0],
        )
        digest = self.write_payload(reordered)
        with self.assertRaisesRegex(ValueError, "reordered, duplicated, or incomplete"):
            kconv.load_workflow(self.manifest, digest)

    def test_tampered_input_and_build_manifest_are_rejected(self) -> None:
        first_input = self.root / self.payload["jobs"][0]["input"]
        first_input.write_text(first_input.read_text() + "# tampered\n")
        digest = self.write_payload(self.payload)
        with self.assertRaisesRegex(ValueError, "input fingerprint"):
            kconv.load_workflow(self.manifest, digest)

        first_input.write_text(
            kconv.render_input(
                self.payload["jobs"][0]["series"],
                self.payload["jobs"][0]["method"],
                self.payload["jobs"][0]["system"],
                self.payload["jobs"][0]["mesh_number"],
            )
        )
        build_manifest = Path(self.binding["campaign_manifest"]["path"])
        build_manifest.write_text('{"campaign_state":"changed"}\n')
        digest = self.write_payload(self.payload)
        with self.assertRaisesRegex(ValueError, "campaign manifest has changed"):
            kconv.load_workflow(self.manifest, digest)

    def test_pre_5582_campaign_and_tampered_command_are_rejected(self) -> None:
        legacy = copy.deepcopy(self.payload)
        legacy["build"]["campaign_manifest"]["campaign_id"] = "legacy-production"
        digest = self.write_payload(legacy)
        with self.assertRaisesRegex(ValueError, "held.*post-5582"):
            kconv.load_workflow(self.manifest, digest)

        record = self.payload["jobs"][0]
        output = self.root / record["output"]
        stamp = self.root / record["stamp"]
        output.write_text("completed output\n")
        stamp.write_text(
            json.dumps(
                {
                    "schema": "expected",
                    "status": "converged",
                    "details": {
                        "returncode": 0,
                        "output": str(output.resolve()),
                        "output_sha256": kconv.sha256(output),
                        "command": ["wrong-cp2k"],
                        "threads": 1,
                    },
                }
            )
            + "\n"
        )
        fake_cp2k = Path(self.temporary.name) / "cp2k.psmp"
        fake_cp2k.write_text("fake executable\n")
        with mock.patch.object(
            kconv, "expected_stamp", return_value={"schema": "expected"}
        ):
            with self.assertRaisesRegex(ValueError, "command/threading"):
                kconv.validate_completed_stamp(
                    self.manifest,
                    "0" * 64,
                    self.payload,
                    record,
                    fake_cp2k,
                )

    def test_finalizer_removes_all_stale_outputs_on_failure(self) -> None:
        output = Path(self.temporary.name) / "publication"
        output.mkdir()
        for suffix in ("csv", "json", "tex"):
            (output / f"{kconv.OUTPUT_STEM}.{suffix}").write_text("stale\n")
        args = SimpleNamespace(
            workflow_manifest=self.manifest,
            workflow_manifest_sha256="0" * 64,
            output_dir=output,
        )
        with mock.patch.object(
            kconv, "load_workflow", side_effect=ValueError("incomplete")
        ):
            with self.assertRaisesRegex(ValueError, "incomplete"):
                kconv.finalize(args)
        for suffix in ("csv", "json", "tex"):
            self.assertFalse((output / f"{kconv.OUTPUT_STEM}.{suffix}").exists())


if __name__ == "__main__":
    unittest.main()
