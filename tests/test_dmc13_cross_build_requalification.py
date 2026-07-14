from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cross_build = load_script(
    "requalify_dmc13_cross_build_test",
    REPOSITORY
    / "DMC-ICE13"
    / "scripts"
    / "requalify_dmc13_cross_build.py",
)
runner = cross_build.load_runner(REPOSITORY)


class DMC13CrossBuildRequalificationTests(unittest.TestCase):
    def test_selection_requires_same_mesh_ih_and_rejects_duplicates(self) -> None:
        self.assertEqual(
            cross_build.parse_selection("k333:VII", runner), ("k333", "VII")
        )
        with self.assertRaises(argparse.ArgumentTypeError):
            cross_build.parse_selection("k333:unknown", runner)
        with self.assertRaisesRegex(ValueError, "same-mesh Ih"):
            cross_build.validate_matrix([("k333", "VII")])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            cross_build.validate_matrix([("k333", "Ih"), ("k333", "Ih")])
        cross_build.validate_matrix([("k333", "Ih"), ("k333", "VII")])

    def test_candidate_manifest_must_bind_concrete_schema_and_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k_source = root / "cp2k-source"
            tblite_source = root / "tblite-source"
            cp2k_source.mkdir()
            tblite_source.mkdir()

            def make_file(name: str, content: str) -> Path:
                path = root / name
                path.write_text(content)
                return path

            cp2k = make_file("cp2k.psmp", "cp2k")
            cp2k_library = make_file("libcp2k.so", "libcp2k")
            cp2k_cache = make_file("cp2k-cache.txt", "cache")
            tblite = make_file("tblite", "cli")
            tblite_library = make_file("libtblite.a", "libtblite")
            tblite_cache = make_file("tblite-cache.txt", "cache")
            oracle_input = make_file("oracle.inp", "input")
            oracle_output = make_file("oracle.out", "output")
            oracle_affinity = make_file("oracle-affinity.json", "affinity")
            oracle_launch = make_file("oracle-launch.json", "launch")
            fd_input = make_file("fd.inp", "input")
            fd_output = make_file("fd.out", "output")
            fd_affinity = make_file("fd-affinity.json", "affinity")

            identity = SimpleNamespace(
                cp2k=cp2k.resolve(),
                cp2k_library=cp2k_library.resolve(),
                tblite_static_library=tblite_library.resolve(),
                cp2k_sha256=runner.sha256(cp2k),
                cp2k_library_sha256=runner.sha256(cp2k_library),
                tblite_static_library_sha256=runner.sha256(tblite_library),
                cp2k_source_revision="4" * 40,
                tblite_source_revision="5" * 40,
            )
            args = SimpleNamespace(
                candidate_build_manifest=root / "build_manifest.json",
                cp2k_source=cp2k_source.resolve(),
                tblite_source=tblite_source.resolve(),
                tblite=tblite.resolve(),
                required_cp2k_ancestor=cross_build.REQUIRED_CP2K_ANCESTOR,
            )
            empty_diff_hash = hashlib.sha256(b"").hexdigest()

            def artifact_fields(path: Path, key: str) -> dict[str, str]:
                return {
                    key: str(path.resolve()),
                    f"{key}_sha256": runner.sha256(path),
                }

            qualification_identity = {
                "cp2k_binary_sha256": identity.cp2k_sha256,
                "cp2k_loaded_library_sha256": identity.cp2k_library_sha256,
                "save_tblite_library_sha256": (
                    identity.tblite_static_library_sha256
                ),
                "source_diff_sha256": empty_diff_hash,
            }
            oracle = {
                "status": "passed",
                **artifact_fields(oracle_input, "input"),
                **artifact_fields(oracle_output, "output"),
                **artifact_fields(oracle_affinity, "affinity_record"),
                **artifact_fields(oracle_launch, "launch_contract"),
                **qualification_identity,
                "exchange_duality_max_relative": 1.0e-10,
                "acp_duality_max_relative": 2.0e-10,
            }
            force_stress = {
                "status": "passed",
                **artifact_fields(fd_input, "input"),
                **artifact_fields(fd_output, "output"),
                **artifact_fields(fd_affinity, "affinity_record"),
                **qualification_identity,
                "force_max_abs_error_hartree_per_bohr": 1.0e-7,
                "force_rms_error_hartree_per_bohr": 1.0e-8,
                "force_max_relative_error": 1.0e-6,
                "virial_max_abs_error_hartree": 1.0e-7,
                "virial_rms_error_hartree": 1.0e-8,
                "virial_max_relative_error": 1.0e-6,
                "debug_displacement_bohr": 1.0e-4,
                "debug_cell_displacement_bohr": 1.0e-4,
            }
            payload = {
                "campaign_id": "post-5582-test",
                "campaign_state": "production_ready",
                "cp2k": {
                    "repository": "DCM-Uni-Paderborn/cp2k",
                    "branch": "g-xTB-pbc",
                    "revision": identity.cp2k_source_revision,
                    "reported_revision": identity.cp2k_source_revision[:12],
                    "required_upstream_ancestor": (
                        cross_build.REQUIRED_CP2K_ANCESTOR
                    ),
                    "source_path": str(cp2k_source.resolve()),
                    "source_diff_sha256": empty_diff_hash,
                    **artifact_fields(cp2k, "binary"),
                    **artifact_fields(cp2k_library, "loaded_library"),
                    **artifact_fields(cp2k_cache, "cmake_cache"),
                },
                "save_tblite": {
                    "repository": "DCM-Uni-Paderborn/save_tblite",
                    "branch": "pbc",
                    "revision": identity.tblite_source_revision,
                    "source_path": str(tblite_source.resolve()),
                    "source_diff_sha256": empty_diff_hash,
                    **artifact_fields(tblite, "cli"),
                    **artifact_fields(tblite_library, "static_library"),
                    **artifact_fields(tblite_cache, "cmake_cache"),
                },
                "qualification": {
                    "oracle": oracle,
                    "force_stress_fd": force_stress,
                    "oracle_binary_is_production_binary": True,
                },
            }
            completed = SimpleNamespace(returncode=0, stdout=b"")
            with mock.patch.object(cross_build.subprocess, "run", return_value=completed):
                cross_build.validate_candidate_manifest(
                    payload, identity, args, runner
                )
                payload["qualification"]["force_stress_fd"][
                    "output_sha256"
                ] = "0" * 64
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    cross_build.validate_candidate_manifest(
                        payload, identity, args, runner
                    )

    def test_full_matrix_is_hash_derived_and_subset_never_authorizes_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index_hash = "a" * 64
            source_identity = {
                "cp2k_source_revision": "1" * 40,
                "tblite_source_revision": "2" * 40,
            }
            summary_payload = {
                "benchmark": "DMC-ICE13",
                "status": "phasewise_kpoint_converged",
                "reference_phase": "Ih",
                "sources": {"validation_index": {"sha256": index_hash}},
                "fixed_k333_same_mesh_comparison": {
                    "mesh": "k333",
                    "not_a_phasewise_converged_result": True,
                },
                "methods": {
                    "GXTB": {
                        "status": "phasewise_kpoint_converged",
                        "provenance": {
                            "cp2k_source_revision": "1" * 40,
                            "provider_source_revision": "2" * 40,
                        },
                        "phases": {
                            phase: {
                                "selected_mesh": "k666",
                                "previous_mesh": "k555",
                            }
                            for phase in runner.PHASES
                            if phase != "Ih"
                        },
                    }
                },
            }
            summary_path = root / "summary.json"
            summary_path.write_text(json.dumps(summary_payload) + "\n")
            args = SimpleNamespace(
                scope=cross_build.FULL_PUBLICATION_SCOPE,
                selection=[],
                final_summary=summary_path,
                final_summary_sha256=runner.sha256(summary_path),
                reference_validation_index_sha256=index_hash,
            )
            selections, _ = cross_build.select_requalification_matrix(
                args, {"source_identity": source_identity}, runner
            )
            self.assertEqual(len(selections), 39)
            self.assertTrue(
                {("k333", phase) for phase in runner.PHASES}.issubset(selections)
            )
            self.assertIn(("k555", "Ih"), selections)
            self.assertIn(("k666", "VII"), selections)

            sentinel = cross_build.qualification_outcome(
                cross_build.SENTINEL_SCOPE
            )
            self.assertEqual(sentinel["status"], "sentinel_passed")
            self.assertFalse(sentinel["old_results_reusable"])
            args.selection = [("k333", "Ih"), ("k333", "VII")]
            with self.assertRaisesRegex(ValueError, "does not accept --selection"):
                cross_build.select_requalification_matrix(
                    args, {"source_identity": source_identity}, runner
                )

    def test_run_removes_stale_report_before_any_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.json"
            report.write_text('{"old_results_reusable": true}\n')
            args = SimpleNamespace(
                report=report,
                required_cp2k_ancestor="0" * 40,
            )
            with self.assertRaisesRegex(ValueError, "protocol-fixed"):
                cross_build.run(args, runner)
            self.assertFalse(report.exists())

    def test_pinned_json_rejects_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            path.write_text('{"status": "ready"}\n')
            digest = runner.sha256(path)
            self.assertEqual(
                cross_build.read_pinned_json(path, digest, runner)["status"],
                "ready",
            )
            path.write_text('{"status": "changed"}\n')
            with self.assertRaisesRegex(ValueError, "SHA256 pin mismatch"):
                cross_build.read_pinned_json(path, digest, runner)


if __name__ == "__main__":
    unittest.main()
