from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


benchmark = load_script(
    "dmc_ice13_kpoint_benchmark_summary_test",
    REPOSITORY / "DMC-ICE13" / "scripts" / "dmc_ice13_kpoint_benchmark.py",
)
summary = load_script(
    "finalize_dmc13_phasewise_summary",
    REPOSITORY
    / "DMC-ICE13"
    / "scripts"
    / "finalize_dmc13_phasewise_summary.py",
)


class DMC13PhasewiseSummaryTests(unittest.TestCase):
    def write_case(self, root: Path) -> None:
        data = root / "data"
        data.mkdir(parents=True)
        geometries = {
            phase: {"counts": {"O": 1}} for phase in benchmark.PHASES
        }
        (data / "geometries.json").write_text(
            json.dumps(geometries, indent=2, sort_keys=True) + "\n"
        )

        results: dict[str, object] = {"results": {}}
        for mesh, delta in (("k111", -0.04), ("k222", 0.0), ("k333", 0.0)):
            mesh_results: dict[str, object] = {}
            for method_index, method in enumerate(summary.METHODS, start=1):
                ih = -10.0 - method_index
                per_h2o = {"Ih": ih}
                relative = {"Ih": 0.0}
                energies = {"Ih": ih}
                for phase in summary.NONREFERENCE_PHASES:
                    dmc = (
                        benchmark.DMC_ABS_KJMOL[phase]
                        - benchmark.DMC_ABS_KJMOL["Ih"]
                    )
                    value = dmc + method_index + delta
                    phase_energy = ih + value / summary.HARTREE_TO_KJMOL
                    per_h2o[phase] = phase_energy
                    energies[phase] = phase_energy
                    relative[phase] = value
                mesh_results[method] = {
                    "complete": True,
                    "energies_hartree": energies,
                    "per_h2o_hartree": per_h2o,
                    "relative_kjmol": relative,
                }
            results["results"][mesh] = mesh_results

        report, rows = benchmark.build_phasewise_kpoint_convergence_report(
            results
        )
        (data / "dmc_ice13_gxtb_spglib_kpoint_results.json").write_text(
            json.dumps(results, indent=2, sort_keys=True) + "\n"
        )
        (data / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        (data / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.csv").write_text(
            "method,phase,status\n"
            + "".join(
                f"{row['method']},{row['phase']},{row['status']}\n"
                for row in rows
            )
        )
        (data / "dmc_ice13_gxtb_spglib_validation_index.json").write_text(
            json.dumps({"schema_version": 2, "records": []}) + "\n"
        )
        (data / "build_provenance.json").write_text(
            json.dumps(
                {
                    "cp2k": {
                        "source_revision": "1" * 40,
                        "executable_sha256": "2" * 64,
                        "library_sha256": "3" * 64,
                    },
                    "tblite": {
                        "main_revision": "4" * 40,
                        "local_merge_revision": "5" * 40,
                        "executable_sha256": "6" * 64,
                        "library_sha256": "7" * 64,
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        (data / "build_provenance_gxtb_spglib.json").write_text(
            json.dumps(
                {
                    "campaign": {
                        "id": "test-campaign",
                        "manifest_sha256": "8" * 64,
                    },
                    "cp2k": {
                        "execution_build_id": "sha256:" + "9" * 64,
                        "source_revision_validated": "a" * 40,
                        "sha256": "b" * 64,
                        "loaded_library_sha256": "c" * 64,
                    },
                    "save_tblite": {
                        "source_revision_validated": "d" * 40,
                        "sha256": "e" * 64,
                        "static_library_sha256": "f" * 64,
                    },
                    "protocol": {
                        "gxtb_protocol_id": "dmc13-gxtb-spglib-reduced-v1"
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    def test_finalizes_three_method_summary_with_raw_energy_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "DMC-ICE13"
            self.write_case(root)
            csv_path, json_path = summary.finalize(root)

            self.assertTrue(csv_path.is_file())
            payload = json.loads(json_path.read_text())
            self.assertEqual(payload["status"], "phasewise_kpoint_converged")
            self.assertEqual(
                payload["publication_qualification"]["status"],
                "gxtb_cross_build_requalification_pending",
            )
            self.assertFalse(
                payload["publication_qualification"]["paper_freeze_authorized"]
            )
            self.assertFalse(
                payload["publication_qualification"]["gxtb_old_results_reusable"]
            )
            self.assertEqual(list(payload["methods"]), list(summary.METHODS))
            self.assertAlmostEqual(
                payload["methods"]["GXTB"]["metrics_kjmol_per_h2o"]["MAE"],
                3.0,
            )
            self.assertEqual(
                payload["methods"]["GXTB"]["selected_mesh_distribution"],
                {"k222": 12},
            )
            fixed = payload["methods"]["GXTB"][
                "fixed_k333_same_mesh_comparison"
            ]
            self.assertEqual(
                fixed["status"],
                "numerically_unconverged_same_mesh_comparator",
            )
            self.assertEqual(
                payload["methods"]["GXTB"]["publication_qualification"][
                    "status"
                ],
                "diagnostic_pre_post_5582_requalification",
            )
            self.assertFalse(
                payload["methods"]["GXTB"]["publication_qualification"][
                    "old_results_reusable"
                ]
            )
            self.assertAlmostEqual(fixed["metrics_kjmol_per_h2o"]["MAE"], 3.0)
            self.assertTrue(
                payload["fixed_k333_same_mesh_comparison"][
                    "not_a_phasewise_converged_result"
                ]
            )
            phase = payload["methods"]["GXTB"]["phases"]["VII"]
            self.assertEqual(phase["selected_mesh"], "k222")
            self.assertAlmostEqual(
                phase["relative_energy_kjmol_per_h2o"]
                - phase["dmc_relative_energy_kjmol_per_h2o"],
                3.0,
            )
            self.assertEqual(
                payload["sources"]["phasewise_json"]["sha256"],
                summary.sha256(
                    root
                    / "data"
                    / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json"
                ),
            )
            self.assertEqual(len(csv_path.read_text().splitlines()), 4)

    def test_incomplete_method_removes_stale_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "DMC-ICE13"
            self.write_case(root)
            data = root / "data"
            report_path = (
                data
                / "dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json"
            )
            report = json.loads(report_path.read_text())
            report["methods"]["GXTB"]["status"] = "unresolved_phases"
            report["methods"]["GXTB"]["phasewise_kpoint_converged"] = False
            report_path.write_text(json.dumps(report, indent=2) + "\n")
            csv_path = data / f"{summary.SUMMARY_STEM}.csv"
            json_path = data / f"{summary.SUMMARY_STEM}.json"
            csv_path.write_text("stale\n")
            json_path.write_text("{}\n")

            with self.assertRaisesRegex(
                ValueError, "GXTB is not phase-wise k-point converged"
            ):
                summary.finalize(root)
            self.assertFalse(csv_path.exists())
            self.assertFalse(json_path.exists())


if __name__ == "__main__":
    unittest.main()
