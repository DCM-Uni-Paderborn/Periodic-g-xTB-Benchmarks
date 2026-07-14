from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "Goldzak12" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import run_goldzak12_benchmark as base  # noqa: E402
import run_goldzak12_eos_benchmark as eos  # noqa: E402
import run_goldzak12_k_convergence as kconv  # noqa: E402


def complete_values() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for method_index, method in enumerate(kconv.METHODS):
        for solid_index, solid in enumerate(kconv.PAPER_SYSTEMS):
            origin_a = 4.0 + 0.01 * solid_index + 0.001 * method_index
            origin_e = 5.0 + 0.02 * solid_index + 0.01 * method_index
            # 3->4 fails both; 4->5 passes both.  k555 is selected.
            for number, da, de in (
                (3, 0.0020, 0.0010),
                (4, 0.0002, 0.0002),
                (5, 0.0001, 0.0001),
            ):
                rows.append(
                    {
                        "solid": solid,
                        "method": method,
                        "mesh": kconv.mesh_name(number),
                        "mesh_n": number,
                        "a0_A": f"{origin_a + da:.10f}",
                        "ecoh_eV_per_atom": f"{origin_e + de:.12f}",
                    }
                )
    return rows


class LC10AdaptiveKConvergenceTests(unittest.TestCase):
    def test_exact_thresholds_and_single_step_take_denser_value(self) -> None:
        rows = complete_values()
        # C/GFN1 passes already at 3->4, exactly on both inclusive limits.
        keyed = {
            (row["method"], row["solid"], row["mesh_n"]): row for row in rows
        }
        coarse = keyed[("GFN1", "C", 3)]
        dense = keyed[("GFN1", "C", 4)]
        coarse["a0_A"] = f"{float(dense['a0_A']) + kconv.LATTICE_THRESHOLD_A:.10f}"
        coarse["ecoh_eV_per_atom"] = (
            f"{float(dense['ecoh_eV_per_atom']) + kconv.ECOH_THRESHOLD_EV_PER_ATOM:.12f}"
        )
        steps, selections, pending = kconv.assess_convergence(rows)
        self.assertEqual(pending, [])
        chosen = next(
            row
            for row in selections
            if row["method"] == "GFN1" and row["solid"] == "C"
        )
        self.assertEqual(chosen["converged_from_mesh"], "k333")
        self.assertEqual(chosen["selected_mesh"], "k444")
        self.assertEqual(chosen["a0_A"], dense["a0_A"])
        self.assertEqual(chosen["ecoh_eV_per_atom"], dense["ecoh_eV_per_atom"])
        selected_step = next(
            row
            for row in steps
            if row["method"] == "GFN1"
            and row["solid"] == "C"
            and row["coarse_mesh"] == "k333"
        )
        self.assertTrue(selected_step["lattice_passed"])
        self.assertTrue(selected_step["ecoh_passed"])
        self.assertTrue(selected_step["both_passed"])
        self.assertEqual(selected_step["decision_rule"], "one_consecutive_step_AND")

    def test_one_failed_criterion_requests_exactly_the_next_mesh(self) -> None:
        rows = complete_values()
        # Make all adjacent C/GXTB a0 steps fail while Ecoh passes.
        keyed = {
            (row["method"], row["solid"], row["mesh_n"]): row for row in rows
        }
        keyed[("GXTB", "C", 3)]["a0_A"] = "4.0060000000"
        keyed[("GXTB", "C", 4)]["a0_A"] = "4.0030000000"
        keyed[("GXTB", "C", 5)]["a0_A"] = "4.0000000000"
        _steps, selections, pending = kconv.assess_convergence(rows)
        self.assertEqual(pending, [("GXTB", "C", 6)])
        self.assertFalse(
            any(
                row.get("method") == "GXTB" and row.get("solid") == "C"
                for row in selections
            )
        )

    def test_k888_without_a_passing_step_is_fail_closed(self) -> None:
        rows = complete_values()
        rows = [
            row
            for row in rows
            if not (row["method"] == "GXTB" and row["solid"] == "C")
        ]
        for number in range(3, 9):
            rows.append(
                {
                    "solid": "C",
                    "method": "GXTB",
                    "mesh": kconv.mesh_name(number),
                    "mesh_n": number,
                    "a0_A": f"{4.0 + 0.003 * (8 - number):.10f}",
                    "ecoh_eV_per_atom": f"{5.0 + 0.001 * (8 - number):.12f}",
                }
            )
        _steps, selections, pending = kconv.assess_convergence(rows)
        self.assertEqual(pending, [])
        chosen = next(
            row
            for row in selections
            if row["method"] == "GXTB" and row["solid"] == "C"
        )
        self.assertEqual(chosen["selection_status"], "unconverged_at_k888")
        self.assertEqual(chosen["selected_mesh"], "")

    def test_manifest_explicitly_forbids_rms_and_two_step_gates(self) -> None:
        rows = complete_values()
        steps, selections, pending = kconv.assess_convergence(rows)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data").mkdir()
            payload = kconv.write_convergence_artifacts(
                root,
                rows,
                steps,
                selections,
                pending,
                campaign={"campaign_id": "test"},
                fits_sha256="f" * 64,
            )
            stored = json.loads((root / "data" / kconv.CONVERGENCE_NAME).read_text())
        self.assertEqual(payload, stored)
        self.assertEqual(stored["algorithm"]["required_consecutive_passing_steps"], 1)
        self.assertFalse(stored["algorithm"]["aggregate_rms_gate"])
        self.assertEqual(stored["algorithm"]["criteria_combination"], "AND")
        self.assertAlmostEqual(
            stored["algorithm"]["cohesive_abs_delta_threshold_eV_per_atom"],
            0.000518213,
            places=9,
        )

    def test_mesh_parser_rejects_non_cubic_meshes(self) -> None:
        self.assertEqual(kconv.mesh_number("k101010"), 10)
        with self.assertRaisesRegex(ValueError, "non-cubic"):
            kconv.mesh_number("k345")

    def test_keyed_merge_preserves_independent_eos_meshes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fits.csv"
            base.write_csv(
                path,
                [
                    {
                        "solid": "C",
                        "method": "GFN1",
                        "eos_mesh": "k333",
                        "a_eos_A": "3.50",
                    }
                ],
            )
            eos.merge_keyed_rows(
                path,
                [
                    {
                        "solid": "C",
                        "method": "GFN1",
                        "eos_mesh": "k444",
                        "a_eos_A": "3.51",
                    }
                ],
                ("solid", "method", "eos_mesh"),
            )
            rows = base.read_csv(path)
        self.assertEqual({row["eos_mesh"] for row in rows}, {"k333", "k444"})

    def test_same_campaign_stamp_is_written_for_gfn_and_gxtb(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k"
            cp2k.write_bytes(b"same-final-binary")
            campaign = base.make_campaign_identity(
                campaign_id="same-binary-test",
                cp2k_executable_sha256=base.sha256(cp2k),
                cp2k_loaded_library_sha256="library",
                cp2k_cmake_cache_sha256="cache",
                cp2k_embedded_source_revision="a" * 10,
                cp2k_source_revision="a" * 40,
                save_tblite_executable_sha256="save-cli",
                save_tblite_source_revision="b" * 40,
                save_tblite_library_sha256="save-library",
                save_tblite_cmake_cache_sha256="save-cache",
                dependency_lock_sha256="dependencies",
            )
            specs: list[tuple[str, Path, Path, bool]] = []
            for method in ("GFN1", "GXTB"):
                input_path = root / f"{method}.inp"
                output_path = root / f"{method}.out"
                input_path.write_text(
                    base.solid_input(
                        base.LC10_PAPER_REFERENCES[0],
                        method,
                        "ENERGY",
                        "k333",
                        base.LC10_PAPER_REFERENCES[0].a_exp,
                        method,
                    )
                )
                specs.append((f"k-eq {method} C k333", input_path, output_path, False))

            def succeed(
                _cp2k: Path, _input_path: Path, output_path: Path, _threads: int
            ) -> int:
                output_path.write_text("PROGRAM ENDED\n")
                return 0

            with patch.object(base, "run_cp2k", side_effect=succeed):
                eos.run_jobs(
                    specs,
                    cp2k,
                    jobs=1,
                    threads=1,
                    force=False,
                    retry_scf=False,
                    campaign_fingerprint=campaign,
                    campaign_bind_all_methods=True,
                )
            stamps = [json.loads(base.job_stamp_path(spec[2]).read_text()) for spec in specs]
        self.assertEqual(stamps[0]["campaign_identity"], campaign)
        self.assertEqual(stamps[1]["campaign_identity"], campaign)
        self.assertEqual(
            stamps[0]["executable_sha256"], stamps[1]["executable_sha256"]
        )


if __name__ == "__main__":
    unittest.main()
