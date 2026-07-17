from __future__ import annotations

import csv
import hashlib
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Goldzak12" / "data"
HARTREE_TO_KJMOL = 2625.4996394798254
HARTREE_TO_EV = 27.211386245988
EV_TO_KJMOL = HARTREE_TO_KJMOL / HARTREE_TO_EV


class LC10GXTBProvisionalTablesTest(unittest.TestCase):
    def test_paper_artifact_sha256_manifest(self) -> None:
        manifest = json.loads(
            (ROOT / "validation" / "paper_artifact_sha256.json").read_text()
        )
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["algorithm"], "SHA-256")
        self.assertEqual(len(manifest["artifacts"]), 6)
        for artifact in manifest["artifacts"]:
            path = ROOT / artifact["path"]
            self.assertTrue(path.is_file(), artifact["path"])
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                artifact["sha256"],
                artifact["path"],
            )

    def test_mixed_k999_recomputes_exactly_from_per_solid_rows(self) -> None:
        with (DATA / "lc10_gxtb_mixed_k999_provisional.csv").open() as handle:
            records = list(csv.DictReader(handle))
        self.assertEqual(len(records), 10)
        self.assertEqual(
            {row["solid"] for row in records},
            {"C", "Si", "SiC", "BN", "BP", "AlN", "AlP", "MgS", "LiF", "LiCl"},
        )
        self.assertEqual(
            {row["mesh"] for row in records},
            {"k777", "k888", "k999"},
        )

        with (DATA / "lc10_gxtb_kmesh_mae_provisional.csv").open() as handle:
            stages = {row["stage"]: row for row in csv.DictReader(handle)}
        aggregate = stages["mixed_k999"]

        a_errors = [float(row["a0_error_A"]) for row in records]
        e_errors = [float(row["ecoh_error_eV_per_atom"]) for row in records]

        def metrics(values: list[float]) -> tuple[float, float, float, float]:
            return (
                sum(values) / len(values),
                sum(abs(value) for value in values) / len(values),
                math.sqrt(sum(value * value for value in values) / len(values)),
                max(abs(value) for value in values),
            )

        for actual, field in zip(
            metrics(a_errors),
            ("a0_me_A", "a0_mae_A", "a0_rmse_A", "a0_maxae_A"),
            strict=True,
        ):
            self.assertAlmostEqual(actual, float(aggregate[field]), places=11)
        for actual, field in zip(
            metrics(e_errors),
            (
                "ecoh_me_eV_per_atom",
                "ecoh_mae_eV_per_atom",
                "ecoh_rmse_eV_per_atom",
                "ecoh_maxae_eV_per_atom",
            ),
            strict=True,
        ):
            self.assertAlmostEqual(actual, float(aggregate[field]), places=11)
        for actual, field in zip(
            metrics([value * EV_TO_KJMOL for value in e_errors]),
            (
                "ecoh_me_kJmol_per_atom",
                "ecoh_mae_kJmol_per_atom",
                "ecoh_rmse_kJmol_per_atom",
                "ecoh_maxae_kJmol_per_atom",
            ),
            strict=True,
        ):
            self.assertAlmostEqual(actual, float(aggregate[field]), places=9)

    def test_final_adaptive_selection_recomputes_exactly(self) -> None:
        with (DATA / "lc10_gxtb_final_selected_values.csv").open() as handle:
            records = list(csv.DictReader(handle))
        self.assertEqual(len(records), 10)
        self.assertEqual(
            {row["solid"] for row in records},
            {"C", "Si", "SiC", "BN", "BP", "AlN", "AlP", "MgS", "LiF", "LiCl"},
        )
        self.assertEqual(
            {row["selected_mesh"] for row in records},
            {"k777", "k888", "k999"},
        )

        aggregate = json.loads(
            (DATA / "lc10_gxtb_final_aggregate.json").read_text()
        )
        self.assertEqual(aggregate["schema_version"], 2)
        self.assertEqual(aggregate["publication_status"], "final_protocol_converged")
        self.assertEqual(
            aggregate["final_mesh_counts"],
            {"k777": 3, "k888": 4, "k999": 3},
        )

        def metrics(values: list[float]) -> dict[str, float]:
            return {
                "ME": sum(values) / len(values),
                "MAE": sum(abs(value) for value in values) / len(values),
                "RMSE": math.sqrt(sum(value * value for value in values) / len(values)),
                "MaxAE": max(abs(value) for value in values),
            }

        a_errors = [float(row["a0_error_A"]) for row in records]
        e_errors = [float(row["ecoh_error_eV_per_atom"]) for row in records]
        for name, value in metrics(a_errors).items():
            self.assertAlmostEqual(value, aggregate["lattice_A"][name], places=12)
        for name, value in metrics(e_errors).items():
            self.assertAlmostEqual(
                value, aggregate["cohesive_eV_per_atom"][name], places=12
            )
        for name, value in metrics(
            [error * EV_TO_KJMOL for error in e_errors]
        ).items():
            self.assertAlmostEqual(
                value, aggregate["cohesive_kJmol_per_atom"][name], places=10
            )

        self.assertEqual(
            aggregate["convergence_rule"]["threshold_abs_delta_a0_A"], 0.025
        )
        self.assertEqual(
            aggregate["convergence_rule"][
                "threshold_abs_delta_ecoh_kJmol_per_atom"
            ],
            0.25,
        )
        decisions = {row["solid"]: row for row in aggregate["selection_decisions"]}
        self.assertEqual(set(decisions), {row["solid"] for row in records})
        for record in records:
            decision = decisions[record["solid"]]
            self.assertEqual(decision["coarse_mesh"], record["coarse_mesh"])
            self.assertEqual(decision["dense_mesh"], record["selected_mesh"])
            self.assertTrue(decision["preceding_interval_failed"])
            self.assertTrue(decision["lattice_passed"])
            self.assertTrue(decision["energy_passed"])
            self.assertAlmostEqual(
                decision["delta_a0_A"], float(record["delta_a0_A"]), places=12
            )
            self.assertAlmostEqual(
                decision["delta_ecoh_kJmol_per_atom"],
                float(record["delta_ecoh_kJmol_per_atom"]),
                places=10,
            )
            self.assertLessEqual(
                abs(decision["delta_a0_A"]),
                aggregate["convergence_rule"]["threshold_abs_delta_a0_A"],
            )
            self.assertLessEqual(
                abs(decision["delta_ecoh_kJmol_per_atom"]),
                aggregate["convergence_rule"][
                    "threshold_abs_delta_ecoh_kJmol_per_atom"
                ],
            )

    def test_final_selected_outputs_and_archive_manifest_are_hash_bound(self) -> None:
        with (DATA / "lc10_gxtb_final_selected_values.csv").open() as handle:
            records = list(csv.DictReader(handle))
        for record in records:
            output = ROOT / record["source_output"]
            self.assertTrue(output.is_file(), record["source_output"])
            self.assertEqual(
                hashlib.sha256(output.read_bytes()).hexdigest(),
                record["output_sha256"],
                record["solid"],
            )

        archive = ROOT / "validation" / "lc10_gxtb_final_adaptive_20260716"
        entries = []
        for line in (archive / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            relative = name.lstrip("*")
            if relative.startswith("./"):
                relative = relative[2:]
            entries.append((digest, relative))
        expected = {
            path.relative_to(archive).as_posix()
            for path in archive.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual({relative for _, relative in entries}, expected)
        for digest, relative in entries:
            path = archive / relative
            self.assertTrue(path.is_file(), relative)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest, relative)

    def test_final_plot_row_and_gfn2_comparison_match_aggregate(self) -> None:
        aggregate = json.loads(
            (DATA / "lc10_gxtb_final_aggregate.json").read_text()
        )
        with (DATA / "lc10_gxtb_adaptive_kmesh_mae.csv").open() as handle:
            stages = {row["stage"]: row for row in csv.DictReader(handle)}
        final = stages["adaptive_final"]
        self.assertAlmostEqual(
            float(final["a0_mae_A"]), aggregate["lattice_A"]["MAE"], places=11
        )
        self.assertAlmostEqual(
            float(final["ecoh_mae_eV_per_atom"]),
            aggregate["cohesive_eV_per_atom"]["MAE"],
            places=11,
        )
        self.assertEqual(final["mesh_counts"], "k999:3;k888:4;k777:3")

        comparison = json.loads(
            (DATA / "lc10_gxtb_vs_gfn2_current.json").read_text()
        )
        self.assertEqual(
            comparison["publication_status"], "final_gxtb_first_passing_gfn2"
        )
        self.assertEqual(comparison["schema_version"], 2)
        self.assertAlmostEqual(
            comparison["gxtb_final_adaptive"]["a0_error_A"]["mae"],
            aggregate["lattice_A"]["MAE"],
            places=12,
        )
        for source in comparison["sources"].values():
            path = ROOT / source["path"]
            self.assertTrue(path.is_file(), source["path"])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])


if __name__ == "__main__":
    unittest.main()
