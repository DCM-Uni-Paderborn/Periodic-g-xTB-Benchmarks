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


if __name__ == "__main__":
    unittest.main()
