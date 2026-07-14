from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]


def load_script():
    path = REPOSITORY / "scripts" / "finalize_paper_benchmark_bundle.py"
    spec = importlib.util.spec_from_file_location("finalize_paper_benchmark_bundle", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


bundle = load_script()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PaperBenchmarkBundleTests(unittest.TestCase):
    def write_dmc(self, root: Path) -> None:
        data = root / "DMC-ICE13" / "data"
        data.mkdir(parents=True)
        methods: dict[str, object] = {}
        csv_rows: list[dict[str, object]] = []
        for index, method in enumerate(bundle.METHODS, start=1):
            metrics = {
                "ME": -float(index),
                "MAE": float(index),
                "RMSE": float(index) + 0.1,
                "MaxAE": float(index) + 0.2,
            }
            methods[method] = {
                "method_label": method,
                "status": "phasewise_kpoint_converged",
                "n_nonreference_phases": 12,
                "metrics_kjmol_per_h2o": metrics,
            }
            csv_rows.append(
                {
                    "method_id": method,
                    "N_nonreference_phases": 12,
                    **{f"{name}_kJmol_per_H2O": value for name, value in metrics.items()},
                }
            )
        (data / "dmc_ice13_gfn_gxtb_phasewise_summary.json").write_text(
            json.dumps(
                {
                    "benchmark": "DMC-ICE13",
                    "status": "phasewise_kpoint_converged",
                    "n_nonreference_phases": 12,
                    "methods": methods,
                },
                indent=2,
            )
            + "\n"
        )
        write_csv(data / "dmc_ice13_gfn_gxtb_phasewise_summary.csv", csv_rows)

    def write_x23b(self, root: Path) -> None:
        data = root / "X23b" / "data"
        data.mkdir(parents=True)
        rows: list[dict[str, object]] = []
        for index, method in enumerate(bundle.METHODS, start=1):
            for quantity, calculation, mesh in (
                ("lattice_energy_kJmol", "cell_opt_single_point", "k333"),
                ("volume_error_percent", "cell_opt", "k222"),
            ):
                rows.append(
                    {
                        "method": method,
                        "method_label": method,
                        "quantity": quantity,
                        "calculation": calculation,
                        "mesh": mesh,
                        "N": 23,
                        "ME": f"{-index:.12f}",
                        "MAE": f"{index:.12f}",
                        "RMSE": f"{index + 0.1:.12f}",
                        "MaxAE": f"{index + 0.2:.12f}",
                    }
                )
        csv_path = data / "x23b_gfn_gxtb_paper_summary.csv"
        write_csv(csv_path, rows)
        (data / "x23b_gfn_gxtb_paper_summary.json").write_text(
            json.dumps(
                {
                    "benchmark": "X23b",
                    "publication_status": "publication_ready",
                    "methods": list(bundle.METHODS),
                    "coverage": {"exact_common_coverage": True, "common": 23},
                    "summary": rows,
                    "publication_csv_sha256": digest(csv_path),
                },
                indent=2,
            )
            + "\n"
        )

    def write_lc12(self, root: Path) -> None:
        data = root / "Goldzak12" / "data"
        data.mkdir(parents=True)
        rows: list[dict[str, object]] = []
        systems = ";".join(f"S{index}" for index in range(10))
        for index, method in enumerate(bundle.METHODS, start=1):
            for scope, n in (
                ("method_available_coverage", 12 if method == "GFN1" else 10),
                ("three_method_common_subset", 10),
            ):
                row: dict[str, object] = {
                    "method_id": method,
                    "method_label": method,
                    "scope": scope,
                    "n_systems": n,
                    "coverage_denominator": 12,
                    "systems": systems if scope == "three_method_common_subset" else systems,
                    "eos_mesh": "k444",
                    "result_mesh": "k555",
                }
                for prefix, suffix in (("lattice", "A"), ("cohesive", "eV_per_atom")):
                    row.update(
                        {
                            f"{prefix}_ME_{suffix}": -float(index),
                            f"{prefix}_MAE_{suffix}": float(index),
                            f"{prefix}_RMSE_{suffix}": float(index) + 0.1,
                            f"{prefix}_MaxAE_{suffix}": float(index) + 0.2,
                        }
                    )
                rows.append(row)
        csv_path = data / "lc12_gfn_gxtb_paper_summary.csv"
        write_csv(csv_path, rows)
        (data / "lc12_gfn_gxtb_paper_summary.json").write_text(
            json.dumps(
                {
                    "benchmark": "LC12 (Goldzak12)",
                    "status": "publication_ready_reduced_coverage",
                    "methods": {method: {} for method in bundle.METHODS},
                    "protocol": {"common_subset_count": 10},
                    "summary_rows": rows,
                    "paper_summary_csv": {"sha256": digest(csv_path)},
                },
                indent=2,
            )
            + "\n"
        )

    def write_complete(self, root: Path) -> None:
        self.write_dmc(root)
        self.write_x23b(root)
        self.write_lc12(root)

    def test_complete_three_benchmark_bundle_is_atomic_and_paper_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            output = root / "paper"
            csv_path, json_path, tex_path = bundle.finalize(root, output)
            payload = json.loads(json_path.read_text())
            self.assertEqual(payload["status"], "publication_ready")
            self.assertEqual(payload["benchmarks"], ["DMC-ICE13", "X23b", "LC12"])
            self.assertEqual(len(payload["rows"]), 21)
            with csv_path.open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 21)
            self.assertIn("\\GXTB", tex_path.read_text())
            self.assertEqual(
                payload["generated_outputs"]["csv_sha256"], digest(csv_path)
            )
            self.assertEqual(
                payload["generated_outputs"]["tex_sha256"], digest(tex_path)
            )

    def test_tampered_child_csv_removes_all_stale_bundle_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            output = root / "paper"
            output.mkdir()
            for suffix in ("csv", "json", "tex"):
                (output / f"{bundle.OUTPUT_STEM}.{suffix}").write_text("stale\n")
            x23_csv = root / "X23b" / "data" / "x23b_gfn_gxtb_paper_summary.csv"
            x23_csv.write_text(x23_csv.read_text() + "tampered\n")
            with self.assertRaisesRegex(ValueError, "X23b publication CSV hash"):
                bundle.finalize(root, output)
            for suffix in ("csv", "json", "tex"):
                self.assertFalse((output / f"{bundle.OUTPUT_STEM}.{suffix}").exists())

    def test_incomplete_dmc_never_emits_a_cross_benchmark_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            dmc = root / "DMC-ICE13" / "data" / "dmc_ice13_gfn_gxtb_phasewise_summary.json"
            payload = json.loads(dmc.read_text())
            payload["methods"].pop("GXTB")
            dmc.write_text(json.dumps(payload, indent=2) + "\n")
            with self.assertRaisesRegex(ValueError, "lacks the three methods"):
                bundle.finalize(root, root / "paper")

    def test_lc12_common_subset_must_be_identical_for_all_methods(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            json_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.json"
            csv_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.csv"
            payload = json.loads(json_path.read_text())
            for row in payload["summary_rows"]:
                if row["method_id"] == "GXTB" and row["scope"] == "three_method_common_subset":
                    row["systems"] = "different"
            json_path.write_text(json.dumps(payload, indent=2) + "\n")
            # The child CSV remains internally hash-valid; the mismatch is between its
            # systems and the JSON/common-method contract and must still be rejected.
            with self.assertRaisesRegex(ValueError, "common-subset systems differ"):
                bundle.finalize(root, root / "paper")


if __name__ == "__main__":
    unittest.main()
