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
            fixed_metrics = {
                "ME": -10.0 * index,
                "MAE": 10.0 * index,
                "RMSE": 10.0 * index + 0.1,
                "MaxAE": 10.0 * index + 0.2,
            }
            fixed_status = (
                "numerically_unconverged_same_mesh_comparator"
                if method == "GXTB"
                else "same_mesh_comparator"
            )
            methods[method] = {
                "method_label": method,
                "status": "phasewise_kpoint_converged",
                "n_nonreference_phases": 12,
                "metrics_kjmol_per_h2o": metrics,
                "fixed_k333_same_mesh_comparison": {
                    "mesh": "k333",
                    "status": fixed_status,
                    "phasewise_kpoint_converged_value": False,
                    "metrics_kjmol_per_h2o": fixed_metrics,
                },
            }
            csv_rows.append(
                {
                    "method_id": method,
                    "N_nonreference_phases": 12,
                    **{f"{name}_kJmol_per_H2O": value for name, value in metrics.items()},
                    "fixed_k333_status": fixed_status,
                    **{
                        f"fixed_k333_{name}_kJmol_per_H2O": value
                        for name, value in fixed_metrics.items()
                    },
                }
            )
        (data / "dmc_ice13_gfn_gxtb_phasewise_summary.json").write_text(
            json.dumps(
                {
                    "benchmark": "DMC-ICE13",
                    "status": "phasewise_kpoint_converged",
                    "n_nonreference_phases": 12,
                    "fixed_k333_same_mesh_comparison": {
                        "mesh": "k333",
                        "not_a_phasewise_converged_result": True,
                    },
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
        common_systems = ";".join(bundle.LC10_SYSTEMS)
        all_systems = ";".join(
            (*bundle.LC10_SYSTEMS, *bundle.LC10_EXCLUDED_SYSTEMS)
        )
        for index, method in enumerate(bundle.METHODS, start=1):
            for scope, n in (
                ("method_available_coverage", 12 if method == "GFN1" else 10),
                ("three_method_common_subset", 10),
            ):
                is_common = scope == bundle.LC10_CHILD_SCOPE
                metric_value = float(index if is_common else index + 100)
                row: dict[str, object] = {
                    "method_id": method,
                    "method_label": method,
                    "scope": scope,
                    "n_systems": n,
                    "coverage_denominator": 12,
                    "systems": (
                        common_systems
                        if is_common or n == len(bundle.LC10_SYSTEMS)
                        else all_systems
                    ),
                    "eos_mesh": "k444",
                    "result_mesh": "k555",
                }
                for prefix, suffix in (("lattice", "A"), ("cohesive", "eV_per_atom")):
                    row.update(
                        {
                            f"{prefix}_ME_{suffix}": -metric_value,
                            f"{prefix}_MAE_{suffix}": metric_value,
                            f"{prefix}_RMSE_{suffix}": metric_value + 0.1,
                            f"{prefix}_MaxAE_{suffix}": metric_value + 0.2,
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
                    "protocol": {
                        "common_subset_count": len(bundle.LC10_SYSTEMS),
                        "common_subset_systems": list(bundle.LC10_SYSTEMS),
                    },
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
            self.assertEqual(payload["schema_version"], 3)
            self.assertEqual(payload["benchmarks"], ["DMC-ICE13", "X23b", "LC12"])
            self.assertEqual(
                payload["lc10_scope"],
                {
                    "scope_id": bundle.LC10_OUTPUT_SCOPE,
                    "systems": list(bundle.LC10_SYSTEMS),
                    "n_systems": len(bundle.LC10_SYSTEMS),
                    "three_method_comparison_only": True,
                    "method_available_coverage_exported": False,
                    "excluded_systems": list(bundle.LC10_EXCLUDED_SYSTEMS),
                    "excluded_systems_note": (
                        "LC10 is the fixed three-method comparison set; LiH and "
                        "MgO are outside its scope."
                    ),
                },
            )
            self.assertEqual(len(payload["rows"]), 18)
            self.assertEqual(len(payload["gxtb_vs_gfn_baseline_comparisons"]), 12)
            lc_rows = [
                row for row in payload["rows"] if row["benchmark"] == "LC12"
            ]
            self.assertEqual(len(lc_rows), 6)
            self.assertEqual(
                {(row["scope"], row["N"]) for row in lc_rows},
                {(bundle.LC10_OUTPUT_SCOPE, len(bundle.LC10_SYSTEMS))},
            )
            self.assertEqual(
                {row["systems"] for row in lc_rows},
                {";".join(bundle.LC10_SYSTEMS)},
            )
            self.assertEqual(
                {(row["method_id"], row["quantity"]) for row in lc_rows},
                {
                    (method, quantity)
                    for method in bundle.METHODS
                    for quantity in ("lattice_constant", "cohesive_energy")
                },
            )
            # The method-available fixtures deliberately carry 100+ errors.
            # Their absence proves that only exact common-10 statistics escape.
            self.assertTrue(all(float(row["MAE"]) < 10.0 for row in lc_rows))
            dmc = next(
                comparison
                for comparison in payload["gxtb_vs_gfn_baseline_comparisons"]
                if comparison["benchmark"] == "DMC-ICE13"
                and comparison["scope"] == "phasewise_kpoint_converged"
                and comparison["baseline_method"] == "GFN2"
            )
            self.assertAlmostEqual(dmc["MAE_delta_GXTB_minus_baseline"], 1.0)
            self.assertAlmostEqual(dmc["MAE_ratio_GXTB_over_baseline"], 1.5)
            self.assertAlmostEqual(dmc["MAE_percent_change_GXTB_vs_baseline"], 50.0)
            with csv_path.open() as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual(len(csv_rows), 18)
            self.assertFalse(
                any(row["scope"] == "method_available_coverage" for row in csv_rows)
            )
            self.assertEqual(
                {
                    row["systems"]
                    for row in csv_rows
                    if row["benchmark"] == "LC12"
                },
                {";".join(bundle.LC10_SYSTEMS)},
            )
            tex = tex_path.read_text()
            self.assertIn("\\newcommand{\\GXTBLC10N}{10}", tex)
            self.assertIn(
                "\\newcommand{\\GXTBLC10Systems}{C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF, LiCl}",
                tex,
            )
            self.assertIn("GXTBvsGFN2MAEPercentChange", tex)
            self.assertIn("GXTBvsGFN1MAEPercentChange", tex)
            self.assertNotIn("MethodAvailableCoverage", tex)
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

    def test_lc10_rejects_an_alternative_ten_system_subset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            json_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.json"
            csv_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.csv"
            payload = json.loads(json_path.read_text())
            with csv_path.open() as handle:
                csv_rows = list(csv.DictReader(handle))
            alternative = ";".join((*bundle.LC10_SYSTEMS[:-1], "MgO"))
            for row in payload["summary_rows"]:
                if row["method_id"] == "GXTB" and row["scope"] == "three_method_common_subset":
                    row["systems"] = alternative
            for row in csv_rows:
                if row["method_id"] == "GXTB" and row["scope"] == "three_method_common_subset":
                    row["systems"] = alternative
            write_csv(csv_path, csv_rows)
            payload["paper_summary_csv"]["sha256"] = digest(csv_path)
            json_path.write_text(json.dumps(payload, indent=2) + "\n")
            # The child JSON/CSV pair agrees and is hash-valid, and still must not
            # redefine LC10 as another arbitrary ten-system subset.
            with self.assertRaisesRegex(ValueError, "exact LC10 set"):
                bundle.finalize(root, root / "paper")

    def test_lc10_protocol_order_is_pinned_and_failure_removes_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            output = root / "paper"
            output.mkdir()
            for suffix in ("csv", "json", "tex"):
                (output / f"{bundle.OUTPUT_STEM}.{suffix}").write_text("stale\n")
            json_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.json"
            payload = json.loads(json_path.read_text())
            payload["protocol"]["common_subset_systems"] = list(
                reversed(bundle.LC10_SYSTEMS)
            )
            json_path.write_text(json.dumps(payload, indent=2) + "\n")
            with self.assertRaisesRegex(ValueError, "pin the exact LC10 systems"):
                bundle.finalize(root, output)
            for suffix in ("csv", "json", "tex"):
                self.assertFalse((output / f"{bundle.OUTPUT_STEM}.{suffix}").exists())

    def test_lc10_rejects_duplicate_child_csv_scope_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repository"
            self.write_complete(root)
            json_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.json"
            csv_path = root / "Goldzak12" / "data" / "lc12_gfn_gxtb_paper_summary.csv"
            with csv_path.open() as handle:
                csv_rows = list(csv.DictReader(handle))
            csv_rows.append(dict(csv_rows[-1]))
            write_csv(csv_path, csv_rows)
            payload = json.loads(json_path.read_text())
            payload["paper_summary_csv"]["sha256"] = digest(csv_path)
            json_path.write_text(json.dumps(payload, indent=2) + "\n")
            with self.assertRaisesRegex(ValueError, "incomplete or duplicate"):
                bundle.finalize(root, root / "paper")

if __name__ == "__main__":
    unittest.main()
