from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "Goldzak12" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_goldzak12_paper_summary as summary  # noqa: E402
import run_goldzak12_benchmark as base  # noqa: E402
import run_goldzak12_eos_benchmark as eos  # noqa: E402
import run_goldzak12_k_convergence as kconv  # noqa: E402


class SyntheticLC10:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.data = root / "data"
        self.data.mkdir(parents=True)
        self.scales = (0.94, 0.98, 1.00, 1.02, 1.06)
        self.manifest_path = (
            root.parent / "campaigns" / "test-lc10" / "build_manifest.json"
        )
        self.manifest_path.parent.mkdir(parents=True)
        self.manifest = {
            "campaign_id": "test-lc10",
            "campaign_state": "production_ready",
            "cp2k": {
                "binary_sha256": "1" * 64,
                "loaded_library_sha256": "2" * 64,
                "cmake_cache_sha256": "3" * 64,
                "reported_revision": "4" * 10,
                "revision": "4" * 40,
            },
            "save_tblite": {
                "cli_sha256": "5" * 64,
                "revision": "6" * 40,
                "static_library_sha256": "7" * 64,
                "cmake_cache_sha256": "8" * 64,
            },
            "fetched_dependencies": {
                "tblite": {"revision": "9" * 40},
            },
        }
        self.manifest_path.write_text(
            json.dumps(self.manifest, indent=2, sort_keys=True) + "\n"
        )
        self.campaign = base.campaign_identity_from_manifest(
            self.manifest, self.manifest_path
        )

    @staticmethod
    def cp2k_text(energy: float) -> str:
        return (
            f" Total energy (extrapolated to T->0) {energy:.15f}\n"
            " PROGRAM ENDED AT synthetic test\n"
        )

    def write_stamp(
        self,
        result: Path,
        input_path: Path,
        executable_role: str,
        command_contract: dict[str, object] | None = None,
    ) -> None:
        executable_field = {
            "cp2k": "cp2k_executable_sha256",
            "save_tblite": "save_tblite_executable_sha256",
        }[executable_role]
        payload = {
            "schema_version": 1,
            "completed": True,
            "return_code": 0,
            "campaign_identity": self.campaign,
            "input": f"/remote/host/{input_path.name}",
            "input_sha256": summary.sha256(input_path),
            "executable": f"/remote/host/{executable_role}",
            "executable_sha256": self.campaign[executable_field],
            "command_contract": command_contract or {"driver": executable_role},
        }
        stamp = base.job_stamp_path(result)
        stamp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    def write_cp2k_job(
        self,
        input_path: Path,
        output_path: Path,
        energy: float,
        *,
        stamped: bool,
        input_text: str = "synthetic CP2K input\n",
    ) -> None:
        input_path.parent.mkdir(parents=True, exist_ok=True)
        input_path.write_text(input_text)
        output_path.write_text(self.cp2k_text(energy))
        if stamped:
            self.write_stamp(output_path, input_path, "cp2k")

    def write_atoms(self) -> dict[tuple[str, str], float]:
        energies: dict[tuple[str, str], float] = {}
        rows: list[dict[str, object]] = []
        elements = base.LC10_PAPER_ELEMENTS
        for method_index, method in enumerate(summary.METHODS, start=1):
            for element_index, element in enumerate(elements, start=1):
                energy = -10.0 * method_index - 0.1 * element_index
                energies[(method, element)] = energy
                run_dir = self.root / "runs" / "atoms_cli" / method / element
                run_dir.mkdir(parents=True)
                xyz = run_dir / f"atom_{element}.xyz"
                result = run_dir / f"atom_{element}_{method}.json"
                stdout = result.with_suffix(".out")
                xyz.write_text(
                    f"1\n{element} atom\n{element} 0.0 0.0 0.0\n"
                )
                result.write_text(json.dumps({"energy": energy}) + "\n")
                stdout.write_text("synthetic atom calculation\n")
                self.write_stamp(
                    result,
                    xyz,
                    "save_tblite",
                    {
                        "driver": "tblite_run",
                        "method": base.method_cli_name(method),
                        "spin_2S": base.ELEMENT_MULTIPLICITY[element] - 1,
                        "accuracy": 0.05,
                        "restart": False,
                    },
                )
                rows.append(
                    {
                        "method": method,
                        "element": element,
                        "energy_hartree": f"{energy:.12f}",
                        "source": "save_tblite_cli",
                        "multiplicity": base.ELEMENT_MULTIPLICITY[element],
                        "spin_2S": base.ELEMENT_MULTIPLICITY[element] - 1,
                    }
                )
        base.write_csv(
            self.data / "atom_energies_tblite_cli.csv",
            [row for row in rows if row["method"] in {"GFN1", "GFN2"}],
        )
        base.write_csv(
            self.data / "atom_energies_save_tblite_cli_gxtb.csv",
            [row for row in rows if row["method"] == "GXTB"],
        )

        check_rows: list[dict[str, object]] = []
        for element in elements:
            cli_energy = energies[("GXTB", element)]
            cp2k_energy = cli_energy + 1.0e-12
            run_dir = self.root / "runs" / "atoms" / "GXTB" / element
            input_path = run_dir / f"atom_{element}_GXTB.inp"
            output_path = run_dir / f"atom_{element}_GXTB.out"
            self.write_cp2k_job(
                input_path, output_path, cp2k_energy, stamped=True
            )
            check_rows.append(
                {
                    "method": "GXTB",
                    "element": element,
                    "multiplicity": base.ELEMENT_MULTIPLICITY[element],
                    "spin_2S": base.ELEMENT_MULTIPLICITY[element] - 1,
                    "cp2k_energy_hartree": f"{cp2k_energy:.15f}",
                    "cli_energy_hartree": f"{cli_energy:.15f}",
                    "delta_cp2k_minus_cli_hartree": f"{cp2k_energy-cli_energy:.15e}",
                    "tolerance_hartree": "1.0e-6",
                    "passed": True,
                    "cli_provider": "save_tblite",
                    "cp2k_scf_contract": "synthetic",
                    "cohesive_energy_atom_reference": "save_tblite_cli_only",
                    "campaign_stamp_issue": "",
                }
            )
        base.write_csv(
            self.data / "atom_reference_cp2k_vs_save_tblite_gxtb.csv",
            check_rows,
        )
        return energies

    def write_eos(
        self,
        atom_energies: dict[tuple[str, str], float],
    ) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        fits: list[dict[str, object]] = []
        points: list[dict[str, object]] = []
        values: list[dict[str, object]] = []
        scale_records: list[dict[str, object]] = []
        mesh_shifts = {
            "k333": (0.0020, 0.0010),
            "k444": (0.0002, 0.0002),
            "k555": (0.0001, 0.0001),
        }
        for solid_index, ref in enumerate(base.LC10_PAPER_REFERENCES):
            for method_index, method in enumerate(summary.METHODS, start=1):
                for mesh, (a_shift, ecoh_shift) in mesh_shifts.items():
                    target = ref.a_exp + 0.01 * method_index + a_shift
                    numeric: list[tuple[float, float, float, bool]] = []
                    scale_records.append(
                        {
                            "method": method,
                            "solid": ref.solid,
                            "mesh": mesh,
                            "requested_scales": list(self.scales),
                        }
                    )
                    for scale in self.scales:
                        a_value = ref.a_exp * scale
                        energy = (
                            -50.0
                            - solid_index
                            - method_index
                            - 0.01 * kconv.mesh_number(mesh)
                            + 0.25 * (a_value - target) ** 2
                        )
                        numeric.append((a_value, scale, energy, True))
                        project = eos.eos_project(ref.solid, method, mesh, scale)
                        run_dir = (
                            self.root
                            / "runs"
                            / "eos"
                            / method
                            / ref.solid
                            / mesh
                            / eos.scale_tag(scale, method)
                        )
                        input_path = run_dir / f"{project}.inp"
                        output_path = run_dir / f"{project}.out"
                        self.write_cp2k_job(
                            input_path,
                            output_path,
                            energy,
                            stamped=True,
                            input_text=base.solid_input(
                                ref,
                                method,
                                "ENERGY",
                                mesh,
                                a_value,
                                project,
                            ),
                        )
                        points.append(
                            {
                                "solid": ref.solid,
                                "method": method,
                                "mesh": mesh,
                                "scale": f"{scale:.5f}",
                                "a_A": f"{a_value:.10f}",
                                "energy_hartree": f"{energy:.12f}",
                                "completed": True,
                                "valid_for_eos": True,
                                "diagnostic": "",
                                "classification_resolution": "",
                                "classification_rationale": "",
                                "scf_strategy": (
                                    "native_gxtb_fdiis"
                                    if method == "GXTB"
                                    else "default_tblite_mixer"
                                ),
                            }
                        )
                    fitted = (
                        eos.fit_gxtb_eos(numeric)
                        if method == "GXTB"
                        else eos.fit_eos(numeric)
                    )
                    fit = {
                        "solid": ref.solid,
                        "structure": ref.structure,
                        "method": method,
                        "eos_mesh": mesh,
                        "a_exp_A": ref.a_exp,
                        "n_requested": len(self.scales),
                        "n_completed": len(self.scales),
                        "n_converged_raw": len(self.scales),
                        "n_charge_collapsed": 0,
                        "n_explicit_excluded": 0,
                        "n_unresolved_branch_candidates": 0,
                        **fitted,
                    }
                    fits.append(fit)

                    counts = base.atom_counts(ref)
                    atom_sum = sum(
                        atom_energies[(method, element)] * count
                        for element, count in counts.items()
                    )
                    cohesive = (
                        ref.ecoh_exp
                        + 0.10 * method_index
                        + 0.01 * solid_index
                        + ecoh_shift
                    )
                    solid_energy = atom_sum - cohesive * len(
                        base.conventional_cell_atoms(ref)
                    ) / base.HARTREE_TO_EV
                    eq_input, eq_output, eq_lineage_path = kconv.equilibrium_paths(
                        self.root, ref.solid, method, mesh
                    )
                    self.write_cp2k_job(
                        eq_input,
                        eq_output,
                        solid_energy,
                        stamped=True,
                        input_text=base.solid_input(
                            ref,
                            method,
                            "ENERGY",
                            mesh,
                            float(fit["a_eos_A"]),
                            kconv.equilibrium_project(ref.solid, method, mesh),
                        ),
                    )
                    eq_lineage = {
                        "schema_version": kconv.EQUILIBRIUM_LINEAGE_SCHEMA,
                        "benchmark": "LC10 independent-EOS k convergence",
                        "solid": ref.solid,
                        "method": method,
                        "eos_mesh": mesh,
                        "energy_mesh": mesh,
                        "a_eos_A": str(fit["a_eos_A"]),
                        "fit_status": "quadratic",
                        "fit_record_sha256": hashlib.sha256(
                            json.dumps(
                                {key: str(value) for key, value in fit.items()},
                                sort_keys=True,
                                separators=(",", ":"),
                            ).encode()
                        ).hexdigest(),
                        "input_sha256": summary.sha256(eq_input),
                        "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
                    }
                    eq_lineage_path.write_text(
                        json.dumps(eq_lineage, indent=2, sort_keys=True) + "\n"
                    )
                    values.append(
                        {
                            "solid": ref.solid,
                            "structure": ref.structure,
                            "method": method,
                            "mesh": mesh,
                            "mesh_n": kconv.mesh_number(mesh),
                            "a0_A": f"{float(fit['a_eos_A']):.10f}",
                            "solid_energy_hartree": f"{solid_energy:.12f}",
                            "ecoh_eV_per_atom": f"{cohesive:.12f}",
                            "fit_energy_hartree": str(fit["energy_fit_hartree"]),
                            "fit_rmse_hartree": str(fit["fit_rmse_hartree"]),
                            "fit_status": "quadratic",
                            "n_eos_points": str(fit["n_completed"]),
                            "input_sha256": summary.sha256(eq_input),
                            "output_sha256": summary.sha256(eq_output),
                            "lineage_sha256": summary.sha256(eq_lineage_path),
                            "campaign_stamp_sha256": summary.sha256(
                                base.job_stamp_path(eq_output)
                            ),
                            "value_source": (
                                "single_point_at_own_independent_eos_minimum"
                            ),
                            "atom_reference_source": "save_tblite_cli",
                        }
                    )
        base.write_csv(self.data / "eos_points.csv", points)
        base.write_csv(self.data / "eos_fits.csv", fits)
        scale_records.sort(
            key=lambda row: (
                summary.METHODS.index(str(row["method"])),
                base.LC10_PAPER_SOLIDS.index(str(row["solid"])),
                kconv.mesh_number(str(row["mesh"])),
            )
        )
        (self.data / kconv.SCALE_MANIFEST_NAME).write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "benchmark": "LC10 independent-EOS k convergence",
                    "methods": list(summary.METHODS),
                    "paper_systems": list(base.LC10_PAPER_SOLIDS),
                    "records": scale_records,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        values.sort(
            key=lambda row: (
                summary.METHODS.index(str(row["method"])),
                base.LC10_PAPER_SOLIDS.index(str(row["solid"])),
                int(row["mesh_n"]),
            )
        )
        steps, selections, pending = kconv.assess_convergence(values)
        self.assert_no_pending(pending)
        kconv.write_convergence_artifacts(
            self.root,
            values,
            steps,
            selections,
            pending,
            campaign=self.campaign,
            fits_sha256=summary.sha256(self.data / "eos_fits.csv"),
        )
        return fits, points

    @staticmethod
    def assert_no_pending(pending: list[tuple[str, str, int]]) -> None:
        if pending:
            raise AssertionError(f"synthetic adaptive fixture did not converge: {pending}")

    def write_results(
        self,
        fits: list[dict[str, object]],
        atom_energies: dict[tuple[str, str], float],
    ) -> None:
        refs = {ref.solid: ref for ref in base.REFERENCES}
        rows: list[dict[str, object]] = []
        for fit in fits:
            if fit.get("eos_mesh") != kconv.FIXED_GEOMETRY_EOS_MESH:
                continue
            method = str(fit["method"])
            solid = str(fit["solid"])
            if not summary.fit_is_valid(fit, method):
                continue
            ref = refs[solid]
            a_calc = float(fit["a_eos_A"])
            atom_sum = sum(
                atom_energies[(method, element)] * count
                for element, count in base.atom_counts(ref).items()
            )
            n_atoms = len(base.conventional_cell_atoms(ref))
            method_index = summary.METHODS.index(method) + 1
            solid_index = list(refs).index(solid)
            cohesive = ref.ecoh_exp + 0.10 * method_index + 0.01 * solid_index
            solid_energy = atom_sum - cohesive * n_atoms / base.HARTREE_TO_EV
            for mesh in kconv.FIXED_GEOMETRY_ENERGY_MESHES:
                project = eos.final_project(solid, method, mesh)
                run_dir = self.root / "runs" / "eos_final_sp" / method / solid / mesh
                input_path = run_dir / f"{project}.inp"
                output_path = run_dir / f"{project}.out"
                self.write_cp2k_job(
                    input_path,
                    output_path,
                    solid_energy,
                    stamped=True,
                    input_text=base.solid_input(
                        ref,
                        method,
                        "ENERGY",
                        mesh,
                        a_calc,
                        project,
                    ),
                )
                if method == "GXTB":
                    lineage = {
                        "schema_version": eos.FINAL_INPUT_LINEAGE_SCHEMA,
                        "benchmark": "LC10 (fixed Goldzak12 subset)",
                        "valid": True,
                        "reason": "synthetic test",
                        "solid": solid,
                        "method": "GXTB",
                        "eos_mesh": kconv.FIXED_GEOMETRY_EOS_MESH,
                        "energy_mesh": mesh,
                        "fit_status": "quadratic",
                        "a_eos_A": str(fit["a_eos_A"]),
                        "input": f"/remote/host/{input_path.name}",
                        "input_sha256": summary.sha256(input_path),
                        "kpoint_mesh_contract": base.KPOINT_MESH_CONTRACT,
                    }
                    eos.final_input_lineage_path(input_path).write_text(
                        json.dumps(lineage, indent=2, sort_keys=True) + "\n"
                    )
                rows.append(
                    {
                        "solid": solid,
                        "structure": ref.structure,
                        "method": method,
                        "eos_mesh": kconv.FIXED_GEOMETRY_EOS_MESH,
                        "energy_mesh": mesh,
                        "fit_status": fit["fit_status"],
                        "sp_completed": True,
                        "sp_scf_strategy": "synthetic",
                        "a_calc_A": f"{a_calc:.8f}",
                        "a_ref_exp_A": f"{ref.a_exp:.8f}",
                        "a_error_A": f"{a_calc-ref.a_exp:.8f}",
                        "a_abs_error_A": f"{abs(a_calc-ref.a_exp):.8f}",
                        "ecoh_calc_eV_per_atom": f"{cohesive:.8f}",
                        "ecoh_ref_exp_eV_per_atom": f"{ref.ecoh_exp:.8f}",
                        "ecoh_error_eV_per_atom": f"{cohesive-ref.ecoh_exp:.8f}",
                        "ecoh_abs_error_eV_per_atom": f"{abs(cohesive-ref.ecoh_exp):.8f}",
                        "solid_energy_hartree": f"{solid_energy:.12f}",
                        "atom_reference_source": "save_tblite_cli",
                    }
                )
        base.write_csv(self.data / "eos_results.csv", rows)

    def write_provenance(self, fits: list[dict[str, object]]) -> None:
        base.write_csv(
            self.data / "reference_goldzak2022.csv", base.reference_rows()
        )
        scale_path = self.data / kconv.SCALE_MANIFEST_NAME
        approved_hash = kconv.fit_fingerprint(fits)
        convergence_path = self.data / kconv.CONVERGENCE_NAME
        protocol = {
            "benchmark": "LC10 (fixed Goldzak12 subset)",
            "methods": list(summary.METHODS),
            "selected_solids": list(base.LC10_PAPER_SOLIDS),
            "paper_systems": list(base.LC10_PAPER_SOLIDS),
            "diagnostic_only_systems": list(base.LC10_DIAGNOSTIC_ONLY_SOLIDS),
            "exact_lc10_scope": True,
            "single_cp2k_binary_for_all_methods": True,
            "single_tblite_provider_for_all_methods": "save_tblite",
            "campaign_manifest_sha256_external_pin": summary.sha256(
                self.manifest_path
            ),
            "adaptive_k_convergence": {
                "initial_meshes": ["k333", "k444", "k555"],
                "scientific_maximum_mesh": None,
                "technical_resource_guard_mesh": None,
                "technical_resource_guard_is_convergence": False,
                "required_consecutive_passing_steps": 1,
                "aggregate_rms_gate": False,
                "criteria_combination": "AND",
                "lattice_abs_delta_threshold_A": kconv.LATTICE_THRESHOLD_A,
                "cohesive_abs_delta_threshold_kJmol_per_atom": (
                    kconv.ECOH_THRESHOLD_KJMOL_PER_ATOM
                ),
                "cohesive_abs_delta_threshold_eV_per_atom": (
                    kconv.ECOH_THRESHOLD_EV_PER_ATOM
                ),
                "selected_value_policy": (
                    "take denser n+1 value from earliest passing step"
                ),
                "equilibrium_energy_protocol": (
                    "single point at each mesh's own independent EOS minimum"
                ),
            },
            "fixed_geometry_single_point_diagnostic": {
                "separate_from_adaptive_selection": True,
                "eos_mesh": kconv.FIXED_GEOMETRY_EOS_MESH,
                "energy_meshes": list(kconv.FIXED_GEOMETRY_ENERGY_MESHES),
            },
            "fit_approval_required": True,
            "fit_approved": True,
            "approved_fit_sha256": approved_hash,
            "current_fit_sha256": approved_hash,
            "k_convergence_artifact_sha256": summary.sha256(convergence_path),
            "k_convergence_scale_manifest_sha256": summary.sha256(scale_path),
        }
        legacy = {
            "cp2k": {
                "sha256": self.campaign["cp2k_executable_sha256"],
                "source": {"revision": "4" * 40},
            },
            "tblite": {
                "sha256": self.campaign["save_tblite_executable_sha256"],
                "source": {"revision": "6" * 40},
            },
            "repository_patches": {},
            "protocol": protocol,
        }
        (self.data / "build_provenance.json").write_text(
            json.dumps(legacy, indent=2, sort_keys=True) + "\n"
        )
        gxtb = {
            "campaign_identity": self.campaign,
            "campaign_manifest": {
                "path": "/remote/host/build_manifest.json",
                "file_sha256": summary.sha256(self.manifest_path),
                "campaign_id": "test-lc10",
                "campaign_state": "production_ready",
            },
            "cp2k": {
                "sha256": self.campaign["cp2k_executable_sha256"]
            },
            "save_tblite": {
                "sha256": self.campaign["save_tblite_executable_sha256"]
            },
            "protocol": protocol,
        }
        (self.data / "build_provenance_gxtb.json").write_text(
            json.dumps(gxtb, indent=2, sort_keys=True) + "\n"
        )

    def complete(self) -> None:
        atom_energies = self.write_atoms()
        fits, _ = self.write_eos(atom_energies)
        self.write_results(fits, atom_energies)
        self.write_provenance(fits)

class Goldzak12PaperSummaryTests(unittest.TestCase):
    def test_mesh_limit_contract_has_no_scientific_cap(self) -> None:
        uncapped = {
            "scientific_maximum_mesh": None,
            "technical_resource_guard_mesh": None,
            "technical_resource_guard_is_convergence": False,
        }
        summary.validate_adaptive_mesh_limit_contract(uncapped, "test")
        guarded = {**uncapped, "technical_resource_guard_mesh": "k121212"}
        summary.validate_adaptive_mesh_limit_contract(guarded, "test")
        with self.assertRaisesRegex(ValueError, "scientific maximum"):
            summary.validate_adaptive_mesh_limit_contract(
                {**uncapped, "maximum_mesh": "k888"}, "test"
            )

    def test_finalizes_exact_three_by_ten_bundle_with_tex_and_raw_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            csv_path, json_path, tex_path = summary.finalize(root)

            with csv_path.open(newline="") as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual(len(csv_rows), 3)
            self.assertEqual({row["n_systems"] for row in csv_rows}, {"10"})
            self.assertEqual(
                {row["systems"] for row in csv_rows},
                {";".join(base.LC10_PAPER_SOLIDS)},
            )
            self.assertIn("lattice_MaxAE_A", csv_rows[0])
            self.assertIn("cohesive_MaxAE_eV_per_atom", csv_rows[0])
            payload = json.loads(json_path.read_text())
            self.assertEqual(payload["status"], "publication_ready")
            self.assertEqual(payload["coverage"]["systems"], list(base.LC10_PAPER_SOLIDS))
            self.assertEqual(payload["coverage"]["common"], 10)
            self.assertEqual(payload["protocol"]["diagnostic_only_systems"], ["LiH", "MgO"])
            self.assertNotIn("LiH", payload["methods"]["GXTB"]["systems"])
            self.assertNotIn("MgO", payload["methods"]["GXTB"]["systems"])
            self.assertIn("GXTB_vs_GFN1", payload["gxtb_vs_gfn_baseline_comparisons"])
            self.assertIn("GXTB_vs_GFN2", payload["gxtb_vs_gfn_baseline_comparisons"])
            self.assertEqual(payload["paper_summary_csv"]["sha256"], summary.sha256(csv_path))
            self.assertEqual(payload["paper_summary_tex"]["sha256"], summary.sha256(tex_path))
            tex = tex_path.read_text()
            self.assertIn("\\providecommand{\\LCtenN}{10}", tex)
            self.assertIn("\\LCtenGxTBvsGfnTwoLatticeMAEPercentChange", tex)
            c_lineage = payload["methods"]["GXTB"]["systems"]["C"]
            raw_output = root / c_lineage["independent_eos_meshes"]["k555"][
                "equilibrium_output"
            ]["path"]
            self.assertEqual(
                c_lineage["independent_eos_meshes"]["k555"][
                    "equilibrium_output"
                ]["sha256"],
                summary.sha256(raw_output),
            )
            self.assertEqual(
                payload["methods"]["GXTB"]["systems"]["C"]["reporting_status"],
                "reported_at_denser_value_of_earliest_single_passing_adjacent_step",
            )
            self.assertEqual(
                payload["methods"]["GXTB"]["systems"]["C"][
                    "adaptive_selection"
                ]["selected_mesh"],
                "k555",
            )
            self.assertFalse(
                payload["protocol"]["adaptive_selection"]["aggregate_rms_gate"]
            )

    def test_unapproved_gxtb_removes_stale_publication_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            provenance_path = root / "data" / "build_provenance_gxtb.json"
            provenance = json.loads(provenance_path.read_text())
            provenance["protocol"]["fit_approved"] = False
            provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
            legacy_path = root / "data" / "build_provenance.json"
            legacy = json.loads(legacy_path.read_text())
            legacy["protocol"]["fit_approved"] = False
            legacy_path.write_text(json.dumps(legacy, indent=2) + "\n")
            csv_path = root / "data" / f"{summary.SUMMARY_STEM}.csv"
            json_path = root / "data" / f"{summary.SUMMARY_STEM}.json"
            tex_path = root / "data" / f"{summary.SUMMARY_STEM}.tex"
            csv_path.write_text("stale\n")
            json_path.write_text("{}\n")
            tex_path.write_text("stale\n")

            with self.assertRaisesRegex(ValueError, "not explicitly approved"):
                summary.finalize(root)
            self.assertFalse(csv_path.exists())
            self.assertFalse(json_path.exists())
            self.assertFalse(tex_path.exists())
            self.assertEqual(list((root / "data").glob("*.tmp.*")), [])

    def test_raw_output_tamper_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            output = (
                root
                / "runs"
                / "eos_final_sp"
                / "GXTB"
                / "C"
                / "k555"
                / "C_GXTB_eos_final_k555.out"
            )
            output.write_text(case.cp2k_text(-999.0))
            with self.assertRaisesRegex(ValueError, "fixed SP energy mismatch"):
                summary.finalize(root)

    def test_atomic_triple_cleanup_if_second_replace_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            csv_path = root / "data" / f"{summary.SUMMARY_STEM}.csv"
            json_path = root / "data" / f"{summary.SUMMARY_STEM}.json"
            tex_path = root / "data" / f"{summary.SUMMARY_STEM}.tex"
            real_replace = summary.os.replace
            calls = 0

            def fail_second(source: Path, destination: Path) -> None:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("synthetic second replace failure")
                real_replace(source, destination)

            with patch.object(summary.os, "replace", side_effect=fail_second):
                with self.assertRaisesRegex(OSError, "second replace failure"):
                    summary.finalize(root)
            self.assertFalse(csv_path.exists())
            self.assertFalse(json_path.exists())
            self.assertFalse(tex_path.exists())
            self.assertEqual(list((root / "data").glob("*.tmp.*")), [])

    def test_missing_one_of_the_fixed_ten_is_fatal_and_removes_stale_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            fit_path = root / "data" / "eos_fits.csv"
            fits = [
                row
                for row in base.read_csv(fit_path)
                if not (row["method"] == "GXTB" and row["solid"] == "LiCl")
            ]
            base.write_csv(fit_path, fits)
            outputs = [
                root / "data" / f"{summary.SUMMARY_STEM}.{suffix}"
                for suffix in ("csv", "json", "tex")
            ]
            for output in outputs:
                output.write_text("stale\n")
            with self.assertRaisesRegex(ValueError, "fingerprint|initial meshes|LC10"):
                summary.finalize(root)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_diagnostic_multistart_artifacts_are_not_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Goldzak12"
            case = SyntheticLC10(root)
            case.complete()
            self.assertFalse((root / "runs" / "eos" / "GXTB" / "LiH").exists())
            self.assertFalse((root / "runs" / "eos" / "GXTB" / "MgO").exists())
            _, json_path, _ = summary.finalize(root)
            payload = json.loads(json_path.read_text())
            self.assertEqual(payload["coverage"]["common"], 10)
            self.assertIn("not publication prerequisites", payload["protocol"]["diagnostic_note"])


if __name__ == "__main__":
    unittest.main()
