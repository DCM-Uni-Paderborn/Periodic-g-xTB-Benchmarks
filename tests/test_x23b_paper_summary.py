from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "X23b" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_x23b_paper_summary as summary  # noqa: E402
import x23b_common as common  # noqa: E402


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("synthetic CSV needs rows")
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class SyntheticX23b:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.data = root / "data"
        self.staging = self.data / "gxtb_staging"
        self.data.mkdir(parents=True)
        metadata = json.loads((REPOSITORY / "X23b" / "data" / "metadata.json").read_text())
        (self.data / "metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n"
        )
        self.metadata = {str(row["id"]): row for row in metadata["systems"]}
        self.systems = tuple(sorted(self.metadata))
        self.preflight_root = root / "runs" / "gxtb_native" / "preflight"
        self.cellopt_root = root / "runs" / "gxtb_native" / "cellopt"
        self.k_roots = {
            "k333": root / "runs" / "gxtb_native" / "k333",
            "k444": root / "runs" / "gxtb_native" / "k444",
        }
        self.final_manifest = self.write_build_manifest("test-direct-acp", "1", "2")
        self.frozen_manifest = self.write_build_manifest("test-frozen-acp", "a", "b")
        self.campaign = summary.campaign_identity_from_manifest(self.final_manifest)
        self.frozen_campaign = summary.campaign_identity_from_manifest(self.frozen_manifest)
        self.gxtb_deltas: list[float] = []

    def write_build_manifest(
        self, campaign_id: str, cp2k_digit: str, save_digit: str
    ) -> Path:
        path = self.root.parent / "campaigns" / campaign_id / "build_manifest.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        manifest = {
            "schema": 1,
            "campaign_id": campaign_id,
            "campaign_state": "production_ready",
            "cp2k": {
                "binary_sha256": cp2k_digit * 64,
                "loaded_library_sha256": chr(ord(cp2k_digit) + 1) * 64,
                "cmake_cache_sha256": chr(ord(cp2k_digit) + 2) * 64,
                "reported_revision": cp2k_digit * 10,
                "revision": cp2k_digit * 40,
            },
            "save_tblite": {
                "cli_sha256": save_digit * 64,
                "revision": save_digit * 40,
                "static_library_sha256": chr(ord(save_digit) + 1) * 64,
                "cmake_cache_sha256": chr(ord(save_digit) + 2) * 64,
            },
            "fetched_dependencies": {"tblite": {"revision": "e" * 40}},
        }
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return path

    @staticmethod
    def output_text(energy: float, *, optimize: bool = False, volume: float | None = None) -> str:
        lines = [f" ENERGY| Total FORCE_EVAL synthetic {energy:.15f}"]
        if volume is not None:
            lines.append(f" CELL| Volume [angstrom^3] {volume:.12f}")
        if optimize:
            lines.append(" GEOMETRY OPTIMIZATION COMPLETED")
        lines.append(" PROGRAM ENDED AT synthetic test")
        return "\n".join(lines) + "\n"

    def job(
        self,
        input_path: Path,
        output_path: Path,
        energy: float,
        campaign: dict[str, object],
        phase: str,
        *,
        optimize: bool = False,
        volume: float | None = None,
        protocol: dict[str, object] | None = None,
        sources: dict[str, Path] | None = None,
        input_text: str = "synthetic CP2K input\n",
    ) -> Path:
        input_path.parent.mkdir(parents=True, exist_ok=True)
        input_path.write_text(input_text)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.output_text(energy, optimize=optimize, volume=volume))
        stamp: dict[str, object] = {
            "schema": common.JOB_STAMP_SCHEMA,
            "method": "GXTB",
            "phase": phase,
            "input": {"path": f"/remote/{input_path.name}", "sha256": summary.sha256(input_path)},
            "cp2k": {"path": "/remote/cp2k.psmp", "sha256": campaign["cp2k_executable_sha256"]},
            "campaign_identity": campaign,
            "status": "converged_synthetic",
            "details": {
                "returncode": 0,
                "output": f"/remote/{output_path.name}",
                "output_sha256": summary.sha256(output_path),
            },
        }
        if protocol is not None:
            stamp["protocol_identity"] = protocol
        if sources is not None:
            stamp["source_artifacts"] = {
                role: {"path": str(path.resolve()), "sha256": summary.sha256(path)}
                for role, path in sorted(sources.items())
            }
        stamp_path = output_path.parent / common.JOB_STAMP_NAME
        stamp_path.write_text(json.dumps(stamp, indent=2, sort_keys=True) + "\n")
        return stamp_path

    def write_legacy(self) -> None:
        final_rows: list[dict[str, object]] = []
        volumes: list[dict[str, object]] = []
        for method_index, method in enumerate(("GFN1", "GFN2"), start=1):
            delta = 0.01 * method_index
            for index, system in enumerate(self.systems):
                meta = self.metadata[system]
                ref = float(meta["ref_energy"])
                error = method_index + 0.01 * index
                k333 = ref + error
                k222 = k333 - 0.05
                k444 = k333 + delta
                final_rows.append(
                    {
                        "method": method,
                        "system": system,
                        "k222_lattice_energy_kJmol": f"{k222:.12f}",
                        "k333_lattice_energy_kJmol": f"{k333:.12f}",
                        "k444_lattice_energy_kJmol": f"{k444:.12f}",
                        "k222_error_kJmol": f"{k222-ref:.12f}",
                        "k333_error_kJmol": f"{error:.12f}",
                        "k444_error_kJmol": f"{k444-ref:.12f}",
                        "delta_k333_minus_k222_kJmol": "0.050000000000",
                        "delta_k444_minus_k333_kJmol": f"{delta:.12f}",
                    }
                )
                ref_volume = float(meta["x23b_same_cell_ref_volume"])
                volume_error = -method_index - 0.01 * index
                volume = ref_volume * (1.0 + volume_error / 100.0)
                volumes.append(
                    {
                        "calculation": "cell_opt",
                        "mesh": "k222",
                        "method": f"{method}-xTB",
                        "system": system,
                        "label": meta["label"],
                        "complete": True,
                        "volume_A3": f"{volume:.12f}",
                        "x23b_same_cell_ref_volume_A3": f"{ref_volume:.12f}",
                        "x23b_reported_ref_volume_A3": f"{ref_volume:.12f}",
                        "volume_error_percent": f"{volume_error:.12f}",
                    }
                )
        write_csv(self.data / "x23b_final_geometry_kpoint_rows.csv", final_rows)
        write_csv(self.data / "x23b_cell_volumes.csv", volumes)
        legacy = {
            "benchmark": "X23b",
            "cp2k": {"source_revision": "c" * 40, "executable_sha256": "d" * 64},
            "tblite": {"main_revision": "e" * 40, "executable_sha256": "f" * 64},
            "repository_patches": {},
        }
        (self.data / "build_provenance.json").write_text(
            json.dumps(legacy, indent=2, sort_keys=True) + "\n"
        )

    def write_production(self) -> None:
        preflight_manifest_rows: list[dict[str, object]] = []
        preflight_csv_rows: list[dict[str, object]] = []
        cell_manifest_rows: list[dict[str, object]] = []
        cell_csv_rows: list[dict[str, object]] = []
        final_manifest_rows: dict[str, list[dict[str, object]]] = {"k333": [], "k444": []}
        final_csv_rows: dict[str, list[dict[str, object]]] = {"k333": [], "k444": []}
        variant = "k222_cellopt_keep_angles_from_experimental_reference"
        lineage = "frozen_x23_reference_structure->experimental_k222_preflight->k222_cellopt_input"
        source_protocol = {
            "source_policy": "experimental_reference",
            "variant": variant,
            "lineage": lineage,
        }
        for index, system in enumerate(self.systems):
            meta = self.metadata[system]
            reference_dir = self.root / "inputs" / "crystal_sp" / "k222" / "GXTB"
            source_input = reference_dir / f"{system}_GXTB_k222_sp.inp"
            structure = self.root / "structures" / "cif" / f"{system}.cif"
            source_input.parent.mkdir(parents=True, exist_ok=True)
            structure.parent.mkdir(parents=True, exist_ok=True)
            source_input.write_text(f"synthetic reference {system}\n")
            structure.write_text(f"synthetic CIF {system}\n")

            pf_dir = self.preflight_root / "GXTB" / system / "experimental_k222_preflight"
            pf_input = pf_dir / f"{system}_GXTB_experimental_k222_preflight.inp"
            pf_output = pf_dir / "cp2k.out"
            pf_sources = {"reference_input": source_input, "reference_structure": structure}
            pf_stamp = self.job(
                pf_input,
                pf_output,
                -5.0 - index,
                self.campaign,
                summary.PREFLIGHT_PHASE,
                protocol=summary.PREFLIGHT_PROTOCOL,
                sources=pf_sources,
            )
            pf_record = {
                "method": "GXTB",
                "system": system,
                "phase": summary.PREFLIGHT_PHASE,
                "variant": "experimental_k222_preflight",
                "source_policy": "experimental_reference",
                "source_input": str(source_input.resolve()),
                "source_input_sha256": summary.sha256(source_input),
                "structure_path": str(structure.resolve()),
                "structure_sha256": summary.sha256(structure),
                "structure_source": meta["structure_source"],
                "input": str(pf_input.resolve()),
                "input_sha256": summary.sha256(pf_input),
                "run_dir": str(pf_dir.resolve()),
                "output": str(pf_output.resolve()),
                "start_volume_A3": meta["input_volume"],
                "atom_count": meta["n_atoms_crystal"],
            }
            preflight_manifest_rows.append(pf_record)
            preflight_csv_rows.append(
                {
                    "method": "GXTB",
                    "system": system,
                    "phase": summary.PREFLIGHT_PHASE,
                    "variant": "experimental_k222_preflight",
                    "source_policy": "experimental_reference",
                    "campaign_fingerprint_sha256": self.campaign["fingerprint_sha256"],
                    "program_ended": True,
                    "scientific_status": "measured_not_approved",
                    "approved": False,
                    "energy_hartree": f"{-5.0-index:.12f}",
                    "max_force_hartree_per_bohr": "0.001",
                    "max_abs_stress_GPa": "0.01",
                    "pressure_GPa": "0.001",
                    "source_input": str(source_input.resolve()),
                    "source_input_sha256": summary.sha256(source_input),
                    "structure_path": str(structure.resolve()),
                    "structure_sha256": summary.sha256(structure),
                    "input": str(pf_input.resolve()),
                    "input_sha256": summary.sha256(pf_input),
                    "output": str(pf_output.resolve()),
                    "output_sha256": summary.sha256(pf_output),
                }
            )

            gas_stem = f"{system}_GXTB_mol_geoopt"
            gas_dir = self.root / "runs" / "molecule_geoopt" / "GXTB" / gas_stem
            gas_input = gas_dir / f"{gas_stem}.inp"
            gas_output = gas_dir / f"{gas_stem}.out"
            gas_energy = -10.0 - 0.1 * index
            self.job(
                gas_input,
                gas_output,
                gas_energy,
                self.campaign,
                "x23b_molecule_geoopt",
                optimize=True,
            )
            (gas_dir / f"{gas_stem}-1.restart").write_text("synthetic gas restart\n")

            ref_energy = float(meta["ref_energy"])
            ref_volume = float(meta["x23b_same_cell_ref_volume"])
            volume_error = 3.0 + 0.01 * index
            volume = ref_volume * (1.0 + volume_error / 100.0)
            lattice222 = ref_energy + 2.95 + 0.01 * index
            lattice333 = ref_energy + 3.00 + 0.01 * index
            lattice444 = lattice333 + 0.01
            self.gxtb_deltas.append(lattice444 - lattice333)
            n_molecules = int(meta["molecules_per_cell"])
            crystal_energies = {
                "k222": n_molecules * (gas_energy - lattice222 / summary.HARTREE_TO_KJMOL),
                "k333": n_molecules * (gas_energy - lattice333 / summary.HARTREE_TO_KJMOL),
                "k444": n_molecules * (gas_energy - lattice444 / summary.HARTREE_TO_KJMOL),
            }
            cell_dir = self.cellopt_root / "GXTB" / system / variant
            cell_input = cell_dir / f"{system}_GXTB_k222.inp"
            cell_output = cell_dir / "cp2k.out"
            cell_sources = {
                "reference_input": source_input,
                "reference_structure": structure,
                "preflight_input": pf_input,
                "preflight_output": pf_output,
                "preflight_stamp": pf_stamp,
            }
            self.job(
                cell_input,
                cell_output,
                crystal_energies["k222"],
                self.campaign,
                summary.CELL_OPT_PHASE,
                optimize=True,
                volume=volume,
                protocol=source_protocol,
                sources=cell_sources,
            )
            cell_restart = cell_dir / f"{system}_GXTB_k222-1.restart"
            cell_restart.write_text("synthetic cell restart\n")
            cell_manifest = {
                "schema": 2,
                "method": "GXTB",
                "system": system,
                "variant": variant,
                "source_policy": "experimental_reference",
                "source_kind": "experimental_reference",
                "lineage": lineage,
                "source_path": str(source_input.resolve()),
                "source_sha256": summary.sha256(source_input),
                "structure_path": str(structure.resolve()),
                "structure_sha256": summary.sha256(structure),
                "preflight_input": str(pf_input.resolve()),
                "preflight_input_sha256": summary.sha256(pf_input),
                "preflight_output": str(pf_output.resolve()),
                "preflight_output_sha256": summary.sha256(pf_output),
                "preflight_stamp": str(pf_stamp.resolve()),
                "preflight_stamp_sha256": summary.sha256(pf_stamp),
                "input": str(cell_input.resolve()),
                "input_sha256": summary.sha256(cell_input),
                "run_dir": str(cell_dir.resolve()),
            }
            cell_manifest_rows.append(cell_manifest)
            cell_csv_rows.append(
                {
                    "method": "GXTB",
                    "system": system,
                    "source_policy": "experimental_reference",
                    "preflight_input_sha256": summary.sha256(pf_input),
                    "preflight_output_sha256": summary.sha256(pf_output),
                    "preflight_stamp_sha256": summary.sha256(pf_stamp),
                    "program_ended": True,
                    "opt_completed": True,
                    "max_iter_reached": False,
                    "energy_hartree": f"{crystal_energies['k222']:.12f}",
                    "gas_energy_hartree": f"{gas_energy:.12f}",
                    "lattice_energy_kJmol": f"{lattice222:.12f}",
                    "volume_A3": f"{volume:.12f}",
                    "volume_error_percent": f"{volume_error:.12f}",
                    "output": str(cell_output.resolve()),
                }
            )
            source_hashes = {role: summary.sha256(path) for role, path in cell_sources.items()}
            for mesh, target_lattice in (("k333", lattice333), ("k444", lattice444)):
                target = int(mesh[1])
                target_dir = self.k_roots[mesh] / "GXTB" / system / f"{mesh}_sp_on_k222"
                target_input = target_dir / f"{system}_GXTB_{mesh}.inp"
                target_output = target_dir / "cp2k.out"
                target_protocol = {
                    "source_policy": "experimental_reference",
                    "source_variant": variant,
                    "source_protocol_identity": source_protocol,
                    "target_mesh": mesh,
                }
                target_sources = {
                    **cell_sources,
                    "cellopt_input": cell_input,
                    "cellopt_output": cell_output,
                    "cellopt_restart": cell_restart,
                }
                self.job(
                    target_input,
                    target_output,
                    crystal_energies[mesh],
                    self.campaign,
                    f"x23b_final_k{target}{target}{target}_on_k222",
                    protocol=target_protocol,
                    sources=target_sources,
                )
                final_manifest_rows[mesh].append(
                    {
                        "method": "GXTB",
                        "system": system,
                        "mesh": target,
                        "source_run_dir": str(cell_dir.resolve()),
                        "source_input": str(cell_input.resolve()),
                        "source_input_sha256": summary.sha256(cell_input),
                        "source_output_sha256": summary.sha256(cell_output),
                        "source_policy": "experimental_reference",
                        "source_variant": variant,
                        "source_protocol_identity": json.dumps(source_protocol, sort_keys=True),
                        "source_artifacts": json.dumps(
                            {role: str(path.resolve()) for role, path in cell_sources.items()},
                            sort_keys=True,
                        ),
                        "source_artifact_hashes": json.dumps(source_hashes, sort_keys=True),
                        "source_restart": str(cell_restart.resolve()),
                        "source_restart_sha256": summary.sha256(cell_restart),
                        "input": str(target_input.resolve()),
                        "input_sha256": summary.sha256(target_input),
                        "run_dir": str(target_dir.resolve()),
                    }
                )
                final_csv_rows[mesh].append(
                    {
                        "method": "GXTB",
                        "system": system,
                        "target_mesh": mesh,
                        "source_mesh": "k222_cellopt",
                        "program_ended": True,
                        "source_energy_hartree": f"{crystal_energies['k222']:.12f}",
                        "target_energy_hartree": f"{crystal_energies[mesh]:.12f}",
                        "source_lattice_energy_kJmol": f"{lattice222:.12f}",
                        "target_lattice_energy_kJmol": f"{target_lattice:.12f}",
                        "delta_target_minus_source_kJmol": f"{target_lattice-lattice222:.12f}",
                        "target_error_kJmol": f"{target_lattice-ref_energy:.12f}",
                        "output": str(target_output.resolve()),
                    }
                )
        preflight_manifest = {
            "schema": 1,
            "phase": summary.PREFLIGHT_PHASE,
            "variant": "experimental_k222_preflight",
            "source_policy": "experimental_reference",
            "campaign_identity": self.campaign,
            "systems": preflight_manifest_rows,
        }
        (self.preflight_root / "experimental_k222_preflight_manifest.json").parent.mkdir(
            parents=True, exist_ok=True
        )
        (self.preflight_root / "experimental_k222_preflight_manifest.json").write_text(
            json.dumps(preflight_manifest, indent=2, sort_keys=True) + "\n"
        )
        write_csv(self.staging / "x23b_experimental_k222_preflight.csv", preflight_csv_rows)
        write_csv(self.cellopt_root.parent / "x23b_k222_cellopt_manifest.csv", cell_manifest_rows)
        write_csv(self.staging / "x23b_k222_cellopt_results.csv", cell_csv_rows)
        for mesh in ("k333", "k444"):
            write_csv(self.k_roots[mesh] / "manifest.csv", final_manifest_rows[mesh])
            write_csv(self.staging / f"x23b_{mesh}_results.csv", final_csv_rows[mesh])

    def write_fd_gate(self) -> None:
        fd_root = self.root / "runs" / "gxtb_native" / "fd"
        cases: list[dict[str, object]] = []
        protocol = {
            "systems": list(summary.FD_SYSTEMS),
            "mesh": "MACDONALD 2 2 2 0.25 0.25 0.25",
            "coordinate_step_bohr": 0.001,
            "coordinate_direction_count": 2,
            "strain_step": 0.0005,
        }
        job_specs = [
            ("baseline", "baseline_energy_force", None, None),
            ("coord_d1_minus", "coordinate_energy", "d1", -1),
            ("coord_d1_plus", "coordinate_energy", "d1", 1),
            ("coord_d2_minus", "coordinate_energy", "d2", -1),
            ("coord_d2_plus", "coordinate_energy", "d2", 1),
            ("strain_s1_minus", "strain_energy", "s1", -1),
            ("strain_s1_plus", "strain_energy", "s1", 1),
            ("strain_s2_minus", "strain_energy", "s2", -1),
            ("strain_s2_plus", "strain_energy", "s2", 1),
        ]
        for system in summary.FD_SYSTEMS:
            source = self.root / "inputs" / "crystal_sp" / "k222" / "GXTB" / f"{system}_GXTB_k222_sp.inp"
            structure = self.root / "structures" / "cif" / f"{system}.cif"
            jobs: list[dict[str, object]] = []
            for job_id, job_type, direction, sign in job_specs:
                run_dir = fd_root / "GXTB" / system / "frozen_reference_shifted_k222_spglib_fd" / job_id
                input_path = run_dir / f"{system}_{job_id}.inp"
                input_path.parent.mkdir(parents=True, exist_ok=True)
                input_path.write_text(f"synthetic FD {system}/{job_id}\n")
                jobs.append(
                    {
                        "job_id": job_id,
                        "job_type": job_type,
                        "direction_id": direction,
                        "direction_sha256": None,
                        "generator": None,
                        "sign": sign,
                        "input": str(input_path.resolve()),
                        "input_sha256": summary.sha256(input_path),
                        "output": str((run_dir / "cp2k.out").resolve()),
                        "run_dir": str(run_dir.resolve()),
                    }
                )
            cases.append(
                {
                    "system": system,
                    "source_input": str(source.resolve()),
                    "source_input_sha256": summary.sha256(source),
                    "structure_path": str(structure.resolve()),
                    "structure_sha256": summary.sha256(structure),
                    "jobs": jobs,
                }
            )
        manifest: dict[str, object] = {
            "schema": 1,
            "phase": summary.FD_PHASE,
            "variant": "frozen_reference_shifted_k222_spglib_fd",
            "source_policy": "experimental_reference",
            "scientific_status": "prepared_not_measured",
            "campaign_identity": self.campaign,
            "protocol": protocol,
            "cases": cases,
        }
        manifest["payload_sha256"] = summary.fingerprint(manifest)
        manifest_path = fd_root / "x23b_k222_fd_gate_manifest.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        for case in cases:
            system = str(case["system"])
            sources = {
                "fd_manifest": manifest_path,
                "reference_input": Path(str(case["source_input"])),
                "reference_structure": Path(str(case["structure_path"])),
            }
            for job_record in case["jobs"]:  # type: ignore[index]
                assert isinstance(job_record, dict)
                job_protocol = {
                    "gate_phase": summary.FD_PHASE,
                    "variant": "frozen_reference_shifted_k222_spglib_fd",
                    "source_policy": "experimental_reference",
                    "manifest_payload_sha256": manifest["payload_sha256"],
                    "mesh": protocol["mesh"],
                    "coordinate_step_bohr": protocol["coordinate_step_bohr"],
                    "strain_step": protocol["strain_step"],
                    "job_id": job_record["job_id"],
                    "job_type": job_record["job_type"],
                    "direction_id": job_record.get("direction_id"),
                    "direction_sha256": job_record.get("direction_sha256"),
                    "generator": job_record.get("generator"),
                    "sign": job_record.get("sign"),
                }
                self.job(
                    Path(str(job_record["input"])),
                    Path(str(job_record["output"])),
                    -20.0,
                    self.campaign,
                    summary.FD_PHASE,
                    protocol=job_protocol,
                    sources=sources,
                    input_text=f"synthetic FD {system}/{job_record['job_id']}\n",
                )
        measured_rows: list[dict[str, object]] = []
        for system in summary.FD_SYSTEMS:
            for measurement, direction in (
                ("coordinate_directional_derivative", "d1"),
                ("coordinate_directional_derivative", "d2"),
                ("symmetric_strain_directional_derivative", "s1"),
                ("symmetric_strain_directional_derivative", "s2"),
            ):
                measured_rows.append(
                    {
                        "method": "GXTB",
                        "system": system,
                        "measurement_type": measurement,
                        "direction_id": direction,
                        "energy_derivative_error_hartree_per_parameter": "0.000001",
                        "stress_conjugation_error_GPa": "0.000001",
                    }
                )
        measured_path = self.staging / "x23b_k222_fd_measured.csv"
        write_csv(measured_path, measured_rows)
        report = {
            "schema": 1,
            "phase": summary.FD_PHASE,
            "variant": "frozen_reference_shifted_k222_spglib_fd",
            "scientific_status": "measured_not_approved",
            "approved": False,
            "campaign_identity": self.campaign,
            "manifest": {"path": str(manifest_path.resolve()), "sha256": summary.sha256(manifest_path)},
            "measured_csv": {"path": str(measured_path.resolve()), "sha256": summary.sha256(measured_path)},
            "systems": list(summary.FD_SYSTEMS),
            "row_count": len(measured_rows),
            "rows": measured_rows,
        }
        report_path = self.staging / "x23b_k222_fd_measured.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        checks = [
            {
                "system": row["system"],
                "direction_id": row["direction_id"],
                "passed": True,
            }
            for row in measured_rows
        ]
        approval = {
            "schema": 1,
            "phase": summary.FD_PHASE,
            "decision": "approved",
            "reviewer": "Synthetic Reviewer",
            "report_json": {"path": str(report_path.resolve()), "sha256": summary.sha256(report_path)},
            "measured_csv": report["measured_csv"],
            "manifest": report["manifest"],
            "thresholds": {
                "coordinate_abs_tolerance_hartree_per_bohr": 1.0e-4,
                "stress_abs_tolerance_gpa": 1.0e-3,
            },
            "checks": checks,
            "passed_count": len(checks),
            "check_count": len(checks),
        }
        (self.staging / "x23b_k222_fd_approval.json").write_text(
            json.dumps(approval, indent=2, sort_keys=True) + "\n"
        )

    def write_cross_build(self) -> None:
        rows: list[dict[str, object]] = []
        deltas: list[float] = []
        for index, system in enumerate(self.systems):
            outputs: dict[str, Path] = {}
            stamps: dict[str, Path] = {}
            inputs: dict[str, Path] = {}
            final_energy = -100.0 - index
            frozen_energy = final_energy - 1.0e-8
            for role, campaign, phase, energy in (
                ("final", self.campaign, summary.CROSS_BUILD_FINAL_JOB_PHASE, final_energy),
                ("frozen", self.frozen_campaign, summary.CROSS_BUILD_FROZEN_JOB_PHASE, frozen_energy),
            ):
                run_dir = self.root / "runs" / "cross_build" / role / system
                inputs[role] = run_dir / f"{system}_{role}.inp"
                outputs[role] = run_dir / "cp2k.out"
                stamps[role] = self.job(
                    inputs[role], outputs[role], energy, campaign, phase
                )
            delta = final_energy - frozen_energy
            deltas.append(delta)
            row: dict[str, object] = {
                "method": "GXTB",
                "system": system,
                "delta_final_minus_frozen_hartree": f"{delta:.15e}",
                "passed": True,
            }
            for role in ("final", "frozen"):
                for name, path in (
                    ("input", inputs[role]),
                    ("output", outputs[role]),
                    ("stamp", stamps[role]),
                ):
                    row[f"{role}_{name}"] = str(path.resolve())
                    row[f"{role}_{name}_sha256"] = summary.sha256(path)
            rows.append(row)
        comparison = self.staging / "x23b_direct_acp_cross_build.csv"
        write_csv(comparison, rows)
        approval = {
            "schema": 1,
            "phase": summary.CROSS_BUILD_PHASE,
            "decision": "approved",
            "reviewer": "Synthetic Reviewer",
            "comparison_csv": {"path": str(comparison.resolve()), "sha256": summary.sha256(comparison)},
            "final_campaign_manifest": {
                "path": str(self.final_manifest.resolve()),
                "sha256": summary.sha256(self.final_manifest),
            },
            "frozen_campaign_manifest": {
                "path": str(self.frozen_manifest.resolve()),
                "sha256": summary.sha256(self.frozen_manifest),
            },
            "absolute_energy_tolerance_hartree": 1.0e-6,
            "max_abs_delta_hartree": max(abs(value) for value in deltas),
            "checks": [{"system": system, "passed": True} for system in self.systems],
        }
        (self.staging / "x23b_direct_acp_cross_build_approval.json").write_text(
            json.dumps(approval, indent=2, sort_keys=True) + "\n"
        )

    def write_kpoint_approval(self) -> None:
        checks = []
        for method, delta in (("GFN1", 0.01), ("GFN2", 0.02), ("GXTB", 0.01)):
            checks.append(
                {
                    "method": method,
                    "N": 23,
                    "mean_abs_change_kJmol": delta,
                    "max_abs_change_kJmol": delta,
                    "mean_abs_tolerance_kJmol": 0.1,
                    "max_abs_tolerance_kJmol": 1.0,
                    "passed": True,
                }
            )
        approval = {
            "schema": 1,
            "phase": summary.KPOINT_APPROVAL_PHASE,
            "decision": "approved",
            "reviewer": "Synthetic Reviewer",
            "checks": checks,
        }
        (self.staging / "x23b_k333_k444_convergence_approval.json").write_text(
            json.dumps(approval, indent=2, sort_keys=True) + "\n"
        )

    def write_gxtb_provenance(self) -> None:
        validation: dict[str, object] = {}
        for name in (
            "gas_optimizations",
            "experimental_k222_preflight",
            "k222_cell_optimizations",
            "k333_single_points",
            "k444_single_points",
        ):
            validation[f"{name}_expected"] = 23
            validation[f"{name}_completed"] = 23
        provenance = {
            "benchmark": "X23b",
            "method": "GXTB",
            "status": "production_complete",
            "campaign_identity": self.campaign,
            "campaign_manifest": {
                "path": "/remote/final/build_manifest.json",
                "file_sha256": summary.sha256(self.final_manifest),
            },
            "cp2k": {"executable_sha256": self.campaign["cp2k_executable_sha256"]},
            "save_tblite": {"executable_sha256": self.campaign["save_tblite_executable_sha256"]},
            "workflow_paths": {
                "experimental_k222_preflight_root": str(self.preflight_root.resolve()),
                "k222_cellopt_root": str(self.cellopt_root.resolve()),
                "k222_source_policy": "experimental_reference",
                "final_single_point_roots": {
                    mesh: str(path.resolve()) for mesh, path in self.k_roots.items()
                },
            },
            "validation": validation,
        }
        (self.data / "build_provenance_gxtb.json").write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n"
        )

    def complete(self) -> None:
        self.write_legacy()
        self.write_production()
        self.write_fd_gate()
        self.write_cross_build()
        self.write_kpoint_approval()
        self.write_gxtb_provenance()

    def args(self) -> argparse.Namespace:
        return argparse.Namespace(
            root=self.root,
            legacy_rows=None,
            legacy_volumes=None,
            preflight_csv=None,
            fd_report=None,
            fd_approval=None,
            cellopt_csv=None,
            k333_csv=None,
            k444_csv=None,
            cross_build_approval=None,
            kpoint_approval=None,
            output_csv=None,
            output_json=None,
        )


class X23bPaperSummaryTests(unittest.TestCase):
    def fixture(self, directory: Path) -> SyntheticX23b:
        fixture = SyntheticX23b(directory / "X23b")
        fixture.complete()
        return fixture

    def test_complete_fixture_writes_one_publication_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            csv_path, json_path, payload = summary.finalize(fixture.args())
            self.assertTrue(csv_path.is_file())
            self.assertTrue(json_path.is_file())
            self.assertEqual(payload["publication_status"], "publication_ready")
            self.assertEqual(payload["coverage"]["common"], 23)
            self.assertEqual(len(read_csv_for_test(csv_path)), 6)
            self.assertEqual(
                payload["gates"]["direct_acp_cross_build"]["coverage_passed"], 23
            )
            self.assertEqual(
                payload["gates"]["experimental_k222_preflight"]["coverage_accepted"], 23
            )

    def test_unapproved_fd_removes_stale_publication_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            out_csv = fixture.data / f"{summary.SUMMARY_STEM}.csv"
            out_json = fixture.data / f"{summary.SUMMARY_STEM}.json"
            out_csv.write_text("stale\n")
            out_json.write_text("stale\n")
            approval_path = fixture.staging / "x23b_k222_fd_approval.json"
            approval = json.loads(approval_path.read_text())
            approval["decision"] = "rejected"
            approval_path.write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "finite-difference pilot gate"):
                summary.finalize(fixture.args())
            self.assertFalse(out_csv.exists())
            self.assertFalse(out_json.exists())
            self.assertEqual(list(fixture.data.glob(f".*{summary.SUMMARY_STEM}*.tmp")), [])

    def test_tampered_k444_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            rows = read_csv_for_test(fixture.staging / "x23b_k444_results.csv")
            Path(rows[0]["output"]).write_text("tampered\n")
            with self.assertRaisesRegex(ValueError, "incomplete CP2K output|fingerprint"):
                summary.finalize(fixture.args())
            self.assertFalse((fixture.data / f"{summary.SUMMARY_STEM}.json").exists())

    def test_missing_preflight_system_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            path = fixture.staging / "x23b_experimental_k222_preflight.csv"
            rows = read_csv_for_test(path)
            write_csv(path, rows[:-1])
            with self.assertRaisesRegex(ValueError, "preflight is not exactly 23/23"):
                summary.finalize(fixture.args())

    def test_cross_build_table_tamper_is_rejected_against_raw_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            approval_path = fixture.staging / "x23b_direct_acp_cross_build_approval.json"
            approval = json.loads(approval_path.read_text())
            comparison = Path(approval["comparison_csv"]["path"])
            rows = read_csv_for_test(comparison)
            rows[0]["delta_final_minus_frozen_hartree"] = "9.0e-7"
            write_csv(comparison, rows)
            # Even a coordinated rewrite of the approval's table pointer must
            # still fail against the independently stamped raw energies.
            approval["comparison_csv"]["sha256"] = summary.sha256(comparison)
            approval_path.write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "cross-build .* delta mismatch"):
                summary.finalize(fixture.args())

    def test_fd_csv_tamper_is_rejected_against_frozen_json_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            approval_path = fixture.staging / "x23b_k222_fd_approval.json"
            approval = json.loads(approval_path.read_text())
            report_path = Path(approval["report_json"]["path"])
            report = json.loads(report_path.read_text())
            measured = Path(report["measured_csv"]["path"])
            rows = read_csv_for_test(measured)
            rows[0]["energy_derivative_error_hartree_per_parameter"] = "0.9"
            write_csv(measured, rows)
            report["measured_csv"]["sha256"] = summary.sha256(measured)
            report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            approval["measured_csv"] = report["measured_csv"]
            approval["report_json"]["sha256"] = summary.sha256(report_path)
            approval_path.write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "FD JSON rows do not exactly reproduce"):
                summary.finalize(fixture.args())

    def test_kpoint_approval_tamper_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            approval_path = fixture.staging / "x23b_k333_k444_convergence_approval.json"
            approval = json.loads(approval_path.read_text())
            approval["checks"][2]["mean_abs_change_kJmol"] = 0.0
            approval_path.write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "approved GXTB k-point mean mismatch"):
                summary.finalize(fixture.args())

    def test_second_atomic_replace_failure_leaves_no_partial_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = self.fixture(Path(temporary))
            original = summary.os.replace
            calls = 0

            def fail_second(source: Path, target: Path) -> None:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("synthetic second replace failure")
                original(source, target)

            with patch.object(summary.os, "replace", side_effect=fail_second):
                with self.assertRaisesRegex(OSError, "synthetic second replace failure"):
                    summary.finalize(fixture.args())
            self.assertFalse((fixture.data / f"{summary.SUMMARY_STEM}.csv").exists())
            self.assertFalse((fixture.data / f"{summary.SUMMARY_STEM}.json").exists())


def read_csv_for_test(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    unittest.main()
