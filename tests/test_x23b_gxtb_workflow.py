from __future__ import annotations

import csv
import fcntl
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "X23b" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import x23b_common as common  # noqa: E402
import x23b_final_kpoint_sp as final_sp  # noqa: E402
import x23b_k222_force_stress_gate as derivative_gate  # noqa: E402
import x23b_experimental_k222_preflight as experimental_preflight  # noqa: E402
import x23b_kpoint_cellopt as cellopt  # noqa: E402
import x23b_pipeline as pipeline  # noqa: E402


def fake_executable(directory: Path, name: str = "cp2k") -> Path:
    path = directory / name
    path.write_text("#!/bin/sh\necho 'fake CP2K 2026.1'\n")
    path.chmod(path.stat().st_mode | 0o111)
    return path


def fake_campaign(cp2k: Path, *, discriminator: str = "1") -> dict[str, object]:
    return common.make_campaign_identity(
        campaign_id="test-campaign",
        cp2k_executable_sha256=common.sha256_file(cp2k),
        cp2k_loaded_library_sha256="2" * 64,
        cp2k_cmake_cache_sha256="3" * 64,
        cp2k_embedded_source_revision="0123456789",
        cp2k_source_revision="0123456789" + "0" * 30,
        save_tblite_executable_sha256=discriminator * 64,
        save_tblite_source_revision=discriminator * 40,
        save_tblite_library_sha256=discriminator * 64,
        save_tblite_cmake_cache_sha256="4" * 64,
        dependency_lock_sha256="5" * 64,
    )


def clean_git_source(directory: Path, name: str) -> tuple[Path, str]:
    source = directory / name
    source.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=source, check=True)
    subprocess.run(["git", "config", "user.name", "Test"], cwd=source, check=True)
    subprocess.run(["git", "config", "user.email", "test@example.invalid"], cwd=source, check=True)
    (source / "source.txt").write_text(name + "\n")
    subprocess.run(["git", "add", "source.txt"], cwd=source, check=True)
    subprocess.run(["git", "commit", "-q", "-m", "source"], cwd=source, check=True)
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True, capture_output=True, check=True
    ).stdout.strip()
    return source, revision


def fake_campaign_artifacts(directory: Path) -> tuple[Path, Path, Path, Path, Path]:
    cp2k_source, revision = clean_git_source(directory, "cp2k-source")
    save_source, save_revision = clean_git_source(directory, "save-source")
    cp2k = directory / "cp2k"
    cp2k.write_text(
        "#!/bin/sh\n"
        "echo 'CP2K version 2026.1'\n"
        f"echo ' Source code revision {revision[:10]}'\n"
    )
    cp2k.chmod(cp2k.stat().st_mode | 0o111)
    tblite = fake_executable(directory, "tblite")
    library = directory / "libtblite.a"
    library.write_bytes(b"static save_tblite test archive")
    libcp2k = directory / "libcp2k.dylib"
    libcp2k.write_bytes(b"test libcp2k")
    manifest = directory / "build_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "campaign_id": "test-campaign",
                "campaign_state": "production_ready",
                "cp2k": {
                    "revision": revision,
                    "source_clean": True,
                    "binary_sha256": common.sha256_file(cp2k),
                    "reported_revision": revision[:10],
                    "loaded_library": str(libcp2k),
                    "loaded_library_sha256": common.sha256_file(libcp2k),
                    "cmake_cache_sha256": "3" * 64,
                },
                "save_tblite": {
                    "revision": save_revision,
                    "source_clean": True,
                    "reported_version": "",
                    "cli_sha256": common.sha256_file(tblite),
                    "static_library": str(library),
                    "static_library_sha256": common.sha256_file(library),
                    "cmake_cache_sha256": "4" * 64,
                },
                "fetched_dependencies": {"dependency": "a" * 40},
            },
            indent=2,
        )
        + "\n"
    )
    return cp2k, cp2k_source, tblite, save_source, manifest


def fake_derivative_output(
    *,
    energy: float = -10.0,
    force_offset: float = 0.0,
    stress_offset: float = 0.0,
    kpoints: int = 4,
    include_kpoint_report: bool = True,
    elements: list[str] | None = None,
) -> str:
    elements = elements or (["N"] * 4 + ["H"] * 12)
    lines: list[str] = []
    if include_kpoint_report:
        lines += [
            f"       Number of Special K-points: {kpoints:45d}",
            "       K-point Mesh:                                               2     2     2",
            "                   Wavevector Basis                   Special Points    Rotation",
        ]
        for index in range(8):
            x, y, z = ((index >> shift) & 1 for shift in (2, 1, 0))
            lines.append(
                f" {index + 1:10d} {0.5 * x:9.5f} {0.5 * y:9.5f} {0.5 * z:9.5f}"
                f" {index + 1:14d} {min(index + 1, kpoints):11d} {1:11d}"
            )
    lines += [
        f" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] {energy:.12f}",
        "",
        " FORCES| Atomic forces [hartree/bohr]",
        " FORCES|   Atom     x                   y                   z                   |f|",
    ]
    for atom, element in enumerate(elements, start=1):
        value = atom * 1.0e-4 + force_offset
        lines.append(
            f" FORCES| {atom:6d} {value:19.12E} {-value:19.12E} "
            f"{0.5 * value:19.12E} {1.5 * value:19.12E}"
        )
    lines += [
        " FORCES| Sum     0.0 0.0 0.0",
        " FORCES| Total atomic force 0.0",
        "",
        " STRESS| Analytical stress tensor [GPa]",
        " STRESS|                         x                   y                   z",
        f" STRESS|      x {1.0 + stress_offset:19.11E} {0.1:19.11E} {0.2:19.11E}",
        f" STRESS|      y {0.1:19.11E} {2.0 + stress_offset:19.11E} {0.3:19.11E}",
        f" STRESS|      z {0.2:19.11E} {0.3:19.11E} {3.0 + stress_offset:19.11E}",
        " PROGRAM ENDED AT 2026-07-14 12:00:00",
        " PROGRAM STOPPED IN /tmp/normal-cp2k-run",
    ]
    return "\n".join(lines) + "\n"


class GXTBK222DerivativeGateTests(unittest.TestCase):
    def test_versioned_pair_is_exact_and_campaign_pinned(self) -> None:
        spec = derivative_gate.load_spec(derivative_gate.DEFAULT_SPEC)
        paths = derivative_gate.validate_gate_inputs(spec)
        self.assertEqual(set(paths), {"full", "spglib"})
        self.assertEqual(
            spec["campaign_identity"]["fingerprint_sha256"],
            "31a22ed28d2cdec9e9b071236a0fa4c1d30cb0d7cfb6a93f652f08758c83223e",
        )
        self.assertEqual(spec["tolerances"]["energy_max_abs_hartree"], 1.0e-9)
        self.assertEqual(spec["tolerances"]["force_max_abs_hartree_per_bohr"], 1.0e-6)
        self.assertEqual(spec["tolerances"]["stress_max_abs_gpa"], 1.0e-5)

    def test_parser_extracts_energy_forces_stress_and_kpoints(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "cp2k.out"
            output.write_text(fake_derivative_output())
            elements = ["N"] * 4 + ["H"] * 12
            parsed = derivative_gate.parse_cp2k_output(output, 16, elements)
            self.assertEqual(parsed["energy_hartree"], -10.0)
            self.assertEqual(parsed["kpoint_count"], 4)
            self.assertEqual(parsed["kpoint_count_source"], "cp2k_output")
            self.assertEqual(parsed["kpoint_mesh"], [2, 2, 2])
            self.assertEqual(parsed["kpoint_mesh_rows"], 8)
            self.assertEqual(len(parsed["forces"]), 16)
            self.assertEqual(parsed["forces"][0]["element"], "N")
            self.assertAlmostEqual(parsed["forces"][0]["vector_hartree_per_bohr"][0], 1.0e-4)
            self.assertEqual(len(parsed["stress_gpa"]), 3)
            self.assertAlmostEqual(parsed["stress_gpa"][2][2], 3.0)

    def test_full_grid_count_can_only_fall_back_to_hashed_input_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "cp2k.out"
            output.write_text(fake_derivative_output(include_kpoint_report=False))
            elements = ["N"] * 4 + ["H"] * 12
            with self.assertRaisesRegex(ValueError, "number of evaluated k-points missing"):
                derivative_gate.parse_cp2k_output(output, 16, elements)
            parsed = derivative_gate.parse_cp2k_output(
                output,
                16,
                elements,
                fallback_kpoint_count=8,
            )
            self.assertEqual(parsed["kpoint_count"], 8)
            self.assertEqual(parsed["kpoint_count_source"], "hashed_input_contract")

    def test_parser_rejects_nonconvergence_even_with_program_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "cp2k.out"
            output.write_text(fake_derivative_output() + "SCF run NOT converged\n")
            with self.assertRaisesRegex(ValueError, "fatal/nonconverged"):
                derivative_gate.parse_cp2k_output(output, 16)

    def test_validation_state_does_not_open_production_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(
                Path(tmp)
            )
            payload = json.loads(manifest.read_text())
            payload["campaign_state"] = "validation_in_progress"
            manifest.write_text(json.dumps(payload, indent=2) + "\n")
            with self.assertRaisesRegex(ValueError, "not production_ready"):
                common.validate_campaign_artifacts(
                    cp2k=cp2k,
                    cp2k_source=cp2k_source,
                    save_tblite=tblite,
                    save_tblite_source=save_source,
                    campaign_manifest=manifest,
                )
            identity, _, _, manifest_record = common.validate_campaign_artifacts(
                cp2k=cp2k,
                cp2k_source=cp2k_source,
                save_tblite=tblite,
                save_tblite_source=save_source,
                campaign_manifest=manifest,
                allowed_campaign_states=("validation_in_progress", "production_ready"),
            )
            self.assertEqual(identity["campaign_id"], "test-campaign")
            self.assertEqual(manifest_record["campaign_state"], "validation_in_progress")


class GXTBExperimentalK222PreflightTests(unittest.TestCase):
    def setUp(self) -> None:
        self.campaign = common.load_campaign_identity(experimental_preflight.ROOT)
        provenance = json.loads(
            (experimental_preflight.ROOT / "data" / common.GXTB_PROVENANCE_NAME).read_text()
        )
        self.cp2k = Path(provenance["cp2k"]["path"])

    def _write_completed_preflight(
        self,
        record: dict[str, object],
        *,
        force_offset: float = 0.0,
        stress_offset: float = 0.0,
    ) -> Path:
        input_path = Path(str(record["input"]))
        output = Path(str(record["output"]))
        elements = cellopt._coord_elements(input_path.read_text())
        output.write_text(
            fake_derivative_output(
                force_offset=force_offset,
                stress_offset=stress_offset,
                elements=elements,
            )
        )
        common.write_job_stamp(
            output.parent,
            input_path,
            self.cp2k,
            "GXTB",
            experimental_preflight.PHASE,
            "converged_measured_not_approved",
            details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
            campaign_identity=self.campaign,
            protocol_identity=experimental_preflight.PROTOCOL_IDENTITY,
            source_artifacts=experimental_preflight.source_artifacts(record),
        )
        return output

    def test_preflight_prepare_covers_all_23_and_preserves_input_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "preflight"
            path = experimental_preflight.prepare(output_root, self.campaign)
            payload = json.loads(path.read_text())
            self.assertEqual(payload["source_policy"], "experimental_reference")
            self.assertEqual(len(payload["systems"]), 23)
            records = experimental_preflight.load_manifest(output_root, self.campaign)
            self.assertEqual(set(records), {str(row["id"]) for row in cellopt.systems()})
            text = Path(str(records["ammonia"]["input"])).read_text()
            self.assertIn("RUN_TYPE ENERGY_FORCE", text)
            self.assertIn("MACDONALD 2 2 2 0.25 0.25 0.25", text)
            self.assertIn("SYMMETRY_BACKEND SPGLIB", text)
            self.assertIn("FULL_GRID F", text)
            self.assertIn("VERBOSE T", text)
            self.assertIn("&FORCES ON", text)
            self.assertIn("NDIGITS 12", text)
            self.assertIn("STRESS_UNIT GPa", text)
            self.assertNotIn("&MOTION", text)
            self.assertFalse(any(output_root.rglob("cp2k.out")))

    def test_large_finite_derivatives_are_measured_not_parse_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "preflight"
            experimental_preflight.prepare(output_root, self.campaign)
            record = experimental_preflight.load_manifest(output_root, self.campaign)["ammonia"]
            output = self._write_completed_preflight(
                record,
                force_offset=0.5,
                stress_offset=100.0,
            )
            summary = experimental_preflight.derivative_summary(output, record)
            self.assertGreater(summary["max_force_hartree_per_bohr"], 0.5)
            self.assertGreater(summary["max_abs_stress_GPa"], 100.0)
            trusted = experimental_preflight.validate_completed_case(
                output_root,
                "ammonia",
                self.campaign,
            )
            self.assertEqual(trusted["kpoint_count"], 4)

    def test_preflight_parser_rejects_incomplete_full_mesh_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "preflight"
            experimental_preflight.prepare(output_root, self.campaign)
            record = experimental_preflight.load_manifest(output_root, self.campaign)["ammonia"]
            output = Path(str(record["output"]))
            text = fake_derivative_output(
                elements=cellopt._coord_elements(Path(str(record["input"])).read_text())
            )
            mapping_row = (
                "          8   0.50000   0.50000   0.50000"
                "              8           4           1\n"
            )
            output.write_text(text.replace(mapping_row, ""))
            with self.assertRaisesRegex(ValueError, "full-mesh mapping"):
                experimental_preflight.derivative_summary(output, record)

    def test_preflight_collector_reports_all_derivatives_and_source_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "preflight"
            experimental_preflight.prepare(output_root, self.campaign)
            records = experimental_preflight.load_manifest(output_root, self.campaign)
            for record in records.values():
                self._write_completed_preflight(record)
            csv_path = base / "preflight.csv"
            with mock.patch.object(experimental_preflight, "_register_provenance_root"):
                experimental_preflight.collect_command(
                    SimpleNamespace(output_root=output_root, csv=csv_path)
                )
            with csv_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 23)
            ammonia = next(row for row in rows if row["system"] == "ammonia")
            self.assertEqual(ammonia["source_policy"], "experimental_reference")
            self.assertEqual(ammonia["scientific_status"], "measured_not_approved")
            self.assertEqual(ammonia["approved"], "False")
            self.assertEqual(float(ammonia["pressure_GPa"]), -2.0)
            self.assertEqual(float(ammonia["pressure_bar"]), -20000.0)
            self.assertGreater(float(ammonia["max_force_hartree_per_bohr"]), 0.0)
            self.assertGreater(float(ammonia["rms_force_hartree_per_bohr"]), 0.0)
            self.assertEqual(len(ammonia["source_input_sha256"]), 64)
            self.assertEqual(len(ammonia["structure_sha256"]), 64)

    def test_experimental_reference_prepare_uses_preflight_not_gamma(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            preflight_root = base / "preflight"
            experimental_preflight.prepare(preflight_root, self.campaign)
            record = experimental_preflight.load_manifest(preflight_root, self.campaign)["ammonia"]
            self._write_completed_preflight(record)
            output_root = base / "cellopt"
            args = SimpleNamespace(
                clean=False,
                method="GXTB",
                source_policy="experimental_reference",
                override=[],
                gamma_root=None,
                preflight_root=preflight_root,
                system=["ammonia"],
            )
            cellopt._prepare(args, output_root)
            with cellopt.manifest_path(output_root).open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            row = cellopt.validate_manifest_row(rows[0], output_root)
            self.assertEqual(row["source_policy"], "experimental_reference")
            self.assertEqual(row["variant"], cellopt.EXPERIMENTAL_VARIANT)
            self.assertEqual(row["source_restart"], "")
            self.assertEqual(row["preflight_output_sha256"], common.sha256_file(Path(row["preflight_output"])))
            text = Path(row["input"]).read_text()
            self.assertIn("RUN_TYPE CELL_OPT", text)
            self.assertIn("KEEP_ANGLES T", text)
            self.assertIn("MACDONALD 2 2 2 0.25 0.25 0.25", text)
            self.assertFalse((output_root / "GXTB" / "ammonia" / cellopt.VARIANT).exists())

    def test_cellopt_collector_reports_experimental_lineage_not_gamma(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            preflight_root = base / "preflight"
            experimental_preflight.prepare(preflight_root, self.campaign)
            preflight_record = experimental_preflight.load_manifest(
                preflight_root, self.campaign
            )["ammonia"]
            self._write_completed_preflight(preflight_record)
            output_root = base / "cellopt"
            cellopt._prepare(
                SimpleNamespace(
                    clean=False,
                    method="GXTB",
                    source_policy="experimental_reference",
                    override=[],
                    gamma_root=None,
                    preflight_root=preflight_root,
                    system=["ammonia"],
                ),
                output_root,
            )
            with cellopt.manifest_path(output_root).open(newline="") as handle:
                manifest_row = cellopt.validate_manifest_row(
                    next(csv.DictReader(handle)), output_root
                )
            run_dir = Path(manifest_row["run_dir"])
            input_path = Path(manifest_row["input"])
            output = run_dir / "cp2k.out"
            output.write_text(
                "ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -10.000000000000\n"
                "CELL| Volume [angstrom^3] 135.045176\n"
                "OPT| Step number 1\n"
                "OPT| Internal pressure [bar] 0.0\n"
                "OPT| Maximum step size 0.0\n"
                "OPT| RMS step size 0.0\n"
                "OPT| Maximum gradient 0.0\n"
                "OPT| RMS gradient 0.0\n"
                "GEOMETRY OPTIMIZATION COMPLETED\n"
                "PROGRAM ENDED\n"
            )
            protocol, sources = cellopt.stamp_context(manifest_row)
            common.write_job_stamp(
                run_dir,
                input_path,
                self.cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=self.campaign,
                protocol_identity=protocol,
                source_artifacts=sources,
            )
            result_csv = base / "cellopt.csv"
            ammonia_metadata = next(
                row for row in cellopt.systems() if str(row["id"]) == "ammonia"
            )
            with (
                mock.patch.object(cellopt, "systems", return_value=[ammonia_metadata]),
                mock.patch.object(cellopt, "update_provenance"),
                mock.patch.object(
                    cellopt,
                    "load_molecule_rows",
                    return_value={
                        ("GXTB", "ammonia"): {"gas_energy_hartree": "-1.0"}
                    },
                ),
            ):
                cellopt.collect(
                    SimpleNamespace(
                        method="GXTB",
                        output_root=output_root,
                        gamma_csv=None,
                        molecule_run_root=base,
                        csv=result_csv,
                        allow_incomplete=False,
                    )
                )
            with result_csv.open(newline="") as handle:
                result = next(csv.DictReader(handle))
            self.assertEqual(result["source_policy"], "experimental_reference")
            self.assertEqual(result["source"], "experimental_reference")
            self.assertEqual(result["source_variant"], "experimental_k222_preflight")
            self.assertEqual(result["source_restart"], "")
            self.assertEqual(result["variant"], cellopt.EXPERIMENTAL_VARIANT)
            self.assertEqual(result["preflight_output_sha256"], common.sha256_file(Path(result["preflight_output"])))

    def test_cellopt_collector_rejects_unstamped_gxtb_gamma_csv(self) -> None:
        with self.assertRaisesRegex(ValueError, "--gamma-csv is not accepted for GXTB"):
            cellopt.collect(
                SimpleNamespace(
                    method="GXTB",
                    output_root=Path("unused"),
                    gamma_csv=Path("untrusted.csv"),
                    molecule_run_root=None,
                    csv=Path("unused.csv"),
                    allow_incomplete=False,
                )
            )

    def test_final_sp_resolves_experimental_cellopt_variant_from_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            preflight_root = base / "preflight"
            experimental_preflight.prepare(preflight_root, self.campaign)
            preflight_record = experimental_preflight.load_manifest(
                preflight_root, self.campaign
            )["ammonia"]
            self._write_completed_preflight(preflight_record)
            cellopt_root = base / "cellopt"
            cellopt._prepare(
                SimpleNamespace(
                    clean=False,
                    method="GXTB",
                    source_policy="experimental_reference",
                    override=[],
                    gamma_root=None,
                    preflight_root=preflight_root,
                    system=["ammonia"],
                ),
                cellopt_root,
            )
            with cellopt.manifest_path(cellopt_root).open(newline="") as handle:
                source_row = cellopt.validate_manifest_row(
                    next(csv.DictReader(handle)), cellopt_root
                )
            run_dir = Path(source_row["run_dir"])
            input_path = Path(source_row["input"])
            output = run_dir / "cp2k.out"
            output.write_text(
                "ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -10.0\n"
                "GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n"
            )
            protocol, sources = cellopt.stamp_context(source_row)
            common.write_job_stamp(
                run_dir,
                input_path,
                self.cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=self.campaign,
                protocol_identity=protocol,
                source_artifacts=sources,
            )
            restart = run_dir / f"{input_path.stem}-1.restart"
            restart.write_text(input_path.read_text())
            ammonia_metadata = next(
                row for row in cellopt.systems() if str(row["id"]) == "ammonia"
            )
            final_root = base / "k333"
            with mock.patch.object(cellopt, "systems", return_value=[ammonia_metadata]):
                final_sp._prepare(
                    SimpleNamespace(
                        clean=False,
                        method="GXTB",
                        system=["ammonia"],
                        mesh=3,
                        cellopt_root=cellopt_root,
                    ),
                    final_root,
                )
            with final_sp.manifest_path(final_root).open(newline="") as handle:
                final_row = next(csv.DictReader(handle))
            self.assertEqual(final_row["source_variant"], cellopt.EXPERIMENTAL_VARIANT)
            self.assertEqual(final_row["source_policy"], "experimental_reference")
            self.assertEqual(Path(final_row["source_input"]).resolve(), input_path.resolve())
            self.assertEqual(Path(final_row["source_run_dir"]).resolve(), run_dir.resolve())
            target_protocol, target_sources = final_sp.target_stamp_context(final_row)
            self.assertEqual(target_protocol["source_policy"], "experimental_reference")
            self.assertEqual(target_protocol["target_mesh"], "k333")
            self.assertIn("preflight_stamp", target_sources)
            self.assertEqual(target_sources["cellopt_restart"].resolve(), restart.resolve())
            self.assertIn("RUN_TYPE ENERGY", Path(final_row["input"]).read_text())
            self.assertIn("MACDONALD 3 3 3 0 0 0", Path(final_row["input"]).read_text())
            Path(final_row["input"]).write_text(Path(final_row["input"]).read_text() + "# drift\n")
            with self.assertRaisesRegex(ValueError, "target input fingerprint differs"):
                final_sp.target_stamp_context(final_row)

    def test_final_sp_collector_rejects_mesh_mismatch_and_duplicates_early(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "final"
            output_root.mkdir()
            row = {
                "method": "GXTB",
                "system": "ammonia",
                "mesh": "4",
                "input": str(base / "input.inp"),
                "run_dir": str(base / "run"),
            }
            with final_sp.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=final_sp.MANIFEST_FIELDS)
                writer.writeheader()
                writer.writerow(row)
            args = SimpleNamespace(
                output_root=output_root,
                method="GXTB",
                mesh=3,
                molecule_run_root=base,
                csv=base / "result.csv",
                allow_incomplete=True,
            )
            with self.assertRaisesRegex(ValueError, "collect --mesh 3 differs"):
                final_sp.collect(args)

            row["mesh"] = "3"
            with final_sp.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=final_sp.MANIFEST_FIELDS)
                writer.writeheader()
                writer.writerow(row)
                writer.writerow(row)
            with self.assertRaisesRegex(ValueError, "duplicate method/system"):
                final_sp.collect(args)

    def test_experimental_reference_is_gxtb_only_and_gxtb_override_is_forbidden(self) -> None:
        base = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(base, ignore_errors=True))
        common_args = {
            "clean": False,
            "gamma_root": None,
            "preflight_root": base / "preflight",
            "system": ["ammonia"],
        }
        with self.assertRaisesRegex(ValueError, "restricted to --method GXTB"):
            cellopt._prepare(
                SimpleNamespace(
                    **common_args,
                    method="GFN2",
                    source_policy="experimental_reference",
                    override=[],
                ),
                base / "cellopt-gfn2",
            )
        with self.assertRaisesRegex(ValueError, "--override is forbidden for GXTB"):
            cellopt._prepare(
                SimpleNamespace(
                    **common_args,
                    method="GXTB",
                    source_policy="experimental_reference",
                    override=["GXTB/ammonia=/tmp/foreign.restart"],
                ),
                base / "cellopt-gxtb",
            )

    def test_source_artifact_mutation_invalidates_completed_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cp2k = fake_executable(base)
            campaign = fake_campaign(cp2k)
            input_path = base / "input.inp"
            output = base / "cp2k.out"
            source = base / "source.inp"
            input_path.write_text("input")
            output.write_text("PROGRAM ENDED")
            source.write_text("source-v1")
            protocol = {"source_policy": "experimental_reference"}
            common.write_job_stamp(
                base,
                input_path,
                cp2k,
                "GXTB",
                experimental_preflight.PHASE,
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=campaign,
                protocol_identity=protocol,
                source_artifacts={"reference_input": source},
            )
            valid, _ = common.recorded_job_stamp_matches(
                base,
                input_path,
                "GXTB",
                experimental_preflight.PHASE,
                output,
                campaign_identity=campaign,
                protocol_identity=protocol,
                source_artifacts={"reference_input": source},
            )
            self.assertTrue(valid)
            source.write_text("source-v2")
            valid, reason = common.recorded_job_stamp_matches(
                base,
                input_path,
                "GXTB",
                experimental_preflight.PHASE,
                output,
                campaign_identity=campaign,
                protocol_identity=protocol,
                source_artifacts={"reference_input": source},
            )
            self.assertFalse(valid)
            self.assertIn("source artifact fingerprint differs", reason)


class GXTBInputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.system = pipeline.SYSTEMS[0]
        self.geometry = {
            "cell": [[8.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 10.0]],
            "atoms": [{"element": "C", "frac": [0.1, 0.2, 0.3]}],
        }

    def test_gxtb_uses_native_mixer_and_spglib_reduction(self) -> None:
        text = pipeline.molecule_input(
            self.system,
            [{"element": "C", "coord": [0.0, 0.0, 0.0]}],
            "GXTB",
        )
        self.assertIn("METHOD GXTB", text)
        self.assertIn("SCC_MIXER TBLITE", text)
        self.assertNotIn("SCC_MIXER CP2K", text)
        self.assertNotIn("&TBLITE_MIXER", text)
        periodic = pipeline.crystal_input(
            self.system,
            self.geometry,
            "GXTB",
            pipeline.MESHES[2],
            "ENERGY",
        )
        self.assertIn("FULL_GRID F", periodic)
        self.assertIn("SYMMETRY T", periodic)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", periodic)
        self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", periodic)
        legacy = pipeline.crystal_input(
            self.system,
            self.geometry,
            "GFN2",
            pipeline.MESHES[2],
            "ENERGY",
        )
        self.assertIn("FULL_GRID F", legacy)
        self.assertIn("SYMMETRY T", legacy)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", legacy)
        common.validate_method_input(text, "GXTB")

    def test_method_mismatch_is_rejected(self) -> None:
        text = pipeline.molecule_input(
            self.system,
            [{"element": "C", "coord": [0.0, 0.0, 0.0]}],
            "GXTB",
        )
        with self.assertRaisesRegex(ValueError, "METHOD GFN2"):
            common.validate_method_input(text, "GFN2")

    def test_common_gate_rejects_legacy_gxtb_full_grid_input(self) -> None:
        text = pipeline.crystal_input(
            self.system,
            self.geometry,
            "GXTB",
            pipeline.MESHES[2],
            "ENERGY",
        )
        legacy = text.replace("SYMMETRY T", "SYMMETRY F").replace("FULL_GRID F", "FULL_GRID T")
        legacy = legacy.replace("      SYMMETRY_BACKEND SPGLIB\n", "")
        legacy = legacy.replace("      SYMMETRY_REDUCTION_METHOD SPGLIB\n", "")
        with self.assertRaisesRegex(ValueError, "SPGLIB-reduced mesh contract"):
            common.validate_method_input(legacy, "GXTB")

    def test_gamma_to_k222_to_k444_keeps_gxtb_provenance(self) -> None:
        gamma = pipeline.crystal_input(
            self.system,
            self.geometry,
            "GXTB",
            pipeline.MESHES[0],
            "CELL_OPT",
        )
        with tempfile.TemporaryDirectory() as tmp:
            gamma_restart = Path(tmp) / "gamma-1.restart"
            gamma_restart.write_text(gamma)
            k222 = cellopt.restart_to_k222_input(gamma_restart, "gxtb_k222", "GXTB")
            self.assertIn("MACDONALD 2 2 2 0.25 0.25 0.25", k222)
            self.assertIn("FULL_GRID F", k222)
            self.assertIn("SYMMETRY T", k222)
            self.assertIn("SYMMETRY_BACKEND SPGLIB", k222)
            self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", k222)
            self.assertIn("KEEP_ANGLES T", k222)
            self.assertIn("SCC_MIXER TBLITE", k222)
            k222_restart = Path(tmp) / "k222-1.restart"
            k222_restart.write_text(k222)
            k444 = final_sp.restart_to_single_point(k222_restart, "gxtb_k444", 4, "GXTB")
            self.assertIn("MACDONALD 4 4 4 0.375 0.375 0.375", k444)
            self.assertIn("RUN_TYPE ENERGY", k444)
            self.assertNotIn("&MOTION", k444)
            common.validate_method_input(k444, "GXTB")

            stale_restart = Path(tmp) / "stale-k222.restart"
            stale_restart.write_text(k222.replace("FULL_GRID F", "FULL_GRID T"))
            with self.assertRaisesRegex(ValueError, "SPGLIB-reduced production mesh"):
                final_sp.restart_to_single_point(stale_restart, "stale_gxtb_k444", 4, "GXTB")

    def test_fixed_reference_inputs_cover_all_si_meshes(self) -> None:
        self.assertNotIn(pipeline.FIXED_REFERENCE_PHASE, pipeline.PRODUCTION_PHASES)
        self.assertIn(pipeline.FIXED_REFERENCE_PHASE, pipeline.RUN_PHASES)
        for mesh in pipeline.MESHES:
            text = pipeline.crystal_input(self.system, self.geometry, "GXTB", mesh, "ENERGY")
            common.validate_method_input(text, "GXTB")
            if mesh["scheme"] is None:
                self.assertNotIn("&KPOINTS", text)
            else:
                self.assertIn(str(mesh["scheme"]), text)
                self.assertIn("SYMMETRY T", text)
                self.assertIn("FULL_GRID F", text)
                self.assertIn("SYMMETRY_BACKEND SPGLIB", text)
                self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", text)
        self.assertTrue(
            pipeline.phase_completed(
                pipeline.FIXED_REFERENCE_PHASE,
                self._completed_sp_output(),
            )
        )

    def _completed_sp_output(self) -> Path:
        temporary = tempfile.NamedTemporaryFile("w", delete=False)
        temporary.write("ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.0\nPROGRAM ENDED\n")
        temporary.close()
        self.addCleanup(Path(temporary.name).unlink, missing_ok=True)
        return Path(temporary.name)


class GXTBPruneTests(unittest.TestCase):
    def test_pruning_is_gxtb_only_and_retains_continuation_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            run_dir = base / "GXTB" / "case"
            run_dir.mkdir(parents=True)
            final_restart = run_dir / "case-1.restart"
            final_restart.write_text("final")
            (run_dir / "case-1_4.restart").write_text("old")
            (run_dir / "case-RESTART.kp").write_text("large")
            (run_dir / "case-RESTART.wfn").write_text("current")
            (run_dir / "case-RESTART.wfn.bak-1").write_text("backup")
            (run_dir / "cp2k.out").write_text("PROGRAM ENDED")

            report = common.prune_gxtb_transients(run_dir)
            self.assertTrue(final_restart.exists())
            self.assertTrue((run_dir / "case-RESTART.wfn").exists())
            self.assertFalse((run_dir / "case-1_4.restart").exists())
            self.assertFalse((run_dir / "case-RESTART.kp").exists())
            self.assertFalse((run_dir / "case-RESTART.wfn.bak-1").exists())
            self.assertGreater(report["bytes_deleted"], 0)

            published = base / "GFN2" / "case"
            published.mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "outside a GXTB directory"):
                common.prune_gxtb_transients(published)


class GXTBCoverageTests(unittest.TestCase):
    def test_cellopt_loader_accepts_exactly_23_gxtb_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cellopt.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("method", "system", "mesh", "opt_completed"),
                    lineterminator="\n",
                )
                writer.writeheader()
                for system in pipeline.SYSTEMS:
                    writer.writerow(
                        {
                            "method": "GXTB",
                            "system": system["id"],
                            "mesh": "k222",
                            "opt_completed": "True",
                        }
                    )
            rows = pipeline.load_cellopt_rows(path, ["GXTB"])
            self.assertEqual(len(rows), 23)

            lines = path.read_text().splitlines()
            path.write_text("\n".join(lines[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "not a complete X23b set"):
                pipeline.load_cellopt_rows(path, ["GXTB"])


class GXTBLockAndResumeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.system = pipeline.SYSTEMS[0]
        self.geometry = {
            "cell": [[8.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 10.0]],
            "atoms": [{"element": "C", "frac": [0.1, 0.2, 0.3]}],
        }

    def test_busy_lock_is_non_success_in_all_three_runners(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cp2k = fake_executable(base)
            molecule = base / "molecule.inp"
            molecule.write_text(
                pipeline.molecule_input(
                    self.system,
                    [{"element": "C", "coord": [0.0, 0.0, 0.0]}],
                    "GXTB",
                )
            )
            old_root = pipeline.ROOT
            pipeline.ROOT = base
            try:
                with molecule.open() as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    _, code, action = pipeline.run_phase_one(
                        molecule, "molecule_geoopt", "GXTB", cp2k, 1, False, False
                    )
                self.assertEqual((code, action), (common.BUSY_RETURN_CODE, "BUSY"))
            finally:
                pipeline.ROOT = old_root

            gamma = base / "gamma.restart"
            gamma.write_text(
                pipeline.crystal_input(self.system, self.geometry, "GXTB", pipeline.MESHES[0], "CELL_OPT")
            )
            primary = base / "GXTB" / "system" / cellopt.VARIANT / "primary.inp"
            primary.parent.mkdir(parents=True)
            primary.write_text(cellopt.restart_to_k222_input(gamma, "primary", "GXTB"))
            with primary.open() as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                _, code, action = cellopt.run_one(primary, cp2k, "GXTB", 1, False, False)
            self.assertEqual((code, action), (common.BUSY_RETURN_CODE, "BUSY"))

            final_input = primary.parent / "final.inp"
            final_input.write_text(final_sp.restart_to_single_point(primary, "final", 3, "GXTB"))
            with final_input.open() as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                _, code, action = final_sp.run_one(
                    final_input, cp2k, 1, 3, "GXTB", False, False
                )
            self.assertEqual((code, action), (common.BUSY_RETURN_CODE, "BUSY"))

    def test_gxtb_skip_requires_matching_input_and_cp2k_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "GXTB" / "system" / final_sp.variant(3)
            run_dir.mkdir(parents=True)
            cp2k = fake_executable(Path(tmp))
            campaign = fake_campaign(cp2k)
            gamma = Path(tmp) / "gamma.restart"
            gamma.write_text(
                pipeline.crystal_input(self.system, self.geometry, "GXTB", pipeline.MESHES[0], "CELL_OPT")
            )
            k222 = Path(tmp) / "k222.restart"
            k222.write_text(cellopt.restart_to_k222_input(gamma, "k222", "GXTB"))
            input_path = run_dir / "primary.inp"
            input_path.write_text(final_sp.restart_to_single_point(k222, "final", 3, "GXTB"))
            output = run_dir / "cp2k.out"
            output.write_text("ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.0\nPROGRAM ENDED\n")

            _, code, action = final_sp.run_one(
                input_path, cp2k, 1, 3, "GXTB", False, False, campaign
            )
            self.assertEqual((code, action), (1, "STALE_OUTPUT"))

            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                "GXTB",
                "x23b_final_k333_on_k222",
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=campaign,
            )
            _, code, action = final_sp.run_one(
                input_path, cp2k, 1, 3, "GXTB", False, False, campaign
            )
            self.assertEqual((code, action), (0, "SKIP"))

            original_output = output.read_text()
            output.write_text(original_output + "# foreign output\n")
            _, code, action = final_sp.run_one(
                input_path, cp2k, 1, 3, "GXTB", False, False, campaign
            )
            self.assertEqual((code, action), (1, "STALE_OUTPUT"))
            output.write_text(original_output)

            other_cp2k = fake_executable(Path(tmp), "cp2k-other")
            other_cp2k.write_text("#!/bin/sh\necho 'different CP2K'\n")
            _, code, action = final_sp.run_one(
                input_path, other_cp2k, 1, 3, "GXTB", False, False, campaign
            )
            self.assertEqual((code, action), (1, "STALE_OUTPUT"))

            input_path.write_text(input_path.read_text() + "\n# changed\n")
            _, code, action = final_sp.run_one(
                input_path, cp2k, 1, 3, "GXTB", False, False, campaign
            )
            self.assertEqual((code, action), (1, "STALE_OUTPUT"))

    def test_continuation_promotes_a_path_bound_collectable_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            run_dir = base / "GXTB" / "system" / cellopt.EXPERIMENTAL_VARIANT
            run_dir.mkdir(parents=True)
            cp2k = fake_executable(base)
            campaign = fake_campaign(cp2k)
            gamma = base / "gamma.restart"
            gamma.write_text(
                pipeline.crystal_input(
                    self.system, self.geometry, "GXTB", pipeline.MESHES[0], "CELL_OPT"
                )
            )
            input_path = run_dir / "primary.inp"
            primary_text = cellopt.restart_to_k222_input(gamma, "primary", "GXTB")
            input_path.write_text(primary_text)
            restart = run_dir / "primary-1_500.restart"
            restart.write_text(
                primary_text.replace("    MAX_ITER 500\n", "    MAX_ITER 500\n    STEP_START_VAL 500\n")
            )
            output = run_dir / "cp2k.out"
            output.write_text("MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED\nPROGRAM ENDED\n")
            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "max_iter",
                details={
                    "output": str(output.resolve()),
                    "output_sha256": common.sha256_file(output),
                },
                campaign_identity=campaign,
            )

            def complete_continuation(command: list[str], *, cwd: Path, **_: object) -> object:
                continued_output = Path(cwd) / command[command.index("-o") + 1]
                continued_output.write_text(
                    "GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n"
                )
                return SimpleNamespace(returncode=0)

            with mock.patch.object(cellopt.subprocess, "run", side_effect=complete_continuation):
                _, code, action = cellopt.continue_one(
                    input_path,
                    cp2k,
                    100,
                    1,
                    "GXTB",
                    1,
                    False,
                    campaign,
                )
            self.assertEqual((code, action), (0, "CONTINUE_1"))
            valid, reason = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                "GXTB",
                "x23b_k222_cellopt",
                output,
                campaign_identity=campaign,
            )
            self.assertTrue(valid, reason)

    def test_bfgs_promotion_promotes_a_path_bound_collectable_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            run_dir = base / "GXTB" / "system" / cellopt.EXPERIMENTAL_VARIANT
            run_dir.mkdir(parents=True)
            cp2k = fake_executable(base)
            campaign = fake_campaign(cp2k)
            gamma = base / "gamma.restart"
            gamma.write_text(
                pipeline.crystal_input(
                    self.system, self.geometry, "GXTB", pipeline.MESHES[0], "CELL_OPT"
                )
            )
            input_path = run_dir / "primary.inp"
            primary_text = cellopt.restart_to_k222_input(gamma, "primary", "GXTB")
            input_path.write_text(primary_text)
            restart = run_dir / "primary-1_500.restart"
            restart.write_text(primary_text)
            output = run_dir / "cp2k.out"
            output.write_text(
                "OPT| Step number 500\n"
                "OPT| Internal pressure [bar] 50.0\n"
                "OPT| Maximum step size 0.002\n"
                "OPT| RMS step size 0.001\n"
                "OPT| Maximum gradient 0.0004\n"
                "OPT| RMS gradient 0.0002\n"
                "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED\n"
                "PROGRAM ENDED\n"
            )
            common.write_job_stamp(
                run_dir,
                input_path,
                cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "max_iter",
                details={
                    "output": str(output.resolve()),
                    "output_sha256": common.sha256_file(output),
                },
                campaign_identity=campaign,
            )

            def complete_polish(command: list[str], *, cwd: Path, **_: object) -> object:
                polish_dir = Path(cwd)
                polish_output = polish_dir / command[command.index("-o") + 1]
                polish_output.write_text(
                    "GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n"
                )
                (polish_dir / "polished-1.restart").write_text(primary_text)
                return SimpleNamespace(returncode=0)

            with mock.patch.object(cellopt.subprocess, "run", side_effect=complete_polish):
                _, code, action = cellopt.polish_one(
                    input_path,
                    cp2k,
                    100,
                    0.05,
                    False,
                    "GXTB",
                    1,
                    False,
                    campaign,
                )
            self.assertEqual((code, action), (0, "POLISHED_FROM_500"))
            valid, reason = common.recorded_job_stamp_matches(
                run_dir,
                input_path,
                "GXTB",
                "x23b_k222_cellopt",
                output,
                campaign_identity=campaign,
            )
            self.assertTrue(valid, reason)


class GXTBManifestTests(unittest.TestCase):
    def test_k222_selects_manifest_primary_and_ignores_continuations(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "k222"
            run_dir = output_root / "GXTB" / "ammonia" / cellopt.VARIANT
            run_dir.mkdir(parents=True)
            primary = run_dir / "ammonia_GXTB_k222_cellopt_keep_angles.inp"
            continuation = run_dir / "ammonia_GXTB_k222_cellopt_keep_angles_continue_500.inp"
            polish = run_dir / "bfgs_polish_from_500" / "polish.inp"
            primary.write_text("primary")
            continuation.write_text("continuation")
            polish.parent.mkdir()
            polish.write_text("polish")
            gamma = base / "gamma.restart"
            gamma.write_text("gamma")
            with cellopt.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=cellopt.MANIFEST_FIELDS,
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "schema": cellopt.MANIFEST_SCHEMA,
                        "method": "GXTB",
                        "system": "ammonia",
                        "variant": cellopt.VARIANT,
                        "source_policy": "gamma_cellopt_restart",
                        "source_kind": "gamma_cellopt_restart",
                        "lineage": "gamma_cellopt_restart->k222_cellopt_input",
                        "source_path": str(gamma),
                        "source_sha256": common.sha256_file(gamma),
                        "source_restart": str(gamma),
                        "input": str(primary),
                        "input_sha256": common.sha256_file(primary),
                        "run_dir": str(run_dir),
                    }
                )
            selected = cellopt.manifest_owned_inputs(output_root, "GXTB", {"ammonia"})
            self.assertEqual(selected, [primary.resolve()])

    def test_policyless_gxtb_manifest_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "k222"
            run_dir = output_root / "GXTB" / "ammonia" / cellopt.VARIANT
            run_dir.mkdir(parents=True)
            primary = run_dir / "primary.inp"
            primary.write_text("primary")
            with cellopt.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("method", "system", "source_restart", "input", "run_dir"),
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "method": "GXTB",
                        "system": "ammonia",
                        "source_restart": str(Path(tmp) / "gamma.restart"),
                        "input": str(primary),
                        "run_dir": str(run_dir),
                    }
                )
            with self.assertRaisesRegex(ValueError, "policyless legacy GXTB"):
                cellopt.manifest_owned_inputs(output_root, "GXTB", {"ammonia"})

    def test_schema2_gamma_manifest_requires_exact_kind_lineage_and_no_structure_label(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "k222"
            run_dir = output_root / "GXTB" / "ammonia" / cellopt.VARIANT
            run_dir.mkdir(parents=True)
            primary = run_dir / "primary.inp"
            gamma = base / "gamma.restart"
            primary.write_text("primary\n")
            gamma.write_text("gamma\n")
            row = {
                "schema": str(cellopt.MANIFEST_SCHEMA),
                "method": "GXTB",
                "system": "ammonia",
                "variant": cellopt.VARIANT,
                "source_policy": "gamma_cellopt_restart",
                "source_kind": "gamma_cellopt_restart",
                "lineage": cellopt.LINEAGE_BY_POLICY["gamma_cellopt_restart"],
                "source_path": str(gamma),
                "source_sha256": common.sha256_file(gamma),
                "source_restart": str(gamma),
                "input": str(primary),
                "input_sha256": common.sha256_file(primary),
                "run_dir": str(run_dir),
            }
            self.assertEqual(
                cellopt.validate_manifest_row(row, output_root)["source_kind"],
                "gamma_cellopt_restart",
            )
            for field, value, message in (
                ("source_kind", "experimental_reference", "source kind/policy mismatch"),
                ("lineage", "invented", "noncanonical source lineage"),
                ("structure_source", "invented", "structure lineage label"),
            ):
                changed = dict(row)
                changed[field] = value
                with self.assertRaisesRegex(ValueError, message):
                    cellopt.validate_manifest_row(changed, output_root)

    def test_final_sp_root_is_single_mesh_and_cannot_be_rebound(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = (base / "final").resolve()
            output_root.mkdir()
            with final_sp.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=final_sp.MANIFEST_FIELDS)
                writer.writeheader()
                writer.writerow({"method": "GXTB", "system": "ammonia", "mesh": "3"})
            with self.assertRaisesRegex(ValueError, "belongs to meshes"):
                final_sp.guard_output_root_mesh(output_root, 4)

            final_sp.manifest_path(output_root).unlink()
            benchmark = base / "X23b"
            (benchmark / "data").mkdir(parents=True)
            (benchmark / "data" / common.GXTB_PROVENANCE_NAME).write_text(
                json.dumps(
                    {
                        "workflow_paths": {
                            "final_single_point_roots": {"k333": str(output_root)}
                        }
                    }
                )
                + "\n"
            )
            with mock.patch.object(cellopt, "ROOT", benchmark):
                with self.assertRaisesRegex(ValueError, "already frozen for k333"):
                    final_sp.guard_output_root_mesh(output_root, 4)

    def test_source_policy_cannot_be_switched_in_place(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "k222"
            run_dir = output_root / "GXTB" / "ammonia" / cellopt.VARIANT
            run_dir.mkdir(parents=True)
            primary = run_dir / "primary.inp"
            gamma = base / "gamma.restart"
            primary.write_text("primary")
            gamma.write_text("gamma")
            with cellopt.manifest_path(output_root).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=cellopt.MANIFEST_FIELDS)
                writer.writeheader()
                writer.writerow(
                    {
                        "schema": cellopt.MANIFEST_SCHEMA,
                        "method": "GXTB",
                        "system": "ammonia",
                        "variant": cellopt.VARIANT,
                        "source_policy": "gamma_cellopt_restart",
                        "source_kind": "gamma_cellopt_restart",
                        "lineage": "gamma_cellopt_restart->k222_cellopt_input",
                        "source_path": str(gamma),
                        "source_sha256": common.sha256_file(gamma),
                        "source_restart": str(gamma),
                        "input": str(primary),
                        "input_sha256": common.sha256_file(primary),
                        "run_dir": str(run_dir),
                    }
                )
            with self.assertRaisesRegex(ValueError, "different GXTB source policy"):
                cellopt._prepare(
                    SimpleNamespace(
                        clean=False,
                        method="GXTB",
                        source_policy="experimental_reference",
                        override=[],
                        gamma_root=None,
                        preflight_root=base / "preflight",
                        system=["ammonia"],
                    ),
                    output_root,
                )

    def test_duplicate_override_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "duplicate override"):
            cellopt.parse_overrides(
                ["GFN2/ammonia=/tmp/a.restart", "GFN2/ammonia=/tmp/b.restart"]
            )


class GXTBProvenanceTests(unittest.TestCase):
    def test_gxtb_build_artifacts_are_mandatory(self) -> None:
        with self.assertRaisesRegex(ValueError, "--cp2k-source"):
            common.require_gxtb_build_artifacts(
                cp2k=Path("cp2k"),
                cp2k_source=None,
                save_tblite=None,
                save_tblite_source=None,
                campaign_manifest=None,
            )

    def test_frozen_campaign_manifest_file_is_immutable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "X23b"
            (root / "data").mkdir(parents=True)
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(base)
            common.update_gxtb_provenance(
                root,
                cp2k=cp2k,
                cp2k_source=cp2k_source,
                save_tblite=tblite,
                save_tblite_source=save_source,
                campaign_manifest=manifest,
            )
            changed = json.loads(manifest.read_text())
            changed["validation_note"] = "same build identity, changed manifest bytes"
            manifest.write_text(json.dumps(changed, indent=2) + "\n")
            with self.assertRaisesRegex(ValueError, "manifest fingerprint differs"):
                common.load_campaign_identity(root)
            with self.assertRaisesRegex(ValueError, "manifest fingerprint differs"):
                common.update_gxtb_provenance(
                    root,
                    cp2k=cp2k,
                    cp2k_source=cp2k_source,
                    save_tblite=tblite,
                    save_tblite_source=save_source,
                    campaign_manifest=manifest,
                )

    def test_recorded_stamp_cannot_be_relocated_to_another_run_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cp2k = fake_executable(base)
            campaign = fake_campaign(cp2k)
            first = base / "first"
            second = base / "second"
            first.mkdir()
            second.mkdir()
            first_input = first / "job.inp"
            first_output = first / "cp2k.out"
            first_input.write_text("METHOD GXTB\n")
            first_output.write_text("PROGRAM ENDED\n")
            common.write_job_stamp(
                first,
                first_input,
                cp2k,
                "GXTB",
                "relocation_test",
                "converged",
                details={
                    "output": str(first_output.resolve()),
                    "output_sha256": common.sha256_file(first_output),
                },
                campaign_identity=campaign,
            )
            second_input = second / "job.inp"
            second_output = second / "cp2k.out"
            shutil.copyfile(first_input, second_input)
            shutil.copyfile(first_output, second_output)
            shutil.copyfile(first / common.JOB_STAMP_NAME, second / common.JOB_STAMP_NAME)
            valid, reason = common.recorded_job_stamp_matches(
                second,
                second_input,
                "GXTB",
                "relocation_test",
                second_output,
                campaign_identity=campaign,
            )
            self.assertFalse(valid)
            self.assertIn("input fingerprint differs", reason)

    def test_provenance_is_gxtb_only_hashed_and_stamp_counted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "X23b"
            (root / "data").mkdir(parents=True)
            (root / "data" / "metadata.json").write_text(
                json.dumps({"systems": [{"id": "ammonia"}]})
            )
            legacy = root / "data" / "build_provenance.json"
            legacy.write_text('{"frozen": true}\n')
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(
                Path(tmp)
            )

            cellopt_root = root / "runs" / "native" / "k222"
            k333_root = root / "runs" / "native" / "k333"
            common.update_gxtb_provenance(
                root,
                cp2k=cp2k,
                cp2k_source=cp2k_source,
                save_tblite=tblite,
                save_tblite_source=save_source,
                campaign_manifest=manifest,
                workflow_paths={
                    "k222_cellopt_root": cellopt_root,
                    "final_single_point_roots": {"k333": k333_root},
                },
            )
            campaign = common.load_campaign_identity(root)

            stem = "ammonia_GXTB_mol_geoopt"
            input_path = root / "inputs" / "molecule_geoopt" / "GXTB" / f"{stem}.inp"
            output = root / "runs" / "molecule_geoopt" / "GXTB" / stem / f"{stem}.out"
            input_path.parent.mkdir(parents=True)
            output.parent.mkdir(parents=True)
            input_path.write_text("METHOD GXTB\n")
            output.write_text("GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n")
            common.write_job_stamp(
                output.parent,
                input_path,
                cp2k,
                "GXTB",
                "x23b_molecule_geoopt",
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=campaign,
            )

            cellopt_run = cellopt_root / "GXTB" / "ammonia" / cellopt.VARIANT
            cellopt_run.mkdir(parents=True)
            cellopt_input = cellopt_run / "ammonia_GXTB_k222_cellopt_keep_angles.inp"
            cellopt_output = cellopt_run / "cp2k.out"
            cellopt_input.write_text("METHOD GXTB\n")
            cellopt_output.write_text("GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n")
            common.write_job_stamp(
                cellopt_run,
                cellopt_input,
                cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "converged",
                details={
                    "output": str(cellopt_output.resolve()),
                    "output_sha256": common.sha256_file(cellopt_output),
                },
                campaign_identity=campaign,
            )

            final_run = k333_root / "GXTB" / "ammonia" / final_sp.variant(3)
            final_run.mkdir(parents=True)
            final_input = final_run / "ammonia_GXTB_k333_sp_on_k222.inp"
            final_output = final_run / "cp2k.out"
            final_input.write_text("METHOD GXTB\n")
            final_output.write_text("ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.0\nPROGRAM ENDED\n")
            common.write_job_stamp(
                final_run,
                final_input,
                cp2k,
                "GXTB",
                "x23b_final_k333_on_k222",
                "converged",
                details={
                    "output": str(final_output.resolve()),
                    "output_sha256": common.sha256_file(final_output),
                },
                campaign_identity=campaign,
            )

            path = common.update_gxtb_provenance(root)
            payload = json.loads(path.read_text())
            self.assertEqual(path.name, common.GXTB_PROVENANCE_NAME)
            self.assertEqual(payload["method"], "GXTB")
            self.assertEqual(payload["validation"]["gas_optimizations_completed"], 1)
            self.assertEqual(payload["validation"]["k222_cell_optimizations_completed"], 1)
            self.assertEqual(payload["validation"]["k333_single_points_completed"], 1)
            self.assertEqual(payload["cp2k"]["executable_sha256"], common.sha256_file(cp2k))
            self.assertTrue(payload["cp2k"]["source_revision"])
            self.assertEqual(payload["campaign_identity"], campaign)
            self.assertEqual(payload["campaign_manifest"]["path"], str(manifest.resolve()))
            self.assertIsInstance(payload["workflow_paths"]["final_single_point_roots"]["k333"], str)
            self.assertIn("SPGLIB-reduced", payload["protocol"]["kpoint_mesh_contract"])
            self.assertIn("diagnostics only", payload["protocol"]["legacy_full_grid_policy"])
            self.assertEqual(legacy.read_text(), '{"frozen": true}\n')

    def test_experimental_completion_uses_preflight_gate_not_gamma(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "X23b"
            (root / "data").mkdir(parents=True)
            (root / "data" / "metadata.json").write_text(
                json.dumps({"systems": [{"id": "case"}]})
            )
            cp2k, cp2k_source, tblite, save_source, campaign_manifest = (
                fake_campaign_artifacts(base)
            )
            preflight_root = root / "runs" / "experimental_preflight"
            cellopt_root = root / "runs" / "k222"
            k333_root = root / "runs" / "k333"
            k444_root = root / "runs" / "k444"
            common.update_gxtb_provenance(
                root,
                cp2k=cp2k,
                cp2k_source=cp2k_source,
                save_tblite=tblite,
                save_tblite_source=save_source,
                campaign_manifest=campaign_manifest,
                workflow_paths={
                    "experimental_k222_preflight_root": preflight_root,
                    "k222_cellopt_root": cellopt_root,
                    "k222_source_policy": "experimental_reference",
                    "final_single_point_roots": {
                        "k333": k333_root,
                        "k444": k444_root,
                    },
                },
            )
            campaign = common.load_campaign_identity(root)

            gas_stem = "case_GXTB_mol_geoopt"
            gas_input = root / "inputs" / "molecule_geoopt" / "GXTB" / f"{gas_stem}.inp"
            gas_output = root / "runs" / "molecule_geoopt" / "GXTB" / gas_stem / f"{gas_stem}.out"
            gas_input.parent.mkdir(parents=True)
            gas_output.parent.mkdir(parents=True)
            gas_input.write_text("METHOD GXTB\n")
            gas_output.write_text("GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n")
            common.write_job_stamp(
                gas_output.parent,
                gas_input,
                cp2k,
                "GXTB",
                "x23b_molecule_geoopt",
                "converged",
                details={
                    "output": str(gas_output.resolve()),
                    "output_sha256": common.sha256_file(gas_output),
                },
                campaign_identity=campaign,
            )

            reference = base / "reference.inp"
            structure = base / "structure.cif"
            reference.write_text("reference\n")
            structure.write_text("structure\n")
            preflight_run = preflight_root / "GXTB" / "case" / experimental_preflight.VARIANT
            preflight_run.mkdir(parents=True)
            preflight_input = preflight_run / "case_GXTB_experimental_k222_preflight.inp"
            preflight_output = preflight_run / "cp2k.out"
            preflight_input.write_text("METHOD GXTB\n")
            preflight_output.write_text(
                "ENERGY| Total FORCE_EVAL -1.0\n"
                "FORCES| Atomic forces [hartree/bohr]\n"
                "STRESS| Analytical stress tensor [GPa]\nPROGRAM ENDED\n"
            )
            preflight_stamp = common.write_job_stamp(
                preflight_run,
                preflight_input,
                cp2k,
                "GXTB",
                experimental_preflight.PHASE,
                "converged_measured_not_approved",
                details={
                    "output": str(preflight_output.resolve()),
                    "output_sha256": common.sha256_file(preflight_output),
                },
                campaign_identity=campaign,
                protocol_identity=experimental_preflight.PROTOCOL_IDENTITY,
                source_artifacts={
                    "reference_input": reference,
                    "reference_structure": structure,
                },
            )
            (preflight_root / experimental_preflight.MANIFEST_NAME).write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "phase": experimental_preflight.PHASE,
                        "campaign_identity": campaign,
                        "systems": [
                            {
                                "system": "case",
                                "input": str(preflight_input),
                                "output": str(preflight_output),
                                "source_input": str(reference),
                                "structure_path": str(structure),
                            }
                        ],
                    }
                )
                + "\n"
            )

            cellopt_run = cellopt_root / "GXTB" / "case" / cellopt.EXPERIMENTAL_VARIANT
            cellopt_run.mkdir(parents=True)
            cellopt_input = cellopt_run / "case_GXTB_experimental.inp"
            cellopt_output = cellopt_run / "cp2k.out"
            cellopt_restart = cellopt_run / "case-1.restart"
            cellopt_input.write_text("METHOD GXTB\n")
            cellopt_output.write_text("GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n")
            cellopt_restart.write_text("METHOD GXTB\n")
            lineage = "frozen_x23_reference_structure->experimental_k222_preflight->k222_cellopt_input"
            source_files = {
                "reference_input": reference,
                "reference_structure": structure,
                "preflight_input": preflight_input,
                "preflight_output": preflight_output,
                "preflight_stamp": preflight_stamp,
            }
            cellopt_protocol = {
                "source_policy": "experimental_reference",
                "variant": cellopt.EXPERIMENTAL_VARIANT,
                "lineage": lineage,
            }
            common.write_job_stamp(
                cellopt_run,
                cellopt_input,
                cp2k,
                "GXTB",
                "x23b_k222_cellopt",
                "converged",
                details={
                    "output": str(cellopt_output.resolve()),
                    "output_sha256": common.sha256_file(cellopt_output),
                },
                campaign_identity=campaign,
                protocol_identity=cellopt_protocol,
                source_artifacts=source_files,
            )
            with cellopt.manifest_path(cellopt_root).open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=cellopt.MANIFEST_FIELDS,
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "schema": cellopt.MANIFEST_SCHEMA,
                        "method": "GXTB",
                        "system": "case",
                        "variant": cellopt.EXPERIMENTAL_VARIANT,
                        "source_policy": "experimental_reference",
                        "source_kind": "experimental_reference",
                        "lineage": lineage,
                        "source_path": str(reference),
                        "source_sha256": common.sha256_file(reference),
                        "structure_path": str(structure),
                        "structure_sha256": common.sha256_file(structure),
                        "preflight_input": str(preflight_input),
                        "preflight_input_sha256": common.sha256_file(preflight_input),
                        "preflight_output": str(preflight_output),
                        "preflight_output_sha256": common.sha256_file(preflight_output),
                        "preflight_stamp": str(preflight_stamp),
                        "preflight_stamp_sha256": common.sha256_file(preflight_stamp),
                        "input": str(cellopt_input),
                        "input_sha256": common.sha256_file(cellopt_input),
                        "run_dir": str(cellopt_run),
                    }
                )

            for mesh, target_root in ((3, k333_root), (4, k444_root)):
                mesh_id = f"k{mesh}{mesh}{mesh}"
                target_run = target_root / "GXTB" / "case" / f"{mesh_id}_sp_on_k222"
                target_run.mkdir(parents=True)
                target_input = target_run / f"case_GXTB_{mesh_id}_sp_on_k222.inp"
                target_output = target_run / "cp2k.out"
                target_input.write_text("METHOD GXTB\n")
                target_output.write_text("ENERGY| Total FORCE_EVAL -1.0\nPROGRAM ENDED\n")
                target_protocol = {
                    "source_policy": "experimental_reference",
                    "source_variant": cellopt.EXPERIMENTAL_VARIANT,
                    "source_protocol_identity": cellopt_protocol,
                    "target_mesh": mesh_id,
                }
                target_sources = dict(source_files)
                target_sources.update(
                    {
                        "cellopt_input": cellopt_input,
                        "cellopt_output": cellopt_output,
                        "cellopt_restart": cellopt_restart,
                    }
                )
                common.write_job_stamp(
                    target_run,
                    target_input,
                    cp2k,
                    "GXTB",
                    f"x23b_final_{mesh_id}_on_k222",
                    "converged",
                    details={
                        "output": str(target_output.resolve()),
                        "output_sha256": common.sha256_file(target_output),
                    },
                    campaign_identity=campaign,
                    protocol_identity=target_protocol,
                    source_artifacts=target_sources,
                )
                with final_sp.manifest_path(target_root).open("w", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=final_sp.MANIFEST_FIELDS,
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow(
                        {
                            "method": "GXTB",
                            "system": "case",
                            "mesh": str(mesh),
                            "source_run_dir": str(cellopt_run),
                            "source_input": str(cellopt_input),
                            "source_input_sha256": common.sha256_file(cellopt_input),
                            "source_output_sha256": common.sha256_file(cellopt_output),
                            "source_policy": "experimental_reference",
                            "source_variant": cellopt.EXPERIMENTAL_VARIANT,
                            "source_protocol_identity": json.dumps(
                                cellopt_protocol, sort_keys=True
                            ),
                            "source_artifacts": json.dumps(
                                {
                                    name: str(path.resolve())
                                    for name, path in source_files.items()
                                },
                                sort_keys=True,
                            ),
                            "source_artifact_hashes": json.dumps(
                                {
                                    name: common.sha256_file(path)
                                    for name, path in source_files.items()
                                },
                                sort_keys=True,
                            ),
                            "source_restart": str(cellopt_restart),
                            "source_restart_sha256": common.sha256_file(cellopt_restart),
                            "input": str(target_input),
                            "input_sha256": common.sha256_file(target_input),
                            "run_dir": str(target_run),
                        }
                    )

            provenance = json.loads(common.update_gxtb_provenance(root).read_text())
            self.assertEqual(provenance["validation"]["gamma_cell_optimizations_completed"], 0)
            self.assertEqual(provenance["validation"]["experimental_k222_preflight_completed"], 1)
            self.assertEqual(provenance["validation"]["k222_cell_optimizations_completed"], 1)
            self.assertEqual(provenance["status"], "production_complete")

    def test_provenance_rejects_cp2k_binary_source_revision_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "X23b"
            (root / "data").mkdir(parents=True)
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(
                Path(tmp)
            )
            cp2k.write_text(
                "#!/bin/sh\n"
                "echo 'CP2K version 2026.1'\n"
                "echo ' Source code revision deadbeef00'\n"
            )
            cp2k.chmod(cp2k.stat().st_mode | 0o111)
            with self.assertRaisesRegex(ValueError, "executable/source revision mismatch"):
                common.update_gxtb_provenance(
                    root,
                    cp2k=cp2k,
                    cp2k_source=cp2k_source,
                    save_tblite=tblite,
                    save_tblite_source=save_source,
                    campaign_manifest=manifest,
                )
            self.assertFalse((root / "data" / common.GXTB_PROVENANCE_NAME).exists())

    def test_mixed_campaign_stamp_is_not_counted_or_collectable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "X23b"
            (root / "data").mkdir(parents=True)
            (root / "data" / "metadata.json").write_text(
                json.dumps({"systems": [{"id": "case"}]})
            )
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(base)
            common.update_gxtb_provenance(
                root,
                cp2k=cp2k,
                cp2k_source=cp2k_source,
                save_tblite=tblite,
                save_tblite_source=save_source,
                campaign_manifest=manifest,
            )
            current = common.load_campaign_identity(root)

            stem = "case_GXTB_mol_geoopt"
            input_path = root / "inputs" / "molecule_geoopt" / "GXTB" / f"{stem}.inp"
            output = root / "runs" / "molecule_geoopt" / "GXTB" / stem / f"{stem}.out"
            input_path.parent.mkdir(parents=True)
            output.parent.mkdir(parents=True)
            input_path.write_text("METHOD GXTB\n")
            output.write_text("GEOMETRY OPTIMIZATION COMPLETED\nPROGRAM ENDED\n")
            foreign_cp2k = fake_executable(base, "foreign-cp2k")
            foreign = fake_campaign(foreign_cp2k, discriminator="a")
            common.write_job_stamp(
                output.parent,
                input_path,
                foreign_cp2k,
                "GXTB",
                "x23b_molecule_geoopt",
                "converged",
                details={"output": str(output.resolve()), "output_sha256": common.sha256_file(output)},
                campaign_identity=foreign,
            )

            path = common.update_gxtb_provenance(root)
            payload = json.loads(path.read_text())
            self.assertEqual(payload["validation"]["gas_optimizations_completed"], 0)
            valid, reason = common.recorded_job_stamp_matches(
                output.parent,
                input_path,
                "GXTB",
                "x23b_molecule_geoopt",
                output,
                campaign_identity=current,
            )
            self.assertFalse(valid)
            self.assertIn("campaign fingerprint differs", reason)

    def test_campaign_manifest_must_be_production_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "X23b"
            (root / "data").mkdir(parents=True)
            cp2k, cp2k_source, tblite, save_source, manifest = fake_campaign_artifacts(base)
            payload = json.loads(manifest.read_text())
            payload["campaign_state"] = "validation_in_progress"
            manifest.write_text(json.dumps(payload) + "\n")
            with self.assertRaisesRegex(ValueError, "not production_ready"):
                common.update_gxtb_provenance(
                    root,
                    cp2k=cp2k,
                    cp2k_source=cp2k_source,
                    save_tblite=tblite,
                    save_tblite_source=save_source,
                    campaign_manifest=manifest,
                )

    def test_thread_environment_pins_veclib_and_passive_wait(self) -> None:
        env = common.thread_environment(3)
        self.assertEqual(env["OMP_NUM_THREADS"], "3")
        self.assertEqual(env["VECLIB_MAXIMUM_THREADS"], "1")
        self.assertEqual(env["OMP_WAIT_POLICY"], "PASSIVE")


if __name__ == "__main__":
    unittest.main()
