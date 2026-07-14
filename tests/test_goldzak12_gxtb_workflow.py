from __future__ import annotations

import csv
import json
import math
import re
import subprocess
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


def fake_campaign(
    *,
    loaded_library_sha: str = "libcp2k-one",
    cp2k_executable_sha: str = "test-launcher",
) -> dict[str, object]:
    return base.make_campaign_identity(
        campaign_id="test-campaign",
        cp2k_executable_sha256=cp2k_executable_sha,
        cp2k_loaded_library_sha256=loaded_library_sha,
        cp2k_cmake_cache_sha256="test-cp2k-cache",
        cp2k_embedded_source_revision="a" * 10,
        cp2k_source_revision="a" * 40,
        save_tblite_executable_sha256="test-save-cli",
        save_tblite_source_revision="b" * 40,
        save_tblite_library_sha256="test-libtblite",
        save_tblite_cmake_cache_sha256="test-save-cache",
        dependency_lock_sha256="test-dependency-lock",
    )


class Goldzak12GXTBInputTests(unittest.TestCase):
    def test_hysteresis_note_metrics_are_derived_from_the_pinned_eos_table(self) -> None:
        table = REPOSITORY / "Goldzak12" / "data" / "eos_fits.csv"
        note = REPOSITORY / "Goldzak12" / "data" / "gxtb_wfn_hysteresis.md"
        with table.open() as handle:
            rows = list(csv.DictReader(handle))
        gxtb = [
            row
            for row in rows
            if row["method"] == "GXTB"
            and row["fit_status"] == "quadratic"
            and row["a_eos_A"]
        ]
        self.assertEqual(
            {row["solid"] for row in gxtb},
            {"C", "Si", "SiC", "BN", "BP", "AlN", "AlP", "MgS", "LiF", "LiCl"},
        )
        errors = [float(row["a_eos_A"]) - float(row["a_exp_A"]) for row in gxtb]
        expected = {
            "ME": sum(errors) / len(errors),
            "MAE": sum(abs(error) for error in errors) / len(errors),
            "RMSE": math.sqrt(sum(error * error for error in errors) / len(errors)),
            "MaxAE": max(abs(error) for error in errors),
        }
        text = note.read_text()
        match = re.search(
            r"SHA256 `(?P<sha>[0-9a-f]{64})`\) are "
            r"ME `(?P<ME>[-+0-9.]+)` angstrom, "
            r"MAE `(?P<MAE>[-+0-9.]+)` angstrom, "
            r"RMSE `(?P<RMSE>[-+0-9.]+)` angstrom, and "
            r"MaxAE `(?P<MaxAE>[-+0-9.]+)` angstrom",
            text,
        )
        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(match.group("sha"), base.sha256(table))
        for label, value in expected.items():
            self.assertAlmostEqual(float(match.group(label)), value, places=10)
        self.assertNotIn("0.16713620093", text)

    def test_solid_input_pins_native_fdiis_and_shared_spglib_contract(self) -> None:
        text = base.solid_input(base.REFERENCES[0], "GXTB", "ENERGY", "k444", 3.553, "C_GXTB")
        self.assertIn("METHOD GXTB", text)
        self.assertIn("SCC_MIXER TBLITE", text)
        self.assertIn("METHOD DIRECT_P_MIXING", text)
        self.assertIn("BACKUP_COPIES 0", text)
        self.assertIn("&RESTART OFF", text)
        self.assertNotIn("SCC_MIXER CP2K", text)
        self.assertNotIn("&TBLITE_MIXER", text)
        self.assertIn("FULL_GRID F", text)
        self.assertIn("SYMMETRY T", text)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", text)
        self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", text)
        self.assertNotIn("STRESS_TENSOR ANALYTICAL", text)

        legacy = base.solid_input(base.REFERENCES[0], "GFN2", "ENERGY", "k444", 3.553, "C_GFN2")
        self.assertIn("FULL_GRID F", legacy)
        self.assertIn("SYMMETRY T", legacy)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", legacy)
        self.assertIn("STRESS_TENSOR ANALYTICAL", legacy)
        stale = text.replace("FULL_GRID F", "FULL_GRID T", 1)
        with self.assertRaisesRegex(ValueError, "SPGLIB-reduced mesh contract"):
            base.validate_method_input(stale, "GXTB")

    def test_gxtb_backup_copies_is_scoped_to_scf_restart_print_key(self) -> None:
        inputs = (
            base.solid_input(
                base.REFERENCES[0], "GXTB", "ENERGY", "k444", 3.553, "C_GXTB"
            ),
            base.atom_input("C", "GXTB"),
        )
        for text in inputs:
            with self.subTest(project=text.split("PROJECT", 1)[1].splitlines()[0].strip()):
                global_section = text.split("&END GLOBAL", 1)[0]
                self.assertNotIn("BACKUP_COPIES", global_section)
                self.assertEqual(text.count("BACKUP_COPIES 0"), 1)
                self.assertRegex(
                    text,
                    r"&RESTART OFF\n\s+BACKUP_COPIES 0\n\s+&END RESTART",
                )

        malformed = inputs[1].replace(
            "&END GLOBAL", "  BACKUP_COPIES 0\n&END GLOBAL", 1
        )
        with self.assertRaisesRegex(ValueError, "not valid in GLOBAL"):
            base.validate_method_input(malformed, "GXTB")

    def test_gxtb_energy_inputs_do_not_request_stress(self) -> None:
        ref = base.REFERENCES[0]
        eos_input = base.solid_input(ref, "GXTB", "ENERGY", "k444", ref.a_exp, "gxtb_eos")
        atom_input = base.atom_input("C", "GXTB")
        cellopt_input = base.solid_input(
            ref, "GXTB", "CELL_OPT", "k444", ref.a_exp, "gxtb_cellopt"
        )
        self.assertNotIn("STRESS_TENSOR", eos_input)
        self.assertNotIn("STRESS_TENSOR", atom_input)
        self.assertIn("STRESS_TENSOR ANALYTICAL", cellopt_input)
        with self.assertRaisesRegex(ValueError, "must not request analytical stress"):
            base.validate_method_input(
                eos_input.replace("  &DFT", "  STRESS_TENSOR ANALYTICAL\n  &DFT", 1),
                "GXTB",
            )

    def test_gxtb_atom_uses_supported_no_smear_ot_exception(self) -> None:
        self.assertEqual(
            base.GXTB_NO_SMEAR_OT_ATOMS,
            frozenset(base.ELEMENT_MULTIPLICITY) - {"Li"},
        )
        text = base.atom_input("H", "GXTB")
        self.assertIn("SCC_MIXER CP2K", text)
        self.assertIn("&OT", text)
        self.assertIn("MINIMIZER DIIS", text)
        self.assertIn("PRECONDITIONER FULL_SINGLE_INVERSE", text)
        self.assertIn("&SMEAR OFF", text)
        self.assertNotIn("ADDED_MOS", text)
        self.assertNotIn("DIRECT_P_MIXING", text)
        self.assertNotIn("&KPOINTS", text)
        base.validate_method_input(text, "GXTB", gxtb_atom_reference=True)
        with self.assertRaisesRegex(ValueError, "SCC_MIXER TBLITE"):
            base.validate_method_input(text, "GXTB")

        lithium = base.atom_input("Li", "GXTB")
        self.assertIn("SCC_MIXER TBLITE", lithium)
        self.assertIn("DIRECT_P_MIXING", lithium)
        self.assertNotIn("&OT", lithium)
        base.validate_method_input(lithium, "GXTB")

        legacy = base.atom_input("H", "GFN2")
        self.assertNotIn("&OT", legacy)
        self.assertIn("DIRECT_P_MIXING", legacy)
        self.assertNotIn("SCC_MIXER CP2K", legacy)

    def test_default_gxtb_plan_has_132_eos_and_36_maximum_sp_jobs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            old_eos_root = eos.ROOT
            base.ROOT = eos.ROOT = Path(tmp)
            try:
                specs = eos.eos_job_specs("k444", eos.DEFAULT_SCALES, ("GXTB",))
                self.assertEqual(len(specs), 12 * len(eos.DEFAULT_SCALES))
                self.assertTrue(all("/GXTB/" in str(spec[1]) for spec in specs))
                common_ten = (
                    "C", "Si", "SiC", "BN", "BP", "AlN", "AlP", "MgS", "LiF", "LiCl"
                )
                selected = eos.eos_job_specs(
                    "k444", eos.DEFAULT_SCALES, ("GXTB",), common_ten
                )
                self.assertEqual(len(selected), 10 * len(eos.DEFAULT_SCALES))
                self.assertFalse(any(" LiH " in spec[0] or " MgO " in spec[0] for spec in selected))
                fits = [
                    {
                        "solid": ref.solid,
                        "method": "GXTB",
                        "eos_mesh": "k444",
                        "a_eos_A": f"{ref.a_exp + 0.01:.8f}",
                        "fit_status": "quadratic",
                    }
                    for ref in base.REFERENCES
                ]
                self.assertEqual(len(eos.final_sp_specs(fits, ["k333", "k444", "k555"])), 36)
            finally:
                base.ROOT = old_root
                eos.ROOT = old_eos_root

    def test_final_input_requires_current_eos_lineage_and_regenerates_stale_reference_cell(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            old_eos_root = eos.ROOT
            base.ROOT = eos.ROOT = Path(tmp)
            try:
                ref = base.REFERENCES[0]
                mesh = "k555"
                input_path = eos.final_input_path(ref.solid, "GXTB", mesh)
                stale_text = base.solid_input(
                    ref,
                    "GXTB",
                    "ENERGY",
                    mesh,
                    ref.a_exp,
                    eos.final_project(ref.solid, "GXTB", mesh),
                )
                base.write_file(input_path, stale_text)
                invalid_fit = {
                    "solid": ref.solid,
                    "method": "GXTB",
                    "eos_mesh": "k444",
                    "a_eos_A": "",
                    "fit_status": "no_local_minimum",
                }
                eos.invalidate_existing_gxtb_final_inputs([invalid_fit], [mesh])
                invalid_lineage = json.loads(eos.final_input_lineage_path(input_path).read_text())
                self.assertFalse(invalid_lineage["valid"])
                self.assertEqual(eos.final_sp_specs([invalid_fit], [mesh]), [])
                self.assertEqual(input_path.read_text(), stale_text)

                a_eos = ref.a_exp + 0.12345678
                valid_fit = {
                    "solid": ref.solid,
                    "method": "GXTB",
                    "eos_mesh": "k444",
                    "a_eos_A": f"{a_eos:.8f}",
                    "fit_status": "quadratic",
                }
                specs = eos.final_sp_specs([valid_fit], [mesh])
                self.assertEqual(len(specs), 1)
                self.assertIn(f"ABC {a_eos:.12f} {a_eos:.12f} {a_eos:.12f}", input_path.read_text())
                lineage = json.loads(eos.final_input_lineage_path(input_path).read_text())
                self.assertTrue(lineage["valid"])
                self.assertEqual(lineage["input_sha256"], base.sha256(input_path))
                self.assertIsNone(eos.final_input_lineage_issue(input_path, valid_fit, mesh))
            finally:
                base.ROOT = old_root
                eos.ROOT = old_eos_root


class Goldzak12GXTBAtomTests(unittest.TestCase):
    def test_save_tblite_cli_uses_13_documented_atomic_spins(self) -> None:
        commands: list[list[str]] = []

        def fake_run(
            command: list[str], cwd: Path, stdout: object, stderr: object, env: dict[str, str]
        ) -> subprocess.CompletedProcess:
            commands.append(command)
            self.assertEqual(env["VECLIB_MAXIMUM_THREADS"], "1")
            self.assertEqual(env["OPENBLAS_NUM_THREADS"], "1")
            self.assertEqual(env["MKL_NUM_THREADS"], "1")
            self.assertEqual(env["OMP_WAIT_POLICY"], "PASSIVE")
            json_name = command[command.index("--json") + 1]
            (Path(cwd) / json_name).write_text(json.dumps({"energy": -1.0}))
            return subprocess.CompletedProcess(command, 0)

        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            base.ROOT = Path(tmp)
            try:
                with patch.object(base.subprocess, "run", side_effect=fake_run):
                    base.run_tblite_atom_jobs(
                        Path("tblite"),
                        jobs=1,
                        force=True,
                        methods=("GXTB",),
                        save_tblite=Path("save_tblite"),
                        campaign_fingerprint=fake_campaign(),
                    )
            finally:
                base.ROOT = old_root
        self.assertEqual(len(commands), 13)
        by_element = {command[-1].split("atom_")[1].split(".")[0]: command for command in commands}
        for element, multiplicity in base.ELEMENT_MULTIPLICITY.items():
            command = by_element[element]
            self.assertEqual(command[0], "save_tblite")
            self.assertEqual(command[command.index("--method") + 1], "gxtb")
            self.assertEqual(int(command[command.index("--spin") + 1]), multiplicity - 1)

    def test_additive_merge_preserves_frozen_gfn_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "energies.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=("method", "element", "energy_hartree"))
                writer.writeheader()
                writer.writerow({"method": "GFN1", "element": "C", "energy_hartree": "-1.23"})
            base.merge_method_rows(
                path,
                [{"method": "GXTB", "element": "C", "energy_hartree": "-4.56"}],
                ("GXTB",),
            )
            rows = base.read_csv(path)
        self.assertEqual(rows[0]["energy_hartree"], "-1.23")
        self.assertEqual({row["method"] for row in rows}, {"GFN1", "GXTB"})

    def test_failed_save_tblite_atom_is_fatal_after_other_atoms_finish(self) -> None:
        calls: list[str] = []

        def fake_run(
            command: list[str], cwd: Path, stdout: object, stderr: object, env: dict[str, str]
        ) -> subprocess.CompletedProcess:
            element = command[-1].split("atom_")[1].split(".")[0]
            calls.append(element)
            if element != "C":
                json_name = command[command.index("--json") + 1]
                (Path(cwd) / json_name).write_text(json.dumps({"energy": -1.0}))
                return subprocess.CompletedProcess(command, 0)
            return subprocess.CompletedProcess(command, 2)

        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            base.ROOT = Path(tmp)
            try:
                with patch.object(base.subprocess, "run", side_effect=fake_run):
                    with self.assertRaises(RuntimeError):
                        base.run_tblite_atom_jobs(
                            Path("tblite"),
                            jobs=4,
                            force=True,
                            methods=("GXTB",),
                            save_tblite=Path("save_tblite"),
                            campaign_fingerprint=fake_campaign(),
                        )
            finally:
                base.ROOT = old_root
        self.assertEqual(set(calls), set(base.ELEMENT_MULTIPLICITY))

    def test_gxtb_atom_resume_requires_matching_input_and_executable_hashes(self) -> None:
        calls: list[list[str]] = []

        def fake_run(
            command: list[str], cwd: Path, stdout: object, stderr: object, env: dict[str, str]
        ) -> subprocess.CompletedProcess:
            calls.append(command)
            json_name = command[command.index("--json") + 1]
            (Path(cwd) / json_name).write_text(json.dumps({"energy": -1.0}))
            return subprocess.CompletedProcess(command, 0)

        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            base.ROOT = Path(tmp)
            executable = Path(tmp) / "save_tblite"
            executable.write_bytes(b"version-one")
            try:
                with patch.object(base.subprocess, "run", side_effect=fake_run):
                    base.run_tblite_atom_jobs(
                        executable,
                        2,
                        True,
                        ("GXTB",),
                        executable,
                        fake_campaign(),
                    )
                    base.run_tblite_atom_jobs(
                        executable,
                        2,
                        False,
                        ("GXTB",),
                        executable,
                        fake_campaign(),
                    )
                    self.assertEqual(len(calls), 13)
                    executable.write_bytes(b"version-two")
                    base.run_tblite_atom_jobs(
                        executable,
                        2,
                        False,
                        ("GXTB",),
                        executable,
                        fake_campaign(),
                    )
            finally:
                base.ROOT = old_root
        self.assertEqual(len(calls), 26)


class Goldzak12GXTBMixerAndPruneTests(unittest.TestCase):
    def test_execution_record_is_additive_and_pool_only_runs_pending_jobs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "cp2k"
            executable.write_bytes(b"cp2k")
            inp = root / "job.inp"
            out = root / "job.out"
            inp.write_text(
                base.solid_input(
                    base.REFERENCES[0], "GXTB", "ENERGY", "k444", 3.553, "job"
                )
            )
            campaign = fake_campaign(
                cp2k_executable_sha=base.sha256(executable)
            )
            signature = base.job_signature(
                executable,
                inp,
                command_contract={"driver": "cp2k", "omp_threads": 1},
                campaign_fingerprint=campaign,
            )

            class FakePool:
                mpi_ranks_per_job = 2

                def __init__(self) -> None:
                    self.complete = False
                    self.run_calls = 0
                    self.write_calls = 0

                def record_issue(self, _output: Path, _stamp: Path) -> str | None:
                    return None if self.complete else "missing execution record"

                def run_cp2k(
                    self, _cp2k: Path, _input: Path, output: Path
                ) -> tuple[int, dict[str, object]]:
                    self.run_calls += 1
                    output.write_text("PROGRAM ENDED\n")
                    return 0, {"separate": True}

                def write_record(
                    self,
                    _output: Path,
                    _observation: dict[str, object],
                    _stamp: Path,
                ) -> dict[str, str]:
                    self.write_calls += 1
                    self.complete = True
                    return {"path": "execution.json", "sha256": "fixture"}

            pool = FakePool()
            spec = [("eos GXTB C k444 s1p000", inp, out, False)]
            eos.run_jobs(
                spec,
                executable,
                1,
                1,
                False,
                campaign_fingerprint=campaign,
                execution_pool=pool,  # type: ignore[arg-type]
            )
            self.assertEqual((pool.run_calls, pool.write_calls), (1, 1))
            stamp = json.loads(base.job_stamp_path(out).read_text())
            self.assertEqual(stamp, {**signature, "completed": True, "return_code": 0})
            self.assertNotIn("execution_provenance", stamp)
            frozen_output = out.read_bytes()
            pool.complete = False
            with self.assertRaisesRegex(RuntimeError, "refusing an implicit destructive rerun"):
                eos.run_jobs(
                    spec,
                    executable,
                    1,
                    1,
                    False,
                    campaign_fingerprint=campaign,
                    execution_pool=pool,  # type: ignore[arg-type]
                )
            self.assertEqual(out.read_bytes(), frozen_output)
            self.assertEqual((pool.run_calls, pool.write_calls), (1, 1))
            pool.complete = True
            eos.run_jobs(
                spec,
                executable,
                1,
                1,
                False,
                campaign_fingerprint=campaign,
                execution_pool=pool,  # type: ignore[arg-type]
            )
            self.assertEqual((pool.run_calls, pool.write_calls), (1, 1))

    def test_failed_gxtb_job_is_not_retried_with_an_alternative_mixer(self) -> None:
        calls: list[Path] = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inp = root / "job.inp"
            out = root / "job.out"
            inp.write_text(
                base.solid_input(base.REFERENCES[0], "GXTB", "ENERGY", "k444", 3.553, "job")
            )

            def fail(_cp2k: Path, input_path: Path, _out: Path, _threads: int) -> int:
                calls.append(input_path)
                return 1

            with patch.object(base, "run_cp2k", side_effect=fail):
                with self.assertRaises(RuntimeError):
                    eos.run_jobs(
                        [("eos GXTB C k444 s1p000", inp, out, False)],
                        Path("cp2k"),
                        1,
                        1,
                        False,
                        campaign_fingerprint=fake_campaign(),
                    )
        self.assertEqual(calls, [inp])

    def test_successful_gxtb_cp2k_resume_requires_matching_hash_stamp(self) -> None:
        calls: list[Path] = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "cp2k"
            executable.write_bytes(b"cp2k-one")
            inp = root / "job.inp"
            out = root / "job.out"
            inp.write_text(
                base.solid_input(base.REFERENCES[0], "GXTB", "ENERGY", "k444", 3.553, "job")
            )

            def succeed(_cp2k: Path, input_path: Path, output: Path, _threads: int) -> int:
                calls.append(input_path)
                output.write_text("PROGRAM ENDED\n")
                return 0

            spec = [("eos GXTB C k444 s1p000", inp, out, False)]
            with patch.object(base, "run_cp2k", side_effect=succeed):
                campaign = fake_campaign(loaded_library_sha="libcp2k-one")
                eos.run_jobs(spec, executable, 1, 1, False, campaign_fingerprint=campaign)
                eos.run_jobs(spec, executable, 1, 1, False, campaign_fingerprint=campaign)
                self.assertEqual(len(calls), 1)
                changed_library = fake_campaign(loaded_library_sha="libcp2k-two")
                eos.run_jobs(
                    spec,
                    executable,
                    1,
                    1,
                    False,
                    campaign_fingerprint=changed_library,
                )
                self.assertEqual(len(calls), 2)
                inp.write_text(inp.read_text() + "! changed input fingerprint\n")
                eos.run_jobs(
                    spec,
                    executable,
                    1,
                    1,
                    False,
                    campaign_fingerprint=changed_library,
                )
        self.assertEqual(len(calls), 3)

    def test_pruning_requires_validated_output_and_never_touches_gfn(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            old_root = base.ROOT
            base.ROOT = Path(tmp)
            try:
                valid = base.ROOT / "runs" / "eos" / "GXTB" / "C"
                invalid = base.ROOT / "runs" / "eos" / "GXTB" / "Si"
                legacy = base.ROOT / "runs" / "eos" / "GFN2" / "C"
                for directory in (valid, invalid, legacy):
                    directory.mkdir(parents=True)
                    (directory / "case-RESTART.kp").write_text("large")
                (valid / "case.out").write_text("PROGRAM ENDED")
                (invalid / "case.out").write_text("SCF run NOT converged")
                (legacy / "case.out").write_text("PROGRAM ENDED")
                count, _ = base.prune_gxtb_transients((base.ROOT / "runs" / "eos" / "GXTB",))
                self.assertEqual(count, 1)
                self.assertFalse((valid / "case-RESTART.kp").exists())
                self.assertTrue((invalid / "case-RESTART.kp").exists())
                self.assertTrue((legacy / "case-RESTART.kp").exists())
            finally:
                base.ROOT = old_root


class Goldzak12GXTBCampaignFingerprintTests(unittest.TestCase):
    def test_campaign_identity_is_path_independent_and_requires_production_gate(self) -> None:
        manifest = {
            "campaign_id": "test",
            "campaign_state": "production_ready",
            "cp2k": {
                "binary_sha256": "1" * 64,
                "loaded_library_sha256": "2" * 64,
                "cmake_cache_sha256": "3" * 64,
                "reported_revision": "a" * 10,
                "revision": "a" * 40,
                "binary": "/host-one/cp2k",
            },
            "save_tblite": {
                "cli_sha256": "4" * 64,
                "revision": "b" * 40,
                "static_library_sha256": "5" * 64,
                "cmake_cache_sha256": "6" * 64,
                "cli": "/host-one/tblite",
            },
            "fetched_dependencies": {"tblite": "c" * 40},
        }
        first = base.campaign_identity_from_manifest(manifest, Path("first.json"))
        relocated = json.loads(json.dumps(manifest))
        relocated["cp2k"]["binary"] = "/host-two/cp2k"
        relocated["save_tblite"]["cli"] = "/host-two/tblite"
        relocated["administrative_note"] = "path and notes are provenance, not identity"
        second = base.campaign_identity_from_manifest(relocated, Path("second.json"))
        self.assertEqual(first, second)
        relocated["campaign_state"] = "validation_in_progress"
        with self.assertRaisesRegex(ValueError, "not production_ready"):
            base.campaign_identity_from_manifest(relocated, Path("second.json"))
        diagnostic = base.campaign_identity_from_manifest(
            relocated,
            Path("second.json"),
            allowed_campaign_states=("validation_in_progress",),
        )
        self.assertEqual(diagnostic, first)
        with self.assertRaisesRegex(ValueError, r"allowed state\(s\): qualification_pending"):
            base.campaign_identity_from_manifest(
                relocated,
                Path("second.json"),
                allowed_campaign_states=("qualification_pending",),
            )

    def test_manifest_rejects_replaced_libtblite_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            artifacts = {
                "cp2k": root / "cp2k",
                "cp2k_library": root / "libcp2k.dylib",
                "save_tblite": root / "tblite",
                "save_tblite_library": root / "libtblite.a",
            }
            for key, path in artifacts.items():
                path.write_bytes(key.encode())
            manifest = root / "build_manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "campaign_id": "test",
                        "campaign_state": "production_ready",
                        "implementation": {"current_version": 1},
                        "cp2k": {
                            "binary": str(artifacts["cp2k"]),
                            "binary_sha256": base.sha256(artifacts["cp2k"]),
                            "loaded_library": str(artifacts["cp2k_library"]),
                            "loaded_library_sha256": base.sha256(artifacts["cp2k_library"]),
                            "cmake_cache_sha256": "cp2k-cache",
                            "revision": "a" * 40,
                            "reported_revision": "a" * 10,
                            "source_clean": True,
                        },
                        "save_tblite": {
                            "cli": str(artifacts["save_tblite"]),
                            "cli_sha256": base.sha256(artifacts["save_tblite"]),
                            "static_library": str(artifacts["save_tblite_library"]),
                            "static_library_sha256": base.sha256(
                                artifacts["save_tblite_library"]
                            ),
                            "cmake_cache_sha256": "save-cache",
                            "revision": "b" * 40,
                            "source_clean": True,
                        },
                        "fetched_dependencies": {"tblite": "c" * 40},
                    }
                )
            )
            source_fingerprint = {
                "schema_version": 1,
                "cp2k": {
                    "source": {"revision": "a" * 40},
                    "embedded_source_revision": "a" * 10,
                },
                "save_tblite": {
                    "source": {"revision": "b" * 40},
                    "version_output": "",
                },
            }
            with patch.object(
                base,
                "validated_gxtb_campaign_fingerprint",
                return_value=source_fingerprint,
            ):
                fingerprint, _ = base.validated_gxtb_campaign_from_manifest(
                    manifest, root, root
                )
                self.assertEqual(fingerprint["campaign_id"], "test")
                base.validate_campaign_identity(fingerprint)
                artifacts["save_tblite_library"].write_bytes(b"replacement")
                with self.assertRaisesRegex(ValueError, "save_tblite_library hash differs"):
                    base.validated_gxtb_campaign_from_manifest(manifest, root, root)

    def test_embedded_cp2k_revision_must_match_source_head(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k"
            cp2k.write_bytes(b"launcher")
            cp2k.chmod(0o755)
            library = root / "libcp2k.dylib"
            library.write_bytes(b"library")
            save_cli = root / "tblite"
            save_cli.write_bytes(b"cli")
            save_cli.chmod(0o755)
            save_library = root / "libtblite.a"
            save_library.write_bytes(b"archive")
            sources = [
                {
                    "available": True,
                    "revision": "a" * 40,
                    "branch": "g-xTB-pbc",
                    "dirty": False,
                },
                {
                    "available": True,
                    "revision": "b" * 40,
                    "branch": "cp2k-integration",
                    "dirty": False,
                },
            ]
            with (
                patch.object(base, "loaded_cp2k_library", return_value=library),
                patch.object(base, "_validated_clean_source", side_effect=sources),
                patch.object(
                    base,
                    "command_output",
                    return_value="CP2K version test\n Source code revision deadbeef00",
                ),
            ):
                with self.assertRaisesRegex(ValueError, "executable/source revision mismatch"):
                    base.validated_gxtb_campaign_fingerprint(
                        cp2k, library, root, save_cli, save_library, root
                    )

    def test_collector_stamp_rejects_mixed_loaded_library_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            executable = root / "cp2k"
            executable.write_bytes(b"launcher")
            inp = root / "case.inp"
            inp.write_text("input")
            result = root / "case.out"
            result.write_text("PROGRAM ENDED")
            campaign = fake_campaign(
                loaded_library_sha="library-one",
                cp2k_executable_sha=base.sha256(executable),
            )
            signature = base.job_signature(
                executable, inp, campaign_fingerprint=campaign
            )
            base.write_job_stamp(result, signature, completed=True, return_code=0)
            self.assertIsNone(
                base.completed_stamp_campaign_issue(
                    result, campaign, executable_role="cp2k"
                )
            )
            mixed = fake_campaign(
                loaded_library_sha="library-two",
                cp2k_executable_sha=base.sha256(executable),
            )
            self.assertIn(
                "campaign identity mismatch",
                base.completed_stamp_campaign_issue(
                    result, mixed, executable_role="cp2k"
                ),
            )


class Goldzak12GXTBProtocolTests(unittest.TestCase):
    def test_result_mesh_default_is_frozen_dense_mesh(self) -> None:
        self.assertEqual(eos.DEFAULT_RESULT_MESH, "k555")

    def test_gxtb_finals_need_explicit_fit_approval(self) -> None:
        self.assertFalse(
            eos.final_stage_is_explicitly_approved(
                ("GXTB",), stop_after_eos=False, fit_only=False, approve_fits=False
            )
        )
        self.assertFalse(
            eos.final_stage_is_explicitly_approved(
                ("GXTB",), stop_after_eos=True, fit_only=False, approve_fits=True
            )
        )
        self.assertTrue(
            eos.final_stage_is_explicitly_approved(
                ("GXTB",), stop_after_eos=False, fit_only=False, approve_fits=True
            )
        )

    def test_adaptive_scale_manifest_persists_exact_requested_grid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            old_eos_root = eos.ROOT
            old_values = dict(eos.ADAPTIVE_SCALES)
            eos.ROOT = Path(tmp)
            try:
                eos.ADAPTIVE_SCALES[("C", "GXTB")] = (0.99, 1.01)
                payload = eos.write_gxtb_scale_manifest("k444", eos.DEFAULT_SCALES, ("GXTB",))
                self.assertIsNotNone(payload)
                record = next(item for item in payload["systems"] if item["solid"] == "C")
                self.assertEqual(record["adaptive_scales"], [0.99, 1.01])
                self.assertIn(0.99, record["requested_scales"])
                self.assertIn(1.01, record["requested_scales"])
                eos.ADAPTIVE_SCALES.pop(("C", "GXTB"))
                eos.restore_gxtb_scale_manifest("k444", ("GXTB",))
                self.assertEqual(eos.ADAPTIVE_SCALES[("C", "GXTB")], (0.99, 1.01))
            finally:
                eos.ADAPTIVE_SCALES.clear()
                eos.ADAPTIVE_SCALES.update(old_values)
                eos.ROOT = old_eos_root

    def test_gxtb_classification_requires_explicit_action_and_rationale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "classifications.json"
            path.write_text(
                json.dumps(
                    {
                        "entries": [
                            {
                                "solid": "C",
                                "method": "GXTB",
                                "mesh": "k444",
                                "scale": 0.94,
                                "classification": "alternate_scc_branch",
                                "action": "retain",
                                "rationale": "Reviewed smooth continuation on the adaptive grid.",
                            }
                        ]
                    }
                )
            )
            loaded = eos.load_gxtb_classifications(path)
            self.assertEqual(loaded[("C", "k444", 0.94)]["action"], "retain")
            path.write_text(json.dumps({"entries": [{"solid": "C", "mesh": "k444", "scale": 0.94}]}))
            with self.assertRaises(ValueError):
                eos.load_gxtb_classifications(path)

    def test_gxtb_topology_gate_accepts_a_single_well(self) -> None:
        scales = (0.82, 0.88, 0.94, 0.98, 1.00, 1.02, 1.06)
        points = [(scale, scale, (scale - 0.98) ** 2, True) for scale in scales]
        self.assertEqual(eos.gxtb_topology_reversals(points), [])
        self.assertEqual(eos.fit_gxtb_eos(points)["fit_status"], "quadratic")

    def test_gxtb_topology_gate_rejects_a_lower_endpoint_branch(self) -> None:
        points = [
            (0.82, 0.82, -0.06710426, True),
            (0.88, 0.88, -0.05030883, True),
            (0.94, 0.94, 0.01189583, True),
            (0.98, 0.98, 0.00000000, True),
            (1.00, 1.00, 0.00947811, True),
            (1.02, 1.02, 0.02661851, True),
        ]
        fit = eos.fit_gxtb_eos(points)
        self.assertEqual(fit["fit_status"], "nonmonotonic_branch")
        self.assertEqual(fit["a_eos_A"], "")
        self.assertEqual(fit["grid_min_scale"], "0.82000")
        self.assertEqual(fit["topology_reversal_count"], 1)
        self.assertAlmostEqual(
            float(fit["topology_max_reversal_hartree"]), 0.01189583, places=8
        )

    def test_local_adaptive_points_cannot_heal_a_nonmonotonic_branch(self) -> None:
        points = [
            (0.82, 0.82, -0.06710426, True),
            (0.88, 0.88, -0.05030883, True),
            (0.94, 0.94, 0.01189583, True),
            (0.96, 0.96, 0.00400000, True),
            (0.97333, 0.97333, 0.00030000, True),
            (0.98, 0.98, 0.00000000, True),
            (0.98667, 0.98667, 0.00100000, True),
            (0.99, 0.99, 0.00300000, True),
            (1.00, 1.00, 0.00947811, True),
            (1.02, 1.02, 0.02661851, True),
        ]
        fit = eos.fit_gxtb_eos(points)
        self.assertEqual(fit["fit_status"], "nonmonotonic_branch")
        self.assertGreaterEqual(fit["topology_reversal_count"], 1)

    def test_reviewed_exclusion_can_leave_a_single_well(self) -> None:
        regular_branch = [
            (0.94, 0.94, 0.01189583, True),
            (0.98, 0.98, 0.00000000, True),
            (1.00, 1.00, 0.00947811, True),
            (1.02, 1.02, 0.02661851, True),
            (1.06, 1.06, 0.07500000, True),
        ]
        self.assertEqual(eos.gxtb_topology_reversals(regular_branch), [])
        self.assertEqual(eos.fit_gxtb_eos(regular_branch)["fit_status"], "quadratic")

    def test_reduced_coverage_needs_meaningful_and_adaptively_investigated_subset(self) -> None:
        fits = [
            {
                "solid": ref.solid,
                "method": "GXTB",
                "fit_status": "quadratic" if index < 8 else "no_local_minimum",
            }
            for index, ref in enumerate(base.REFERENCES)
        ]
        followup = [
            {"solid": ref.solid, "adaptive_investigated": True}
            for ref in base.REFERENCES[8:]
        ]
        with tempfile.TemporaryDirectory() as tmp:
            old_root = eos.ROOT
            eos.ROOT = Path(tmp)
            try:
                (eos.ROOT / "data").mkdir(parents=True)
                base.write_csv(eos.ROOT / "data" / "gxtb_eos_branch_diagnostics.csv", [])
                eos.enforce_gxtb_coverage(
                    fits,
                    followup,
                    allow_reduced_coverage=True,
                    minimum_valid_fits=8,
                )
                with self.assertRaises(RuntimeError):
                    eos.enforce_gxtb_coverage(
                        fits,
                        followup,
                        allow_reduced_coverage=True,
                        minimum_valid_fits=9,
                    )
            finally:
                eos.ROOT = old_root


if __name__ == "__main__":
    unittest.main()
