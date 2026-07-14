from __future__ import annotations

import importlib.util
import json
import math
import shutil
import sys
import tempfile
import threading
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock


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
    "dmc_ice13_kpoint_benchmark",
    REPOSITORY / "DMC-ICE13" / "scripts" / "dmc_ice13_kpoint_benchmark.py",
)
runner = load_script(
    "run_dmc13_kpoint_jobs",
    REPOSITORY / "scripts" / "run_dmc13_kpoint_jobs.py",
)
gamma_cli = load_script(
    "dmc_gxtb_gamma_cli_check",
    REPOSITORY / "DMC-ICE13" / "scripts" / "dmc_gxtb_gamma_cli_check.py",
)


def fake_production_identity(root: Path, cp2k: Path) -> runner.ProductionIdentity:
    cp2k_library = root / "libcp2k.test.dylib"
    cp2k_library.write_text("CP2K implementation library\n")
    tblite_static_library = root / "libtblite.a"
    tblite_static_library.write_text("static save_tblite archive\n")
    return runner.ProductionIdentity(
        campaign_id="test-gxtb-pbc-v1",
        cp2k=cp2k,
        cp2k_sha256=runner.sha256(cp2k),
        cp2k_library=cp2k_library,
        cp2k_library_sha256=runner.sha256(cp2k_library),
        tblite_static_library=tblite_static_library,
        tblite_static_library_sha256=runner.sha256(tblite_static_library),
        cp2k_source_revision="1" * 40,
        tblite_source_revision="2" * 40,
    )


def write_fake_campaign_manifest(
    path: Path,
    identity: runner.ProductionIdentity,
    *,
    frozen_identity: dict[str, str] | None = None,
) -> dict[str, object]:
    frozen = frozen_identity or runner.execution_build_identity(identity)
    payload = {
        "campaign_id": identity.campaign_id,
        "cp2k": {
            "revision": frozen["cp2k_source_revision"],
            "binary_sha256": frozen["cp2k_sha256"],
            "loaded_library_sha256": frozen["cp2k_library_sha256"],
        },
        "save_tblite": {
            "revision": frozen["tblite_source_revision"],
            "static_library_sha256": frozen["tblite_static_library_sha256"],
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def fake_gamma_identity(
    root: Path,
    cp2k: Path,
    tblite: Path,
) -> gamma_cli.BuildIdentity:
    cp2k_library = root / "libcp2k.test.dylib"
    cp2k_library.write_text("CP2K implementation library\n")
    tblite_static_library = root / "libtblite.a"
    tblite_static_library.write_text("static save_tblite archive\n")
    return gamma_cli.BuildIdentity(
        campaign_id="test-gxtb-pbc-v1",
        cp2k=cp2k,
        cp2k_sha256=gamma_cli.sha256(cp2k),
        cp2k_library=cp2k_library,
        cp2k_library_sha256=gamma_cli.sha256(cp2k_library),
        tblite=tblite,
        tblite_sha256=gamma_cli.sha256(tblite),
        tblite_static_library=tblite_static_library,
        tblite_static_library_sha256=gamma_cli.sha256(tblite_static_library),
        cp2k_source_revision="1" * 40,
        tblite_source_revision="2" * 40,
    )


def complete_job(
    job: runner.Job,
    identity: runner.ProductionIdentity,
    energy: float,
) -> None:
    job.run_dir.mkdir(parents=True, exist_ok=True)
    project = f"ice_{job.phase}_GXTB_{job.mesh}"
    (job.run_dir / job.output_name).write_text(
        "PROGRAM STARTED AT 2026-07-14 09:56:29.000\n"
        "PROGRAM STARTED ON Mac.test\n"
        f"CP2K| source code revision number: {identity.cp2k_source_revision[:10]}\n"
        "CP2K| cp2kflags: tblite tblite_gxtb\n"
        f"CP2K| Input file name {project}.inp\n"
        f"GLOBAL| Project name {project}\n"
        f"tblite source revision: {identity.tblite_source_revision}\n"
        "CP2K| Program compiled on \n"
        "CP2K| Program compiled for arm64\n"
        "SCF run converged\n"
        f"ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] {energy}\n"
        "PROGRAM ENDED\n"
    )
    shutil.copy2(job.input_path, runner.frozen_input_path(job))
    runner.write_stamp(job, identity)


def write_qualified_execution_manifest(
    root: Path,
    identity: runner.ProductionIdentity,
    tblite: Path,
    campaign_path: Path,
) -> tuple[
    Path,
    dict[str, object],
    dict[tuple[str, str], dict[str, object]],
]:
    evidence_dir = root / "qualification_evidence"
    evidence_dir.mkdir(exist_ok=True)
    canonical_input_root = root / runner.GXTB_INPUT_DIRECTORY
    benchmark.prepare_inputs(
        ["GXTB"], canonical_input_root, ["k666"], ["VII", "Ih"]
    )
    build_identity = {
        **runner.execution_build_identity(identity),
        "tblite_cli_sha256": runner.sha256(tblite),
    }
    campaign = json.loads(campaign_path.read_text())
    frozen_identity = runner.frozen_build_identity_from_manifest(campaign)
    remote_build_id = runner.build_id(build_identity)
    reference_build_id = runner.build_id(frozen_identity)
    outputs: dict[str, Path] = {}
    stamps: dict[str, Path] = {}
    evidence_inputs: dict[str, Path] = {}
    reference_records: dict[tuple[str, str], dict[str, object]] = {}
    for system, phase, energy in (("phase", "VII", -0.99), ("ih", "Ih", -1.0)):
        project = f"ice_{phase}_GXTB_k666"
        input_name = f"{project}.inp"
        output_name = f"{project}.out"
        canonical_input = canonical_input_root / "k666" / input_name

        remote_dir = evidence_dir / "remote" / "k666" / phase
        remote_dir.mkdir(parents=True, exist_ok=True)
        remote_input = remote_dir / input_name
        shutil.copy2(canonical_input, remote_input)
        remote_output = remote_dir / output_name
        remote_output.write_text(
            "PROGRAM STARTED AT 2026-07-14 13:07:01.000\n"
            "PROGRAM STARTED ON terok\n"
            f"CP2K| source code revision number: {identity.cp2k_source_revision[:10]}\n"
            "CP2K| cp2kflags: tblite tblite_gxtb\n"
            f"CP2K| Input file name {input_name}\n"
            f"GLOBAL| Project name {project}\n"
            f"tblite source revision: {identity.tblite_source_revision}\n"
            "CP2K| Program compiled on terok\n"
            "CP2K| Program compiled for x86_64\n"
            "SCF run converged\n"
            f"ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] {energy}\n"
            "PROGRAM ENDED\n"
        )
        remote_stamp = remote_output.with_suffix(".run.json")
        remote_stamp.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "campaign_id": identity.campaign_id,
                    "method": "GXTB",
                    "mesh": "k666",
                    "phase": phase,
                    "input": input_name,
                    "input_sha256": runner.sha256(remote_input),
                    "frozen_input": input_name,
                    "frozen_input_sha256": runner.sha256(remote_input),
                    "output": output_name,
                    "output_sha256": runner.sha256(remote_output),
                    "gxtb_protocol_id": runner.GXTB_PROTOCOL_ID,
                    "input_contract_valid": True,
                    "adopted_existing_output": False,
                    "build_id": remote_build_id,
                    **{
                        field: build_identity[field]
                        for field in runner.BUILD_IDENTITY_FIELDS
                    },
                },
                sort_keys=True,
            )
        )

        reference_dir = root / runner.GXTB_RUN_DIRECTORY / "k666" / phase
        reference_dir.mkdir(parents=True, exist_ok=True)
        reference_input = reference_dir / input_name
        shutil.copy2(canonical_input, reference_input)
        reference_output = reference_dir / output_name
        if not reference_output.exists():
            reference_output.write_text(
                "PROGRAM STARTED AT 2026-07-14 09:56:29.000\n"
                "PROGRAM STARTED ON Mac.test\n"
                f"CP2K| source code revision number: {str(frozen_identity['cp2k_source_revision'])[:10]}\n"
                "CP2K| cp2kflags: tblite tblite_gxtb\n"
                f"CP2K| Input file name {input_name}\n"
                f"GLOBAL| Project name {project}\n"
                "tblite source revision: unknown\n"
                "CP2K| Program compiled on \n"
                "CP2K| Program compiled for arm64\n"
                "SCF run converged\n"
                f"ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] {energy}\n"
                "PROGRAM ENDED\n"
            )
        reference_stamp = reference_output.with_suffix(".run.json")
        if not reference_stamp.exists():
            reference_stamp.write_text(
                json.dumps(
                    {
                        "campaign_id": identity.campaign_id,
                        "method": "GXTB",
                        "mesh": "k666",
                        "phase": phase,
                        "input": input_name,
                        "input_sha256": runner.sha256(canonical_input),
                        "output": output_name,
                        "output_sha256": runner.sha256(reference_output),
                        "gxtb_protocol_id": runner.GXTB_PROTOCOL_ID,
                        "input_contract_valid": True,
                        "adopted_existing_output": False,
                        **{
                            field: frozen_identity[field]
                            for field in runner.BUILD_IDENTITY_FIELDS
                        },
                    },
                    sort_keys=True,
                )
            )
        reference_records[("k666", phase)] = {
            "campaign_id": identity.campaign_id,
            "gxtb_protocol_id": runner.GXTB_PROTOCOL_ID,
            "build_id": reference_build_id,
            "mesh": "k666",
            "phase": phase,
            "input": str(canonical_input.relative_to(root)),
            "input_sha256": runner.sha256(canonical_input),
            "frozen_input": str(reference_input.relative_to(root)),
            "frozen_input_sha256": runner.sha256(reference_input),
            "output": str(reference_output.relative_to(root)),
            "output_sha256": runner.sha256(reference_output),
            "stamp": str(reference_stamp.relative_to(root)),
            "stamp_sha256": runner.sha256(reference_stamp),
        }
        evidence_inputs[system] = remote_input
        outputs[f"remote_{system}_output"] = remote_output
        outputs[f"reference_{system}_output"] = reference_output
        stamps[f"remote_{system}_stamp"] = remote_stamp
        stamps[f"reference_{system}_stamp"] = reference_stamp

    phase_input = evidence_inputs["phase"]
    ih_input = evidence_inputs["ih"]
    relative_energy = 0.01 / 12 * runner.HARTREE_TO_KJMOL
    sentinel = {
        "kind": "same_mesh_dense_relative_energy",
        "mesh": "k666",
        "phase": "VII",
        "remote_build_id": remote_build_id,
        "reference_build_id": reference_build_id,
        "phase_input": str(phase_input.relative_to(root)),
        "phase_input_sha256": runner.sha256(phase_input),
        "ih_input": str(ih_input.relative_to(root)),
        "ih_input_sha256": runner.sha256(ih_input),
        **{
            label: str(output.relative_to(root))
            for label, output in outputs.items()
        },
        **{
            f"{label}_sha256": runner.sha256(output)
            for label, output in outputs.items()
        },
        **{
            label: str(stamp.relative_to(root))
            for label, stamp in stamps.items()
        },
        **{
            f"{label}_sha256": runner.sha256(stamp)
            for label, stamp in stamps.items()
        },
        "phase_water_count": 12,
        "ih_water_count": 12,
        "hartree_to_kjmol": runner.HARTREE_TO_KJMOL,
        "phase_total_energy_delta_hartree": 0.0,
        "ih_total_energy_delta_hartree": 0.0,
        "remote_relative_energy_kjmol_per_h2o": relative_energy,
        "reference_relative_energy_kjmol_per_h2o": relative_energy,
        "relative_energy_delta_kjmol_per_h2o": 0.0,
    }
    payload = {
        "schema_version": 1,
        "campaign_id": identity.campaign_id,
        "campaign_manifest_sha256": runner.sha256(campaign_path),
        "gxtb_protocol_id": runner.GXTB_PROTOCOL_ID,
        "build_id": runner.build_id(build_identity),
        "build_identity": build_identity,
        "qualification": {
            "status": "passed",
            "evidence_schema_version": 3,
            "remote_execution_environment": {
                "program_started_on": "terok",
                "program_compiled_on": "terok",
                "program_compiled_for": "x86_64",
            },
            "total_energy_tolerance_hartree": 1.0e-10,
            "relative_energy_tolerance_kjmol_per_h2o": 1.0e-3,
            "observed_max_abs_total_energy_delta_hartree": 0.0,
            "observed_max_abs_relative_energy_delta_kjmol_per_h2o": 0.0,
            "same_mesh_dense_relative_sentinel_count": 1,
            "same_mesh_dense_relative_sentinels": [sentinel],
        },
    }
    manifest = root / "execution-build-qualified.json"
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return manifest, payload, reference_records


class DMC13GXTBInputTests(unittest.TestCase):
    def test_analysis_and_runner_share_the_production_protocol_id(self) -> None:
        self.assertEqual(benchmark.GXTB_PROTOCOL_ID, runner.GXTB_PROTOCOL_ID)

    def test_tblite_revision_metadata_accepts_portable_strings_layouts(self) -> None:
        revision = "1" * 40
        self.assertEqual(
            runner.embedded_tblite_revision(
                f"metadata\ntblite source revision: {revision}\ntrailer\n"
            ),
            revision,
        )
        self.assertEqual(
            runner.embedded_tblite_revision(
                "metadata\ntblite source revision: \n(16x,a)\ng-xTB\n"
                + "pooled-data\n" * 100
                + revision
                + "\ntrailer\n",
                revision,
            ),
            revision,
        )
        with self.assertRaisesRegex(ValueError, "cannot find"):
            runner.embedded_tblite_revision(
                "tblite source revision:\nunrelated\n" + "2" * 40 + "\n",
                revision,
            )
        with self.assertRaisesRegex(ValueError, "cannot find"):
            runner.embedded_tblite_revision(
                "unlabelled pooled data\n" + revision + "\n",
                revision,
            )
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            runner.embedded_tblite_revision(
                "tblite source revision: "
                + "1" * 40
                + "\ntblite source revision: "
                + "2" * 40
                + "\n"
            )

    def test_dense_extensions_do_not_change_the_frozen_core_default(self) -> None:
        self.assertEqual(
            runner.MESHES,
            ["gamma", "k111", "k222", "k333", "k444", "k555"],
        )
        self.assertEqual(
            runner.DENSE_EXTENSION_MESHES,
            [
                "k666",
                "k777",
                "k888",
                "k999",
                "k101010",
                "k111111",
                "k121212",
                "k131313",
            ],
        )
        self.assertEqual(
            len(runner.jobs(Path("/tmp/dmc13-core"), ["GXTB"])),
            78,
        )
        self.assertEqual(
            [mesh["id"] for mesh in benchmark.MESHES], runner.MESHES
        )
        self.assertEqual(
            [mesh["id"] for mesh in benchmark.DENSE_EXTENSION_MESHES],
            runner.DENSE_EXTENSION_MESHES,
        )
        for mesh in benchmark.SUPPORTED_MESHES:
            if mesh["id"] == "gamma":
                continue
            self.assertEqual(runner.MESH_SCHEMES[mesh["id"]], mesh["scheme"])

    def test_prepare_can_select_only_dense_extensions_with_parity_shifts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "inputs"
            benchmark.prepare_inputs(
                ["GXTB"],
                output,
                runner.DENSE_EXTENSION_MESHES,
            )
            self.assertFalse((output / "k555").exists())
            paths = {
                mesh: output / mesh / f"ice_Ih_GXTB_{mesh}.inp"
                for mesh in runner.DENSE_EXTENSION_MESHES
            }
            for mesh, path in paths.items():
                self.assertIn(
                    f"SCHEME {runner.MESH_SCHEMES[mesh]}", path.read_text()
                )
                job = runner.Job(mesh, "GXTB", "Ih", path, output / "run", "x.out")
                self.assertEqual(runner.gxtb_input_contract_errors(job), [])

    def test_prefixed_mae_plot_extends_through_latest_dense_mesh(self) -> None:
        rows: list[dict[str, object]] = []
        for mesh in benchmark.MESHES:
            for method in benchmark.METHODS:
                rows.append(
                    {
                        "mesh": mesh["id"],
                        "method": benchmark.method_label(method),
                        "MAE": "1.0",
                    }
                )
        for mesh in benchmark.DENSE_EXTENSION_MESHES:
            rows.append(
                {
                    "mesh": mesh["id"],
                    "method": benchmark.method_label("GXTB"),
                    "MAE": "2.0",
                }
            )

        scripts: list[str] = []

        def fake_run(command, **kwargs):
            if command == ["gnuplot"]:
                script = kwargs["input"].decode()
                scripts.append(script)
                output_line = next(
                    line
                    for line in script.splitlines()
                    if line.startswith("set output ")
                )
                svg = Path(output_line.removeprefix("set output ").strip("'"))
                svg.write_text("<svg/>\n")
            return mock.Mock(returncode=0)

        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(
            benchmark, "DATA", Path(tmp) / "data"
        ), mock.patch.object(
            benchmark, "FIGURES", Path(tmp) / "figures"
        ), mock.patch.object(
            benchmark.shutil, "which", return_value="/test/tool"
        ), mock.patch.object(
            benchmark.subprocess, "run", side_effect=fake_run
        ):
            benchmark.DATA.mkdir()
            benchmark.make_plots(rows, output_prefix="dense_test")
            data = (
                benchmark.DATA
                / "dmc_ice13_dense_test_kpoint_stats_for_plot.dat"
            ).read_text()

        self.assertIn('10 "9x9x9" NaN NaN 2.0', data)
        self.assertIn("set xrange [0.75:14.25]", scripts[0])
        self.assertIn("'9x9x9' 10", scripts[0])

    def test_duplicate_job_axes_are_rejected_and_paths_are_unique(self) -> None:
        root = Path("/tmp/dmc13-test")
        cases = (
            ((["GXTB", "GXTB"], ["k666"], ["Ih"]), "method"),
            ((["GXTB"], ["k666", "k666"], ["Ih"]), "mesh"),
            ((["GXTB"], ["k666"], ["Ih", "Ih"]), "phase"),
        )
        for selections, label in cases:
            with self.subTest(label=label), self.assertRaisesRegex(ValueError, label):
                runner.jobs(root, *selections)
        selected = runner.jobs(
            root,
            ["GXTB"],
            ["k666", "k777", "k888", "k999"],
            ["Ih", "II"],
            gxtb_input_root=root / "inputs",
            gxtb_run_root=root / "runs",
        )
        targets = [(job.run_dir, job.output_name) for job in selected]
        self.assertEqual(len(targets), len(set(targets)))

    def test_duplicate_cli_axes_fail_before_lock_or_campaign_reads(self) -> None:
        base = [
            "run_dmc13_kpoint_jobs.py",
            "--root",
            "/tmp/dmc13-cli-test",
            "--cp2k",
            "/tmp/cp2k",
            "--tblite",
            "/tmp/save/bin/tblite",
            "--tblite-static-library",
            "/tmp/save/lib/libtblite.a",
            "--cp2k-source",
            "/tmp/cp2k-source",
            "--tblite-source",
            "/tmp/save-source",
        ]
        duplicates = (
            ["--method", "GXTB", "--method", "GXTB"],
            ["--method", "GXTB", "--mesh", "k666", "--mesh", "k666"],
            ["--method", "GXTB", "--phase", "Ih", "--phase", "Ih"],
        )
        for selection in duplicates:
            with (
                self.subTest(selection=selection),
                mock.patch.object(sys, "argv", [*base, *selection]),
                mock.patch.object(runner, "read_campaign_manifest") as manifest,
                mock.patch.object(runner, "acquire_runner_lock") as lock,
                mock.patch.object(
                    runner.argparse.ArgumentParser,
                    "error",
                    side_effect=ValueError("CLI rejected"),
                ),
                self.assertRaisesRegex(ValueError, "CLI rejected"),
            ):
                runner.main()
            manifest.assert_not_called()
            lock.assert_not_called()

    def test_alternate_main_validates_pinned_base_references_before_work(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("remote build\n")
            identity = fake_production_identity(root, cp2k)
            tblite_prefix = root / "install"
            tblite = tblite_prefix / "bin" / "tblite"
            tblite_archive = tblite_prefix / "lib" / "libtblite.a"
            tblite.parent.mkdir(parents=True)
            tblite_archive.parent.mkdir(parents=True)
            tblite.write_text("cli\n")
            tblite_archive.write_text("archive\n")
            campaign_path = root / "campaign.json"
            campaign = write_fake_campaign_manifest(campaign_path, identity)
            execution_path = root / "execution.json"
            execution_path.write_text("{}\n")
            base_path = root / "base.json"
            base_path.write_text("{}\n")
            base_digest = "a" * 64
            frozen_id = runner.build_id(
                runner.frozen_build_identity_from_manifest(campaign)
            )
            base_payload = {
                "records": [
                    {"mesh": "k666", "phase": phase, "build_id": frozen_id}
                    for phase in ("Ih", "VII")
                ]
            }
            argv = [
                "run_dmc13_kpoint_jobs.py",
                "--campaign-manifest", str(campaign_path),
                "--execution-build-manifest", str(execution_path),
                "--base-validation-index", str(base_path),
                "--base-validation-index-sha256", base_digest,
                "--root", str(root),
                "--cp2k", str(cp2k),
                "--tblite", str(tblite),
                "--tblite-static-library", str(tblite_archive),
                "--cp2k-source", str(root / "cp2k-source"),
                "--tblite-source", str(root / "tblite-source"),
                "--method", "GXTB",
                "--mesh", "k666",
                "--phase", "Ih",
            ]

            def reject_execution_manifest(*args):
                references = args[-1]
                self.assertEqual(set(references), {("k666", "Ih"), ("k666", "VII")})
                raise ValueError("sentinel rejected")

            def digest(path):
                return base_digest if Path(path).resolve() == base_path.resolve() else "b" * 64

            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(runner, "acquire_runner_lock", return_value=object()),
                mock.patch.object(runner, "release_runner_lock"),
                mock.patch.object(runner.atexit, "register"),
                mock.patch.object(runner, "read_campaign_manifest", return_value=campaign),
                mock.patch.object(runner, "production_identity", return_value=identity),
                mock.patch.object(runner, "read_validation_index", return_value=base_payload),
                mock.patch.object(runner, "sha256", side_effect=digest),
                mock.patch.object(
                    runner,
                    "validate_execution_build_manifest",
                    side_effect=reject_execution_manifest,
                ) as execution_gate,
                mock.patch.object(runner.subprocess, "run") as process_run,
                mock.patch.object(runner.subprocess, "Popen") as process_open,
                mock.patch.object(
                    runner.argparse.ArgumentParser,
                    "error",
                    side_effect=ValueError("CLI rejected"),
                ),
                self.assertRaisesRegex(ValueError, "CLI rejected"),
            ):
                runner.main()
            execution_gate.assert_called_once()
            process_run.assert_not_called()
            process_open.assert_not_called()

    def test_runner_lock_is_nonblocking_and_released(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lock_path = Path(tmp) / ".dmc13-runner.lock"
            first_owner = {"pid": 1, "argv": ["first"]}
            first = runner.acquire_runner_lock(lock_path, first_owner)
            first_bytes = lock_path.read_bytes()
            self.assertEqual(json.loads(first_bytes), first_owner)
            with self.assertRaisesRegex(ValueError, "already active"):
                runner.acquire_runner_lock(lock_path, {"pid": 2})
            self.assertEqual(lock_path.read_bytes(), first_bytes)
            runner.release_runner_lock(first)

            second_owner = {"pid": 2, "argv": ["second"]}
            second = runner.acquire_runner_lock(lock_path, second_owner)
            self.assertEqual(json.loads(lock_path.read_text()), second_owner)
            runner.release_runner_lock(second)

    def test_native_mixer_has_explicit_iteration_limit(self) -> None:
        source = (REPOSITORY / "DMC-ICE13" / "inputs" / "ice_Ih_GFN2.inp").read_text()
        text = benchmark.gxtb_from_gfn2_template(source, "Ih")
        self.assertIn(
            f"# DMC13_GXTB_PROTOCOL {benchmark.GXTB_PROTOCOL_ID}", text
        )
        self.assertIn("METHOD GXTB", text)
        self.assertIn("SCC_MIXER TBLITE", text)
        self.assertIn("&TBLITE_MIXER\n          ITERATIONS 300", text)

    def test_non_gamma_mesh_uses_spglib_reduction_for_every_method(self) -> None:
        mesh = next(item for item in benchmark.MESHES if item["id"] == "k333")
        base = "    &QS\n    &END QS\n"
        text = benchmark.insert_kpoints(base, mesh)
        self.assertIn("SCHEME MACDONALD 3 3 3 0.0 0.0 0.0", text)
        self.assertIn("FULL_GRID F", text)
        self.assertIn("SYMMETRY T", text)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", text)
        self.assertIn("SYMMETRY_REDUCTION_METHOD SPGLIB", text)
        self.assertNotIn("FULL_GRID T", text)
        self.assertNotIn("SYMMETRY F", text)

    def test_implicit_gamma_and_explicit_k111_remain_distinct(self) -> None:
        base = "    &QS\n    &END QS\n"
        gamma = next(item for item in benchmark.MESHES if item["id"] == "gamma")
        k111 = next(item for item in benchmark.MESHES if item["id"] == "k111")
        implicit = benchmark.insert_kpoints(base, gamma, implicit_gamma=True)
        historical = benchmark.insert_kpoints(base, gamma)
        explicit = benchmark.insert_kpoints(base, k111)
        self.assertNotIn("&KPOINTS", implicit)
        self.assertIn("SCHEME GAMMA", historical)
        self.assertIn("SCHEME MACDONALD 1 1 1 0.0 0.0 0.0", explicit)
        self.assertIn("SYMMETRY_BACKEND SPGLIB", explicit)

    def test_production_contract_rejects_legacy_full_grid_input(self) -> None:
        common = f"""# DMC13_GXTB_PROTOCOL {runner.GXTB_PROTOCOL_ID}
METHOD XTB
METHOD GXTB
ACCURACY 0.1
SCC_MIXER TBLITE
ITERATIONS 300
EPS_SCF 1.0E-9
METHOD DIRECT_P_MIXING
ALPHA 0.2
CANONICALIZE TRUE
&KPOINTS
SCHEME MACDONALD 3 3 3 0.0 0.0 0.0
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inp = root / "ice_Ih_GXTB_k333.inp"
            inp.write_text(common + "SYMMETRY F\nFULL_GRID T\n&END KPOINTS\n")
            job = runner.Job("k333", "GXTB", "Ih", inp, root / "run", "result.out")
            errors = runner.gxtb_input_contract_errors(job)
            self.assertIn("missing SYMMETRY T", errors)
            self.assertIn("missing FULL_GRID F", errors)
            self.assertIn("forbidden legacy setting SYMMETRY F", errors)
            self.assertIn("forbidden legacy setting FULL_GRID T", errors)
            self.assertTrue(runner.is_legacy_full_grid_input(inp))

            inp.write_text(
                common
                + "SYMMETRY T\nFULL_GRID F\n"
                + "SYMMETRY_BACKEND SPGLIB\n"
                + "SYMMETRY_REDUCTION_METHOD SPGLIB\n&END KPOINTS\n"
            )
            self.assertEqual(runner.gxtb_input_contract_errors(job), [])
            self.assertFalse(runner.is_legacy_full_grid_input(inp))

            inp.write_text(inp.read_text() + "METHOD GFN2\n")
            self.assertIn(
                "conflicting tblite method METHOD GFN2",
                runner.gxtb_input_contract_errors(job),
            )
            self.assertIn(
                "conflicting tblite method METHOD GFN2",
                benchmark._gxtb_evidence_input_contract_errors(
                    inp.read_text(), "k333"
                ),
            )

            inp.write_text(inp.read_text().replace("METHOD GFN2\n", "METHOD GXTB\n"))
            self.assertIn(
                "duplicate critical setting METHOD GXTB",
                runner.gxtb_input_contract_errors(job),
            )

            inp.write_text(inp.read_text().replace("METHOD XTB\n", ""))
            self.assertIn(
                "missing METHOD XTB", runner.gxtb_input_contract_errors(job)
            )


class DMC13RunnerSafetyTests(unittest.TestCase):
    def test_job_selection_can_limit_a_dense_mesh_to_pilot_phases(self) -> None:
        selected = runner.jobs(
            REPOSITORY / "DMC-ICE13",
            ["GXTB"],
            ["k555"],
            ["Ih", "VII", "XVII"],
        )
        self.assertEqual(len(selected), 3)
        self.assertEqual([job.phase for job in selected], ["Ih", "VII", "XVII"])
        self.assertTrue(all(job.mesh == "k555" and job.method == "GXTB" for job in selected))

        dense = runner.jobs(
            REPOSITORY / "DMC-ICE13",
            ["GXTB"],
            ["k666"],
            ["Ih", "VII", "XIV", "XI", "VIII", "VI", "XV", "XVII"],
        )
        self.assertEqual(len(dense), 8)
        self.assertTrue(all(job.mesh == "k666" for job in dense))

    def test_production_roots_are_separate_from_legacy_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / "gxtb_spglib_inputs"
            run_root = root / "runs_gxtb_spglib"
            selected = runner.jobs(
                root,
                ["GXTB"],
                ["gamma", "k333"],
                ["Ih"],
                gxtb_input_root=input_root,
                gxtb_run_root=run_root,
            )
            self.assertEqual(selected[0].input_path, input_root / "gamma" / "ice_Ih_GXTB_gamma.inp")
            self.assertEqual(selected[0].run_dir, run_root / "gamma" / "Ih")
            self.assertEqual(selected[0].output_name, "ice_Ih_GXTB_gamma.out")
            self.assertEqual(selected[1].run_dir, run_root / "k333" / "Ih")

    def test_invalid_output_is_archived_before_clean_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_path = root / "source.inp"
            input_path.write_text("&GLOBAL\n&END GLOBAL\n")
            run_dir = root / "run"
            run_dir.mkdir()
            output = run_dir / "result.out"
            output.write_text("PROGRAM ENDED\nstale output without an energy\n")
            (run_dir / "run.log").write_text("old diagnostic\n")
            cp2k = root / "fake-cp2k"
            cp2k.write_text(
                "#!/bin/sh\n"
                "while [ \"$#\" -gt 0 ]; do\n"
                "  if [ \"$1\" = \"-o\" ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf '%s\\n' 'SCF run converged' "
                "'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] -1.0' "
                "'PROGRAM ENDED' > \"$out\"\n"
            )
            cp2k.chmod(cp2k.stat().st_mode | 0o111)
            job = runner.Job("gamma", "GFN1", "Ih", input_path, run_dir, output.name)
            identity = fake_production_identity(root, cp2k)

            _, returncode = runner.run_job(
                identity,
                job,
                False,
                threading.Event(),
            )

            self.assertEqual(returncode, 0)
            self.assertIn("-1.0", output.read_text())
            self.assertEqual(len(list(run_dir.glob("result.out.stale-*"))), 1)
            self.assertEqual(len(list(run_dir.glob("run.log.stale-*"))), 1)
            self.assertTrue(runner.stamp_path(job).is_file())
            self.assertTrue(runner.stamp_valid(job, identity))

            source_bytes = input_path.read_bytes()
            frozen = runner.frozen_input_path(job)
            frozen_bytes = frozen.read_bytes()
            stamp = json.loads(runner.stamp_path(job).read_text())
            self.assertEqual(stamp["frozen_input"], frozen.name)
            self.assertEqual(stamp["input_sha256"], stamp["frozen_input_sha256"])

            stamp.pop("frozen_input")
            stamp.pop("frozen_input_sha256")
            runner.stamp_path(job).write_text(json.dumps(stamp))
            self.assertTrue(runner.stamp_valid(job, identity))

            for field in (
                "mpi_ranks_per_job",
                "threads_per_job",
                "omp_num_threads",
                "omp_schedule",
                "omp_dynamic",
                "omp_wait_policy",
                "omp_proc_bind",
                "omp_places",
                "blas_threads",
            ):
                stamp.pop(field)
            runner.stamp_path(job).write_text(json.dumps(stamp))
            self.assertTrue(runner.stamp_valid(job, identity))

            input_path.write_bytes(source_bytes + b"central mutation\n")
            self.assertFalse(runner.stamp_valid(job, identity))
            input_path.write_bytes(source_bytes)
            self.assertTrue(runner.stamp_valid(job, identity))
            frozen.write_bytes(frozen_bytes + b"local mutation\n")
            self.assertFalse(runner.stamp_valid(job, identity))
            frozen.write_bytes(frozen_bytes)
            self.assertTrue(runner.stamp_valid(job, identity))

            identity.cp2k_library.write_text("replacement CP2K library\n")
            replaced_cp2k = replace(
                identity,
                cp2k_library_sha256=runner.sha256(identity.cp2k_library),
            )
            self.assertFalse(runner.stamp_valid(job, replaced_cp2k))

            identity.tblite_static_library.write_text(
                "different static save_tblite archive\n"
            )
            mixed_tblite = replace(
                identity,
                tblite_static_library_sha256=runner.sha256(
                    identity.tblite_static_library
                ),
            )
            self.assertFalse(runner.stamp_valid(job, mixed_tblite))

    def test_threaded_job_records_and_exports_bounded_threading(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_path = root / "source.inp"
            input_path.write_text("&GLOBAL\n&END GLOBAL\n")
            run_dir = root / "run"
            cp2k = root / "fake-cp2k"
            cp2k.write_text(
                "#!/bin/sh\n"
                "printf '%s\\n' \"$OMP_NUM_THREADS\" \"$OMP_SCHEDULE\" "
                "\"$OMP_DYNAMIC\" \"$OMP_WAIT_POLICY\" "
                "\"$OPENBLAS_NUM_THREADS\" \"$MKL_NUM_THREADS\" "
                "\"$VECLIB_MAXIMUM_THREADS\" "
                "\"${OMP_PROC_BIND-unset}\" \"${OMP_PLACES-unset}\" "
                "> thread-env.txt\n"
                "while [ \"$#\" -gt 0 ]; do\n"
                "  if [ \"$1\" = \"-o\" ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf '%s\\n' 'SCF run converged' "
                "'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] -1.0' "
                "'PROGRAM ENDED' > \"$out\"\n"
            )
            cp2k.chmod(cp2k.stat().st_mode | 0o111)
            job = runner.Job("gamma", "GFN1", "Ih", input_path, run_dir, "result.out")
            identity = fake_production_identity(root, cp2k)

            with mock.patch.dict(
                runner.os.environ,
                {"OMP_PROC_BIND": "true", "OMP_PLACES": "cores"},
                clear=True,
            ):
                _, returncode = runner.run_job(
                    identity,
                    job,
                    False,
                    threading.Event(),
                    threads_per_job=3,
                )

            self.assertEqual(returncode, 0)
            self.assertEqual(
                (run_dir / "thread-env.txt").read_text().splitlines(),
                ["3", "static", "FALSE", "PASSIVE", "1", "1", "1", "unset", "unset"],
            )
            stamp = json.loads(runner.stamp_path(job).read_text())
            self.assertEqual(stamp["mpi_ranks_per_job"], 1)
            self.assertEqual(stamp["threads_per_job"], 3)
            self.assertEqual(stamp["omp_num_threads"], 3)
            self.assertEqual(stamp["omp_schedule"], "static")
            self.assertIs(stamp["omp_dynamic"], False)
            self.assertEqual(stamp["omp_wait_policy"], "PASSIVE")
            self.assertIsNone(stamp["omp_proc_bind"])
            self.assertIsNone(stamp["omp_places"])
            self.assertEqual(stamp["blas_threads"], 1)
            self.assertTrue(runner.stamp_valid(job, identity))

    def test_positive_int_rejects_nonpositive_thread_counts(self) -> None:
        self.assertEqual(runner.positive_int("4"), 4)
        for value in ("0", "-1", "not-an-integer"):
            with self.subTest(value=value):
                with self.assertRaises(runner.argparse.ArgumentTypeError):
                    runner.positive_int(value)

    def test_execution_parallelism_records_nominal_core_budget(self) -> None:
        parallelism = runner.execution_parallelism(5, 3)
        self.assertEqual(parallelism["jobs"], 5)
        self.assertEqual(parallelism["threads_per_job"], 3)
        self.assertEqual(parallelism["nominal_cores"], 15)
        self.assertEqual(parallelism["mpi_ranks_per_job"], 1)
        self.assertEqual(parallelism["blas_threads"], 1)

    def test_write_stamp_refuses_source_frozen_input_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.inp"
            source.write_text("source\n")
            run_dir = root / "run"
            run_dir.mkdir()
            job = runner.Job("gamma", "GFN1", "Ih", source, run_dir, "result.out")
            (run_dir / job.output_name).write_text(
                "SCF run converged\nENERGY| Total FORCE_EVAL -1.0\nPROGRAM ENDED\n"
            )
            runner.frozen_input_path(job).write_text("different\n")
            cp2k = root / "cp2k"
            cp2k.write_text("binary\n")
            identity = fake_production_identity(root, cp2k)
            with self.assertRaisesRegex(ValueError, "source/frozen input mismatch"):
                runner.write_stamp(job, identity)
            self.assertFalse(runner.stamp_path(job).exists())

    def test_pre_cancelled_job_does_not_touch_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_path = root / "source.inp"
            input_path.write_text("input\n")
            run_dir = root / "run"
            run_dir.mkdir()
            output = run_dir / "result.out"
            output.write_text("untouched\n")
            cp2k = root / "unused"
            cp2k.write_text("unused\n")
            job = runner.Job("k333", "GXTB", "III", input_path, run_dir, output.name)
            identity = fake_production_identity(root, cp2k)
            stop_event = threading.Event()
            stop_event.set()

            _, returncode = runner.run_job(
                identity,
                job,
                False,
                stop_event,
            )

            self.assertEqual(returncode, 130)
            self.assertEqual(output.read_text(), "untouched\n")
            self.assertFalse(list(run_dir.glob("*.stale-*")))

    def test_linux_libcp2k_resolution_is_unique_and_resolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("launcher\n")
            library = root / "libcp2k.so.2026.1"
            library.write_text("library\n")
            completed = runner.subprocess.CompletedProcess(
                ["ldd", str(cp2k)],
                0,
                stdout=(
                    "linux-vdso.so.1 (0x0000)\n"
                    f"libcp2k.so.2026.1 => {library} (0x0001)\n"
                    "libm.so.6 => /lib/libm.so.6 (0x0002)\n"
                ),
                stderr="",
            )
            with mock.patch.object(runner.subprocess, "run", return_value=completed):
                self.assertEqual(
                    runner.resolve_cp2k_library_linux(cp2k), library.resolve()
                )
                with mock.patch.object(runner.sys, "platform", "linux"):
                    self.assertEqual(
                        runner.resolve_cp2k_library(cp2k), library.resolve()
                    )

            unresolved = runner.subprocess.CompletedProcess(
                ["ldd", str(cp2k)],
                0,
                stdout="libcp2k.so.2026.1 => not found\n",
                stderr="",
            )
            with mock.patch.object(runner.subprocess, "run", return_value=unresolved):
                with self.assertRaisesRegex(ValueError, "could not resolve"):
                    runner.resolve_cp2k_library_linux(cp2k)

            direct = runner.subprocess.CompletedProcess(
                ["ldd", str(cp2k)],
                0,
                stdout=f"{library} (0x0001)\n",
                stderr="",
            )
            with mock.patch.object(runner.subprocess, "run", return_value=direct):
                self.assertEqual(
                    runner.resolve_cp2k_library_linux(cp2k), library.resolve()
                )

    def test_production_paths_reject_prefix_and_root_escapes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            runner.validate_production_paths(
                root,
                runner.GXTB_ANALYSIS_PREFIX,
                root / runner.GXTB_INPUT_DIRECTORY,
                root / runner.GXTB_RUN_DIRECTORY,
            )
            for prefix in ("", "../gxtb_spglib", "nested/gxtb_spglib"):
                with self.subTest(prefix=prefix):
                    with self.assertRaisesRegex(ValueError, "analysis-prefix"):
                        runner.validate_production_paths(
                            root,
                            prefix,
                            root / runner.GXTB_INPUT_DIRECTORY,
                            root / runner.GXTB_RUN_DIRECTORY,
                        )
            with self.assertRaisesRegex(ValueError, "outside the campaign root"):
                runner.validate_production_paths(
                    root,
                    runner.GXTB_ANALYSIS_PREFIX,
                    root.parent / "outside-inputs",
                    root / runner.GXTB_RUN_DIRECTORY,
                )

    def test_qualified_build_rejects_unknown_embedded_tblite_revision(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            library = root / "libcp2k.so"
            archive = root / "libtblite.a"
            for path in (cp2k, library, archive):
                path.write_text(path.name)
            cp2k_source = root / "cp2k-source"
            tblite_source = root / "tblite-source"
            cp2k_source.mkdir()
            tblite_source.mkdir()
            cp2k_revision = "1" * 40
            tblite_revision = "2" * 40

            def output(command, cwd=None):
                if len(command) == 2 and Path(command[0]).name == "cp2k.psmp" and command[1] == "--version":
                    return f"Source code revision {cp2k_revision}"
                if command[0] == "strings":
                    return "tblite source revision: unknown"
                if command[:2] == ["git", "rev-parse"]:
                    return cp2k_revision if cwd == cp2k_source else tblite_revision
                raise AssertionError(command)

            with mock.patch.object(
                runner, "resolve_cp2k_library", return_value=library.resolve()
            ):
                with mock.patch.object(runner, "command_output", side_effect=output):
                    with self.assertRaisesRegex(ValueError, "must embed"):
                        runner.production_identity(
                            "campaign",
                            cp2k,
                            library,
                            archive,
                            cp2k_source,
                            tblite_source,
                            tblite_revision,
                            require_embedded_tblite_revision=True,
                        )

    def test_equivalent_execution_build_requires_explicit_qualification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("remote build\n")
            tblite = root / "tblite"
            tblite.write_text("remote cli\n")
            identity = fake_production_identity(root, cp2k)
            campaign_path = root / "campaign.json"
            frozen = runner.execution_build_identity(identity)
            frozen.update(
                {
                    "cp2k_sha256": "a" * 64,
                    "cp2k_library_sha256": "b" * 64,
                    "tblite_static_library_sha256": "c" * 64,
                }
            )
            campaign = write_fake_campaign_manifest(
                campaign_path, identity, frozen_identity=frozen
            )
            execution_path, execution, reference_records = write_qualified_execution_manifest(
                root, identity, tblite, campaign_path
            )
            cp2k_source = root / "cp2k-source"
            tblite_source = root / "tblite-source"
            cp2k_source.mkdir()
            tblite_source.mkdir()
            with mock.patch.object(runner, "command_output", return_value=""):
                runner.validate_execution_build_manifest(
                    identity,
                    tblite,
                    campaign_path,
                    campaign,
                    execution_path,
                    cp2k_source,
                    tblite_source,
                    root,
                    reference_records,
                )
            execution["qualification"]["status"] = "unqualified"
            execution_path.write_text(json.dumps(execution))
            with mock.patch.object(runner, "command_output", return_value=""):
                with self.assertRaisesRegex(ValueError, "qualification.status"):
                    runner.validate_execution_build_manifest(
                        identity,
                        tblite,
                        campaign_path,
                        campaign,
                        execution_path,
                        cp2k_source,
                        tblite_source,
                        root,
                        reference_records,
                    )
class DMC13GammaCLISafetyTests(unittest.TestCase):
    def test_gamma_cli_uses_the_same_production_protocol(self) -> None:
        self.assertEqual(gamma_cli.GXTB_PROTOCOL_ID, runner.GXTB_PROTOCOL_ID)
        self.assertEqual(gamma_cli.GXTB_RUN_DIRECTORY, runner.GXTB_RUN_DIRECTORY)

    def test_cp2k_reference_requires_current_binary_and_hash_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("current executable\n")
            tblite = root / "tblite"
            tblite.write_text("current CLI\n")
            identity = fake_gamma_identity(root, cp2k, tblite)
            input_path = root / "ice_Ih_GXTB_gamma.inp"
            input_path.write_text(
                f"""# DMC13_GXTB_PROTOCOL {gamma_cli.GXTB_PROTOCOL_ID}
METHOD GXTB
ACCURACY 0.1
SCC_MIXER TBLITE
ITERATIONS 300
EPS_SCF 1.0E-9
METHOD DIRECT_P_MIXING
ALPHA 0.2
CANONICALIZE TRUE
"""
            )
            output = root / "ice_Ih_GXTB_gamma.out"
            output.write_text(
                "SCF run converged\n"
                "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] -1.0\n"
                "PROGRAM ENDED\n"
            )
            gamma_cli.cp2k_stamp_path(output).write_text(
                json.dumps(
                    {
                        "campaign_id": identity.campaign_id,
                        "method": "GXTB",
                        "mesh": "gamma",
                        "phase": "Ih",
                        "gxtb_protocol_id": gamma_cli.GXTB_PROTOCOL_ID,
                        "cp2k_sha256": identity.cp2k_sha256,
                        "cp2k_library_sha256": identity.cp2k_library_sha256,
                        "tblite_static_library_sha256": (
                            identity.tblite_static_library_sha256
                        ),
                        "cp2k_source_revision": identity.cp2k_source_revision,
                        "tblite_source_revision": identity.tblite_source_revision,
                        "input_sha256": gamma_cli.sha256(input_path),
                        "output_sha256": gamma_cli.sha256(output),
                        "input_contract_valid": True,
                        "adopted_existing_output": False,
                    }
                )
            )

            self.assertEqual(
                gamma_cli.cp2k_validation_errors(
                    "Ih", input_path, output, identity
                ),
                [],
            )
            different_binary = replace(
                identity,
                cp2k_sha256="different-current-binary",
            )
            errors = gamma_cli.cp2k_validation_errors(
                "Ih", input_path, output, different_binary
            )
            self.assertIn("CP2K stamp cp2k_sha256 mismatch", errors)
            different_library = replace(
                identity,
                cp2k_library_sha256="replacement-libcp2k",
            )
            errors = gamma_cli.cp2k_validation_errors(
                "Ih", input_path, output, different_library
            )
            self.assertIn("CP2K stamp cp2k_library_sha256 mismatch", errors)
            mixed_tblite = replace(
                identity,
                tblite_static_library_sha256="mixed-libtblite",
            )
            errors = gamma_cli.cp2k_validation_errors(
                "Ih", input_path, output, mixed_tblite
            )
            self.assertIn(
                "CP2K stamp tblite_static_library_sha256 mismatch",
                errors,
            )

            stamp_path = gamma_cli.cp2k_stamp_path(output)
            stamp = json.loads(stamp_path.read_text())
            stamp["adopted_existing_output"] = True
            stamp_path.write_text(json.dumps(stamp))
            errors = gamma_cli.cp2k_validation_errors(
                "Ih", input_path, output, identity
            )
            self.assertIn(
                "CP2K output was adopted rather than produced by this executable",
                errors,
            )
            stamp["adopted_existing_output"] = False
            stamp_path.write_text(json.dumps(stamp))

            input_path.write_text(input_path.read_text() + "&KPOINTS\n&END KPOINTS\n")
            errors = gamma_cli.cp2k_validation_errors(
                "Ih", input_path, output, identity
            )
            self.assertIn(
                "Gamma production input must be implicit without &KPOINTS",
                errors,
            )
            self.assertIn("CP2K stamp input_sha256 mismatch", errors)

    def test_cli_result_requires_current_tblite_and_artifact_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            tblite = run_dir / "tblite"
            tblite.write_text("current save_tblite CLI\n")
            cp2k = run_dir / "cp2k.psmp"
            cp2k.write_text("current CP2K launcher\n")
            identity = fake_gamma_identity(run_dir, cp2k, tblite)
            poscar = run_dir / "POSCAR"
            poscar.write_text("validated structure\n")
            result = run_dir / "tblite.json"
            result.write_text(json.dumps({"energy": -1.0, "energies": [-0.4, -0.6]}))
            (run_dir / "tblite.out").write_text("completed\n")
            gamma_cli.write_cli_stamp(
                "Ih",
                run_dir,
                identity,
                [str(tblite), "run"],
            )
            poscar_digest = gamma_cli.sha256(poscar)

            self.assertTrue(
                gamma_cli.cli_stamp_valid(
                    "Ih", run_dir, poscar_digest, identity
                )
            )
            old_tblite = replace(identity, tblite_sha256="old-tblite-binary")
            self.assertFalse(
                gamma_cli.cli_stamp_valid(
                    "Ih", run_dir, poscar_digest, old_tblite
                )
            )
            result.write_text(
                json.dumps({"energy": -2.0, "energies": [-0.8, -1.2]})
            )
            self.assertFalse(
                gamma_cli.cli_stamp_valid(
                    "Ih", run_dir, poscar_digest, identity
                )
            )


class DMC13AnalysisValidationTests(unittest.TestCase):
    @staticmethod
    def synthetic_results(
        meshes: dict[str, dict[str, float]],
        coverage: dict[str, list[str]],
    ) -> dict[str, object]:
        return {
            "validated_gxtb_phases": coverage,
            "results": {
                mesh: {"GXTB": {"relative_kjmol": relative}}
                for mesh, relative in meshes.items()
            },
        }

    @staticmethod
    def offset_relative(offset: float) -> dict[str, float]:
        return {
            "Ih": 0.0,
            **{
                phase: offset
                for phase in benchmark.PHASES
                if phase != "Ih"
            },
        }

    @staticmethod
    def phasewise_results(
        method: str,
        relative_by_mesh: dict[str, dict[str, float]],
        ih_by_mesh: dict[str, float] | None = None,
    ) -> dict[str, object]:
        result_meshes: dict[str, object] = {}
        for mesh, relative in relative_by_mesh.items():
            ih = (ih_by_mesh or {}).get(mesh, -10.0)
            per_h2o = {"Ih": ih}
            stored_relative = {"Ih": 0.0}
            for phase, value in relative.items():
                per_h2o[phase] = ih + value / benchmark.HARTREE_TO_KJMOL
                stored_relative[phase] = value
            result_meshes[mesh] = {
                method: {
                    "per_h2o_hartree": per_h2o,
                    "relative_kjmol": stored_relative,
                }
            }
        return {"results": result_meshes}

    @staticmethod
    def phase_values(offset: float) -> dict[str, float]:
        return {
            phase: float(index) + offset
            for index, phase in enumerate(
                (phase for phase in benchmark.PHASES if phase != "Ih"),
                start=1,
            )
        }

    def test_unprefixed_analysis_remains_gfn1_gfn2_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / "data"
            figures = root / "figures"
            data.mkdir()
            (data / "geometries.json").write_text(
                json.dumps(
                    {phase: {"counts": {"O": 1}} for phase in benchmark.PHASES}
                )
            )
            complete = {
                "complete": True,
                "energies_hartree": {phase: -1.0 for phase in benchmark.PHASES},
                "per_h2o_hartree": {phase: -1.0 for phase in benchmark.PHASES},
                "relative_kjmol": {phase: 0.0 for phase in benchmark.PHASES},
            }
            prior = {
                mesh["id"]: {
                    method: complete
                    for method in ("GFN1", "GFN2", "GXTB")
                }
                for mesh in benchmark.MESHES
            }
            (data / "kpoint_results.json").write_text(
                json.dumps({"results": prior})
            )
            with mock.patch.object(benchmark, "ROOT", root):
                with mock.patch.object(benchmark, "DATA", data):
                    with mock.patch.object(benchmark, "FIGURES", figures):
                        with mock.patch.object(benchmark.shutil, "which", return_value=None):
                            result = benchmark.analyse()
            self.assertEqual(result["methods"], ["GFN1", "GFN2"])
            for mesh in benchmark.MESHES:
                self.assertEqual(
                    set(result["results"][mesh["id"]]), {"GFN1", "GFN2"}
                )
            written = json.loads((data / "kpoint_results.json").read_text())
            self.assertEqual(written["methods"], ["GFN1", "GFN2"])

    def test_convergence_requires_two_full_pairs_and_resists_cancellation(self) -> None:
        phases = [phase for phase in benchmark.PHASES if phase != "Ih"]
        base = {"Ih": 0.0, **{phase: float(i) for i, phase in enumerate(phases)}}
        k444 = {
            phase: value + (0.02 if i % 2 else -0.02)
            for i, (phase, value) in enumerate(base.items())
        }
        k555 = {
            phase: value + (0.01 if i % 2 else -0.01)
            for i, (phase, value) in enumerate(k444.items())
        }
        full = {mesh: list(benchmark.PHASES) for mesh in ("k333", "k444", "k555")}
        report, rows = benchmark.build_convergence_report(
            self.synthetic_results(
                {"k333": base, "k444": k444, "k555": k555}, full
            )
        )
        self.assertEqual(report["stopping_assessment"]["status"], "converged")
        self.assertEqual(
            report["stopping_assessment"]["acceptance"]["earliest_reportable_mesh"],
            "k444",
        )
        self.assertEqual(len(rows), 24)

        bad = dict(k444)
        for i, phase in enumerate(phases):
            bad[phase] += 1.0 if i % 2 else -1.0
        failed, _ = benchmark.build_convergence_report(
            self.synthetic_results(
                {"k444": k444, "k555": bad},
                {mesh: list(benchmark.PHASES) for mesh in ("k444", "k555")},
            )
        )
        pair = failed["comparisons"][0]
        self.assertAlmostEqual(pair["rms_delta_kjmol_per_h2o"], 1.0)
        self.assertFalse(pair["passed_numeric_thresholds"])

    def test_convergence_uses_only_the_trailing_contiguous_full_sequence(self) -> None:
        cases = (
            (
                "pass_pass_fail",
                [0.0, 0.01, 0.02, 0.22],
                "not_converged",
                None,
            ),
            (
                "pass_pass_pass",
                [0.0, 0.01, 0.02, 0.03],
                "converged",
                ("k444", "k555"),
            ),
            (
                "fail_pass_pass",
                [0.0, 0.20, 0.21, 0.22],
                "converged",
                ("k555", "k666"),
            ),
        )
        meshes = ("k333", "k444", "k555", "k666")
        for label, offsets, expected_status, expected_acceptance in cases:
            with self.subTest(label=label):
                report, _ = benchmark.build_convergence_report(
                    self.synthetic_results(
                        {
                            mesh: self.offset_relative(offset)
                            for mesh, offset in zip(meshes, offsets)
                        },
                        {mesh: list(benchmark.PHASES) for mesh in meshes},
                    )
                )
                assessment = report["stopping_assessment"]
                self.assertEqual(assessment["status"], expected_status)
                if expected_acceptance is None:
                    self.assertIsNone(assessment["acceptance"])
                else:
                    self.assertEqual(
                        assessment["acceptance"]["earliest_reportable_mesh"],
                        expected_acceptance[0],
                    )
                    self.assertEqual(
                        assessment["acceptance"]["validation_mesh"],
                        expected_acceptance[1],
                    )

    def test_partial_gap_resets_formal_convergence_sequence(self) -> None:
        meshes = {
            "k333": self.offset_relative(0.00),
            "k444": self.offset_relative(0.01),
            "k555": self.offset_relative(0.02),
            "k666": self.offset_relative(0.03),
            "k777": self.offset_relative(0.04),
        }
        pilot_phases = ["Ih", "II"]
        coverage = {
            "k333": list(benchmark.PHASES),
            "k444": list(benchmark.PHASES),
            "k555": pilot_phases,
            "k666": list(benchmark.PHASES),
            "k777": list(benchmark.PHASES),
        }
        report, _ = benchmark.build_convergence_report(
            self.synthetic_results(meshes, coverage)
        )
        assessment = report["stopping_assessment"]
        self.assertEqual(assessment["status"], "insufficient_coverage")
        self.assertIsNone(assessment["acceptance"])
        self.assertEqual(
            assessment["trailing_consecutive_passing_fully_covered_pairs"], 1
        )

    def test_partial_dense_pilot_is_never_formal_stopping_evidence(self) -> None:
        pilot = ["Ih", "VII", "XIV", "XI", "VIII", "VI", "XV", "XVII"]
        full_relative = {phase: 0.0 for phase in benchmark.PHASES}
        pilot_relative = {phase: 0.001 for phase in pilot}
        report, _ = benchmark.build_convergence_report(
            self.synthetic_results(
                {"k555": full_relative, "k666": pilot_relative},
                {"k555": list(benchmark.PHASES), "k666": pilot},
            )
        )
        pair = report["comparisons"][0]
        self.assertEqual(pair["coverage"], "pilot")
        self.assertIsNone(pair["passed_numeric_thresholds"])
        self.assertFalse(pair["eligible_for_stopping"])
        self.assertTrue(pair["pilot_inconclusive"])
        self.assertFalse(pair["pilot_definitively_rejects_candidate"])
        self.assertEqual(
            report["stopping_assessment"]["status"], "insufficient_coverage"
        )

    def test_later_definitive_pilot_revokes_but_inconclusive_pilot_preserves(self) -> None:
        full_meshes = {
            "k444": self.offset_relative(0.00),
            "k555": self.offset_relative(0.01),
            "k666": self.offset_relative(0.02),
        }
        full_coverage = {
            mesh: list(benchmark.PHASES) for mesh in full_meshes
        }
        for label, pilot_delta, expected_status in (
            ("definitive", 0.20, "not_converged"),
            ("inconclusive", 0.09, "converged"),
        ):
            with self.subTest(label=label):
                pilot = {"Ih": 0.0, "II": 0.02 + pilot_delta}
                report, _ = benchmark.build_convergence_report(
                    self.synthetic_results(
                        {**full_meshes, "k777": pilot},
                        {**full_coverage, "k777": ["Ih", "II"]},
                    )
                )
                assessment = report["stopping_assessment"]
                self.assertEqual(assessment["status"], expected_status)
                if label == "definitive":
                    self.assertIsNone(assessment["acceptance"])
                    self.assertEqual(
                        assessment["revoked_acceptance"]["validation_mesh"],
                        "k666",
                    )
                    self.assertEqual(
                        assessment["revoking_dense_pilot_pairs"],
                        ["k666->k777"],
                    )
                else:
                    self.assertEqual(
                        assessment["acceptance"]["validation_mesh"], "k666"
                    )
                    self.assertIsNone(assessment["revoked_acceptance"])

    def test_dense_stopping_sequence_extends_through_k999_and_respects_gaps(self) -> None:
        meshes = ("k666", "k777", "k888", "k999")
        offsets = (0.00, 0.20, 0.21, 0.22)
        values = {
            mesh: self.offset_relative(offset)
            for mesh, offset in zip(meshes, offsets)
        }
        full = {mesh: list(benchmark.PHASES) for mesh in meshes}
        report, _ = benchmark.build_convergence_report(
            self.synthetic_results(values, full)
        )
        assessment = report["stopping_assessment"]
        self.assertEqual(assessment["status"], "converged")
        self.assertEqual(
            assessment["acceptance"]["earliest_reportable_mesh"], "k888"
        )
        self.assertEqual(assessment["acceptance"]["validation_mesh"], "k999")

        failed_values = dict(values)
        failed_values["k999"] = self.offset_relative(0.42)
        failed, _ = benchmark.build_convergence_report(
            self.synthetic_results(failed_values, full)
        )
        self.assertEqual(failed["stopping_assessment"]["status"], "not_converged")
        self.assertIsNone(failed["stopping_assessment"]["acceptance"])

        gap_coverage = dict(full)
        gap_coverage["k888"] = ["Ih", "II"]
        gap, _ = benchmark.build_convergence_report(
            self.synthetic_results(values, gap_coverage)
        )
        self.assertEqual(
            gap["stopping_assessment"]["status"], "insufficient_coverage"
        )
        self.assertIsNone(gap["stopping_assessment"]["acceptance"])

    def test_dense_mesh_adjacency_and_phasewise_selection_extend_through_k13(
        self,
    ) -> None:
        self.assertEqual(
            benchmark.CONVERGENCE_MESH_IDS[-5:],
            ["k999", "k101010", "k111111", "k121212", "k131313"],
        )
        self.assertEqual(
            [
                benchmark.mesh_parity_direction(left, right)
                for left, right in zip(
                    benchmark.CONVERGENCE_MESH_IDS[-5:],
                    benchmark.CONVERGENCE_MESH_IDS[-4:],
                )
            ],
            ["odd->even", "even->odd", "odd->even", "even->odd"],
        )
        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results(
                "GXTB",
                {
                    "k999": self.phase_values(0.00),
                    "k101010": self.phase_values(0.04),
                },
            )
        )
        method = report["methods"]["GXTB"]
        self.assertTrue(method["phasewise_kpoint_converged"])
        for phase in (value for value in benchmark.PHASES if value != "Ih"):
            convergence = method["phase_convergence"][phase]
            self.assertEqual(convergence["previous_mesh"], "k999")
            self.assertEqual(convergence["smallest_required_mesh"], "k101010")
            self.assertEqual(convergence["mesh_n"], 10)
            self.assertEqual(convergence["nk_total"], 1000)

    def test_sparse_phasewise_frontier_requests_earliest_missing_adjacent_mesh(
        self,
    ) -> None:
        relative = {
            "k999": self.phase_values(0.00),
            "k101010": self.phase_values(0.20),
            "k131313": self.phase_values(0.40),
        }
        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", relative)
        )
        method = report["methods"]["GXTB"]
        self.assertFalse(method["phasewise_kpoint_converged"])
        self.assertEqual(
            method["next_required_by_mesh"],
            {
                "k111111": [
                    phase for phase in benchmark.PHASES if phase != "Ih"
                ]
            },
        )
        self.assertFalse(method["exhausted_mesh_sequence_phases"])

    def test_pilot_rejection_uses_max_or_full_set_rms_lower_bound(self) -> None:
        cases = (
            ("max", [0.11], True, False),
            ("rms_lower_bound", [0.09] * 4, True, False),
            ("inconclusive", [0.09], False, True),
        )
        nonreference = [phase for phase in benchmark.PHASES if phase != "Ih"]
        for label, deltas, rejects, inconclusive in cases:
            with self.subTest(label=label):
                pilot_phases = ["Ih", *nonreference[: len(deltas)]]
                pilot_relative = {
                    "Ih": 0.0,
                    **{
                        phase: delta
                        for phase, delta in zip(nonreference, deltas)
                    },
                }
                report, _ = benchmark.build_convergence_report(
                    self.synthetic_results(
                        {
                            "k555": self.offset_relative(0.0),
                            "k666": pilot_relative,
                        },
                        {
                            "k555": list(benchmark.PHASES),
                            "k666": pilot_phases,
                        },
                    )
                )
                pair = report["comparisons"][0]
                expected_observed_rms = math.sqrt(
                    sum(delta * delta for delta in deltas) / len(deltas)
                )
                expected_lower_bound = math.sqrt(
                    sum(delta * delta for delta in deltas) / 12
                )
                self.assertAlmostEqual(
                    pair["observed_subset_rms_delta_kjmol_per_h2o"],
                    expected_observed_rms,
                )
                self.assertAlmostEqual(
                    pair["full_set_rms_lower_bound_kjmol_per_h2o"],
                    expected_lower_bound,
                )
                self.assertEqual(
                    pair["pilot_definitively_rejects_candidate"], rejects
                )
                self.assertEqual(pair["pilot_inconclusive"], inconclusive)
                self.assertIsNone(pair["passed_numeric_thresholds"])
                self.assertFalse(pair["passed_formal_stopping_pair"])

    def test_phasewise_convergence_reports_first_pair_and_same_mesh_ih(
        self,
    ) -> None:
        relative = {
            "k111": self.phase_values(0.00),
            "k222": self.phase_values(0.04),
            "k333": self.phase_values(0.07),
        }
        results = self.phasewise_results(
            "GXTB",
            relative,
            {"k111": -10.0, "k222": -20.0, "k333": -30.0},
        )
        report, rows = (
            benchmark.build_phasewise_kpoint_convergence_report(results)
        )
        method = report["methods"]["GXTB"]
        self.assertEqual(
            report["protocol_id"], "phasewise-kpoint-convergence-v1"
        )
        self.assertEqual(
            method["result_label"], "phase-wise k-point-converged MAE"
        )
        self.assertEqual(method["status"], "phasewise_kpoint_converged")
        self.assertTrue(method["phasewise_kpoint_converged"])
        self.assertEqual(method["converged_phase_count"], 12)
        self.assertEqual(len(rows), 36)
        for phase in (phase for phase in benchmark.PHASES if phase != "Ih"):
            convergence = method["phase_convergence"][phase]
            self.assertEqual(
                (
                    convergence["previous_mesh"],
                    convergence["smallest_required_mesh"],
                    convergence["mesh_label"],
                    convergence["mesh_n"],
                    convergence["nk_total"],
                ),
                ("k111", "k222", "2x2x2", 2, 8),
            )
            self.assertAlmostEqual(
                convergence["relative_energy_kjmol_per_h2o"],
                relative["k222"][phase],
            )
            self.assertAlmostEqual(
                convergence["error_kjmol_per_h2o"],
                relative["k222"][phase]
                - (
                    benchmark.DMC_ABS_KJMOL[phase]
                    - benchmark.DMC_ABS_KJMOL["Ih"]
                ),
            )
            self.assertAlmostEqual(
                convergence["last_delta_kjmol_per_h2o"], 0.04
            )
        diagnostic = method["last_delta_statistics_diagnostic"]
        self.assertTrue(diagnostic["diagnostic_only"])
        self.assertAlmostEqual(
            diagnostic["rms_delta_kjmol_per_h2o"], 0.04
        )
        self.assertIsNotNone(
            method["phasewise_kpoint_converged_stats_nonreference"]
        )
        self.assertIsNotNone(
            method["previous_mesh_comparison_stats_nonreference"]
        )
        bound = method["mae_mesh_difference_bound"]
        self.assertAlmostEqual(
            bound[
                "provable_absolute_mae_difference_upper_bound_kjmol_per_h2o"
            ],
            0.04,
        )
        self.assertTrue(bound["bound_satisfied"])

    def test_phasewise_partial_coverage_never_emits_a_full_mae(self) -> None:
        relative = {
            "k111": self.phase_values(0.00),
            "k222": self.phase_values(0.04),
        }
        del relative["k222"]["XVII"]
        report, rows = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", relative)
        )
        method = report["methods"]["GXTB"]
        self.assertEqual(method["status"], "unresolved_phases")
        self.assertFalse(method["phasewise_kpoint_converged"])
        self.assertEqual(method["unresolved_phases"], ["XVII"])
        self.assertEqual(method["next_required_phases"], ["XVII"])
        self.assertEqual(method["next_required_by_mesh"], {"k222": ["XVII"]})
        self.assertEqual(
            method["required_same_mesh_reference_by_mesh"],
            {"k222": "Ih"},
        )
        self.assertIsNone(
            method["phasewise_kpoint_converged_stats_nonreference"]
        )
        self.assertIsNone(method["mae_mesh_difference_bound"])
        unresolved_row = next(
            row
            for row in rows
            if row["method"] == "GXTB" and row["phase"] == "XVII"
        )
        for field in (
            "smallest_required_mesh",
            "mesh_label",
            "mesh_n",
            "nk_total",
            "relative_energy_kJmol_per_H2O",
            "DMC_relative_kJmol_per_H2O",
            "error_kJmol_per_H2O",
            "last_delta_kJmol_per_H2O",
        ):
            self.assertIsNone(unresolved_row[field])

    def test_phasewise_uses_earliest_available_adjacent_pair_across_gap(
        self,
    ) -> None:
        relative = {
            "k111": self.phase_values(0.00),
            "k222": self.phase_values(0.04),
            "k333": self.phase_values(0.07),
            "k444": self.phase_values(0.08),
            "k555": self.phase_values(0.09),
        }
        del relative["k222"]["II"]
        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", relative)
        )
        method = report["methods"]["GXTB"]
        self.assertTrue(method["phasewise_kpoint_converged"])
        self.assertEqual(
            (
                method["phase_convergence"]["II"]["previous_mesh"],
                method["phase_convergence"]["II"][
                    "smallest_required_mesh"
                ],
            ),
            ("k333", "k444"),
        )

    def test_phasewise_delta_statistics_are_diagnostic_only(self) -> None:
        relative = {
            "k111": self.phase_values(0.00),
            "k222": self.phase_values(0.049),
            "k333": self.phase_values(0.098),
        }
        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", relative)
        )
        method = report["methods"]["GXTB"]
        self.assertTrue(method["phasewise_kpoint_converged"])
        self.assertEqual(method["status"], "phasewise_kpoint_converged")
        diagnostic = method["last_delta_statistics_diagnostic"]
        self.assertAlmostEqual(
            diagnostic["rms_delta_kjmol_per_h2o"], 0.049
        )
        self.assertAlmostEqual(
            diagnostic["mean_absolute_delta_kjmol_per_h2o"], 0.049
        )
        self.assertAlmostEqual(
            diagnostic["max_absolute_delta_kjmol_per_h2o"], 0.049
        )
        self.assertNotIn("passes_overall_rms_limit", diagnostic)
        self.assertIsNotNone(
            method["phasewise_kpoint_converged_stats_nonreference"]
        )

    def test_phasewise_later_evidence_revokes_and_later_pair_recovers(
        self,
    ) -> None:
        phases = [phase for phase in benchmark.PHASES if phase != "Ih"]
        early = {
            "k111": self.phase_values(0.00),
            "k222": self.phase_values(0.01),
            "k333": self.phase_values(0.02),
            "k444": self.phase_values(0.30),
        }
        revoked, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", early)
        )
        revoked_method = revoked["methods"]["GXTB"]
        self.assertEqual(revoked_method["converged_phase_count"], 0)
        self.assertEqual(
            revoked_method["next_required_by_mesh"], {"k555": phases}
        )
        for phase in phases:
            candidate = revoked_method["unresolved"][phase][
                "revoked_candidates"
            ][0]
            self.assertTrue(candidate["revoked_by_later_evidence"])
            self.assertEqual(
                candidate["later_available_contradictions"][0]["from_mesh"],
                "k333",
            )

        recovered_relative = {
            **early,
            "k555": self.phase_values(0.31),
        }
        recovered, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", recovered_relative)
        )
        recovered_method = recovered["methods"]["GXTB"]
        self.assertTrue(recovered_method["phasewise_kpoint_converged"])
        for convergence in recovered_method["phase_convergence"].values():
            self.assertEqual(
                (
                    convergence["previous_mesh"],
                    convergence["smallest_required_mesh"],
                ),
                ("k444", "k555"),
            )

        gap_relative = dict(early)
        gap_relative.pop("k444")
        gap, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", gap_relative)
        )
        self.assertTrue(
            gap["methods"]["GXTB"]["phasewise_kpoint_converged"]
        )
        for convergence in gap["methods"]["GXTB"][
            "phase_convergence"
        ].values():
            self.assertFalse(convergence["revoked_by_later_evidence"])

    def test_phasewise_gfn_retrospective_matches_frozen_baselines(self) -> None:
        frozen = json.loads(
            (
                REPOSITORY
                / "DMC-ICE13"
                / "data"
                / "kpoint_results.json"
            ).read_text()
        )
        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            frozen
        )
        expected = {
            "GFN1": {
                "mae": 8.0064220272827,
                "previous_mae": 8.003066937205277,
                "rms": 0.011073969859886084,
                "maximum": 0.03375683355920267,
                "mean_absolute": 0.005838396106517836,
            },
            "GFN2": {
                "mae": 3.4614288230295376,
                "previous_mae": 3.459322043260181,
                "rms": 0.01769033950512505,
                "maximum": 0.040888492065901616,
                "mean_absolute": 0.011900562418979197,
            },
        }
        for method_name, values in expected.items():
            method = report["methods"][method_name]
            self.assertTrue(method["phasewise_kpoint_converged"])
            self.assertEqual(method["unresolved_phases"], [])
            self.assertEqual(method["next_required_phases"], [])
            self.assertAlmostEqual(
                method[
                    "phasewise_kpoint_converged_stats_nonreference"
                ]["MAE"],
                values["mae"],
            )
            self.assertAlmostEqual(
                method["previous_mesh_comparison_stats_nonreference"]["MAE"],
                values["previous_mae"],
            )
            diagnostic = method["last_delta_statistics_diagnostic"]
            self.assertAlmostEqual(
                diagnostic["rms_delta_kjmol_per_h2o"], values["rms"]
            )
            self.assertAlmostEqual(
                diagnostic["max_absolute_delta_kjmol_per_h2o"],
                values["maximum"],
            )
            self.assertAlmostEqual(
                diagnostic["mean_absolute_delta_kjmol_per_h2o"],
                values["mean_absolute"],
            )
            self.assertTrue(method["mae_mesh_difference_bound"]["bound_satisfied"])
        self.assertEqual(
            report["methods"]["GFN1"]["phase_convergence"]["VIII"][
                "smallest_required_mesh"
            ],
            "k444",
        )

    def test_phasewise_gxtb_k666_frontier_has_expected_four_phases(
        self,
    ) -> None:
        nonreference = [phase for phase in benchmark.PHASES if phase != "Ih"]
        mesh_ids = ("k111", "k222", "k333", "k444", "k555", "k666")
        relative: dict[str, dict[str, float]] = {
            mesh: {
                phase: float(phase_index) + 0.25 * mesh_index
                for phase_index, phase in enumerate(nonreference)
            }
            for mesh_index, mesh in enumerate(mesh_ids)
        }
        relative["k555"]["XIII"] = relative["k444"]["XIII"] + 0.06
        relative["k666"]["XIII"] = relative["k555"]["XIII"] + 0.01
        for phase in ("III", "IV", "IX"):
            relative["k666"][phase] = relative["k555"][phase] + 0.04

        report, _ = benchmark.build_phasewise_kpoint_convergence_report(
            self.phasewise_results("GXTB", relative)
        )
        method = report["methods"]["GXTB"]
        expected_converged = {"III", "IV", "IX", "XIII"}
        self.assertEqual(method["converged_phase_count"], 4)
        self.assertEqual(
            set(method["phase_convergence"]), expected_converged
        )
        self.assertEqual(
            (
                method["phase_convergence"]["XIII"]["previous_mesh"],
                method["phase_convergence"]["XIII"][
                    "smallest_required_mesh"
                ],
            ),
            ("k555", "k666"),
        )
        for phase in ("III", "IV", "IX"):
            self.assertEqual(
                (
                    method["phase_convergence"][phase]["previous_mesh"],
                    method["phase_convergence"][phase][
                        "smallest_required_mesh"
                    ],
                ),
                ("k555", "k666"),
            )
        expected_unresolved = [
            phase for phase in nonreference if phase not in expected_converged
        ]
        self.assertEqual(method["unresolved_phases"], expected_unresolved)
        self.assertEqual(
            method["next_required_by_mesh"],
            {"k777": expected_unresolved},
        )
        self.assertEqual(
            method["required_same_mesh_reference_by_mesh"],
            {"k777": "Ih"},
        )
        self.assertIsNone(
            method["phasewise_kpoint_converged_stats_nonreference"]
        )


    def test_same_mesh_relative_energy_rejects_inconsistent_cached_value(self) -> None:
        results = self.phasewise_results(
            "GXTB", {"k333": {"II": 1.0}}, {"k333": -123.0}
        )
        results["results"]["k333"]["GXTB"]["relative_kjmol"]["II"] = 2.0
        with self.assertRaisesRegex(ValueError, "same-mesh Ih"):
            benchmark.same_mesh_relative_energy(
                results["results"], "k333", "GXTB", "II"
            )

    def test_validation_index_is_hash_bound_and_rejects_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / runner.GXTB_INPUT_DIRECTORY
            run_root = root / runner.GXTB_RUN_DIRECTORY
            benchmark.prepare_inputs(["GXTB"], input_root, ["k666"])
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("test launcher\n")
            identity = fake_production_identity(root, cp2k)
            campaign_path = root / "campaign.json"
            write_fake_campaign_manifest(campaign_path, identity)
            job = runner.jobs(
                root,
                ["GXTB"],
                ["k666"],
                ["Ih"],
                gxtb_input_root=input_root,
                gxtb_run_root=run_root,
            )[0]
            output = job.run_dir / job.output_name
            complete_job(job, identity, -1.0)
            legacy_stamp = json.loads(runner.stamp_path(job).read_text())
            for field in (
                "schema_version",
                "build_id",
                "frozen_input",
                "frozen_input_sha256",
            ):
                legacy_stamp.pop(field)
            runner.stamp_path(job).write_text(json.dumps(legacy_stamp))
            self.assertTrue(runner.stamp_valid(job, identity))
            args = type(
                "Args",
                (),
                {
                    "root": root,
                    "analysis_prefix": runner.GXTB_ANALYSIS_PREFIX,
                    "gxtb_input_root": input_root,
                    "gxtb_run_root": run_root,
                    "campaign_manifest": campaign_path,
                },
            )()
            index = runner.write_convergence_validation_index(args, identity)
            first_snapshot_bytes = index.read_bytes()
            first_snapshot_digest = runner.sha256(index)
            self.assertIn(first_snapshot_digest, index.name)
            current_index = runner.convergence_validation_index_path(args)
            self.assertEqual(current_index.read_bytes(), first_snapshot_bytes)
            self.assertEqual(
                runner.write_convergence_validation_index(args, identity), index
            )
            with mock.patch.object(benchmark, "ROOT", root):
                payload = benchmark.read_validation_index(
                    index, campaign_manifest_path=campaign_path
                )
                self.assertEqual(payload["validated_phase_coverage"], {"k666": ["Ih"]})
                with self.assertRaisesRegex(ValueError, "output path mismatch"):
                    benchmark.validated_gxtb_output_paths(
                        payload, root / "unvalidated-derived-run-root"
                    )
                record = payload["records"][0]
                self.assertEqual(
                    record["input_sha256"], record["frozen_input_sha256"]
                )

                def assert_record_index_rejected(candidate, pattern, name):
                    candidate_path = root / f"record-{name}.json"
                    candidate_path.write_text(json.dumps(candidate))
                    with self.assertRaisesRegex(ValueError, pattern):
                        runner.read_validation_index(
                            candidate_path,
                            root,
                            campaign_manifest_path=campaign_path,
                        )
                    with self.assertRaisesRegex(ValueError, pattern):
                        benchmark.read_validation_index(
                            candidate_path,
                            campaign_manifest_path=campaign_path,
                        )

                relabeled = json.loads(index.read_text())
                relabeled["records"][0]["phase"] = "VII"
                assert_record_index_rejected(
                    relabeled, "input path is not canonical", "relabel-v2"
                )

                shared_input = json.loads(index.read_text())
                shared_input["records"][0]["frozen_input"] = shared_input[
                    "records"
                ][0]["input"]
                assert_record_index_rejected(
                    shared_input, "frozen input path mismatch", "shared-input"
                )

                original_output_bytes = output.read_bytes()
                original_stamp_bytes = runner.stamp_path(job).read_bytes()
                original_output_text = original_output_bytes.decode()
                energy_mutations = (
                    ("nan", "energy [a.u.] -1.0", "energy [a.u.] nan", "finite energy"),
                    ("inf", "energy [a.u.] -1.0", "energy [a.u.] inf", "finite energy"),
                    (
                        "malformed",
                        "energy [a.u.] -1.0",
                        "energy [a.u.] not-a-number",
                        "invalid energy",
                    ),
                    (
                        "multiple",
                        "PROGRAM ENDED",
                        "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.] -1.0\n"
                        "PROGRAM ENDED",
                        "exactly one total energy",
                    ),
                )
                for name, old, new, pattern in energy_mutations:
                    try:
                        output.write_text(original_output_text.replace(old, new))
                        mutated_stamp = json.loads(original_stamp_bytes)
                        mutated_stamp["output_sha256"] = runner.sha256(output)
                        runner.stamp_path(job).write_text(
                            json.dumps(mutated_stamp, sort_keys=True)
                        )
                        candidate = json.loads(index.read_text())
                        candidate["records"][0]["output_sha256"] = runner.sha256(
                            output
                        )
                        candidate["records"][0]["stamp_sha256"] = runner.sha256(
                            runner.stamp_path(job)
                        )
                        assert_record_index_rejected(
                            candidate, pattern, f"energy-{name}"
                        )
                    finally:
                        output.write_bytes(original_output_bytes)
                        runner.stamp_path(job).write_bytes(original_stamp_bytes)
                for missing_field in ("frozen_input", "frozen_input_sha256"):
                    missing_payload = json.loads(index.read_text())
                    missing_payload["records"][0].pop(missing_field)
                    missing_path = root / f"missing-{missing_field}.json"
                    missing_path.write_text(json.dumps(missing_payload))
                    with self.assertRaisesRegex(ValueError, "frozen"):
                        benchmark.read_validation_index(
                            missing_path, campaign_manifest_path=campaign_path
                        )
                    with self.assertRaisesRegex(ValueError, "frozen"):
                        runner.read_validation_index(
                            missing_path,
                            root,
                            campaign_manifest_path=campaign_path,
                        )
                legacy_payload = json.loads(json.dumps(payload))
                legacy_payload["schema_version"] = 1
                legacy_payload["build_identity"] = next(
                    iter(legacy_payload["build_identities"].values())
                )
                legacy_payload.pop("build_identities")
                legacy_payload.pop("source_identity")
                legacy_payload.pop("parents")
                legacy_payload.pop("campaign_manifest_sha256", None)
                legacy_payload["supported_dense_extensions"] = ["k666", "k777"]
                legacy_payload["records"][0].pop("frozen_input")
                legacy_payload["records"][0].pop("frozen_input_sha256")
                legacy_payload["records"][0].pop("campaign_id")
                legacy_payload["records"][0].pop("gxtb_protocol_id")
                legacy_payload["records"][0].pop("build_id")
                legacy_index = root / "legacy-validation-index.json"
                legacy_index.write_text(json.dumps(legacy_payload))
                loaded_legacy = benchmark.read_validation_index(legacy_index)
                self.assertEqual(
                    loaded_legacy["supported_dense_extensions"],
                    ["k666", "k777"],
                )
                relabeled_legacy = json.loads(legacy_index.read_text())
                relabeled_legacy["records"][0]["phase"] = "VII"
                assert_record_index_rejected(
                    relabeled_legacy,
                    "input path is not canonical",
                    "relabel-v1",
                )
                output.write_text(output.read_text() + "tampered\n")
                with self.assertRaisesRegex(ValueError, "invalid (?:output )?hash"):
                    benchmark.read_validation_index(
                        index, campaign_manifest_path=campaign_path
                    )

            complete_job(job, identity, -2.0)
            second_index = runner.write_convergence_validation_index(args, identity)
            self.assertNotEqual(second_index, index)
            self.assertEqual(index.read_bytes(), first_snapshot_bytes)
            self.assertEqual(current_index.read_bytes(), second_index.read_bytes())

    def test_schema_v1_is_read_immutably_and_mixed_builds_merge_per_record(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / runner.GXTB_INPUT_DIRECTORY
            run_root = root / runner.GXTB_RUN_DIRECTORY
            benchmark.prepare_inputs(["GXTB"], input_root, ["k666"])
            cp2k_arm = root / "cp2k-arm"
            cp2k_arm.write_text("arm launcher\n")
            arm = fake_production_identity(root, cp2k_arm)
            campaign_path = root / "campaign.json"
            write_fake_campaign_manifest(campaign_path, arm)
            ih, ii, vii = runner.jobs(
                root,
                ["GXTB"],
                ["k666"],
                ["Ih", "II", "VII"],
                gxtb_input_root=input_root,
                gxtb_run_root=run_root,
            )
            complete_job(ih, arm, -1.0)
            complete_job(vii, arm, -0.99)
            args = type(
                "Args",
                (),
                {
                    "root": root,
                    "analysis_prefix": runner.GXTB_ANALYSIS_PREFIX,
                    "gxtb_input_root": input_root,
                    "gxtb_run_root": run_root,
                    "base_validation_index": None,
                    "base_validation_index_payload": None,
                    "execution_build_manifest": None,
                    "campaign_manifest": campaign_path,
                },
            )()
            initial_v2 = runner.write_convergence_validation_index(args, arm)
            initial_payload = json.loads(initial_v2.read_text())
            legacy_records = []
            for value in initial_payload["records"]:
                record = dict(value)
                record.pop("campaign_id")
                record.pop("gxtb_protocol_id")
                record.pop("build_id")
                legacy_records.append(record)
            legacy_payload = {
                "schema_version": 1,
                "benchmark": "DMC-ICE13",
                "method": "GXTB",
                "campaign_id": arm.campaign_id,
                "gxtb_protocol_id": runner.GXTB_PROTOCOL_ID,
                "core_meshes": runner.MESHES,
                "supported_dense_extensions": runner.DENSE_EXTENSION_MESHES,
                "build_identity": runner.execution_build_identity(arm),
                "validated_phase_coverage": {"k666": ["Ih", "VII"]},
                "records": legacy_records,
            }
            legacy = root / "immutable-schema-v1.json"
            legacy.write_text(json.dumps(legacy_payload, indent=2, sort_keys=True) + "\n")
            legacy_bytes = legacy.read_bytes()
            base = runner.read_validation_index(
                legacy,
                root,
                expected_campaign_id=arm.campaign_id,
                expected_source_identity={
                    "cp2k_source_revision": arm.cp2k_source_revision,
                    "tblite_source_revision": arm.tblite_source_revision,
                },
            )
            with mock.patch.object(benchmark, "ROOT", root):
                self.assertEqual(
                    benchmark.read_validation_index(legacy)["source_schema_version"],
                    1,
                )
            cp2k_remote = root / "cp2k-x86_64"
            cp2k_remote.write_text("remote launcher\n")
            remote = fake_production_identity(root, cp2k_remote)
            tblite = root / "tblite"
            tblite.write_text("remote cli\n")
            execution_path, _, _ = write_qualified_execution_manifest(
                root, remote, tblite, campaign_path
            )
            complete_job(ii, remote, -9.0)
            args.base_validation_index = legacy
            args.base_validation_index_payload = base
            args.execution_build_manifest = execution_path
            merged = runner.write_convergence_validation_index(args, remote)
            self.assertEqual(legacy.read_bytes(), legacy_bytes)
            merged_payload = runner.read_validation_index(
                merged,
                root,
                expected_campaign_id=arm.campaign_id,
                expected_source_identity={
                    "cp2k_source_revision": arm.cp2k_source_revision,
                    "tblite_source_revision": arm.tblite_source_revision,
                },
                campaign_manifest_path=campaign_path,
            )
            self.assertEqual(
                merged_payload["validated_phase_coverage"],
                {"k666": ["Ih", "II", "VII"]},
            )
            self.assertEqual(len(merged_payload["build_identities"]), 2)
            builds_by_phase = {
                record["phase"]: record["build_id"]
                for record in merged_payload["records"]
            }
            self.assertEqual(
                builds_by_phase,
                {"Ih": runner.build_id(runner.execution_build_identity(arm)),
                 "II": runner.build_id(runner.execution_build_identity(remote)),
                 "VII": runner.build_id(runner.execution_build_identity(arm))},
            )
            with mock.patch.object(benchmark, "ROOT", root):
                analysed = benchmark.read_validation_index(
                    merged, campaign_manifest_path=campaign_path
                )
                self.assertEqual(len(analysed["build_identities"]), 2)

            tampered = json.loads(merged.read_text())
            first_identity = next(iter(tampered["build_identities"].values()))
            first_identity["cp2k_sha256"] = "0" * 64
            tampered_path = root / "tampered-build.json"
            tampered_path.write_text(json.dumps(tampered))
            with self.assertRaisesRegex(ValueError, "identity digest mismatch"):
                runner.read_validation_index(
                    tampered_path, root, campaign_manifest_path=campaign_path
                )
            record_tampered = json.loads(merged.read_text())
            remote_build_id = builds_by_phase["II"]
            next(
                record
                for record in record_tampered["records"]
                if record["phase"] == "Ih"
            )["build_id"] = remote_build_id
            record_tampered_path = root / "tampered-record-build.json"
            record_tampered_path.write_text(json.dumps(record_tampered))
            with self.assertRaisesRegex(ValueError, "stamp mismatch"):
                runner.read_validation_index(
                    record_tampered_path,
                    root,
                    campaign_manifest_path=campaign_path,
                )
            with self.assertRaisesRegex(ValueError, "source revision mismatch"):
                runner.read_validation_index(
                    merged,
                    root,
                    expected_source_identity={
                        "cp2k_source_revision": "f" * 40,
                        "tblite_source_revision": arm.tblite_source_revision,
                    },
                    campaign_manifest_path=campaign_path,
                )

    def test_v2_revalidates_execution_manifest_and_numeric_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / runner.GXTB_INPUT_DIRECTORY
            run_root = root / runner.GXTB_RUN_DIRECTORY
            benchmark.prepare_inputs(
                ["GXTB"], input_root, ["k666"], ["Ih", "II"]
            )
            cp2k = root / "cp2k-remote"
            cp2k.write_text("remote launcher\n")
            tblite = root / "tblite"
            tblite.write_text("remote cli\n")
            identity = fake_production_identity(root, cp2k)
            campaign_path = root / "campaign.json"
            frozen = runner.execution_build_identity(identity)
            frozen.update(
                {
                    "cp2k_sha256": "a" * 64,
                    "cp2k_library_sha256": "b" * 64,
                    "tblite_static_library_sha256": "c" * 64,
                }
            )
            write_fake_campaign_manifest(
                campaign_path, identity, frozen_identity=frozen
            )
            execution_path, execution_payload, _ = write_qualified_execution_manifest(
                root, identity, tblite, campaign_path
            )
            frozen_production = replace(
                identity,
                cp2k_sha256=str(frozen["cp2k_sha256"]),
                cp2k_library_sha256=str(frozen["cp2k_library_sha256"]),
                tblite_static_library_sha256=str(
                    frozen["tblite_static_library_sha256"]
                ),
            )
            args = type(
                "Args",
                (),
                {
                    "root": root,
                    "analysis_prefix": runner.GXTB_ANALYSIS_PREFIX,
                    "gxtb_input_root": input_root,
                    "gxtb_run_root": run_root,
                    "campaign_manifest": campaign_path,
                    "base_validation_index": None,
                    "base_validation_index_payload": None,
                    "base_validation_index_sha256": None,
                    "execution_build_manifest": None,
                },
            )()
            base_index = runner.write_convergence_validation_index(
                args, frozen_production
            )
            base_payload = runner.read_validation_index(
                base_index, root, campaign_manifest_path=campaign_path
            )
            job = runner.jobs(
                root,
                ["GXTB"],
                ["k666"],
                ["II"],
                gxtb_input_root=input_root,
                gxtb_run_root=run_root,
            )[0]
            complete_job(job, identity, -10.0)
            args.base_validation_index = base_index
            args.base_validation_index_payload = base_payload
            args.base_validation_index_sha256 = runner.sha256(base_index)
            args.execution_build_manifest = execution_path
            index = runner.write_convergence_validation_index(args, identity)
            runner.read_validation_index(
                index, root, campaign_manifest_path=campaign_path
            )
            with mock.patch.object(benchmark, "ROOT", root):
                benchmark.read_validation_index(
                    index, campaign_manifest_path=campaign_path
                )
            original = json.loads(index.read_text())
            identity_id = runner.build_id(runner.execution_build_identity(identity))

            def assert_both_reject(payload, pattern, name):
                candidate = root / f"invalid-{name}.json"
                candidate.write_text(json.dumps(payload))
                with self.assertRaisesRegex(ValueError, pattern):
                    runner.read_validation_index(
                        candidate, root, campaign_manifest_path=campaign_path
                    )
                with mock.patch.object(benchmark, "ROOT", root):
                    with self.assertRaisesRegex(ValueError, pattern):
                        benchmark.read_validation_index(
                            candidate, campaign_manifest_path=campaign_path
                        )

            def assert_bad_manifest_rejected(bad_manifest, pattern, name):
                bad_manifest_path = root / f"execution-{name}.json"
                bad_manifest_path.write_text(json.dumps(bad_manifest))
                bad_index = json.loads(json.dumps(original))
                bad_identity = bad_index["build_identities"][identity_id]
                bad_identity["execution_build_manifest"] = bad_manifest_path.name
                bad_identity[
                    "execution_build_manifest_sha256"
                ] = runner.sha256(bad_manifest_path)
                assert_both_reject(bad_index, pattern, name)

            alternate_record = next(
                record for record in original["records"]
                if record["build_id"] == identity_id
            )
            alternate_stamp_path = root / alternate_record["stamp"]
            alternate_stamp_bytes = alternate_stamp_path.read_bytes()
            try:
                downgraded_stamp = json.loads(alternate_stamp_bytes)
                for field in (
                    "schema_version",
                    "build_id",
                    "frozen_input",
                    "frozen_input_sha256",
                ):
                    downgraded_stamp.pop(field)
                alternate_stamp_path.write_text(json.dumps(downgraded_stamp))
                downgraded_index = json.loads(json.dumps(original))
                next(
                    record for record in downgraded_index["records"]
                    if record["build_id"] == identity_id
                )["stamp_sha256"] = runner.sha256(alternate_stamp_path)
                assert_both_reject(
                    downgraded_index, "schema_version", "stamp-downgrade"
                )
            finally:
                alternate_stamp_path.write_bytes(alternate_stamp_bytes)

            def assert_same_bytes_survive_race(target, replacement, name):
                original_bytes = target.read_bytes()

                def run_one(reader, reader_name):
                    replaced = False
                    original_read_bytes = Path.read_bytes

                    def racing_read_bytes(path):
                        nonlocal replaced
                        content = original_read_bytes(path)
                        if path.resolve() == target.resolve() and not replaced:
                            path.write_bytes(replacement)
                            replaced = True
                        return content

                    try:
                        with mock.patch.object(
                            Path, "read_bytes", racing_read_bytes
                        ):
                            reader()
                    finally:
                        target.write_bytes(original_bytes)
                    self.assertTrue(replaced, f"{name}/{reader_name} race not triggered")

                run_one(
                    lambda: runner.read_validation_index(
                        index, root, campaign_manifest_path=campaign_path
                    ),
                    "runner",
                )
                with mock.patch.object(benchmark, "ROOT", root):
                    run_one(
                        lambda: benchmark.read_validation_index(
                            index, campaign_manifest_path=campaign_path
                        ),
                        "analyzer",
                    )

            remote_phase_output = root / execution_payload["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]["remote_phase_output"]
            assert_same_bytes_survive_race(
                remote_phase_output,
                b"SCF run converged\nPROGRAM ENDED\n",
                "evidence-output",
            )
            remote_record = next(
                record for record in original["records"]
                if record["build_id"] == identity_id
            )
            assert_same_bytes_survive_race(
                root / remote_record["output"],
                b"SCF run converged\nPROGRAM ENDED\n",
                "record-output",
            )
            assert_same_bytes_survive_race(
                execution_path,
                b"{ invalid execution manifest",
                "execution-manifest",
            )

            incomplete = json.loads(json.dumps(original))
            incomplete["build_identities"][identity_id].pop(
                "execution_build_manifest_sha256"
            )
            assert_both_reject(incomplete, "incomplete execution manifest", "missing")

            unqualified = json.loads(json.dumps(original))
            unqualified_identity = unqualified["build_identities"][identity_id]
            unqualified_identity.pop("execution_build_manifest")
            unqualified_identity.pop("execution_build_manifest_sha256")
            assert_both_reject(
                unqualified,
                "alternate validation build identity lacks execution manifest",
                "unqualified-alternate",
            )

            hash_tampered = json.loads(json.dumps(original))
            hash_tampered["build_identities"][identity_id][
                "execution_build_manifest_sha256"
            ] = "0" * 64
            assert_both_reject(hash_tampered, "manifest hash mismatch", "hash")

            unused = json.loads(json.dumps(original))
            extra_identity = dict(frozen)
            extra_identity["cp2k_sha256"] = "d" * 64
            unused["build_identities"][runner.build_id(extra_identity)] = extra_identity
            assert_both_reject(
                unused,
                "alternate validation build identity lacks execution manifest",
                "unused",
            )

            campaign_tampered = json.loads(json.dumps(original))
            campaign_tampered["campaign_manifest_sha256"] = "4" * 64
            assert_both_reject(
                campaign_tampered,
                "campaign manifest hash mismatch",
                "campaign-hash",
            )

            for field, replacement in (
                ("benchmark", "X23b"),
                ("method", "GFN2"),
            ):
                wrong_scope = json.loads(json.dumps(original))
                wrong_scope[field] = replacement
                assert_both_reject(
                    wrong_scope, "benchmark/method mismatch", f"wrong-{field}"
                )

            missing_campaign = root / "does-not-exist-build-manifest.json"
            with self.assertRaisesRegex(ValueError, "trusted campaign manifest"):
                runner.read_validation_index(
                    index, root, campaign_manifest_path=missing_campaign
                )
            with mock.patch.object(benchmark, "ROOT", root):
                with self.assertRaisesRegex(
                    ValueError, "trusted campaign manifest"
                ):
                    benchmark.read_validation_index(
                        index, campaign_manifest_path=missing_campaign
                    )

            traversal = json.loads(json.dumps(original))
            traversal["build_identities"][identity_id][
                "execution_build_manifest"
            ] = "../execution-build-qualified.json"
            assert_both_reject(traversal, "not relative and safe", "traversal")

            sentinel_template = execution_payload["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]

            def replace_manifest_artifact(
                bad_manifest, artifact, replacement_path
            ):
                sentinel = bad_manifest["qualification"][
                    "same_mesh_dense_relative_sentinels"
                ][0]
                sentinel[artifact] = str(replacement_path.relative_to(root))
                sentinel[f"{artifact}_sha256"] = runner.sha256(replacement_path)

            def replace_remote_output(
                bad_manifest, system, content, case_name
            ):
                phase = "VII" if system == "phase" else "Ih"
                project = f"ice_{phase}_GXTB_k666"
                replacement_dir = (
                    root
                    / "qualification_evidence"
                    / case_name
                    / "k666"
                    / phase
                )
                replacement_dir.mkdir(parents=True, exist_ok=True)
                replacement_input = replacement_dir / f"{project}.inp"
                good_input = root / sentinel_template[f"{system}_input"]
                shutil.copy2(good_input, replacement_input)
                replacement_output = replacement_dir / f"{project}.out"
                replacement_output.write_text(content)
                good_stamp = json.loads(
                    (root / sentinel_template[f"remote_{system}_stamp"]).read_text()
                )
                good_stamp["output_sha256"] = runner.sha256(replacement_output)
                replacement_stamp = replacement_output.with_suffix(".run.json")
                replacement_stamp.write_text(json.dumps(good_stamp, sort_keys=True))
                sentinel = bad_manifest["qualification"][
                    "same_mesh_dense_relative_sentinels"
                ][0]
                for artifact, path in (
                    (f"{system}_input", replacement_input),
                    (f"remote_{system}_output", replacement_output),
                    (f"remote_{system}_stamp", replacement_stamp),
                ):
                    sentinel[artifact] = str(path.relative_to(root))
                    sentinel[f"{artifact}_sha256"] = runner.sha256(path)

            bad_input_dir = root / "qualification_evidence" / "bad-input"
            bad_input_dir.mkdir()
            bad_phase_input = bad_input_dir / "ice_VII_GXTB_k666.inp"
            good_phase_input = root / sentinel_template["phase_input"]
            bad_phase_input.write_text(
                good_phase_input.read_text().replace("METHOD GXTB", "METHOD GFN2")
            )
            bad_input_manifest = json.loads(json.dumps(execution_payload))
            replace_manifest_artifact(
                bad_input_manifest, "phase_input", bad_phase_input
            )
            assert_bad_manifest_rejected(
                bad_input_manifest, "violates input contract", "bad-input-contract"
            )

            for field, replacement in (("mesh", "k777"), ("phase", "II")):
                wrong_semantics_manifest = json.loads(
                    json.dumps(execution_payload)
                )
                wrong_semantics_manifest["qualification"][
                    "same_mesh_dense_relative_sentinels"
                ][0][field] = replacement
                assert_bad_manifest_rejected(
                    wrong_semantics_manifest,
                    "input filename mismatch",
                    f"wrong-sentinel-{field}",
                )

            wrong_count_manifest = json.loads(json.dumps(execution_payload))
            wrong_count_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]["phase_water_count"] = 11
            assert_bad_manifest_rejected(
                wrong_count_manifest, "water_count does not match input", "water-count"
            )

            good_remote_phase = root / sentinel_template["remote_phase_output"]
            output_mutations = (
                (
                    "output-input-header",
                    "CP2K| Input file name ice_VII_GXTB_k666.inp",
                    "CP2K| Input file name ice_II_GXTB_k666.inp",
                    "input header mismatch",
                ),
                (
                    "output-project-header",
                    "GLOBAL| Project name ice_VII_GXTB_k666",
                    "GLOBAL| Project name ice_II_GXTB_k666",
                    "project header mismatch",
                ),
                (
                    "output-revision",
                    identity.cp2k_source_revision[:10],
                    "f" * 10,
                    "source revision mismatch",
                ),
                (
                    "output-tblite-revision",
                    identity.tblite_source_revision,
                    "f" * 40,
                    "tblite source revision mismatch",
                ),
            )
            for name, old, new, pattern in output_mutations:
                bad_output_manifest = json.loads(json.dumps(execution_payload))
                replace_remote_output(
                    bad_output_manifest,
                    "phase",
                    good_remote_phase.read_text().replace(old, new),
                    name,
                )
                assert_bad_manifest_rejected(
                    bad_output_manifest, pattern, name
                )

            copied_reference_manifest = json.loads(json.dumps(execution_payload))
            reference_phase_output = root / sentinel_template[
                "reference_phase_output"
            ]
            replace_remote_output(
                copied_reference_manifest,
                "phase",
                reference_phase_output.read_text(),
                "copied-reference-output",
            )
            assert_bad_manifest_rejected(
                copied_reference_manifest,
                "copied reference phase output",
                "copied-reference-output",
            )

            renamed_manifest = json.loads(json.dumps(execution_payload))
            renamed_sentinel = renamed_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]
            original_remote_output = root / sentinel_template[
                "remote_phase_output"
            ]
            original_remote_stamp = root / sentinel_template[
                "remote_phase_stamp"
            ]
            renamed_output = original_remote_output.with_name("renamed.out")
            renamed_stamp = original_remote_stamp.with_name("renamed.run.json")
            shutil.copy2(original_remote_output, renamed_output)
            shutil.copy2(original_remote_stamp, renamed_stamp)
            for artifact, replacement in (
                ("remote_phase_output", renamed_output),
                ("remote_phase_stamp", renamed_stamp),
            ):
                renamed_sentinel[artifact] = str(replacement.relative_to(root))
                renamed_sentinel[f"{artifact}_sha256"] = runner.sha256(replacement)
            assert_bad_manifest_rejected(
                renamed_manifest,
                "remote phase output filename mismatch",
                "renamed-remote-output",
            )

            missing_stamp_manifest = json.loads(json.dumps(execution_payload))
            missing_stamp_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0].pop("remote_phase_stamp")
            assert_bad_manifest_rejected(
                missing_stamp_manifest,
                "empty artifact path",
                "missing-remote-stamp",
            )

            bad_stamp_manifest = json.loads(json.dumps(execution_payload))
            replace_remote_output(
                bad_stamp_manifest,
                "phase",
                original_remote_output.read_text(),
                "bad-remote-stamp",
            )
            bad_stamp_sentinel = bad_stamp_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]
            bad_stamp_path = root / bad_stamp_sentinel["remote_phase_stamp"]
            bad_stamp = json.loads(bad_stamp_path.read_text())
            bad_stamp["tblite_source_revision"] = "f" * 40
            bad_stamp_path.write_text(json.dumps(bad_stamp, sort_keys=True))
            bad_stamp_sentinel["remote_phase_stamp_sha256"] = runner.sha256(
                bad_stamp_path
            )
            assert_bad_manifest_rejected(
                bad_stamp_manifest,
                "stamp mismatch: tblite_source_revision",
                "bad-remote-stamp",
            )

            mislabeled_geometry_manifest = json.loads(
                json.dumps(execution_payload)
            )
            replace_remote_output(
                mislabeled_geometry_manifest,
                "phase",
                original_remote_output.read_text(),
                "mislabeled-phase-geometry",
            )
            mislabeled_sentinel = mislabeled_geometry_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]
            mislabeled_input = root / mislabeled_sentinel["phase_input"]
            mislabeled_input.write_text(
                (root / sentinel_template["ih_input"])
                .read_text()
                .replace("ice_Ih_GXTB_k666", "ice_VII_GXTB_k666")
            )
            mislabeled_sentinel["phase_input_sha256"] = runner.sha256(
                mislabeled_input
            )
            mislabeled_stamp_path = root / mislabeled_sentinel[
                "remote_phase_stamp"
            ]
            mislabeled_stamp = json.loads(mislabeled_stamp_path.read_text())
            mislabeled_stamp["input_sha256"] = runner.sha256(mislabeled_input)
            mislabeled_stamp["frozen_input_sha256"] = runner.sha256(
                mislabeled_input
            )
            mislabeled_stamp_path.write_text(
                json.dumps(mislabeled_stamp, sort_keys=True)
            )
            mislabeled_sentinel["remote_phase_stamp_sha256"] = runner.sha256(
                mislabeled_stamp_path
            )
            assert_bad_manifest_rejected(
                mislabeled_geometry_manifest,
                "input is not canonical reference bytes",
                "mislabeled-phase-geometry",
            )

            duplicate_output_manifest = json.loads(json.dumps(execution_payload))
            duplicate_sentinel = duplicate_output_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]
            duplicate_sentinel["remote_ih_output"] = duplicate_sentinel[
                "remote_phase_output"
            ]
            duplicate_sentinel[
                "remote_ih_output_sha256"
            ] = duplicate_sentinel["remote_phase_output_sha256"]
            assert_bad_manifest_rejected(
                duplicate_output_manifest,
                "output paths are not distinct",
                "duplicate-phase-Ih-output",
            )

            for build_field in ("remote_build_id", "reference_build_id"):
                wrong_build_manifest = json.loads(json.dumps(execution_payload))
                wrong_build_manifest["qualification"][
                    "same_mesh_dense_relative_sentinels"
                ][0][build_field] = "sha256:" + "0" * 64
                assert_bad_manifest_rejected(
                    wrong_build_manifest,
                    f"{build_field.removesuffix('_id').replace('_', ' ')} mismatch",
                    f"wrong-{build_field}",
                )

            for name, mutate, pattern in (
                (
                    "loose",
                    lambda value: value["qualification"].__setitem__(
                        "total_energy_tolerance_hartree", 1.0e-9
                    ),
                    "looser than 1e-10",
                ),
                (
                    "numeric",
                    lambda value: value["qualification"][
                        "same_mesh_dense_relative_sentinels"
                    ][0].__setitem__(
                        "relative_energy_delta_kjmol_per_h2o", 5.0e-4
                    ),
                    "relative-energy mismatch",
                ),
                (
                    "fabricated-relative",
                    lambda value: value["qualification"][
                        "same_mesh_dense_relative_sentinels"
                    ][0].__setitem__(
                        "remote_relative_energy_kjmol_per_h2o", 1.0e9
                    ),
                    "remote relative-energy mismatch",
                ),
                (
                    "missing-Ih",
                    lambda value: value["qualification"][
                        "same_mesh_dense_relative_sentinels"
                    ][0].pop("remote_ih_output"),
                    "empty artifact path",
                ),
                (
                    "core",
                    lambda value: value["build_identity"].__setitem__(
                        "cp2k_sha256", "0" * 64
                    ),
                    "execution-build manifest mismatch",
                ),
            ):
                bad_manifest = json.loads(json.dumps(execution_payload))
                mutate(bad_manifest)
                bad_manifest_path = root / f"execution-{name}.json"
                bad_manifest_path.write_text(json.dumps(bad_manifest))
                bad_index = json.loads(json.dumps(original))
                bad_identity = bad_index["build_identities"][identity_id]
                bad_identity["execution_build_manifest"] = bad_manifest_path.name
                bad_identity["execution_build_manifest_sha256"] = runner.sha256(
                    bad_manifest_path
                )
                assert_both_reject(bad_index, pattern, name)

            bad_maximum_manifest = json.loads(json.dumps(execution_payload))
            good_remote_ih = root / sentinel_template["remote_ih_output"]
            replace_remote_output(
                bad_maximum_manifest,
                "ih",
                good_remote_ih.read_text().replace(
                    "energy [a.u.] -1.0", "energy [a.u.] -0.99999999995"
                ),
                "observed-Ih-maximum",
            )
            maximum_sentinel = bad_maximum_manifest["qualification"][
                "same_mesh_dense_relative_sentinels"
            ][0]
            maximum_sentinel["ih_total_energy_delta_hartree"] = 5.0e-11
            remote_relative = (
                ((-0.99) - (-0.99999999995)) / 12
            ) * runner.HARTREE_TO_KJMOL
            reference_relative = 0.01 / 12 * runner.HARTREE_TO_KJMOL
            maximum_sentinel[
                "remote_relative_energy_kjmol_per_h2o"
            ] = remote_relative
            maximum_sentinel[
                "reference_relative_energy_kjmol_per_h2o"
            ] = reference_relative
            maximum_sentinel[
                "relative_energy_delta_kjmol_per_h2o"
            ] = abs(remote_relative - reference_relative)
            bad_maximum_manifest["qualification"][
                "observed_max_abs_relative_energy_delta_kjmol_per_h2o"
            ] = abs(remote_relative - reference_relative)
            bad_maximum_path = root / "execution-bad-total-maximum.json"
            bad_maximum_path.write_text(json.dumps(bad_maximum_manifest))
            bad_maximum_index = json.loads(json.dumps(original))
            bad_maximum_identity = bad_maximum_index["build_identities"][
                identity_id
            ]
            bad_maximum_identity[
                "execution_build_manifest"
            ] = bad_maximum_path.name
            bad_maximum_identity[
                "execution_build_manifest_sha256"
            ] = runner.sha256(bad_maximum_path)
            assert_both_reject(
                bad_maximum_index,
                "observed total-energy maximum mismatch",
                "observed-Ih-maximum",
            )

    def test_base_index_resume_preserves_files_and_refuses_record_collision(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / runner.GXTB_INPUT_DIRECTORY
            run_root = root / runner.GXTB_RUN_DIRECTORY
            benchmark.prepare_inputs(["GXTB"], input_root, ["k666"])
            cp2k_a = root / "cp2k-a"
            cp2k_a.write_text("a\n")
            identity_a = fake_production_identity(root, cp2k_a)
            campaign_path = root / "campaign.json"
            write_fake_campaign_manifest(campaign_path, identity_a)
            job = runner.jobs(
                root,
                ["GXTB"],
                ["k666"],
                ["Ih"],
                gxtb_input_root=input_root,
                gxtb_run_root=run_root,
            )[0]
            complete_job(job, identity_a, -1.0)
            args = type(
                "Args",
                (),
                {
                    "root": root,
                    "analysis_prefix": runner.GXTB_ANALYSIS_PREFIX,
                    "gxtb_input_root": input_root,
                    "gxtb_run_root": run_root,
                    "base_validation_index": None,
                    "base_validation_index_payload": None,
                    "execution_build_manifest": None,
                    "campaign_manifest": campaign_path,
                },
            )()
            base_path = runner.write_convergence_validation_index(args, identity_a)
            base_payload = runner.read_validation_index(
                base_path, root, campaign_manifest_path=campaign_path
            )
            with self.assertRaisesRegex(
                ValueError, "validation index SHA256 pin mismatch"
            ):
                runner.read_validation_index(
                    base_path,
                    root,
                    campaign_manifest_path=campaign_path,
                    expected_index_sha256="0" * 64,
                )
            mtimes = {
                path: path.stat().st_mtime_ns
                for path in (
                    job.input_path,
                    runner.frozen_input_path(job),
                    job.run_dir / job.output_name,
                    runner.stamp_path(job),
                )
            }
            cp2k_b = root / "cp2k-b"
            cp2k_b.write_text("b\n")
            identity_b = fake_production_identity(root, cp2k_b)
            args.base_validation_index = base_path
            args.base_validation_index_payload = base_payload
            args.base_validation_index_sha256 = runner.sha256(base_path)
            preserved = runner.write_convergence_validation_index(args, identity_b)
            self.assertEqual(
                runner.read_validation_index(
                    preserved, root, campaign_manifest_path=campaign_path
                )["validated_phase_coverage"],
                {"k666": ["Ih"]},
            )
            self.assertEqual(
                mtimes,
                {path: path.stat().st_mtime_ns for path in mtimes},
            )

            base_bytes = base_path.read_bytes()
            base_path.write_bytes(base_bytes + b" ")
            with self.assertRaisesRegex(ValueError, "changed after initial"):
                runner.write_convergence_validation_index(args, identity_b)
            base_path.write_bytes(base_bytes)

            args.base_validation_index_sha256 = None
            complete_job(job, identity_b, -2.0)
            with self.assertRaisesRegex(ValueError, "validation record collision"):
                runner.write_convergence_validation_index(args, identity_b)

    def test_restricted_analysis_rejects_unvalidated_gxtb_prior_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / "data"
            figures = root / "figures"
            data.mkdir()
            geometries = {
                phase: {"counts": {"O": 1}}
                for phase in benchmark.PHASES
            }
            (data / "geometries.json").write_text(json.dumps(geometries))
            complete_prior = {
                "complete": True,
                "energies_hartree": {phase: -1.0 for phase in benchmark.PHASES},
                "per_h2o_hartree": {phase: -1.0 for phase in benchmark.PHASES},
                "relative_kjmol": {phase: 0.0 for phase in benchmark.PHASES},
            }
            (data / "kpoint_results.json").write_text(
                json.dumps(
                    {
                        "results": {
                            "gamma": {"GFN2": complete_prior},
                            "k444": {"GXTB": complete_prior},
                        }
                    }
                )
            )
            baseline = (data / "kpoint_results.json").read_text()

            with (
                mock.patch.object(benchmark, "ROOT", root),
                mock.patch.object(benchmark, "DATA", data),
                mock.patch.object(benchmark, "FIGURES", figures),
            ):
                result = benchmark.analyse(
                    validated_gxtb_meshes=set(),
                    output_prefix=benchmark.GXTB_PRODUCTION_PREFIX,
                )

                with (
                    mock.patch.object(
                        benchmark,
                        "read_validation_index",
                        return_value={"source_index_sha256": "a" * 64},
                    ),
                    self.assertRaisesRegex(
                        ValueError, "validation index SHA256 pin mismatch"
                    ),
                ):
                    benchmark.analyse(
                        gxtb_run_root=root / runner.GXTB_RUN_DIRECTORY,
                        output_prefix=benchmark.GXTB_PRODUCTION_PREFIX,
                        validation_index_path=root / "validation-index.json",
                        validation_index_sha256="b" * 64,
                    )

            self.assertTrue(result["results"]["gamma"]["GFN2"]["complete"])
            self.assertFalse(result["results"]["k444"]["GXTB"]["complete"])
            self.assertEqual(result["validated_gxtb_meshes"], [])
            self.assertEqual((data / "kpoint_results.json").read_text(), baseline)
            self.assertTrue(
                (
                    data
                    / "dmc_ice13_gxtb_spglib_kpoint_results.json"
                ).is_file()
            )

    def test_validation_index_analysis_respects_sparse_phase_coverage(self) -> None:
        cases = {
            "k777-Ih-only": [
                ("k777", "Ih", -7.0),
            ],
            "mixed-dense-coverage": [
                ("k666", "Ih", -6.0),
                ("k666", "II", -5.9),
                ("k777", "Ih", -7.0),
                ("k888", "Ih", -8.0),
                ("k888", "VII", -7.9),
            ],
        }
        for label, entries in cases.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                data = root / "data"
                figures = root / "figures"
                data.mkdir()
                (data / "geometries.json").write_text(
                    json.dumps(
                        {
                            phase: {"counts": {"O": 1}}
                            for phase in benchmark.PHASES
                        }
                    )
                )
                index_path = root / "validation-index.json"
                index_path.write_text("{}\n")
                index_hash = benchmark.sha256(index_path)
                records = [
                    {
                        "mesh": mesh,
                        "phase": phase,
                        "validated_energy_hartree": energy,
                    }
                    for mesh, phase, energy in entries
                ]
                validation_index = {
                    "records": records,
                    "source_index_sha256": index_hash,
                }
                verified_paths = {
                    (mesh, phase): root / f"{mesh}-{phase}.out"
                    for mesh, phase, _ in entries
                }
                with (
                    mock.patch.object(benchmark, "ROOT", root),
                    mock.patch.object(benchmark, "DATA", data),
                    mock.patch.object(benchmark, "FIGURES", figures),
                    mock.patch.object(
                        benchmark,
                        "read_validation_index",
                        return_value=validation_index,
                    ),
                    mock.patch.object(
                        benchmark,
                        "validated_gxtb_output_paths",
                        return_value=verified_paths,
                    ),
                ):
                    result = benchmark.analyse(
                        gxtb_run_root=root / runner.GXTB_RUN_DIRECTORY,
                        output_prefix=benchmark.GXTB_PRODUCTION_PREFIX,
                        validation_index_path=index_path,
                        validation_index_sha256=index_hash,
                    )

                expected_coverage: dict[str, list[str]] = {}
                for mesh, phase, _ in entries:
                    expected_coverage.setdefault(mesh, []).append(phase)
                expected_coverage = {
                    mesh: [
                        phase for phase in benchmark.PHASES
                        if phase in phases
                    ]
                    for mesh, phases in expected_coverage.items()
                }
                self.assertEqual(
                    result["validated_gxtb_phases"], expected_coverage
                )
                self.assertEqual(result["validated_gxtb_meshes"], [])
                for mesh, phase, energy in entries:
                    method = result["results"][mesh]["GXTB"]
                    self.assertFalse(method["complete"])
                    self.assertEqual(method["energies_hartree"][phase], energy)
                    for missing_phase in benchmark.PHASES:
                        if missing_phase not in expected_coverage[mesh]:
                            self.assertIsNone(
                                method["energies_hartree"][missing_phase]
                            )


if __name__ == "__main__":
    unittest.main()
