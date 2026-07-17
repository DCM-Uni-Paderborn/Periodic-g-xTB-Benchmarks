from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parent
TRUE = Path(shutil.which("true") or "/usr/bin/true").resolve()
os.environ.setdefault("CP2K_EXE", str(TRUE))
os.environ.setdefault("CP2K_LIB", str(TRUE))
os.environ.setdefault("MPIEXEC_EXE", str(TRUE))


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load("exact_binding_runner", ROOT / "run_test_matrix.py")
verifier = load("exact_binding_verifier", ROOT / "verify_test_matrix.py")


class ExactBindingTests(unittest.TestCase):
    @staticmethod
    def rank_sample(
        *,
        pid: int = 100,
        rank: int | None = 0,
        cpu: str = "96",
        starttime: int = 424242,
        observation_status: str = "explicit",
        identity_status: str = "stable",
        state: str = "R (running)",
        stat_state: str = "R",
        snapshot_status: str | None = None,
    ) -> dict:
        if snapshot_status is None:
            snapshot_status = {
                "stable": "consistent",
                "terminal_state": "consistent",
                "disappeared_after_sample": "process_disappeared",
                "pid_reused_during_sample": "pid_reused",
                "identity_unreadable_after_sample": "final_identity_unreadable",
                "executable_changed_during_sample": "executable_changed",
                "cpu_mask_changed_during_sample": "cpu_mask_changed",
            }.get(identity_status, "consistent")
        return {
            "pid": pid,
            "rank": rank,
            "raw_rank": rank,
            "cpus_allowed_list": cpu,
            "process_starttime": starttime,
            "process_identity_status": identity_status,
            "snapshot_consistency_status": snapshot_status,
            "rank_observation_status": observation_status,
            "state": state,
            "stat_state": stat_state,
            "executable": str(runner.CP2K),
            "arguments": [str(runner.CP2K)],
            "is_cp2k_rank": True,
        }

    @staticmethod
    def metadata(*children: dict) -> dict:
        return {
            "cp2k": str(runner.CP2K),
            "observed_child_processes": list(children),
            "concurrent_duplicate_rank_samples": [],
            "concurrent_duplicate_rank_ids_ever": [],
            "concurrent_duplicate_rank_processes_ever": False,
            "live_compute_overlap_runtime_gate": True,
            "live_compute_overlap_runtime_samples": [],
        }

    def test_fixed_matrix_has_48_runs(self) -> None:
        self.assertEqual(len(runner.jobs()), 48)
        self.assertEqual(runner.INPUT_ROOT, ROOT / "test_inputs")
        self.assertTrue(
            all((runner.INPUT_ROOT / case["input"]).is_file() for case in runner.MATRIX["cases"])
        )
        self.assertEqual(verifier.SUMMARY_PATH.parent, ROOT)
        self.assertNotEqual(verifier.SUMMARY_PATH.parent, verifier.RUN_ROOT)

    def test_ordered_reservation_is_literal_and_unique(self) -> None:
        self.assertEqual(runner.parse_ordered_pe_list("96,97"), (96, 97))
        for invalid in ("96-97", "96,96", "96,,97"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    runner.parse_ordered_pe_list(invalid)

    def test_bad_then_correct_rank_sample_remains_failed(self) -> None:
        bad = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=10, cpu="96-97"),
            (96, 97),
        )
        corrected = runner.accumulate_rank_snapshot(
            bad,
            self.rank_sample(pid=10, cpu="96"),
            (96, 97),
        )
        self.assertTrue(corrected["current_sample_matches_assigned_singleton"])
        self.assertTrue(corrected["affinity_violation_ever"])
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof({10: corrected}, 1, (96,))

    def test_openmpi_binding_report_is_rank_complete(self) -> None:
        text = (
            "[terok:1] Rank 0 bound to package 1[core 0]\n"
            "[terok:1] Rank 1 bound to package 1[core 1]\n"
        )
        self.assertEqual(runner.reported_binding_rank_ids(text), [0, 1])
        self.assertEqual(verifier.reported_binding_rank_ids(text), [0, 1])

    def test_cross_process_cpu_lock_rejects_overlap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lock_root = Path(tmp)
            first = runner.acquire_cpu_locks((96,), lock_root)
            try:
                with self.assertRaisesRegex(RuntimeError, "already reserved"):
                    runner.acquire_cpu_locks((96,), lock_root)
            finally:
                for handle in first:
                    handle.close()

    def test_lock_metadata_baseexception_releases_current_handle(self) -> None:
        class InjectedMetadataFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            lock_root = Path(tmp)
            retained: list[BaseException] = []
            with mock.patch.object(
                runner.json,
                "dump",
                side_effect=InjectedMetadataFailure("injected json.dump"),
            ):
                try:
                    runner.acquire_cpu_locks((1000001,), lock_root)
                except InjectedMetadataFailure as error:
                    retained.append(error)
                else:
                    self.fail("injected lock-metadata BaseException was swallowed")
            self.assertIsNotNone(retained[0].__traceback__)
            handles = runner.acquire_cpu_locks((1000001,), lock_root)
            for handle in handles:
                handle.close()
            retained.clear()

    def test_live_overlap_preflight_catches_cp2k_without_rank_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
            (process / "stat").write_text(
                "4242 (cp2k.psmp) "
                + " ".join(["R", *("0" for _ in range(18)), "424242"])
            )
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "Cpus_allowed_list:\t96\n"
            )
            (process / "environ").write_bytes(b"")
            with self.assertRaisesRegex(RuntimeError, "PID 4242.*overlaps"):
                runner.require_no_live_compute_overlap((96,), proc_root)
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tZ (zombie)\n"
                "Cpus_allowed_list:\t96\n"
            )
            runner.require_no_live_compute_overlap((96,), proc_root)

    def test_rank_snapshot_records_executable_and_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
            (process / "stat").write_text(
                "4242 (cp2k.psmp) "
                + " ".join(["R", *("0" for _ in range(18)), "424242"])
            )
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "Cpus_allowed_list:\t96\n"
            )
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=0\0")
            arguments = [str(runner.CP2K), "-i", "/tmp/example.inp"]
            (process / "cmdline").write_bytes(
                b"\0".join(argument.encode() for argument in arguments) + b"\0"
            )
            (process / "exe").symlink_to(runner.CP2K)
            snapshot = runner.rank_snapshot(4242, proc_root)
            self.assertIsNotNone(snapshot)
            self.assertEqual(snapshot["executable"], str(runner.CP2K))
            self.assertEqual(snapshot["arguments"], arguments)

            (process / "cmdline").unlink()
            self.assertIsNone(runner.rank_snapshot(4242, proc_root))

    def test_runtime_owner_exclusion_revalidates_pid_starttime(self) -> None:
        def stat(starttime: int) -> str:
            return "4242 (cp2k.psmp) " + " ".join(
                ["R", *("0" for _ in range(18)), str(starttime)]
            )

        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
            stat_path = process / "stat"
            stat_path.write_text(stat(111))
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "Cpus_allowed_list:\t96\n"
            )
            (process / "environ").write_bytes(b"")
            self.assertEqual(
                runner.live_compute_cpu_owners(
                    (96,),
                    proc_root,
                    ignore_process_identities={4242: 111},
                ),
                [],
            )
            self.assertEqual(
                runner.live_compute_cpu_owners(
                    (96,),
                    proc_root,
                    ignore_process_identities={4242: 7},
                )[0]["process_identity_status"],
                "stable",
            )

            real_read_text = Path.read_text
            stat_reads = iter((stat(111), stat(222)))

            def reused_during_scan(path: Path, *args, **kwargs):
                if path == stat_path:
                    return next(stat_reads)
                return real_read_text(path, *args, **kwargs)

            with mock.patch.object(Path, "read_text", reused_during_scan):
                owners = runner.live_compute_cpu_owners(
                    (96,),
                    proc_root,
                    ignore_process_identities={4242: 111},
                )
            self.assertEqual(
                owners[0]["process_identity_status"],
                "pid_reused_during_scan",
            )
            self.assertEqual(owners[0]["overlap"], [96])

    def test_only_one_rank_pid_generation_is_scaling_eligible(self) -> None:
        parent = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=10, cpu="96", starttime=1000),
            (96, 97),
        )
        runner.resolve_rank_process_lifetime(parent, "process_disappeared")
        successor = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=11, cpu="96", starttime=1001),
            (96, 97),
        )
        runner.resolve_rank_process_lifetime(successor, "process_disappeared")
        single = runner.ordered_rank_proof({10: parent}, 1, (96,))
        self.assertEqual(single[0]["pid_generations"], [10])
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof({10: parent, 11: successor}, 1, (96,))

        wrong = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=12, cpu="97", starttime=1002),
            (96, 97),
        )
        runner.resolve_rank_process_lifetime(wrong, "process_disappeared")
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof({10: parent, 12: wrong}, 1, (96,))
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof(
                {10: parent, 11: successor}, 1, (96,), {0}
            )

    def test_verifier_reconstructs_generation_and_duplicate_summaries(self) -> None:
        child0 = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=10, rank=0, cpu="96", starttime=1000),
            (96, 97),
        )
        child1 = runner.accumulate_rank_snapshot(
            None,
            self.rank_sample(pid=11, rank=1, cpu="97", starttime=1001),
            (96, 97),
        )
        runner.resolve_rank_process_lifetime(child0, "process_disappeared")
        runner.resolve_rank_process_lifetime(child1, "process_disappeared")
        metadata = self.metadata(child0, child1)
        proof = verifier.revalidated_rank_evidence(metadata, 2, (96, 97))
        self.assertEqual(
            [item["pid_generations"] for item in proof], [[10], [11]]
        )

        broken = dict(metadata)
        broken["concurrent_duplicate_rank_ids_ever"] = [0]
        with self.assertRaisesRegex(RuntimeError, "concurrent-rank summary"):
            verifier.revalidated_rank_evidence(broken, 2, (96, 97))

        broken = {
            **metadata,
            "observed_child_processes": [
                {**child0, "observed_cpu_masks": ["97", "96"]},
                child1,
            ],
        }
        with self.assertRaisesRegex(RuntimeError, "child history"):
            verifier.revalidated_rank_evidence(broken, 2, (96, 97))

    def test_terminal_environment_loss_is_deferred_and_revalidated(self) -> None:
        expected = (96,)
        first = runner.accumulate_rank_snapshot(
            None, self.rank_sample(), expected
        )
        pending = runner.accumulate_rank_snapshot(
            first,
            self.rank_sample(
                rank=None, observation_status="environment_empty"
            ),
            expected,
        )
        self.assertIs(pending["rank_environment_unavailable_pending"], True)
        self.assertIs(pending["affinity_violation_ever"], False)
        runner.ordered_rank_proof(
            {100: pending}, 1, expected, final=False
        )
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof(
                {100: pending}, 1, expected, final=True
            )

        second = runner.accumulate_rank_snapshot(
            pending,
            self.rank_sample(
                rank=None, observation_status="environment_unreadable"
            ),
            expected,
        )
        runner.resolve_pending_rank_environment(second, "process_disappeared")
        proof = runner.ordered_rank_proof(
            {100: second}, 1, expected, final=True
        )
        metadata = self.metadata(second)
        self.assertEqual(
            verifier.revalidated_rank_evidence(metadata, 1, expected), proof
        )

        duplicate = runner.accumulate_rank_snapshot(
            None, self.rank_sample(pid=101), expected
        )
        self.assertEqual(
            runner.concurrent_rank_pid_groups([second, duplicate]),
            {0: [100, 101]},
        )

    def test_terminal_environment_loss_rejects_anomalies_and_tampering(self) -> None:
        expected = (96,)
        first = runner.accumulate_rank_snapshot(
            None, self.rank_sample(), expected
        )
        cases = {
            "initial-loss": runner.accumulate_rank_snapshot(
                None,
                self.rank_sample(
                    rank=None, observation_status="environment_empty"
                ),
                expected,
            ),
            "initial-pid-reuse": runner.accumulate_rank_snapshot(
                None,
                self.rank_sample(identity_status="pid_reused_during_sample"),
                expected,
            ),
            "initial-identity-unreadable": runner.accumulate_rank_snapshot(
                None,
                self.rank_sample(
                    identity_status="identity_unreadable_after_sample"
                ),
                expected,
            ),
            "changed-starttime": runner.accumulate_rank_snapshot(
                first,
                self.rank_sample(
                    rank=None,
                    starttime=424243,
                    observation_status="environment_empty",
                ),
                expected,
            ),
            "changed-mask": runner.accumulate_rank_snapshot(
                first,
                self.rank_sample(
                    rank=None,
                    cpu="96-97",
                    observation_status="environment_empty",
                ),
                expected,
            ),
            "missing-rank": runner.accumulate_rank_snapshot(
                first,
                self.rank_sample(rank=None, observation_status="explicit_missing"),
                expected,
            ),
            "invalid-rank": runner.accumulate_rank_snapshot(
                first,
                self.rank_sample(rank=None, observation_status="explicit_invalid"),
                expected,
            ),
        }
        for name, record in cases.items():
            with self.subTest(name=name):
                self.assertIs(record["affinity_violation_ever"], True)
                self.assertIs(
                    record.get("rank_environment_unavailable_pending"), False
                )

        pending = runner.accumulate_rank_snapshot(
            first,
            self.rank_sample(rank=None, observation_status="environment_empty"),
            expected,
        )
        reappeared = runner.accumulate_rank_snapshot(
            pending, self.rank_sample(), expected
        )
        self.assertIs(reappeared["rank_identity_changed_ever"], True)
        self.assertIs(reappeared["affinity_violation_ever"], True)

        runner.resolve_pending_rank_environment(pending, "process_disappeared")
        base_metadata = self.metadata(pending)
        verifier.revalidated_rank_evidence(base_metadata, 1, expected)
        for name, mutation in (
            ("starttime", {"process_starttime": 7}),
            (
                "status",
                {"observed_rank_observation_statuses": ["explicit"]},
            ),
            (
                "resolution",
                {"rank_environment_terminal_confirmation": "terminal_state_Z"},
            ),
        ):
            with self.subTest(tamper=name):
                broken_child = {**pending, **mutation}
                with self.assertRaisesRegex(
                    RuntimeError, "rank-environment|process provenance"
                ):
                    verifier.revalidated_rank_evidence(
                        {**base_metadata, "observed_child_processes": [broken_child]},
                        1,
                        expected,
                    )
        broken_event = {
            **pending,
            "rank_environment_events": [
                {**pending["rank_environment_events"][0], "sample_index": 1}
            ],
        }
        with self.assertRaisesRegex(RuntimeError, "rank-environment"):
            verifier.revalidated_rank_evidence(
                {**base_metadata, "observed_child_processes": [broken_event]},
                1,
                expected,
            )

    def test_terminal_resolver_does_not_accept_live_or_reused_pid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp)
            self.assertEqual(
                runner.process_terminal_resolution(100, 424242, proc_root),
                "process_disappeared",
            )
            process = proc_root / "100"
            process.mkdir()

            def process_stat(state: str, starttime: int) -> str:
                return "100 (cp2k) " + " ".join(
                    [state, *("0" for _ in range(18)), str(starttime)]
                )

            (process / "stat").write_text(process_stat("R", 424242))
            self.assertIsNone(
                runner.process_terminal_resolution(100, 424242, proc_root)
            )
            tracked = runner.accumulate_rank_snapshot(
                None, self.rank_sample(), (96,)
            )
            self.assertIs(
                runner.observed_rank_process_is_still_live(
                    100, tracked, proc_root
                ),
                True,
            )
            self.assertIs(tracked["snapshot_unavailable_ever"], True)
            self.assertIs(tracked["affinity_violation_ever"], True)
            (process / "stat").write_text(process_stat("R", 424243))
            self.assertEqual(
                runner.process_terminal_resolution(100, 424242, proc_root),
                "pid_reused",
            )
            self.assertIs(
                runner.observed_rank_process_is_still_live(
                    100, tracked, proc_root
                ),
                False,
            )
            self.assertIs(tracked["process_starttime_changed_ever"], True)
            (process / "stat").write_text(process_stat("Z", 424242))
            self.assertEqual(
                runner.process_terminal_resolution(100, 424242, proc_root),
                "terminal_state_Z",
            )

    def test_verifier_rejects_missing_bool_aliased_and_mismatched_provenance(
        self,
    ) -> None:
        explicit = runner.accumulate_rank_snapshot(
            None, self.rank_sample(), (96,)
        )
        runner.resolve_rank_process_lifetime(explicit, "process_disappeared")
        metadata = self.metadata(explicit)
        self.assertEqual(
            verifier.revalidated_rank_evidence(metadata, 1, (96,))[0][
                "pid_generations"
            ],
            [100],
        )

        mandatory_process_fields = (
            "raw_rank",
            "rank_identity_source",
            "executable",
            "arguments",
            "process_starttime",
            "observed_process_starttimes",
            "process_identity_status",
            "snapshot_consistency_status",
            "state",
            "observed_process_states",
            "process_terminally_confirmed",
            "process_terminal_confirmation",
        )
        mandatory_environment_fields = (
            "rank_environment_unavailable_ever",
            "rank_environment_unavailable_sample_count",
            "rank_environment_unavailable_pending",
            "rank_environment_terminally_confirmed",
            "rank_environment_terminal_confirmation",
            "rank_environment_events",
        )
        for field in (*mandatory_process_fields, *mandatory_environment_fields):
            with self.subTest(missing=field):
                broken = dict(explicit)
                del broken[field]
                with self.assertRaises(RuntimeError):
                    verifier.revalidated_rank_evidence(
                        self.metadata(broken), 1, (96,)
                    )

        for name, mutation in (
            ("raw-rank-bool", {"raw_rank": False}),
            ("rank-history-bool", {"observed_rank_ids": [False]}),
            (
                "starttime-history-bool",
                {
                    "process_starttime": 1,
                    "observed_process_starttimes": [True],
                },
            ),
            ("arguments-bool", {"arguments": [True]}),
            ("arguments-mismatch", {"arguments": ["/tmp/not-cp2k"]}),
            ("wrong-executable", {"executable": "/tmp/not-cp2k"}),
        ):
            with self.subTest(fuzz=name):
                with self.assertRaises(RuntimeError):
                    verifier.revalidated_rank_evidence(
                        self.metadata({**explicit, **mutation}), 1, (96,)
                    )

        first = runner.accumulate_rank_snapshot(
            None, self.rank_sample(), (96,)
        )
        pending = runner.accumulate_rank_snapshot(
            first,
            self.rank_sample(
                rank=None, observation_status="environment_empty"
            ),
            (96,),
        )
        runner.resolve_pending_rank_environment(pending, "process_disappeared")
        mismatched = dict(pending)
        mismatched["rank_environment_terminal_confirmation"] = "terminal_state_Z"
        mismatched["rank_environment_events"] = [
            {**event, "terminal_resolution": "terminal_state_Z"}
            for event in pending["rank_environment_events"]
        ]
        with self.assertRaisesRegex(
            RuntimeError, "rank-environment|process provenance"
        ):
            verifier.revalidated_rank_evidence(
                self.metadata(mismatched), 1, (96,)
            )

        bad_event = dict(pending)
        bad_event["rank_environment_events"] = [
            {**event, "pid": True, "process_starttime": True}
            for event in pending["rank_environment_events"]
        ]
        with self.assertRaisesRegex(RuntimeError, "rank-environment"):
            verifier.revalidated_rank_evidence(
                self.metadata(bad_event), 1, (96,)
            )

    def test_checked_run_rejects_runtime_overlap_gate_and_samples(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_root = root / "inputs"
            run_root = root / "runs"
            input_root.mkdir()
            case = {
                "name": "tiny",
                "ranks": 1,
                "input": "tiny.inp",
                "nfull": 1,
            }
            input_path = (input_root / case["input"]).resolve()
            input_path.write_text("&GLOBAL\n&END GLOBAL\n")
            run_dir = run_root / "tiny_p1_dense"
            run_dir.mkdir(parents=True)
            output_path = run_dir / "cp2k.out"
            output_path.write_text(
                "GXTB-KGROUP-PARTIAL-ROOT groups=1, nred=1, "
                "nfull=1, batch=1;\nPROGRAM ENDED\n"
            )
            (run_dir / "returncode.txt").write_text("0\n")
            launcher_log = run_dir / "launcher.log"
            launcher_log.write_text("Rank 0 bound\n")

            child = runner.accumulate_rank_snapshot(
                None, self.rank_sample(), (96,)
            )
            runner.resolve_rank_process_lifetime(
                child, "process_disappeared"
            )
            affinity = runner.ordered_rank_proof(
                {int(child["pid"]): child}, 1, (96,)
            )
            launcher = TRUE
            cp2k = TRUE
            command = [
                str(launcher),
                "--map-by",
                "pe-list=96:ordered",
                "--bind-to",
                "core",
                "--report-bindings",
                "-np",
                "1",
                str(cp2k),
                "-i",
                str(input_path),
            ]
            metadata = {
                **self.metadata(child),
                "schema_version": 2,
                "case": case["name"],
                "ranks": case["ranks"],
                "variant": "DENSE",
                "returncode": 0,
                "input": str(input_path),
                "input_sha256": verifier.sha256(input_path),
                "working_directory": str(run_dir.resolve()),
                "cp2k": str(cp2k),
                "cp2k_sha256": verifier.sha256(cp2k),
                "cp2k_lib": str(cp2k),
                "cp2k_lib_sha256": verifier.sha256(cp2k),
                "ordered_pe_list": "96",
                "affinity_proof": affinity,
                "all_observed_rank_samples_match_ordered_pe_list": True,
                "runtime_affinity_gate": True,
                "cross_process_cpu_reservation_gate": True,
                "live_compute_overlap_preflight_gate": True,
                "live_compute_overlap_preflight_owners": [],
                "observed_cp2k_rank_pid_generations": [[100]],
                "observed_cp2k_process_generation_count": 1,
                "mpi_launcher": str(launcher),
                "mpi_launcher_sha256": verifier.sha256(launcher),
                "command": command,
                "launcher_log": launcher_log.name,
                "launcher_log_sha256": verifier.sha256(launcher_log),
                "reported_binding_rank_ids": [0],
                "timing_classification": "production_scaling_eligible",
                "output_sha256": verifier.sha256(output_path),
            }
            record_path = run_dir / "run.json"

            def write_record(record: dict) -> None:
                record_path.write_text(json.dumps(record, sort_keys=True))

            with (
                mock.patch.object(verifier, "INPUT_ROOT", input_root),
                mock.patch.object(verifier, "RUN_ROOT", run_root),
            ):
                write_record(metadata)
                verifier.checked_run(case, "DENSE")
                for name, mutation in (
                    (
                        "false-gate",
                        {"live_compute_overlap_runtime_gate": False},
                    ),
                    (
                        "recorded-sample",
                        {
                            "live_compute_overlap_runtime_samples": [
                                {
                                    "sample_index": 1,
                                    "owners": [{"pid": 1234}],
                                }
                            ]
                        },
                    ),
                ):
                    with self.subTest(mutation=name):
                        write_record({**metadata, **mutation})
                        with self.assertRaisesRegex(
                            RuntimeError, "runtime live CPU-overlap"
                        ):
                            verifier.checked_run(case, "DENSE")

    def test_all_inherited_mpi_mca_controls_are_scrubbed(self) -> None:
        self.assertEqual(
            runner.mpi_control_environment_keys(
                {
                    "OMPI_MCA_mca_base_param_files": "/tmp/indirect.conf",
                    "OMPI_MCA_pml": "ucx",
                    "PRTE_MCA_rmaps_default_mapping_policy": "core",
                    "PATH": "/bin",
                }
            ),
            [
                "OMPI_MCA_mca_base_param_files",
                "OMPI_MCA_pml",
                "PRTE_MCA_rmaps_default_mapping_policy",
            ],
        )

    @unittest.skipUnless(Path("/proc").is_dir(), "Linux /proc is required")
    def test_term_timeout_kills_and_reaps_process_group(self) -> None:
        process = subprocess.Popen(
            [
                sys.executable,
                "-c",
                (
                    "import signal,subprocess,time;"
                    "signal.signal(signal.SIGTERM,signal.SIG_IGN);"
                    "subprocess.Popen(['sleep','60']);"
                    "time.sleep(60)"
                ),
            ],
            start_new_session=True,
        )
        process_group = process.pid
        time.sleep(0.05)
        runner.terminate_and_reap_process_group(process, term_timeout=0.05)
        self.assertIsNotNone(process.returncode)
        self.assertEqual(runner.live_process_group_members(process_group), set())

    @unittest.skipUnless(Path("/proc").is_dir(), "Linux /proc is required")
    def test_cleanup_kills_tracked_rank_that_escaped_launcher_group(self) -> None:
        launcher = subprocess.Popen(["sleep", "60"], start_new_session=True)
        rank = subprocess.Popen(["sleep", "60"], start_new_session=True)
        try:
            rank_starttime = runner.linux_process_starttime(rank.pid)
            self.assertIsNotNone(rank_starttime)
            runner.terminate_and_reap_process_group(
                launcher,
                term_timeout=0.05,
                tracked_rank_starttimes={rank.pid: int(rank_starttime)},
            )
            rank.wait(timeout=1.0)
            self.assertIsNotNone(launcher.returncode)
            self.assertIsNotNone(rank.returncode)
            self.assertEqual(
                runner.live_process_group_members(launcher.pid), set()
            )
            self.assertEqual(runner.live_process_group_members(rank.pid), set())
        finally:
            for process in (launcher, rank):
                if process.poll() is None:
                    process.kill()
                process.wait()

    @unittest.skipUnless(Path("/proc").is_dir(), "Linux /proc is required")
    def test_run_one_baseexception_drains_launcher_and_escaped_rank(self) -> None:
        class InjectedMonitorFailure(BaseException):
            pass

        spawned: list[subprocess.Popen] = []

        def fail_after_spawn(job, slot, lifecycle):
            launcher = subprocess.Popen(
                ["sleep", "60"], start_new_session=True
            )
            rank = subprocess.Popen(["sleep", "60"], start_new_session=True)
            spawned.extend((launcher, rank))
            rank_starttime = runner.linux_process_starttime(rank.pid)
            self.assertIsNotNone(rank_starttime)
            lifecycle["process"] = launcher
            lifecycle["observed"] = {
                rank.pid: {
                    "is_cp2k_rank": True,
                    "process_starttime": int(rank_starttime),
                }
            }
            raise InjectedMonitorFailure("after Popen")

        try:
            with mock.patch.object(
                runner, "_run_one_inner", side_effect=fail_after_spawn
            ):
                with self.assertRaisesRegex(
                    InjectedMonitorFailure, "after Popen"
                ):
                    runner.run_one(({}, 1, "DENSE"), 0)
            for process in spawned:
                process.wait(timeout=1.0)
                self.assertIsNotNone(process.returncode)
                self.assertEqual(
                    runner.live_process_group_members(process.pid), set()
                )
        finally:
            for process in spawned:
                if process.poll() is None:
                    process.kill()
                process.wait()

    def test_new_launcher_contains_only_exact_ordered_binding(self) -> None:
        source = (ROOT / "run_test_matrix.py").read_text()
        self.assertNotIn('"taskset", "-c"', source)
        self.assertNotIn('"--bind-to", "none"', source)
        self.assertIn('f"pe-list={pe_list}:ordered"', source)
        self.assertIn('"--bind-to", "core", "--report-bindings"', source)
        self.assertIn("cwd=run_dir", source)


if __name__ == "__main__":
    unittest.main()
