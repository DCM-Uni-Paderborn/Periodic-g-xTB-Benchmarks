from __future__ import annotations

import inspect
import json
import os
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))

import benchmark_execution as execution  # noqa: E402


class BenchmarkExecutionValidationTests(unittest.TestCase):
    @staticmethod
    def rank_sample(
        *,
        pid: int = 100,
        rank: int | None = 1,
        cpu: str = "197",
        starttime: int = 424242,
        observation_status: str = "explicit",
        identity_status: str = "stable",
        state: str = "R (running)",
        stat_state: str = "R",
    ) -> dict[str, object]:
        thread_process_status = {
            "terminal_state": "terminal_process",
            "disappeared_after_sample": "process_disappeared_after_sample",
        }.get(identity_status, "live")
        return {
            "pid": pid,
            "is_cp2k_rank": True,
            "ompi_comm_world_rank": rank,
            "cpus_allowed_list": cpu,
            "process_starttime": starttime,
            "process_identity_status": identity_status,
            "snapshot_consistency_status": "consistent",
            "rank_observation_status": observation_status,
            "state": state,
            "stat_state": stat_state,
            "executable": "/tmp/cp2k.psmp",
            "arguments": ["/tmp/cp2k.psmp"],
            "thread_affinity_evidence_source": (
                execution.THREAD_AFFINITY_EVIDENCE_SOURCE
            ),
            "thread_affinity_scan_status": "consistent",
            "thread_affinity_scan_issues": [],
            "thread_affinity_process_status": thread_process_status,
            "live_thread_affinity": [
                {
                    "tid": pid,
                    "thread_starttime": starttime,
                    "observed_thread_starttimes": [starttime],
                    "state": state,
                    "stat_state": stat_state,
                    "cpus_allowed_list": cpu,
                    "observed_cpu_masks": [cpu],
                    "live": stat_state not in {"Z", "X"},
                    "identity_status": "stable",
                }
            ],
        }

    def valid_record_fixture(
        self, root: Path
    ) -> tuple[
        Path,
        Path,
        Path,
        Path,
        dict[str, object],
        dict[str, object],
    ]:
        cp2k = root / "bin" / "cp2k.psmp"
        inp = root / "inputs" / "job.inp"
        out = root / "outputs" / "job.out"
        stamp = root / "outputs" / "job.out.job.json"
        for parent in (cp2k.parent, inp.parent, out.parent):
            parent.mkdir(parents=True, exist_ok=True)
        cp2k.write_text("#!/bin/sh\nexit 0\n")
        cp2k.chmod(0o755)
        inp.write_text("input\n")
        out.write_text("output\n")
        signature = {
            "schema_version": 1,
            "executable": str(cp2k.resolve()),
            "executable_sha256": execution.sha256(cp2k),
            "input": str(inp.resolve()),
            "input_sha256": execution.sha256(inp),
            "command_contract": {"driver": "cp2k", "omp_threads": 1},
            "completed": True,
            "return_code": 0,
        }
        stamp.write_text(json.dumps(signature, indent=2, sort_keys=True) + "\n")
        contract: dict[str, object] = {
            "schema_version": 1,
            "mode": "taskset_mpi",
            "mpi_ranks_per_job": 2,
            "mpi_launcher": "/fixture/mpiexec",
            "mpi_launcher_args": ["--bind-to", "none"],
            "taskset": "/fixture/taskset",
            "cpu_sets": ["48,49,50,51"],
        }
        record: dict[str, object] = {
            "schema_version": 1,
            "contract": contract,
            "contract_sha256": execution.canonical_sha256(contract),
            "assigned_cpu_set": "48,49,50,51",
            "command": execution.cp2k_command(
                taskset="/fixture/taskset",
                cpu_set="48,49,50,51",
                mpi_launcher="/fixture/mpiexec",
                mpi_launcher_args=["--bind-to", "none"],
                mpi_ranks_per_job=2,
                cp2k=cp2k,
                inp=inp,
                out=out,
            ),
            "runtime_affinity_gate": True,
            "mpiexec_internal_rebinding_detected": False,
            "observed_cp2k_rank_pids": [101, 102],
            "observed_cp2k_rank_masks": ["48-51", "48-51"],
            "cp2k": str(cp2k.resolve()),
            "cp2k_sha256_at_launch": execution.sha256(cp2k),
            "input": str(inp.resolve()),
            "input_sha256_at_launch": execution.sha256(inp),
            "output": str(out.resolve()),
            "output_sha256": execution.sha256(out),
            "scientific_job_stamp": str(stamp.resolve()),
            "scientific_job_stamp_sha256": execution.sha256(stamp),
        }
        path = execution.execution_record_path(out)
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return cp2k, inp, out, stamp, contract, record

    def test_ordered_pe_lists_must_be_literal_exact_and_disjoint(self) -> None:
        self.assertEqual(execution.parse_ordered_pe_list("96,97,98,99"), (96, 97, 98, 99))
        for value, message in (
            ("96-99", "comma-separated"),
            ("96,97,96", "duplicate"),
            ("96,,97", "comma-separated"),
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, message):
                    execution.parse_ordered_pe_list(value)
        with self.assertRaisesRegex(ValueError, "overlap at"):
            execution.validate_pe_lists(
                ["96,97", "97,98"],
                concurrent_jobs=2,
                mpi_ranks_per_job=2,
                threads_per_rank=1,
            )

    def test_pe_list_length_and_thread_count_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected exactly 4"):
            execution.validate_pe_lists(
                ["96,97,98"],
                concurrent_jobs=1,
                mpi_ranks_per_job=4,
                threads_per_rank=1,
            )
        with self.assertRaisesRegex(ValueError, "exactly one OpenMP thread"):
            execution.validate_pe_lists(
                ["96,97"],
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=2,
            )
        with self.assertRaisesRegex(ValueError, "unavailable CPUs"):
            execution.validate_pe_lists(
                ["96,97"],
                concurrent_jobs=1,
                mpi_ranks_per_job=2,
                threads_per_rank=1,
                available_cpus={96},
            )

    def test_user_launcher_arguments_cannot_override_binding(self) -> None:
        execution.validate_mpi_launcher_args([])
        for arguments in (
            ["--mca", "pml", "ucx"],
            ["--bind-to", "none"],
            ["--bind-to=core"],
            ["--bind-to-socket"],
            ["--map-by", "core"],
            ["--cpu-list", "96,97"],
            ["--cpu-set", "96,97"],
            ["--rank-by", "core"],
            ["--cpus-per-rank", "1"],
            ["--use-hwthread-cpus"],
            ["--report-bindings"],
            ["-np", "4"],
            ["-c", "4"],
            ["--n", "4"],
            ["--rankfile", "ranks.txt"],
            ["--oversubscribe"],
            ["--prtemca", "mca_base_param_files", "unsafe.conf"],
            [":"],
            ["--"],
            ["/usr/bin/taskset"],
        ):
            with self.subTest(arguments=arguments):
                with self.assertRaises(ValueError):
                    execution.validate_mpi_launcher_args(arguments)

    def test_record_accepts_external_output_and_rejects_artifact_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, inp, out, stamp, contract, record = self.valid_record_fixture(root)
            path = execution.execution_record_path(out)
            self.assertNotEqual(inp.parent, out.parent)
            command = record["command"]
            self.assertIsInstance(command, list)
            assert isinstance(command, list)
            self.assertEqual(command[-1], str(out.resolve()))
            self.assertIsNone(
                execution.recorded_execution_issue(path, contract, out, stamp)
            )
            rebound = dict(record)
            rebound["observed_cp2k_rank_masks"] = ["48", "49"]
            rebound["mpiexec_internal_rebinding_detected"] = True
            rebound["runtime_affinity_gate"] = False
            path.write_text(json.dumps(rebound, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "runtime MPI/affinity gate failed",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            out.write_text("tampered\n")
            self.assertIn(
                "output hash mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )
            out.write_text("output\n")
            stamp.write_text('{"tampered": true}\n')
            self.assertIn(
                "job-stamp hash mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )

    def test_all_inherited_ompi_prte_mca_environment_is_removed(self) -> None:
        environment = {
            "OMPI_MCA_hwloc_base_binding_policy": "none",
            "OMPI_MCA_rmaps_future_override": "unsafe",
            "PRTE_MCA_hwloc_default_binding_policy": "none",
            "OMPI_MCA_mca_base_param_files": "/tmp/indirect-override.conf",
            "PRTE_MCA_mca_base_param_files": "/tmp/indirect-override.conf",
            "OMPI_MCA_pml": "ucx",
            "PATH": "/bin",
        }
        self.assertEqual(
            execution.mpi_control_environment_keys(environment),
            [
                "OMPI_MCA_hwloc_base_binding_policy",
                "OMPI_MCA_mca_base_param_files",
                "OMPI_MCA_pml",
                "OMPI_MCA_rmaps_future_override",
                "PRTE_MCA_hwloc_default_binding_policy",
                "PRTE_MCA_mca_base_param_files",
            ],
        )

    def test_poisoned_binding_and_preload_environment_is_hermetic(self) -> None:
        poisoned = {
            "OMPI_MCA_hwloc_base_binding_policy": "none",
            "PRTE_MCA_rmaps_default_mapping_policy": "slot",
            "PMIX_MCA_gds": "hash",
            "HWLOC_XMLFILE": "/tmp/foreign-topology.xml",
            "GOMP_CPU_AFFINITY": "0-255",
            "KMP_AFFINITY": "disabled",
            "OMP_NUM_THREADS": "99",
            "OMP_PLACES": "threads",
            "LD_PRELOAD": "/tmp/repin.so",
            "LD_AUDIT": "/tmp/audit.so",
            "I_MPI_PIN_PROCESSOR_LIST": "0",
            "MV2_ENABLE_AFFINITY": "0",
            "SLURM_CPU_BIND": "none",
            "LD_LIBRARY_PATH": "/pinned/cp2k/provider",
            "PATH": "/usr/bin:/bin",
        }
        sanitized, removed, library_path_preserved = (
            execution.sanitized_launch_environment(poisoned)
        )
        expected_removed = set(poisoned) - {"LD_LIBRARY_PATH", "PATH"}
        self.assertEqual(set(removed), expected_removed)
        self.assertTrue(library_path_preserved)
        self.assertEqual(
            sanitized["LD_LIBRARY_PATH"], "/pinned/cp2k/provider"
        )
        self.assertEqual(sanitized["PATH"], "/usr/bin:/bin")
        self.assertEqual(
            {key: sanitized[key] for key in execution.SANCTIONED_THREAD_ENVIRONMENT},
            execution.SANCTIONED_THREAD_ENVIRONMENT,
        )
        self.assertEqual(
            [
                key
                for key in execution.mpi_control_environment_keys(sanitized)
                if key not in execution.SANCTIONED_THREAD_ENVIRONMENT
            ],
            [],
        )
        self.assertNotIn("LD_PRELOAD", sanitized)
        self.assertNotIn("LD_AUDIT", sanitized)

    def test_openmpi5_binding_report_is_recognized_case_insensitively(self) -> None:
        text = (
            "[terok:48321] Rank 0 bound to package 1[core 36[hwt 0]]\n"
            "[terok:48321] Rank 1 bound to package 1[core 37[hwt 0]]\n"
        )
        self.assertEqual(execution._reported_binding_rank_ids(text), [0, 1])

    def test_rank_processes_are_ordered_by_ompi_rank_not_pid(self) -> None:
        observed = {
            101: {
                "pid": 101,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 1,
            },
            999: {
                "pid": 999,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
            },
        }
        ordered = execution._ordered_cp2k_rank_processes(observed, 2)
        self.assertEqual([item["pid"] for item in ordered], [999, 101])

    def test_affinity_violation_is_sticky_across_rank_samples(self) -> None:
        bad = execution._accumulate_process_snapshot(
            None,
            self.rank_sample(pid=10, rank=0, cpu="48-49"),
            (48, 49),
        )
        corrected = execution._accumulate_process_snapshot(
            bad,
            self.rank_sample(pid=10, rank=0, cpu="48"),
            (48, 49),
        )
        self.assertTrue(bad["affinity_violation_ever"])
        self.assertTrue(corrected["current_sample_matches_assigned_singleton"])
        self.assertTrue(corrected["affinity_violation_ever"])
        self.assertEqual(corrected["observed_cpu_masks"], ["48-49", "48"])

    def test_pool_contract_injects_ordered_core_binding_without_taskset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=2,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["96,97"],
                check_current_affinity=False,
                cpu_reservation_lock_root=root / "cpu-locks",
            )
            self.assertEqual(pool.pe_lists, ("96,97",))
            self.assertEqual(pool.contract["mpi_bind_to"], "core")
            self.assertEqual(pool.contract["outer_taskset"], False)
            self.assertEqual(pool.contract["exact_cpus_per_rank"], 1)
            self.assertEqual(
                pool.contract["cross_process_cpu_reservation"],
                "flock_per_logical_cpu",
            )
            self.assertEqual(
                pool.contract["pool_close_policy"],
                "reject_while_run_admitted_or_active",
            )
            self.assertEqual(
                pool.contract["binding_environment_scrub_contract"],
                execution.binding_environment_scrub_contract(),
            )
            pool.close()

    def test_cross_process_cpu_reservation_fails_closed_on_overlap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = execution.acquire_cpu_reservation_locks((48,), root)
            try:
                with self.assertRaisesRegex(ValueError, "already reserved"):
                    execution.acquire_cpu_reservation_locks((48,), root)
            finally:
                for handle in first:
                    handle.close()
            second = execution.acquire_cpu_reservation_locks((48,), root)
            for handle in second:
                handle.close()

    def test_close_rejects_while_run_is_active_and_preserves_reservations(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            started = root / "launcher.started"
            release = root / "launcher.release"
            launcher = root / "mpirun"
            launcher.write_text(
                f"#!{sys.executable}\n"
                "import time\n"
                "from pathlib import Path\n"
                f"Path({str(started)!r}).touch()\n"
                f"release = Path({str(release)!r})\n"
                "while not release.exists():\n"
                "    time.sleep(0.01)\n"
                "raise SystemExit(1)\n"
            )
            launcher.chmod(0o755)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            inp = root / "job.inp"
            inp.write_text("&GLOBAL\n&END GLOBAL\n")
            out = root / "job.out"
            lock_root = root / "cpu-locks"
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["1000031"],
                check_current_affinity=False,
                cpu_reservation_lock_root=lock_root,
            )
            outcome: list[object] = []

            def run() -> None:
                try:
                    outcome.append(pool.run_cp2k(cp2k, inp, out))
                except BaseException as error:
                    outcome.append(error)

            worker = threading.Thread(target=run)
            worker.start()
            deadline = time.monotonic() + 5.0
            while not started.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            try:
                self.assertTrue(started.exists())
                with self.assertRaisesRegex(RuntimeError, "cannot release"):
                    pool.close()
                self.assertFalse(pool._closed)
                self.assertTrue(
                    all(
                        not handle.closed for handle in pool._reservation_handles
                    )
                )
                with self.assertRaisesRegex(ValueError, "already reserved"):
                    execution.acquire_cpu_reservation_locks(
                        (1000031,), lock_root
                    )
            finally:
                release.touch()
                worker.join(timeout=5.0)
                if worker.is_alive():
                    self.fail("active run did not finish after launcher release")
                pool.close()
            self.assertEqual(len(outcome), 1)
            self.assertIsInstance(outcome[0], tuple)
            reacquired = execution.acquire_cpu_reservation_locks(
                (1000031,), lock_root
            )
            for handle in reacquired:
                handle.close()

    def test_queue_to_active_lock_baseexception_restores_exact_ownership(
        self,
    ) -> None:
        class InjectedTransitionFailure(BaseException):
            pass

        class FailSecondEntryLock:
            def __init__(self, lock) -> None:
                self.lock = lock
                self.entries = 0

            def __enter__(self):
                self.entries += 1
                if self.entries == 2:
                    raise InjectedTransitionFailure("injected transition lock")
                self.lock.acquire()
                return self

            def __exit__(self, exc_type, exc_value, traceback) -> None:
                self.lock.release()

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            inp = root / "job.inp"
            inp.write_text("&GLOBAL\n&END GLOBAL\n")
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["1000032"],
                check_current_affinity=False,
                cpu_reservation_lock_root=root / "cpu-locks",
            )
            pool._lifecycle_lock = FailSecondEntryLock(pool._lifecycle_lock)
            try:
                with self.assertRaisesRegex(
                    InjectedTransitionFailure, "transition lock"
                ):
                    pool.run_cp2k(cp2k, inp, root / "job.out")
                self.assertEqual(pool._admitted_runs, 0)
                self.assertEqual(pool._active, set())
                self.assertEqual(pool._available.qsize(), 1)
                restored = pool._available.get_nowait()
                self.assertEqual(restored, "1000032")
                pool._available.put_nowait(restored)
                self.assertTrue(
                    all(not handle.closed for handle in pool._reservation_handles)
                )
            finally:
                pool.close()

    def test_active_add_baseexception_restores_exact_ownership(self) -> None:
        class InjectedActiveFailure(BaseException):
            pass

        class AddThenFailSet(set):
            def add(self, item) -> None:
                super().add(item)
                raise InjectedActiveFailure("injected active.add")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            inp = root / "job.inp"
            inp.write_text("&GLOBAL\n&END GLOBAL\n")
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["1000033"],
                check_current_affinity=False,
                cpu_reservation_lock_root=root / "cpu-locks",
            )
            pool._active = AddThenFailSet()
            try:
                with self.assertRaisesRegex(
                    InjectedActiveFailure, "active.add"
                ):
                    pool.run_cp2k(cp2k, inp, root / "job.out")
                self.assertEqual(pool._admitted_runs, 0)
                self.assertEqual(pool._active, set())
                self.assertEqual(pool._available.qsize(), 1)
                restored = pool._available.get_nowait()
                self.assertEqual(restored, "1000033")
                pool._available.put_nowait(restored)
                self.assertTrue(
                    all(not handle.closed for handle in pool._reservation_handles)
                )
            finally:
                pool.close()

    def test_trace_baseexception_after_queue_get_does_not_leak_pe_list(
        self,
    ) -> None:
        class InjectedPostGetFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 1\n")
            launcher.chmod(0o755)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            inp = root / "job.inp"
            inp.write_text("&GLOBAL\n&END GLOBAL\n")
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["1000034"],
                check_current_affinity=False,
                cpu_reservation_lock_root=root / "cpu-locks",
            )
            source, first_line = inspect.getsourcelines(
                execution.ExecutionPool.run_cp2k
            )
            get_offset = next(
                offset
                for offset, line in enumerate(source)
                if "pe_list = self._available.get()" in line
            )
            first_post_get_line = first_line + get_offset + 1
            target_code = execution.ExecutionPool.run_cp2k.__code__
            injected = False

            def inject_after_get(frame, event, argument):
                nonlocal injected
                if (
                    not injected
                    and frame.f_code is target_code
                    and event == "line"
                    and frame.f_lineno == first_post_get_line
                    and frame.f_locals.get("pe_list") == "1000034"
                ):
                    injected = True
                    raise InjectedPostGetFailure(
                        "first line event after Queue.get return"
                    )
                return inject_after_get

            previous_trace = sys.gettrace()
            try:
                sys.settrace(inject_after_get)
                with self.assertRaisesRegex(
                    InjectedPostGetFailure, "after Queue.get"
                ):
                    pool.run_cp2k(cp2k, inp, root / "faulted.out")
            finally:
                sys.settrace(previous_trace)
            try:
                self.assertTrue(injected)
                self.assertEqual(pool._admitted_runs, 0)
                self.assertEqual(pool._available.qsize(), 1)
                self.assertEqual(pool._active, set())
                outcome: list[tuple[int, dict[str, object]]] = []

                def subsequent_run() -> None:
                    outcome.append(
                        pool.run_cp2k(cp2k, inp, root / "subsequent.out")
                    )

                worker = threading.Thread(target=subsequent_run)
                worker.start()
                worker.join(timeout=5.0)
                self.assertFalse(worker.is_alive(), "PE-list token was leaked")
                self.assertEqual(len(outcome), 1)
                self.assertEqual(
                    outcome[0][1]["assigned_ordered_pe_list"], "1000034"
                )
                self.assertEqual(pool._admitted_runs, 0)
                self.assertEqual(pool._available.qsize(), 1)
                self.assertEqual(pool._active, set())
            finally:
                pool.close()

    def test_pool_constructor_releases_reservations_after_late_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            lock_root = root / "cpu-locks"
            with mock.patch.object(
                execution,
                "canonical_sha256",
                side_effect=RuntimeError("injected late constructor failure"),
            ):
                with self.assertRaisesRegex(RuntimeError, "injected late"):
                    execution.ExecutionPool(
                        concurrent_jobs=1,
                        mpi_ranks_per_job=1,
                        threads_per_rank=1,
                        mpi_launcher=launcher,
                        mpi_launcher_args=[],
                        pe_lists=["1000000"],
                        check_current_affinity=False,
                        cpu_reservation_lock_root=lock_root,
                    )
            reacquired = execution.acquire_cpu_reservation_locks((1000000,), lock_root)
            for handle in reacquired:
                handle.close()

    def test_pool_releases_current_lock_after_metadata_baseexception(self) -> None:
        class InjectedMetadataFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            lock_root = root / "cpu-locks"
            retained_errors: list[BaseException] = []
            with mock.patch.object(
                execution.json,
                "dump",
                side_effect=InjectedMetadataFailure("injected json.dump"),
            ):
                try:
                    execution.ExecutionPool(
                        concurrent_jobs=1,
                        mpi_ranks_per_job=1,
                        threads_per_rank=1,
                        mpi_launcher=launcher,
                        mpi_launcher_args=[],
                        pe_lists=["1000001"],
                        check_current_affinity=False,
                        cpu_reservation_lock_root=lock_root,
                    )
                except InjectedMetadataFailure as error:
                    retained_errors.append(error)
                    frame_names: list[str] = []
                    traceback = error.__traceback__
                    while traceback is not None:
                        frame_names.append(traceback.tb_frame.f_code.co_name)
                        traceback = traceback.tb_next
                    self.assertIn("acquire_cpu_reservation_locks", frame_names)
                    self.assertIn("__init__", frame_names)
                else:
                    self.fail("injected lock-metadata BaseException was swallowed")
            self.assertIsNotNone(retained_errors[0].__traceback__)
            reacquired = execution.acquire_cpu_reservation_locks((1000001,), lock_root)
            for handle in reacquired:
                handle.close()
            retained_errors.clear()

    def test_pool_releases_locks_after_final_initialization_baseexception(self) -> None:
        class InjectedInitializationFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nexit 0\n")
            launcher.chmod(0o755)
            lock_root = root / "cpu-locks"
            real_lock = execution.threading.Lock
            calls = 0
            retained_errors: list[BaseException] = []

            def fail_second_lock():
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise InjectedInitializationFailure("injected second Lock")
                return real_lock()

            with mock.patch.object(execution.threading, "Lock", fail_second_lock):
                try:
                    execution.ExecutionPool(
                        concurrent_jobs=1,
                        mpi_ranks_per_job=1,
                        threads_per_rank=1,
                        mpi_launcher=launcher,
                        mpi_launcher_args=[],
                        pe_lists=["1000011"],
                        check_current_affinity=False,
                        cpu_reservation_lock_root=lock_root,
                    )
                except InjectedInitializationFailure as error:
                    retained_errors.append(error)
                    frame_names: list[str] = []
                    traceback = error.__traceback__
                    while traceback is not None:
                        frame_names.append(traceback.tb_frame.f_code.co_name)
                        traceback = traceback.tb_next
                    self.assertIn("__init__", frame_names)
                    self.assertIn("fail_second_lock", frame_names)
                else:
                    self.fail("injected final-initialization BaseException was swallowed")
            self.assertIsNotNone(retained_errors[0].__traceback__)
            reacquired = execution.acquire_cpu_reservation_locks((1000011,), lock_root)
            for handle in reacquired:
                handle.close()
            retained_errors.clear()

    def test_reservation_metadata_flush_and_fsync_baseexceptions_release_lock(
        self,
    ) -> None:
        class InjectedMetadataFailure(BaseException):
            pass

        class FlushFailingHandle:
            def __init__(self, handle):
                self._handle = handle

            def flush(self) -> None:
                raise InjectedMetadataFailure("injected flush")

            def __getattr__(self, name):
                return getattr(self._handle, name)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original_open = Path.open

            def flush_failing_open(path, *args, **kwargs):
                return FlushFailingHandle(original_open(path, *args, **kwargs))

            with mock.patch.object(Path, "open", flush_failing_open):
                with self.assertRaisesRegex(InjectedMetadataFailure, "flush"):
                    execution.acquire_cpu_reservation_locks((1000021,), root)
            reacquired = execution.acquire_cpu_reservation_locks((1000021,), root)
            for handle in reacquired:
                handle.close()

            with mock.patch.object(
                execution.os,
                "fsync",
                side_effect=InjectedMetadataFailure("injected fsync"),
            ):
                with self.assertRaises(InjectedMetadataFailure):
                    execution.acquire_cpu_reservation_locks((1000022,), root)
            reacquired = execution.acquire_cpu_reservation_locks((1000022,), root)
            for handle in reacquired:
                handle.close()

    def test_procfs_preflight_rejects_live_cp2k_overlap_and_ignores_zombie(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
            (process / "stat").write_text(
                "4242 (cp2k.psmp) "
                + " ".join(["R", *("0" for _ in range(18)), "424242"])
            )
            status = process / "status"
            status.write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "Cpus_allowed_list:\t48-49\n"
            )
            (process / "environ").write_bytes(b"")
            owners = execution.live_compute_cpu_owners((48,), proc_root)
            self.assertEqual(owners[0]["pid"], 4242)
            self.assertEqual(
                execution.live_compute_cpu_owners(
                    (48,),
                    proc_root,
                    ignore_process_identities={4242: 424242},
                ),
                [],
            )
            with self.assertRaisesRegex(ValueError, "PID 4242.*overlaps"):
                execution.require_no_live_compute_overlap((48,), proc_root)

            status.write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tZ (zombie)\n"
                "Cpus_allowed_list:\t48-49\n"
            )
            self.assertEqual(execution.live_compute_cpu_owners((48,), proc_root), [])

            status.write_text(
                "Name:\tother-solver\n"
                "State:\tS (sleeping)\n"
                "Cpus_allowed_list:\t48\n"
            )
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=3\0")
            self.assertEqual(
                execution.live_compute_cpu_owners((48,), proc_root)[0][
                    "mpi_rank_process"
                ],
                True,
            )

    def test_post_popen_baseexception_drains_group_before_reusing_pe_list(
        self,
    ) -> None:
        class InjectedMonitorFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            launcher = root / "mpirun"
            launcher.write_text("#!/bin/sh\nsleep 60 &\nwait\n")
            launcher.chmod(0o755)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            inp = root / "job.inp"
            out = root / "job.out"
            inp.write_text("&GLOBAL\n&END GLOBAL\n")
            pool = execution.ExecutionPool(
                concurrent_jobs=1,
                mpi_ranks_per_job=1,
                threads_per_rank=1,
                mpi_launcher=launcher,
                mpi_launcher_args=[],
                pe_lists=["1000041"],
                check_current_affinity=False,
                cpu_reservation_lock_root=root / "cpu-locks",
            )
            scan_count = 0
            cleaned_groups: list[int] = []
            real_cleanup = execution._terminate_and_reap_process_group

            def fail_second_owner_scan(*args, **kwargs):
                nonlocal scan_count
                scan_count += 1
                if scan_count == 1:
                    return []
                raise InjectedMonitorFailure("injected runtime monitor failure")

            def recording_cleanup(
                process,
                term_timeout=30.0,
                tracked_rank_starttimes=None,
            ):
                cleaned_groups.append(process.pid)
                return real_cleanup(
                    process,
                    term_timeout=0.2,
                    tracked_rank_starttimes=tracked_rank_starttimes,
                )

            try:
                with mock.patch.object(
                    execution,
                    "live_compute_cpu_owners",
                    side_effect=fail_second_owner_scan,
                ), mock.patch.object(
                    execution,
                    "_terminate_and_reap_process_group",
                    side_effect=recording_cleanup,
                ):
                    with self.assertRaisesRegex(
                        InjectedMonitorFailure, "runtime monitor"
                    ):
                        pool.run_cp2k(cp2k, inp, out)
                self.assertEqual(len(cleaned_groups), 1)
                self.assertEqual(
                    execution._live_process_group_members(cleaned_groups[0]), set()
                )
                self.assertEqual(pool._active, set())
                self.assertEqual(pool._available.qsize(), 1)
            finally:
                pool.close()

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
                "Cpus_allowed_list:\t48\n"
            )
            (process / "environ").write_bytes(b"")
            self.assertEqual(
                execution.live_compute_cpu_owners(
                    (48,),
                    proc_root,
                    ignore_process_identities={4242: 111},
                ),
                [],
            )
            self.assertEqual(
                execution.live_compute_cpu_owners(
                    (48,),
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
                owners = execution.live_compute_cpu_owners(
                    (48,),
                    proc_root,
                    ignore_process_identities={4242: 111},
                )
            self.assertEqual(
                owners[0]["process_identity_status"], "pid_reused_during_scan"
            )
            self.assertEqual(owners[0]["overlap"], [48])

    @unittest.skipUnless(Path("/proc").is_dir(), "Linux /proc is required")
    def test_cleanup_also_drains_tracked_rank_outside_launcher_group(self) -> None:
        launcher = execution.subprocess.Popen(["sleep", "60"], start_new_session=True)
        rank = execution.subprocess.Popen(["sleep", "60"], start_new_session=True)
        try:
            _, rank_starttime = execution._linux_proc_stat_identity(
                (Path("/proc") / str(rank.pid) / "stat").read_text()
            )
            execution._terminate_and_reap_process_group(
                launcher,
                term_timeout=0.2,
                tracked_rank_starttimes={rank.pid: rank_starttime},
            )
            rank.wait(timeout=1.0)
            self.assertIsNotNone(launcher.returncode)
            self.assertIsNotNone(rank.returncode)
            self.assertEqual(
                execution._live_process_group_members(launcher.pid), set()
            )
        finally:
            for process in (launcher, rank):
                if process.poll() is None:
                    process.kill()
                process.wait()

    def test_rank_generation_aggregation_rejects_successors_and_mask_change(
        self,
    ) -> None:
        parent = execution._accumulate_process_snapshot(
            None,
            self.rank_sample(pid=100, rank=0, cpu="48", starttime=1000),
            (48, 49),
        )
        execution._resolve_rank_process_lifetime(parent, "process_disappeared")
        successor = execution._accumulate_process_snapshot(
            None,
            self.rank_sample(pid=101, rank=0, cpu="48", starttime=1001),
            (48, 49),
        )
        execution._resolve_rank_process_lifetime(successor, "process_disappeared")
        exact = execution._aggregate_cp2k_rank_generations(
            {100: parent, 101: successor}, (48, 49)
        )
        self.assertEqual(exact[0]["pid_generations"], [100, 101])
        self.assertIs(exact[0]["all_samples_match_assigned_singleton"], False)
        single = execution._aggregate_cp2k_rank_generations({100: parent}, (48, 49))
        self.assertIs(single[0]["all_samples_match_assigned_singleton"], True)

        wrong_successor = execution._accumulate_process_snapshot(
            None,
            self.rank_sample(pid=102, rank=0, cpu="49", starttime=1002),
            (48, 49),
        )
        execution._resolve_rank_process_lifetime(
            wrong_successor, "process_disappeared"
        )
        changed = execution._aggregate_cp2k_rank_generations(
            {100: parent, 102: wrong_successor}, (48, 49)
        )
        self.assertIs(changed[0]["all_samples_match_assigned_singleton"], False)

    def test_concurrently_live_duplicate_rank_and_rank_migration_fail_sticky(
        self,
    ) -> None:
        snapshots = [
            {
                "pid": pid,
                "state": "R (running)",
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
            }
            for pid in (100, 101)
        ]
        self.assertEqual(
            execution._concurrent_live_duplicate_rank_ids(snapshots), {0}
        )
        first = execution._accumulate_process_snapshot(
            None,
            self.rank_sample(pid=100, rank=0, cpu="48"),
            (48, 49),
        )
        migrated = execution._accumulate_process_snapshot(
            first,
            self.rank_sample(pid=100, rank=1, cpu="49"),
            (48, 49),
        )
        self.assertIs(migrated["rank_identity_changed_ever"], True)
        self.assertIs(migrated["affinity_violation_ever"], True)

    def test_proc_stat_parser_and_stable_rank_snapshot(self) -> None:
        stat = "123 (cp2k) worker ) name) " + " ".join(
            ["R", *("0" for _ in range(18)), "424242"]
        )
        self.assertEqual(execution._linux_proc_stat_identity(stat), ("R", 424242))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            process = root / "proc" / "123"
            process.mkdir(parents=True)
            (process / "stat").write_text(stat)
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "PPid:\t1\n"
                "Cpus_allowed_list:\t197\n"
            )
            (process / "cmdline").write_bytes(str(cp2k).encode() + b"\0")
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=1\0")
            (process / "exe").symlink_to(cp2k)
            task = process / "task" / "123"
            task.mkdir(parents=True)
            (task / "stat").write_text(stat)
            (task / "status").write_text((process / "status").read_text())
            snapshot = execution._linux_process_snapshot(
                123, cp2k, root / "proc"
            )
            assert snapshot is not None
            self.assertEqual(snapshot["ompi_comm_world_rank"], 1)
            self.assertEqual(snapshot["rank_observation_status"], "explicit")
            self.assertEqual(snapshot["process_starttime"], 424242)
            self.assertEqual(snapshot["process_identity_status"], "stable")
            self.assertEqual(snapshot["snapshot_consistency_status"], "consistent")
            self.assertEqual(snapshot["thread_affinity_scan_status"], "consistent")
            self.assertTrue(execution._thread_affinity_sample_matches(snapshot, 197))

            initial_status = (process / "status").read_text()
            changed_status = initial_status.replace(
                "Cpus_allowed_list:\t197", "Cpus_allowed_list:\t198"
            )
            real_read_text = Path.read_text
            status_reads = iter((initial_status, changed_status))

            def changing_status(path: Path, *args, **kwargs):
                if path == process / "status":
                    return next(status_reads)
                return real_read_text(path, *args, **kwargs)

            with mock.patch.object(Path, "read_text", changing_status):
                changed_mask = execution._linux_process_snapshot(
                    123, cp2k, root / "proc"
                )
            assert changed_mask is not None
            self.assertEqual(
                changed_mask["process_identity_status"],
                "cpu_mask_changed_during_sample",
            )
            self.assertEqual(
                changed_mask["snapshot_consistency_status"], "cpu_mask_changed"
            )

            other = root / "other"
            other.write_text("#!/bin/sh\nexit 0\n")
            other.chmod(0o755)
            real_resolve = Path.resolve
            executable_reads = iter((cp2k.resolve(), other.resolve()))

            def changing_executable(path: Path, *args, **kwargs):
                if path == process / "exe":
                    return next(executable_reads)
                return real_resolve(path, *args, **kwargs)

            with mock.patch.object(Path, "resolve", changing_executable):
                changed_executable = execution._linux_process_snapshot(
                    123, cp2k, root / "proc"
                )
            assert changed_executable is not None
            self.assertEqual(
                changed_executable["process_identity_status"],
                "executable_changed_during_sample",
            )
            self.assertEqual(
                changed_executable["snapshot_consistency_status"],
                "executable_changed",
            )
            self.assertIs(changed_executable["is_cp2k_rank"], True)

    def test_repinned_child_thread_is_detected_and_remains_sticky(self) -> None:
        def stat(tid: int, starttime: int) -> str:
            return f"{tid} (cp2k-thread) " + " ".join(
                ["R", *("0" for _ in range(18)), str(starttime)]
            )

        def status(cpu: int) -> str:
            return (
                "Name:\tcp2k-thread\n"
                "State:\tR (running)\n"
                f"Cpus_allowed_list:\t{cpu}\n"
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            process = root / "proc" / "123"
            process.mkdir(parents=True)
            (process / "stat").write_text(stat(123, 424242))
            (process / "status").write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "PPid:\t1\n"
                "Cpus_allowed_list:\t197\n"
            )
            (process / "cmdline").write_bytes(str(cp2k).encode() + b"\0")
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=1\0")
            (process / "exe").symlink_to(cp2k)
            for tid, starttime, cpu in (
                (123, 424242, 197),
                (124, 424250, 198),
            ):
                task = process / "task" / str(tid)
                task.mkdir(parents=True)
                (task / "stat").write_text(stat(tid, starttime))
                (task / "status").write_text(status(cpu))

            poisoned = execution._linux_process_snapshot(
                123, cp2k, root / "proc"
            )
            assert poisoned is not None
            self.assertEqual(poisoned["thread_affinity_scan_status"], "consistent")
            accumulated = execution._accumulate_process_snapshot(
                None, poisoned, (196, 197)
            )
            self.assertIs(
                accumulated["current_sample_matches_assigned_singleton"], True
            )
            self.assertIs(accumulated["current_thread_affinity_sample_exact"], False)
            self.assertIs(accumulated["thread_affinity_violation_ever"], True)
            self.assertIs(accumulated["affinity_violation_ever"], True)
            self.assertEqual(
                accumulated["observed_thread_cpu_masks"], ["197", "198"]
            )

            (process / "task" / "124" / "status").write_text(status(197))
            corrected = execution._linux_process_snapshot(
                123, cp2k, root / "proc"
            )
            assert corrected is not None
            accumulated = execution._accumulate_process_snapshot(
                accumulated, corrected, (196, 197)
            )
            self.assertIs(accumulated["current_thread_affinity_sample_exact"], True)
            self.assertIs(accumulated["thread_affinity_violation_ever"], True)
            self.assertIs(accumulated["all_thread_affinity_samples_exact"], False)
            self.assertIs(accumulated["affinity_violation_ever"], True)

    def test_terminal_process_does_not_launder_tid_reuse(self) -> None:
        def stat(tid: int, state: str, starttime: int) -> str:
            return f"{tid} (cp2k-thread) " + " ".join(
                [state, *("0" for _ in range(18)), str(starttime)]
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            process = root / "proc" / "123"
            task = process / "task" / "123"
            task.mkdir(parents=True)
            process_status = (
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "PPid:\t1\n"
                "Cpus_allowed_list:\t197\n"
            )
            (process / "stat").write_text(stat(123, "R", 424242))
            (process / "status").write_text(process_status)
            (process / "cmdline").write_bytes(str(cp2k).encode() + b"\0")
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=0\0")
            (process / "exe").symlink_to(cp2k)
            (task / "stat").write_text(stat(123, "R", 424242))
            (task / "status").write_text(process_status)
            real_read_text = Path.read_text
            process_stats = iter(
                (stat(123, "R", 424242), stat(123, "Z", 424242))
            )
            thread_stats = iter(
                (stat(123, "R", 424242), stat(123, "Z", 424243))
            )

            def reused_tid(path: Path, *args, **kwargs):
                if path == process / "stat":
                    return next(process_stats)
                if path == task / "stat":
                    return next(thread_stats)
                return real_read_text(path, *args, **kwargs)

            with mock.patch.object(Path, "read_text", reused_tid):
                snapshot = execution._linux_process_snapshot(
                    123, cp2k, root / "proc"
                )
            assert snapshot is not None
            self.assertEqual(
                snapshot["thread_affinity_process_status"], "terminal_process"
            )
            self.assertEqual(
                snapshot["thread_affinity_scan_status"],
                "tid_reused_during_sample",
            )
            self.assertEqual(
                snapshot["thread_affinity_scan_issues"],
                ["tid_reused_during_sample"],
            )
            self.assertFalse(
                execution._thread_affinity_sample_matches(snapshot, 197)
            )
            accumulated = execution._accumulate_process_snapshot(
                None, snapshot, (197,)
            )
            self.assertEqual(
                accumulated["thread_affinity_scan_issues_ever"],
                ["tid_reused_during_sample"],
            )
            self.assertIs(accumulated["thread_affinity_violation_ever"], True)

    def test_process_disappearance_does_not_launder_thread_set_race(self) -> None:
        def stat(tid: int, starttime: int) -> str:
            return f"{tid} (cp2k-thread) " + " ".join(
                ["R", *("0" for _ in range(18)), str(starttime)]
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k = root / "cp2k.psmp"
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            process = root / "proc" / "123"
            leader = process / "task" / "123"
            newcomer = process / "task" / "124"
            leader.mkdir(parents=True)
            newcomer.mkdir(parents=True)
            process_status = (
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "PPid:\t1\n"
                "Cpus_allowed_list:\t197\n"
            )
            (process / "stat").write_text(stat(123, 424242))
            (process / "status").write_text(process_status)
            (process / "cmdline").write_bytes(str(cp2k).encode() + b"\0")
            (process / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=0\0")
            (process / "exe").symlink_to(cp2k)
            for task, tid, starttime in (
                (leader, 123, 424242),
                (newcomer, 124, 424250),
            ):
                (task / "stat").write_text(stat(tid, starttime))
                (task / "status").write_text(process_status)
            real_iterdir = Path.iterdir
            task_scans = iter(([leader], [leader, newcomer]))

            def changing_task_set(path: Path):
                if path == process / "task":
                    return iter(next(task_scans))
                return real_iterdir(path)

            real_read_text = Path.read_text
            process_stat_reads = 0

            def disappearing_process(path: Path, *args, **kwargs):
                nonlocal process_stat_reads
                if path == process / "stat":
                    process_stat_reads += 1
                    if process_stat_reads == 2:
                        raise FileNotFoundError(path)
                return real_read_text(path, *args, **kwargs)

            with mock.patch.object(
                Path, "iterdir", changing_task_set
            ), mock.patch.object(Path, "read_text", disappearing_process):
                snapshot = execution._linux_process_snapshot(
                    123, cp2k, root / "proc"
                )
            assert snapshot is not None
            self.assertEqual(
                snapshot["thread_affinity_process_status"],
                "process_disappeared_after_sample",
            )
            self.assertEqual(
                snapshot["thread_affinity_scan_status"],
                "thread_set_changed_during_sample",
            )
            self.assertEqual(
                snapshot["thread_affinity_scan_issues"],
                ["thread_set_changed_during_sample"],
            )
            self.assertFalse(
                execution._thread_affinity_sample_matches(snapshot, 197)
            )
            accumulated = execution._accumulate_process_snapshot(
                None, snapshot, (197,)
            )
            self.assertEqual(
                accumulated["thread_affinity_scan_issues_ever"],
                ["thread_set_changed_during_sample"],
            )
            self.assertIs(accumulated["thread_affinity_violation_ever"], True)

    def test_terminal_environment_loss_is_pending_then_strictly_resolved(self) -> None:
        assigned = (196, 197, 198, 199)
        first = execution._accumulate_process_snapshot(
            None, self.rank_sample(), assigned
        )
        lost = execution._accumulate_process_snapshot(
            first,
            self.rank_sample(
                rank=None,
                observation_status="environment_empty",
            ),
            assigned,
        )
        self.assertIs(lost["rank_environment_unavailable_pending"], True)
        self.assertIs(lost["affinity_violation_ever"], False)
        self.assertEqual(lost["ompi_comm_world_rank"], 1)
        self.assertIsNone(lost["raw_ompi_comm_world_rank"])
        self.assertEqual(lost["observed_rank_ids"], [1])

        second_lost = execution._accumulate_process_snapshot(
            lost,
            self.rank_sample(
                rank=None,
                observation_status="environment_unreadable",
            ),
            assigned,
        )
        self.assertEqual(second_lost["rank_environment_unavailable_sample_count"], 2)
        self.assertIs(second_lost["rank_environment_unavailable_pending"], True)
        execution._resolve_pending_rank_environment(
            second_lost, "process_disappeared"
        )
        self.assertIs(second_lost["rank_environment_unavailable_pending"], False)
        self.assertIs(second_lost["rank_environment_terminally_confirmed"], True)
        self.assertIs(second_lost["affinity_violation_ever"], False)
        self.assertIsNone(
            execution._rank_environment_evidence_issue(
                second_lost, "197", Path("record.json")
            )
        )
        self.assertIsNone(
            execution._rank_process_provenance_issue(
                second_lost,
                1,
                Path("record.json"),
                Path("/tmp/cp2k.psmp"),
            )
        )
        aggregate = execution._aggregate_cp2k_rank_generations(
            {100: second_lost}, assigned
        )
        self.assertIs(aggregate[0]["all_samples_match_assigned_singleton"], True)

        tampering = (
            ("rank_observation_status", "explicit"),
            ("process_identity_status", "pid_reused_during_sample"),
            ("stat_state", "Q"),
            ("state", "bogus"),
            ("observed_process_states", ["bogus"]),
            ("sample_count", 99),
        )
        for field, value in tampering:
            with self.subTest(tamper=field):
                broken = {**second_lost, field: value}
                issue = execution._rank_process_provenance_issue(
                    broken,
                    1,
                    Path("record.json"),
                    Path("/tmp/cp2k.psmp"),
                ) or execution._rank_environment_evidence_issue(
                    broken, "197", Path("record.json")
                )
                self.assertIsNotNone(issue)
        broken = {
            **second_lost,
            "rank_environment_events": [
                {
                    **event,
                    "state": "bogus",
                }
                for event in second_lost["rank_environment_events"]
            ],
        }
        self.assertIsNotNone(
            execution._rank_environment_evidence_issue(
                broken, "197", Path("record.json")
            )
        )
        mismatched_terminal_proof = {
            **second_lost,
            "rank_environment_terminal_confirmation": "terminal_state_Z",
            "rank_environment_events": [
                {**event, "terminal_resolution": "terminal_state_Z"}
                for event in second_lost["rank_environment_events"]
            ],
        }
        self.assertIsNotNone(
            execution._rank_environment_evidence_issue(
                mismatched_terminal_proof, "197", Path("record.json")
            )
        )

    def test_terminal_environment_loss_preserves_duplicate_rank_detection(self) -> None:
        assigned = (196, 197, 198, 199)
        first = execution._accumulate_process_snapshot(
            None, self.rank_sample(), assigned
        )
        pending = execution._accumulate_process_snapshot(
            first,
            self.rank_sample(
                rank=None,
                observation_status="environment_empty",
            ),
            assigned,
        )
        duplicate = execution._accumulate_process_snapshot(
            None, self.rank_sample(pid=101), assigned
        )
        self.assertEqual(
            execution._concurrent_live_rank_pid_groups([pending, duplicate]),
            {1: [100, 101]},
        )

    def test_rank_environment_loss_rejects_nonterminal_and_identity_anomalies(
        self,
    ) -> None:
        assigned = (196, 197, 198, 199)
        first = execution._accumulate_process_snapshot(
            None, self.rank_sample(), assigned
        )
        cases = {
            "initial-pid-reuse": execution._accumulate_process_snapshot(
                None,
                self.rank_sample(identity_status="pid_reused_during_sample"),
                assigned,
            ),
            "initial-identity-unreadable": execution._accumulate_process_snapshot(
                None,
                self.rank_sample(identity_status="identity_unreadable_after_sample"),
                assigned,
            ),
            "initial-executable-change": execution._accumulate_process_snapshot(
                None,
                self.rank_sample(
                    identity_status="executable_changed_during_sample"
                ),
                assigned,
            ),
            "initial-mask-change": execution._accumulate_process_snapshot(
                None,
                self.rank_sample(identity_status="cpu_mask_changed_during_sample"),
                assigned,
            ),
            "initial-loss": execution._accumulate_process_snapshot(
                None,
                self.rank_sample(
                    rank=None, observation_status="environment_empty"
                ),
                assigned,
            ),
            "pid-reuse": execution._accumulate_process_snapshot(
                first,
                self.rank_sample(
                    rank=None,
                    starttime=424243,
                    observation_status="environment_empty",
                ),
                assigned,
            ),
            "changed-mask": execution._accumulate_process_snapshot(
                first,
                self.rank_sample(
                    rank=None,
                    cpu="197-198",
                    observation_status="environment_empty",
                ),
                assigned,
            ),
            "missing-rank": execution._accumulate_process_snapshot(
                first,
                self.rank_sample(rank=None, observation_status="explicit_missing"),
                assigned,
            ),
            "invalid-rank": execution._accumulate_process_snapshot(
                first,
                self.rank_sample(rank=None, observation_status="explicit_invalid"),
                assigned,
            ),
        }
        for name, record in cases.items():
            with self.subTest(name=name):
                self.assertIs(record["affinity_violation_ever"], True)
                self.assertIs(
                    record.get("rank_environment_unavailable_pending"), False
                )

        pending = execution._accumulate_process_snapshot(
            first,
            self.rank_sample(rank=None, observation_status="environment_empty"),
            assigned,
        )
        for rank, cpu in ((1, "197"), (2, "198")):
            with self.subTest(reappearing_rank=rank):
                reappeared = execution._accumulate_process_snapshot(
                    pending, self.rank_sample(rank=rank, cpu=cpu), assigned
                )
                self.assertIs(reappeared["rank_identity_changed_ever"], True)
                self.assertIs(reappeared["affinity_violation_ever"], True)
                self.assertIs(
                    reappeared["rank_environment_unavailable_pending"], False
                )

        still_live = dict(pending)
        execution._resolve_pending_rank_environment(
            still_live, "launcher_ended_while_process_live"
        )
        self.assertIs(still_live["rank_environment_terminally_confirmed"], False)
        self.assertIs(still_live["affinity_violation_ever"], True)
        self.assertIn(
            "invalid terminal rank-environment evidence",
            execution._rank_environment_evidence_issue(
                still_live, "197", Path("record.json")
            )
            or "",
        )

        executable_changed = execution._accumulate_process_snapshot(
            first,
            {
                **self.rank_sample(rank=None),
                "is_cp2k_rank": False,
                "executable": "/tmp/not-cp2k",
                "arguments": [],
            },
            assigned,
        )
        self.assertIs(executable_changed["is_cp2k_rank"], True)
        self.assertIs(executable_changed["executable_changed_ever"], True)
        self.assertIs(executable_changed["affinity_violation_ever"], True)

    def test_terminal_resolver_requires_same_process_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp)
            self.assertEqual(
                execution._linux_process_terminal_resolution(
                    100, 424242, proc_root
                ),
                "process_disappeared",
            )
            process = proc_root / "100"
            process.mkdir()
            def process_stat(state: str, starttime: int) -> str:
                return "100 (cp2k) " + " ".join(
                    [state, *("0" for _ in range(18)), str(starttime)]
                )

            (process / "stat").write_text(process_stat("R", 424243))
            self.assertEqual(
                execution._linux_process_terminal_resolution(
                    100, 424242, proc_root
                ),
                "pid_reused",
            )
            (process / "stat").write_text(process_stat("R", 424242))
            self.assertIsNone(
                execution._linux_process_terminal_resolution(
                    100, 424242, proc_root
                )
            )
            pending = execution._accumulate_process_snapshot(
                execution._accumulate_process_snapshot(
                    None, self.rank_sample(), (196, 197, 198, 199)
                ),
                self.rank_sample(
                    rank=None, observation_status="environment_empty"
                ),
                (196, 197, 198, 199),
            )
            self.assertIs(
                execution._observed_rank_process_is_still_live(
                    100, pending, proc_root
                ),
                True,
            )
            self.assertIs(pending["snapshot_unavailable_ever"], True)
            self.assertIs(pending["affinity_violation_ever"], True)
            duplicate = execution._accumulate_process_snapshot(
                None, self.rank_sample(pid=101), (196, 197, 198, 199)
            )
            self.assertEqual(
                execution._concurrent_live_rank_pid_groups([pending, duplicate]),
                {1: [100, 101]},
            )
            (process / "stat").write_text(process_stat("Z", 424242))
            self.assertEqual(
                execution._linux_process_terminal_resolution(
                    100, 424242, proc_root
                ),
                "terminal_state_Z",
            )

    def test_smt_core_topology_fails_before_launch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            topology = Path(tmp)
            siblings = topology / "cpu96" / "topology" / "thread_siblings_list"
            siblings.parent.mkdir(parents=True)
            siblings.write_text("96,192\n")
            with self.assertRaisesRegex(ValueError, "SMT siblings"):
                execution.require_single_pu_cores([(96,)], topology)
            siblings.write_text("96\n")
            execution.require_single_pu_cores([(96,)], topology)

    def test_schema_v2_record_requires_per_rank_singletons_and_hashed_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inp = root / "job.inp"
            out = root / "job.out"
            stamp = root / "job.out.job.json"
            log = execution.launcher_log_path(out)
            cp2k = root / "cp2k"
            launcher = root / "mpiexec"
            inp.write_text("input\n")
            out.write_text("output\n")
            cp2k.write_text("#!/bin/sh\nexit 0\n")
            cp2k.chmod(0o755)
            stamp.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "executable": str(cp2k.resolve()),
                        "executable_sha256": execution.sha256(cp2k),
                        "input": str(inp.resolve()),
                        "input_sha256": execution.sha256(inp),
                        "completed": True,
                        "return_code": 0,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
            launcher.write_text("launcher\n")
            log.write_text(
                "[node:1] MCW rank 0 bound to package 0[core 0]\n"
                "[node:2] MCW rank 1 bound to package 0[core 1]\n"
            )
            contract = {
                "schema_version": 2,
                "mode": "openmpi_ordered_pe_list",
                "mpi_ranks_per_job": 2,
                "mpi_launcher": str(launcher),
                "mpi_launcher_sha256": execution.sha256(launcher),
                "mpi_launcher_args": [],
                "ordered_pe_lists": ["48,49"],
                "binding_environment_scrub_contract": (
                    execution.binding_environment_scrub_contract()
                ),
                "pool_close_policy": "reject_while_run_admitted_or_active",
                "rank_affinity_observation": "linux_proc_per_task_tid_starttime",
            }
            command = [
                str(launcher),
                "--map-by",
                "pe-list=48,49:ordered",
                "--bind-to",
                "core",
                "--report-bindings",
                "-np",
                "2",
                str(cp2k),
                "-i",
                str(inp.resolve()),
                "-o",
                str(out.resolve()),
            ]
            record = {
                "schema_version": 2,
                "contract": contract,
                "contract_sha256": execution.canonical_sha256(contract),
                "assigned_ordered_pe_list": "48,49",
                "assigned_cpu_count": 2,
                "command": command,
                "return_code": 0,
                "runtime_affinity_gate": True,
                "cross_process_cpu_reservation_gate": True,
                "binding_environment_scrub_gate": True,
                "residual_binding_environment_keys": [],
                "ld_library_path_preserved": True,
                "thread_environment": dict(
                    execution.SANCTIONED_THREAD_ENVIRONMENT
                ),
                "removed_mpi_binding_environment_keys": [
                    "LD_PRELOAD",
                    "OMP_NUM_THREADS",
                ],
                "live_compute_overlap_preflight_gate": True,
                "mpi_bind_to": "core",
                "timing_classification": "production_scaling_eligible",
                "observed_cp2k_rank_pids": [900, 100],
                "observed_cp2k_rank_ids": [0, 1],
                "observed_cp2k_rank_masks": ["48", "49"],
                "launcher_log": str(log.resolve()),
                "launcher_log_sha256": execution.sha256(log),
                "reported_binding_rank_ids": [0, 1],
                "binding_report_complete": True,
                "all_observed_rank_samples_match_ordered_pe_list": True,
                "live_compute_overlap_preflight_owners": [],
                "live_compute_overlap_runtime_gate": True,
                "live_compute_overlap_runtime_samples": [],
                "local_affinity_violation_gate": True,
                "local_affinity_violation_samples": [],
                "observed_child_processes": [
                    {
                        "pid": 900,
                        "is_cp2k_rank": True,
                        "executable": str(cp2k.resolve()),
                        "arguments": [str(cp2k.resolve())],
                        "ompi_comm_world_rank": 0,
                        "cpus_allowed_list": "48",
                        "sample_count": 2,
                        "observed_rank_ids": [0],
                        "observed_cpu_masks": ["48"],
                        "observed_rank_observation_statuses": ["explicit"],
                        "rank_observation_status": "explicit",
                        "raw_ompi_comm_world_rank": 0,
                        "rank_identity_source": "explicit_environment",
                        "process_starttime": 9000,
                        "observed_process_starttimes": [9000],
                        "process_starttime_changed_ever": False,
                        "process_identity_status": "stable",
                        "snapshot_consistency_status": "consistent",
                        "state": "R (running)",
                        "stat_state": "R",
                        "observed_process_states": ["R (running)"],
                        "process_terminally_confirmed": True,
                        "process_terminal_confirmation": "process_disappeared",
                        "process_reappeared_after_terminal_ever": False,
                        "executable_changed_ever": False,
                        "cpu_mask_changed_during_sample_ever": False,
                        "snapshot_unavailable_ever": False,
                        "rank_environment_unavailable_ever": False,
                        "rank_environment_unavailable_sample_count": 0,
                        "rank_environment_unavailable_pending": False,
                        "rank_environment_terminally_confirmed": False,
                        "rank_environment_terminal_confirmation": None,
                        "rank_environment_events": [],
                        "current_sample_matches_assigned_singleton": True,
                        "rank_identity_changed_ever": False,
                        "affinity_violation_ever": False,
                        "thread_affinity_evidence_source": (
                            execution.THREAD_AFFINITY_EVIDENCE_SOURCE
                        ),
                        "thread_affinity_scan_status": "consistent",
                        "thread_affinity_scan_issues": [],
                        "thread_affinity_process_status": "live",
                        "live_thread_affinity": [
                            {
                                "tid": 900,
                                "thread_starttime": 9000,
                                "observed_thread_starttimes": [9000],
                                "state": "R (running)",
                                "stat_state": "R",
                                "cpus_allowed_list": "48",
                                "observed_cpu_masks": ["48"],
                                "live": True,
                                "identity_status": "stable",
                            }
                        ],
                        "thread_affinity_sample_count": 2,
                        "thread_affinity_scan_statuses": ["consistent"],
                        "thread_affinity_evidence_sources": [
                            execution.THREAD_AFFINITY_EVIDENCE_SOURCE
                        ],
                        "thread_affinity_scan_issues_ever": [],
                        "thread_affinity_process_statuses": ["live"],
                        "observed_thread_cpu_masks": ["48"],
                        "observed_thread_identities": ["900:9000"],
                        "current_thread_affinity_sample_exact": True,
                        "all_thread_affinity_samples_exact": True,
                        "thread_affinity_violation_ever": False,
                    },
                    {
                        "pid": 100,
                        "is_cp2k_rank": True,
                        "executable": str(cp2k.resolve()),
                        "arguments": [str(cp2k.resolve())],
                        "ompi_comm_world_rank": 1,
                        "cpus_allowed_list": "49",
                        "sample_count": 2,
                        "observed_rank_ids": [1],
                        "observed_cpu_masks": ["49"],
                        "observed_rank_observation_statuses": ["explicit"],
                        "rank_observation_status": "explicit",
                        "raw_ompi_comm_world_rank": 1,
                        "rank_identity_source": "explicit_environment",
                        "process_starttime": 1000,
                        "observed_process_starttimes": [1000],
                        "process_starttime_changed_ever": False,
                        "process_identity_status": "stable",
                        "snapshot_consistency_status": "consistent",
                        "state": "S (sleeping)",
                        "stat_state": "S",
                        "observed_process_states": ["S (sleeping)"],
                        "process_terminally_confirmed": True,
                        "process_terminal_confirmation": "process_disappeared",
                        "process_reappeared_after_terminal_ever": False,
                        "executable_changed_ever": False,
                        "cpu_mask_changed_during_sample_ever": False,
                        "snapshot_unavailable_ever": False,
                        "rank_environment_unavailable_ever": False,
                        "rank_environment_unavailable_sample_count": 0,
                        "rank_environment_unavailable_pending": False,
                        "rank_environment_terminally_confirmed": False,
                        "rank_environment_terminal_confirmation": None,
                        "rank_environment_events": [],
                        "current_sample_matches_assigned_singleton": True,
                        "rank_identity_changed_ever": False,
                        "affinity_violation_ever": False,
                        "thread_affinity_evidence_source": (
                            execution.THREAD_AFFINITY_EVIDENCE_SOURCE
                        ),
                        "thread_affinity_scan_status": "consistent",
                        "thread_affinity_scan_issues": [],
                        "thread_affinity_process_status": "live",
                        "live_thread_affinity": [
                            {
                                "tid": 100,
                                "thread_starttime": 1000,
                                "observed_thread_starttimes": [1000],
                                "state": "S (sleeping)",
                                "stat_state": "S",
                                "cpus_allowed_list": "49",
                                "observed_cpu_masks": ["49"],
                                "live": True,
                                "identity_status": "stable",
                            }
                        ],
                        "thread_affinity_sample_count": 2,
                        "thread_affinity_scan_statuses": ["consistent"],
                        "thread_affinity_evidence_sources": [
                            execution.THREAD_AFFINITY_EVIDENCE_SOURCE
                        ],
                        "thread_affinity_scan_issues_ever": [],
                        "thread_affinity_process_statuses": ["live"],
                        "observed_thread_cpu_masks": ["49"],
                        "observed_thread_identities": ["100:1000"],
                        "current_thread_affinity_sample_exact": True,
                        "all_thread_affinity_samples_exact": True,
                        "thread_affinity_violation_ever": False,
                    },
                ],
                "observed_cp2k_rank_pid_generations": [[900], [100]],
                "observed_cp2k_rank_evidence": [
                    {
                        "ompi_comm_world_rank": 0,
                        "canonical_pid": 900,
                        "pid_generations": [900],
                        "observed_cpu_masks": ["48"],
                        "all_samples_match_assigned_singleton": True,
                        "concurrent_duplicate_pid_ever": False,
                    },
                    {
                        "ompi_comm_world_rank": 1,
                        "canonical_pid": 100,
                        "pid_generations": [100],
                        "observed_cpu_masks": ["49"],
                        "all_samples_match_assigned_singleton": True,
                        "concurrent_duplicate_pid_ever": False,
                    },
                ],
                "observed_cp2k_rank_count": 2,
                "observed_cp2k_process_generation_count": 2,
                "concurrent_duplicate_rank_ids_ever": [],
                "concurrent_duplicate_rank_samples": [],
                "concurrent_duplicate_rank_processes_ever": False,
                "unranked_cp2k_process_seen": False,
                "expected_cp2k_rank_count": 2,
                "rank_count_matches": True,
                "rank_ids_exactly_0_to_n_minus_1": True,
                "rank_masks_complete": True,
                "rank_masks_exactly_match_ordered_pe_list": True,
                "mpi_launcher_sha256_at_launch": execution.sha256(launcher),
                "cp2k": str(cp2k),
                "cp2k_sha256_at_launch": execution.sha256(cp2k),
                "input": str(inp.resolve()),
                "input_sha256_at_launch": execution.sha256(inp),
                "working_directory": str(inp.parent.resolve()),
                "output": str(out.resolve()),
                "output_sha256": execution.sha256(out),
                "scientific_job_stamp": str(stamp.resolve()),
                "scientific_job_stamp_sha256": execution.sha256(stamp),
            }
            path = execution.execution_record_path(out)
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertIsNone(execution.recorded_execution_issue(path, contract, out, stamp))
            self.assertEqual(
                execution.execution_record_timing_classification(
                    path, contract, out, stamp
                ),
                "production_scaling_eligible",
            )
            prior_v2_contract = json.loads(json.dumps(contract))
            for field in (
                "binding_environment_scrub_contract",
                "pool_close_policy",
                "rank_affinity_observation",
            ):
                prior_v2_contract.pop(field)
            prior_v2_record = json.loads(json.dumps(record))
            prior_v2_record["contract"] = prior_v2_contract
            prior_v2_record["contract_sha256"] = execution.canonical_sha256(
                prior_v2_contract
            )
            for field in (
                "binding_environment_scrub_gate",
                "residual_binding_environment_keys",
                "ld_library_path_preserved",
                "thread_environment",
                "removed_mpi_binding_environment_keys",
                "local_affinity_violation_gate",
                "local_affinity_violation_samples",
            ):
                prior_v2_record.pop(field)
            thread_fields = (
                "thread_affinity_evidence_source",
                "thread_affinity_scan_status",
                "thread_affinity_scan_issues",
                "thread_affinity_process_status",
                "live_thread_affinity",
                "thread_affinity_sample_count",
                "thread_affinity_scan_statuses",
                "thread_affinity_evidence_sources",
                "thread_affinity_scan_issues_ever",
                "thread_affinity_process_statuses",
                "observed_thread_cpu_masks",
                "observed_thread_identities",
                "current_thread_affinity_sample_exact",
                "all_thread_affinity_samples_exact",
                "thread_affinity_violation_ever",
            )
            for child in prior_v2_record["observed_child_processes"]:
                for field in thread_fields:
                    child.pop(field)
            path.write_text(
                json.dumps(prior_v2_record, indent=2, sort_keys=True) + "\n"
            )
            self.assertIsNone(
                execution.recorded_execution_issue(
                    path, prior_v2_contract, out, stamp
                )
            )
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            for field, invalid, message in (
                ("observed_cp2k_rank_ids", [1, 0], "rank ordering"),
                ("observed_cp2k_rank_masks", ["48", "48"], "CPU-mask"),
                ("reported_binding_rank_ids", [0], "binding report"),
                (
                    "live_compute_overlap_runtime_gate",
                    False,
                    "runtime live CP2K/MPI overlap gate",
                ),
                (
                    "live_compute_overlap_runtime_samples",
                    [{"sample_index": 2, "owners": [{"pid": 77}]}],
                    "runtime live CP2K/MPI overlap",
                ),
                (
                    "concurrent_duplicate_rank_processes_ever",
                    True,
                    "duplicate-rank evidence",
                ),
                (
                    "binding_environment_scrub_gate",
                    False,
                    "environment scrub gate",
                ),
                (
                    "local_affinity_violation_gate",
                    False,
                    "local rank/thread affinity gate",
                ),
            ):
                with self.subTest(field=field):
                    broken = dict(record)
                    broken[field] = invalid
                    path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
                    self.assertIn(
                        message,
                        execution.recorded_execution_issue(path, contract, out, stamp) or "",
                    )
            for name, mutate, message in (
                (
                    "thread-environment",
                    lambda item: item["thread_environment"].update(
                        {"OMP_NUM_THREADS": "2"}
                    ),
                    "thread environment",
                ),
                (
                    "foreign-removed-key",
                    lambda item: item.update(
                        {
                            "removed_mpi_binding_environment_keys": [
                                "LD_PRELOAD",
                                "PATH",
                            ]
                        }
                    ),
                    "removed binding-environment key evidence",
                ),
                (
                    "duplicate-removed-key",
                    lambda item: item.update(
                        {
                            "removed_mpi_binding_environment_keys": [
                                "LD_PRELOAD",
                                "LD_PRELOAD",
                            ]
                        }
                    ),
                    "removed binding-environment key evidence",
                ),
            ):
                with self.subTest(environment_evidence=name):
                    broken = json.loads(json.dumps(record))
                    mutate(broken)
                    path.write_text(
                        json.dumps(broken, indent=2, sort_keys=True) + "\n"
                    )
                    self.assertIn(
                        message,
                        execution.recorded_execution_issue(
                            path, contract, out, stamp
                        )
                        or "",
                    )
            broken = dict(record)
            broken["command"] = [*command, "--bind-to", "none"]
            path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "command/affinity mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            log.write_text("tampered\n")
            self.assertIn(
                "launcher-log hash mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )
            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            broken = dict(record)
            broken["return_code"] = 1
            path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "return code",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )
            broken = dict(record)
            broken["all_observed_rank_samples_match_ordered_pe_list"] = False
            path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "rank-sample gate mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )

            for name, mutate, message in (
                (
                    "duplicate-rank-list",
                    lambda item: item.update(
                        {"concurrent_duplicate_rank_ids_ever": [0]}
                    ),
                    "duplicate-rank evidence",
                ),
                (
                    "generation-count",
                    lambda item: item.update(
                        {"observed_cp2k_process_generation_count": 999}
                    ),
                    "process-generation count",
                ),
                (
                    "rank-count-gate",
                    lambda item: item.update({"rank_count_matches": False}),
                    "rank-count gate",
                ),
                (
                    "mask-history",
                    lambda item: item["observed_child_processes"][0].update(
                        {"observed_cpu_masks": ["49", "48"]}
                    ),
                    "affinity history",
                ),
                (
                    "rank-history",
                    lambda item: item["observed_child_processes"][0].update(
                        {"observed_rank_ids": [1, 0]}
                    ),
                    "affinity history",
                ),
                (
                    "rank-starttime",
                    lambda item: item["observed_child_processes"][0].update(
                        {"process_starttime": 7}
                    ),
                    "process provenance",
                ),
                (
                    "rank-not-terminal",
                    lambda item: item["observed_child_processes"][0].update(
                        {"process_terminally_confirmed": False}
                    ),
                    "process provenance",
                ),
                (
                    "rank-raw-identity",
                    lambda item: item["observed_child_processes"][0].update(
                        {"raw_ompi_comm_world_rank": 1}
                    ),
                    "process provenance",
                ),
                (
                    "rank-snapshot-loss",
                    lambda item: item["observed_child_processes"][0].update(
                        {"snapshot_unavailable_ever": True}
                    ),
                    "process provenance",
                ),
                (
                    "child-thread-repinned",
                    lambda item: item["observed_child_processes"][0][
                        "live_thread_affinity"
                    ][0].update(
                        {
                            "cpus_allowed_list": "49",
                            "observed_cpu_masks": ["48", "49"],
                        }
                    ),
                    "thread-affinity evidence",
                ),
                (
                    "sticky-thread-scan-issue",
                    lambda item: item["observed_child_processes"][0].update(
                        {
                            "thread_affinity_scan_issues_ever": [
                                "thread_set_changed_during_sample"
                            ]
                        }
                    ),
                    "thread-affinity evidence",
                ),
            ):
                with self.subTest(tamper=name):
                    broken = json.loads(json.dumps(record))
                    mutate(broken)
                    path.write_text(
                        json.dumps(broken, indent=2, sort_keys=True) + "\n"
                    )
                    self.assertIn(
                        message,
                        execution.recorded_execution_issue(
                            path, contract, out, stamp
                        )
                        or "",
                    )

            broken = json.loads(json.dumps(record))
            successor = json.loads(
                json.dumps(broken["observed_child_processes"][0])
            )
            successor.update(
                {
                    "pid": 901,
                    "process_starttime": 9010,
                    "observed_process_starttimes": [9010],
                }
            )
            successor["live_thread_affinity"][0].update(
                {
                    "tid": 901,
                    "thread_starttime": 9010,
                    "observed_thread_starttimes": [9010],
                }
            )
            successor["observed_thread_identities"] = ["901:9010"]
            broken["observed_child_processes"].append(successor)
            path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "multiple CP2K PID generations",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )

            for field in (
                "rank_environment_unavailable_ever",
                "rank_environment_unavailable_sample_count",
                "rank_environment_unavailable_pending",
                "rank_environment_terminally_confirmed",
                "rank_environment_terminal_confirmation",
                "rank_environment_events",
            ):
                with self.subTest(missing_environment_field=field):
                    broken = json.loads(json.dumps(record))
                    del broken["observed_child_processes"][0][field]
                    path.write_text(
                        json.dumps(broken, indent=2, sort_keys=True) + "\n"
                    )
                    self.assertIn(
                        "rank-environment evidence",
                        execution.recorded_execution_issue(
                            path, contract, out, stamp
                        )
                        or "",
                    )

            broken = json.loads(json.dumps(record))
            del broken["observed_child_processes"][0][
                "raw_ompi_comm_world_rank"
            ]
            path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "missing CP2K rank process provenance",
                execution.recorded_execution_issue(path, contract, out, stamp)
                or "",
            )

            for field, invalid in (
                ("raw_ompi_comm_world_rank", None),
                ("executable", "/tmp/not-cp2k"),
                ("rank_environment_unavailable_ever", 0),
                ("rank_environment_unavailable_sample_count", False),
            ):
                with self.subTest(invalid_provenance_field=field):
                    broken = json.loads(json.dumps(record))
                    broken["observed_child_processes"][0][field] = invalid
                    path.write_text(
                        json.dumps(broken, indent=2, sort_keys=True) + "\n"
                    )
                    self.assertIsNotNone(
                        execution.recorded_execution_issue(
                            path, contract, out, stamp
                        )
                    )

            for name, mutation in (
                (
                    "raw-rank-bool-alias",
                    {"raw_ompi_comm_world_rank": True},
                ),
                (
                    "rank-history-bool-alias",
                    {"observed_rank_ids": [True]},
                ),
                (
                    "starttime-history-bool-alias",
                    {
                        "process_starttime": 1,
                        "observed_process_starttimes": [True],
                    },
                ),
                (
                    "unrelated-executable-with-cp2k-second-argument",
                    {
                        "executable": "/bin/echo",
                        "arguments": ["/bin/echo", str(cp2k.resolve())],
                    },
                ),
            ):
                with self.subTest(exact_type_or_classifier_fuzz=name):
                    broken = json.loads(json.dumps(record))
                    broken["observed_child_processes"][1].update(mutation)
                    path.write_text(
                        json.dumps(broken, indent=2, sort_keys=True) + "\n"
                    )
                    self.assertIsNotNone(
                        execution.recorded_execution_issue(
                            path, contract, out, stamp
                        )
                    )

    def test_schema_v1_remains_readable_but_timing_is_non_scaling(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out, stamp, contract, _ = self.valid_record_fixture(root)
            path = execution.execution_record_path(out)
            self.assertIsNone(execution.recorded_execution_issue(path, contract, out, stamp))
            self.assertEqual(
                execution.execution_record_timing_classification(path),
                "legacy_timing_non_scaling",
            )

    def test_record_rejects_command_and_cp2k_binary_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp2k, inp, out, stamp, contract, record = self.valid_record_fixture(root)
            path = execution.execution_record_path(out)

            original_command = record["command"]
            self.assertIsInstance(original_command, list)
            assert isinstance(original_command, list)
            bad_command = dict(record)
            bad_command["command"] = [*original_command[:-1], str(root / "wrong.out")]
            path.write_text(json.dumps(bad_command, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "full execution command/affinity mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )

            alternate = cp2k.with_name("cp2k-alternate.psmp")
            alternate.write_bytes(cp2k.read_bytes())
            alternate.chmod(0o755)
            wrong_binary = dict(record)
            wrong_binary["cp2k"] = str(alternate.resolve())
            wrong_binary["cp2k_sha256_at_launch"] = execution.sha256(alternate)
            wrong_binary["command"] = execution.cp2k_command(
                taskset="/fixture/taskset",
                cpu_set="48,49,50,51",
                mpi_launcher="/fixture/mpiexec",
                mpi_launcher_args=["--bind-to", "none"],
                mpi_ranks_per_job=2,
                cp2k=alternate,
                inp=inp,
                out=out,
            )
            path.write_text(json.dumps(wrong_binary, indent=2, sort_keys=True) + "\n")
            self.assertIn(
                "differs from scientific job stamp",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )

            path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            cp2k.write_text("#!/bin/sh\nexit 7\n")
            self.assertIn(
                "CP2K executable hash mismatch",
                execution.recorded_execution_issue(path, contract, out, stamp) or "",
            )


if __name__ == "__main__":
    unittest.main()
