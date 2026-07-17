from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))

import benchmark_execution as execution  # noqa: E402


class BenchmarkExecutionValidationTests(unittest.TestCase):
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
            {
                "pid": 10,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "48-49",
            },
            (48, 49),
        )
        corrected = execution._accumulate_process_snapshot(
            bad,
            {
                "pid": 10,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "48",
            },
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

    def test_procfs_preflight_rejects_live_cp2k_overlap_and_ignores_zombie(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
            status = process / "status"
            status.write_text(
                "Name:\tcp2k.psmp\n"
                "State:\tR (running)\n"
                "Cpus_allowed_list:\t48-49\n"
            )
            (process / "environ").write_bytes(b"")
            owners = execution.live_compute_cpu_owners((48,), proc_root)
            self.assertEqual(owners[0]["pid"], 4242)
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

    def test_rank_generation_aggregation_accepts_sequential_exact_and_rejects_mask_change(
        self,
    ) -> None:
        parent = execution._accumulate_process_snapshot(
            None,
            {
                "pid": 100,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "48",
            },
            (48, 49),
        )
        successor = execution._accumulate_process_snapshot(
            None,
            {
                "pid": 101,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "48",
            },
            (48, 49),
        )
        exact = execution._aggregate_cp2k_rank_generations(
            {100: parent, 101: successor}, (48, 49)
        )
        self.assertEqual(exact[0]["pid_generations"], [100, 101])
        self.assertIs(exact[0]["all_samples_match_assigned_singleton"], True)

        wrong_successor = execution._accumulate_process_snapshot(
            None,
            {
                "pid": 102,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "49",
            },
            (48, 49),
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
            {
                "pid": 100,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 0,
                "cpus_allowed_list": "48",
            },
            (48, 49),
        )
        migrated = execution._accumulate_process_snapshot(
            first,
            {
                "pid": 100,
                "is_cp2k_rank": True,
                "ompi_comm_world_rank": 1,
                "cpus_allowed_list": "49",
            },
            (48, 49),
        )
        self.assertIs(migrated["rank_identity_changed_ever"], True)
        self.assertIs(migrated["affinity_violation_ever"], True)

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
                "observed_child_processes": [
                    {
                        "pid": 900,
                        "is_cp2k_rank": True,
                        "ompi_comm_world_rank": 0,
                        "cpus_allowed_list": "48",
                        "sample_count": 2,
                        "observed_rank_ids": [0],
                        "observed_cpu_masks": ["48"],
                        "current_sample_matches_assigned_singleton": True,
                        "rank_identity_changed_ever": False,
                        "affinity_violation_ever": False,
                    },
                    {
                        "pid": 100,
                        "is_cp2k_rank": True,
                        "ompi_comm_world_rank": 1,
                        "cpus_allowed_list": "49",
                        "sample_count": 2,
                        "observed_rank_ids": [1],
                        "observed_cpu_masks": ["49"],
                        "current_sample_matches_assigned_singleton": True,
                        "rank_identity_changed_ever": False,
                        "affinity_violation_ever": False,
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
            for field, invalid, message in (
                ("observed_cp2k_rank_ids", [1, 0], "rank ordering"),
                ("observed_cp2k_rank_masks", ["48", "48"], "CPU-mask"),
                ("reported_binding_rank_ids", [0], "binding report"),
                (
                    "concurrent_duplicate_rank_processes_ever",
                    True,
                    "duplicate-rank evidence",
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
