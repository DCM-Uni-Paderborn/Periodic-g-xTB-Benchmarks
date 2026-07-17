from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


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
            {"pid": 10, "rank": 0, "cpus_allowed_list": "96-97"},
            (96, 97),
        )
        corrected = runner.accumulate_rank_snapshot(
            bad,
            {"pid": 10, "rank": 0, "cpus_allowed_list": "96"},
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

    def test_live_overlap_preflight_catches_cp2k_without_rank_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc_root = Path(tmp) / "proc"
            process = proc_root / "4242"
            process.mkdir(parents=True)
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

    def test_sequential_rank_generation_is_exact_but_changed_mask_and_concurrency_fail(
        self,
    ) -> None:
        parent = runner.accumulate_rank_snapshot(
            None,
            {"pid": 10, "rank": 0, "cpus_allowed_list": "96"},
            (96, 97),
        )
        successor = runner.accumulate_rank_snapshot(
            None,
            {"pid": 11, "rank": 0, "cpus_allowed_list": "96"},
            (96, 97),
        )
        exact = runner.ordered_rank_proof(
            {10: parent, 11: successor}, 1, (96,)
        )
        self.assertEqual(exact[0]["pid_generations"], [10, 11])

        wrong = runner.accumulate_rank_snapshot(
            None,
            {"pid": 12, "rank": 0, "cpus_allowed_list": "97"},
            (96, 97),
        )
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof({10: parent, 12: wrong}, 1, (96,))
        with self.assertRaisesRegex(RuntimeError, "singleton-exact"):
            runner.ordered_rank_proof(
                {10: parent, 11: successor}, 1, (96,), {0}
            )

    def test_verifier_reconstructs_generation_and_duplicate_summaries(self) -> None:
        child0 = runner.accumulate_rank_snapshot(
            None,
            {"pid": 10, "rank": 0, "cpus_allowed_list": "96"},
            (96, 97),
        )
        child1 = runner.accumulate_rank_snapshot(
            None,
            {"pid": 11, "rank": 1, "cpus_allowed_list": "97"},
            (96, 97),
        )
        metadata = {
            "observed_child_processes": [child0, child1],
            "concurrent_duplicate_rank_samples": [],
            "concurrent_duplicate_rank_ids_ever": [],
            "concurrent_duplicate_rank_processes_ever": False,
        }
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

    def test_new_launcher_contains_only_exact_ordered_binding(self) -> None:
        source = (ROOT / "run_test_matrix.py").read_text()
        self.assertNotIn('"taskset", "-c"', source)
        self.assertNotIn('"--bind-to", "none"', source)
        self.assertIn('f"pe-list={pe_list}:ordered"', source)
        self.assertIn('"--bind-to", "core", "--report-bindings"', source)
        self.assertIn("cwd=run_dir", source)


if __name__ == "__main__":
    unittest.main()
