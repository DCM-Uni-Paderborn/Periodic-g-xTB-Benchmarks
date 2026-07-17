from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
HELPER_PATH = (
    REPOSITORY
    / "validation"
    / "accelerated_exchange"
    / "cp2k_streamed_reverse_consumer"
    / "scripts"
    / "run_with_rss.py"
)
SCRIPTS = HELPER_PATH.parent
SPEC = importlib.util.spec_from_file_location("streamed_reverse_run_with_rss", HELPER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {HELPER_PATH}")
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)


class AcceleratedExchangeAffinityTests(unittest.TestCase):
    def test_shell_launchers_lock_evidence_and_each_reserved_cpu(self) -> None:
        for name in ("run_linux_matrix.sh", "run_linux_mode_rss.sh"):
            with self.subTest(name=name):
                text = (SCRIPTS / name).read_text()
                self.assertIn("flock -n 9", text)
                self.assertIn('cpu-${cpu}.lock', text)
                self.assertIn('flock -n "$lock_fd"', text)
                self.assertLess(text.index("flock -n 9"), text.index("rm -f PASS FAIL"))

    def test_affinity_violation_is_sticky_after_a_corrected_sample(self) -> None:
        bad = helper.accumulate_rank_snapshot(
            None,
            {"pid": 100, "rank": 0, "cpus_allowed_list": "48-49"},
            (48, 49),
        )
        corrected = helper.accumulate_rank_snapshot(
            bad,
            {"pid": 100, "rank": 0, "cpus_allowed_list": "48"},
            (48, 49),
        )
        self.assertTrue(corrected["current_sample_matches_assigned_singleton"])
        self.assertTrue(corrected["affinity_violation_ever"])
        self.assertEqual(corrected["observed_cpu_masks"], ["48-49", "48"])

    def test_ordered_pe_list_is_literal_and_unique(self) -> None:
        self.assertEqual(helper.parse_ordered_pe_list("96,97"), (96, 97))
        for invalid in ("96-97", "96,96", "96,,97"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    helper.parse_ordered_pe_list(invalid)

    def test_openmpi5_capitalized_binding_report_is_parsed(self) -> None:
        self.assertEqual(
            helper.reported_binding_ranks(
                "[terok:123] Rank 0 bound to package 1[core 36]\n"
            ),
            [0],
        )

    def test_unknown_future_hwloc_rmaps_environment_is_removed(self) -> None:
        self.assertEqual(
            helper.mpi_control_environment_keys(
                {
                    "OMPI_MCA_rmaps_future_override": "unsafe",
                    "PRTE_MCA_hwloc_default_binding_policy": "none",
                    "OMPI_MCA_mca_base_param_files": "/tmp/indirect.conf",
                    "OMPI_MCA_pml": "ucx",
                }
            ),
            [
                "OMPI_MCA_mca_base_param_files",
                "OMPI_MCA_pml",
                "OMPI_MCA_rmaps_future_override",
                "PRTE_MCA_hwloc_default_binding_policy",
            ],
        )

    def test_smt_core_is_rejected_before_helper_launch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            topology = Path(tmp)
            siblings = topology / "cpu7" / "topology" / "thread_siblings_list"
            siblings.parent.mkdir(parents=True)
            siblings.write_text("7,31\n")
            with self.assertRaisesRegex(ValueError, "SMT siblings"):
                helper.require_single_pu_cores((7,), topology)

    def test_missing_rank_to_singleton_proof_changes_success_to_97(self) -> None:
        executable = shutil.which("true")
        if executable is None:
            self.skipTest("true executable is unavailable")
        cpu = min(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 0
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "launcher.log"
            log.write_text("[terok:123] Rank 0 bound to package 0[core 0]\n")
            result = root / "result.json"
            process = subprocess.run(
                [
                    sys.executable,
                    str(HELPER_PATH),
                    "--mpi-ranks",
                    "1",
                    "--ordered-pe-list",
                    str(cpu),
                    "--cp2k",
                    executable,
                    "--launcher-log",
                    str(log),
                    str(result),
                    executable,
                ],
                check=False,
            )
            self.assertEqual(process.returncode, 97)
            payload = json.loads(result.read_text())
            self.assertIs(payload["affinity"]["runtime_affinity_gate"], False)
            self.assertEqual(
                payload["affinity"]["timing_classification"],
                "timing_non_scaling",
            )


if __name__ == "__main__":
    unittest.main()
