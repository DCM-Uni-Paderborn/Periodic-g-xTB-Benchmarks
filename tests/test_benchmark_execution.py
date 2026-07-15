from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


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

    def test_cpu_sets_must_be_disjoint(self) -> None:
        with self.assertRaisesRegex(ValueError, "overlap at"):
            execution.validate_cpu_sets(
                ["0-3", "3-6"],
                concurrent_jobs=2,
                mpi_ranks_per_job=2,
                threads_per_rank=1,
            )

    def test_cpu_set_must_have_enough_rank_thread_slots(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least 4"):
            execution.validate_cpu_sets(
                ["0-2"],
                concurrent_jobs=1,
                mpi_ranks_per_job=2,
                threads_per_rank=2,
            )

    def test_mpi_launcher_cannot_rebind_inside_taskset(self) -> None:
        execution.require_bind_to_none(["--bind-to", "none"])
        for arguments in (
            ["--bind-to", "core"],
            ["--bind-to=socket"],
            ["--bind-to", "none", "--bind-to", "core"],
            [],
        ):
            with self.subTest(arguments=arguments):
                with self.assertRaisesRegex(ValueError, "exactly one '--bind-to none'"):
                    execution.require_bind_to_none(arguments)

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
