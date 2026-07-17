#!/usr/bin/env python3
"""One-job Terok smoke for the production exact-binding ExecutionPool."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import benchmark_execution as execution


def main() -> int:
    if len(sys.argv) != 6:
        raise SystemExit(
            "usage: exact_binding_smoke_driver.py MPIEXEC CP2K INPUT OUTPUT PE_LIST"
        )
    launcher, cp2k, input_path, output_path, pe_list = sys.argv[1:]
    result_path = Path(output_path).with_suffix(".execution-smoke.json")
    pool = execution.ExecutionPool(
        concurrent_jobs=1,
        mpi_ranks_per_job=2,
        threads_per_rank=1,
        mpi_launcher=launcher,
        mpi_launcher_args=[],
        pe_lists=[pe_list],
    )
    try:
        returncode, observation = pool.run_cp2k(
            Path(cp2k), Path(input_path), Path(output_path)
        )
    finally:
        pool.close()
        if "observation" in locals():
            result_path.write_text(
                json.dumps(observation, indent=2, sort_keys=True) + "\n"
            )

    # Unreachable, but retained for type checkers.
    return returncode


if __name__ == "__main__":
    rc = main()
    payload = json.loads(
        Path(sys.argv[4]).with_suffix(".execution-smoke.json").read_text()
    )
    print(
        json.dumps(
            {
                "cp2k_return_code": rc,
                "runtime_affinity_gate": payload.get("runtime_affinity_gate"),
                "rank_ids": payload.get("observed_cp2k_rank_ids"),
                "rank_masks": payload.get("observed_cp2k_rank_masks"),
                "binding_report_complete": payload.get("binding_report_complete"),
                "reservation_gate": payload.get(
                    "cross_process_cpu_reservation_gate"
                ),
            },
            sort_keys=True,
        )
    )
    raise SystemExit(
        0
        if rc == 0
        and payload.get("runtime_affinity_gate") is True
        and payload.get("observed_cp2k_rank_ids") == [0, 1]
        and payload.get("observed_cp2k_rank_masks") == [
            sys.argv[5].split(",")[0],
            sys.argv[5].split(",")[1],
        ]
        else 97
    )
