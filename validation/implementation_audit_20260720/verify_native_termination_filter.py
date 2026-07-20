#!/usr/bin/env python3
"""Exercise positive and negative native-endpoint termination filters."""

from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
ASSEMBLER = PACKAGE / "scripts/assemble_comparison_tables.py"
ABSOLUTE_TABLE = PACKAGE / "tables/cp2k_native_absolute_energies_by_mesh.csv"
QUALIFIED_CP2K_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)


def load_assembler():
    spec = importlib.util.spec_from_file_location("seidler_assembler", ASSEMBLER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import assembler: {ASSEMBLER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    assembler = load_assembler()
    with tempfile.TemporaryDirectory(prefix="gxtb-native-termination-") as temporary:
        work = Path(temporary)
        status_cases = {
            "zero": ("0\n", 0),
            "positive": ("1\n", 1),
            "negative": ("-9\n", -9),
            "malformed": ("success\n", None),
            "empty": ("", None),
        }
        parsed_statuses = {}
        for label, (text, expected) in status_cases.items():
            path = work / f"{label}.status"
            path.write_text(text, encoding="utf-8")
            actual = assembler.cp2k_exit_status(path)
            if actual != expected:
                raise AssertionError(
                    f"status parser mismatch for {label}: {actual!r} != {expected!r}"
                )
            parsed_statuses[label] = actual
        if assembler.cp2k_exit_status(work / "missing.status") is not None:
            raise AssertionError("missing exit status was not rejected")

        ended = work / "ended.out"
        ended.write_text(
            " **** PROGRAM STARTED AT 2026-07-20 00:00:00.000\n"
            " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.250000000000000\n"
            " **** PROGRAM ENDED AT 2026-07-20 00:00:00.000\n",
            encoding="utf-8",
        )
        incomplete = work / "incomplete.out"
        incomplete.write_text(
            " **** PROGRAM STARTED AT 2026-07-20 00:00:00.000\n"
            " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.250000000000000\n",
            encoding="utf-8",
        )
        appended_incomplete = work / "appended-incomplete.out"
        appended_incomplete.write_text(
            ended.read_text(encoding="utf-8")
            + " **** PROGRAM STARTED AT 2026-07-20 01:00:00.000\n"
            + " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -2.500000000000000\n",
            encoding="utf-8",
        )
        ended_energy = assembler.cp2k_energy(ended)
        incomplete_energy = assembler.cp2k_energy(incomplete)
        appended_incomplete_energy = assembler.cp2k_energy(appended_incomplete)
        if (
            ended_energy != -1.25
            or incomplete_energy is not None
            or appended_incomplete_energy is not None
        ):
            raise AssertionError("normal-end-marker filter did not fail closed")

    with ABSOLUTE_TABLE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    admitted = {(int(row["mesh_n"]), row["phase"]) for row in rows}
    admitted_rows_pass = bool(rows) and all(
        row["qualification"] == "PASS"
        and row["cp2k_binary_sha256"] == QUALIFIED_CP2K_SHA256
        and row["exit_status"] == "0"
        and row["normal_termination_qualification"] == "PASS"
        for row in rows
    )
    if not admitted_rows_pass:
        raise AssertionError("an admitted native row lacks exact build/termination proof")

    missing_status_run = PACKAGE / "raw/cp2k_native/k777-reduced/II"
    nonzero_status_run = PACKAGE / "raw/cp2k_native/k666-reduced/XIII"
    real_negative_cases_pass = (
        (7, "II") not in admitted
        and assembler.cp2k_energy(missing_status_run / "cp2k.out") is not None
        and assembler.cp2k_exit_status(missing_status_run / "exit_status") is None
        and (6, "XIII") not in admitted
        and assembler.cp2k_energy(nonzero_status_run / "cp2k.out") is None
        and assembler.cp2k_exit_status(nonzero_status_run / "exit_status") == 137
    )
    if not real_negative_cases_pass:
        raise AssertionError("archived negative native endpoints were not rejected")

    print(json.dumps({
        "schema": "periodic-gxtb-native-termination-filter-v1",
        "status": "PASS",
        "synthetic_status_cases": parsed_statuses,
        "synthetic_normal_end_marker_filter": "PASS",
        "synthetic_appended_incomplete_segment_rejected": "PASS",
        "admitted_endpoint_count": len(rows),
        "admitted_rows_exact_build_and_exit_zero": admitted_rows_pass,
        "archived_missing_status_endpoint_excluded": "PASS",
        "archived_nonzero_status_endpoint_excluded": "PASS",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
