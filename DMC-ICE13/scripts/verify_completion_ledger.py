#!/usr/bin/env python3
"""Verify the deduplicated DMC-ICE13 qualified-completion ledger."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LEDGER = ROOT / "data/dmc_ice13_completion_ledger.csv"
RAW = ROOT / "reproduction/seidler_dmc13_recalculation/raw/cp2k_native"
PHASES = {"II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"}
QUALIFIED_CP2K_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sidecar(path: Path) -> str:
    fields = path.read_text(encoding="utf-8").split()
    return fields[0] if fields else ""


def main() -> None:
    with LEDGER.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise AssertionError("completion ledger is empty")

    timestamps = [datetime.fromisoformat(row["result_completed_at"]) for row in rows]
    keys = [row["notification_key"] for row in rows]
    endpoints = [(row["phase"], int(row["mesh_n"])) for row in rows]
    if timestamps != sorted(timestamps) or len(timestamps) != len(set(timestamps)):
        raise AssertionError("completion timestamps are not unique and chronological")
    if len(keys) != len(set(keys)) or any(not key for key in keys):
        raise AssertionError("completion notification keys are not unique and populated")
    if len(endpoints) != len(set(endpoints)):
        raise AssertionError("a qualified phase/mesh completion appears more than once")

    raw_proofs = []
    for index, row in enumerate(rows):
        phase = row["phase"]
        mesh = int(row["mesh_n"])
        if phase not in PHASES or mesh < 1:
            raise AssertionError(f"invalid endpoint in ledger row {index + 1}")
        for field in ("binary_sha256", "input_sha256", "output_sha256"):
            if SHA256_PATTERN.fullmatch(row[field]) is None:
                raise AssertionError(f"invalid {field} in ledger row {index + 1}")
        if row["binary_sha256"] != QUALIFIED_CP2K_SHA256:
            raise AssertionError(f"wrong CP2K build in ledger row {index + 1}")
        if row["output_sha256"][:12] not in row["notification_key"]:
            raise AssertionError(f"notification key lacks output identity in row {index + 1}")

        current = float(row["mixed_mae_kj_mol_per_water"])
        previous = float(row["previous_mixed_mae_kj_mol_per_water"])
        change = float(row["mixed_mae_change_kj_mol_per_water"])
        if not all(math.isfinite(value) for value in (current, previous, change)):
            raise AssertionError(f"non-finite aggregate in ledger row {index + 1}")
        if not math.isclose(current - previous, change, rel_tol=0.0, abs_tol=1.0e-12):
            raise AssertionError(f"incorrect MAE change in ledger row {index + 1}")
        if index and not math.isclose(
            previous,
            float(rows[index - 1]["mixed_mae_kj_mol_per_water"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise AssertionError(f"broken aggregate continuity in ledger row {index + 1}")

        unresolved = [item for item in row["unresolved_phases"].split(";") if item]
        converged = int(row["converged_phase_count"])
        if set(unresolved) - PHASES or converged + len(unresolved) != len(PHASES):
            raise AssertionError(f"invalid convergence partition in ledger row {index + 1}")

        run = RAW / f"k{mesh}{mesh}{mesh}-reduced" / phase
        output = run / "cp2k.out"
        if (
            not output.is_file()
            or sha256(output) != row["output_sha256"]
            or read_sidecar(run / "binary.sha256") != row["binary_sha256"]
            or read_sidecar(run / "input.sha256") != row["input_sha256"]
            or (run / "exit_status").read_text(encoding="utf-8").strip() != "0"
        ):
            raise AssertionError(f"raw proof mismatch in ledger row {index + 1}")
        text = output.read_text(encoding="utf-8", errors="replace")
        if "PROGRAM STARTED AT" not in text or "PROGRAM ENDED AT" not in text:
            raise AssertionError(f"missing CP2K lifecycle marker in ledger row {index + 1}")
        raw_proofs.append(run.relative_to(ROOT).as_posix())

    print(json.dumps({
        "schema": "periodic-gxtb-dmc13-completion-ledger-v1",
        "status": "PASS",
        "row_count": len(rows),
        "unique_notification_key_count": len(set(keys)),
        "unique_endpoint_count": len(set(endpoints)),
        "chronological": True,
        "aggregate_continuity": True,
        "raw_endpoint_proofs": raw_proofs,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
