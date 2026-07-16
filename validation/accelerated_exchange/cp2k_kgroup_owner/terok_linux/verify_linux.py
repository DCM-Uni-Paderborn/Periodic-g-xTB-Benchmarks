#!/usr/bin/env python3
"""Re-evaluate the frozen Terok rank-count qualification."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from compare import maximum_delta, parse  # noqa: E402

ORACLE_RE = re.compile(
    r"KGROUP-ORACLE iter=1 dE=\s*([-+0-9.Ee]+) "
    r"dVsh=\s*([-+0-9.Ee]+) dFfold=\s*([-+0-9.Ee]+)"
)
GROUP_RE = re.compile(r"KGROUP-OWNER groups=(\d+), nred=(\d+), nfull=(\d+)")
MEMORY_RE = re.compile(r"Estimated peak process memory \[MiB\]\s+(\d+)")


def close(actual: float, expected: float, *, printed: bool = False) -> bool:
    if printed:
        return math.isclose(actual, expected, rel_tol=5.0e-7, abs_tol=5.0e-18)
    return math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-18)


def metrics(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text().splitlines():
        key, value = line.split("=", 1)
        values[key] = float(value)
    return values


def main() -> None:
    with (ROOT / "summary.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    reference = parse(ROOT / "runs/ch4_p1/cp2k.out")

    for row in rows:
        run = ROOT / "runs" / row["case"].lower().replace("ch4_", "ch4_")
        output = run / "cp2k.out"
        text = output.read_text(errors="replace")
        if text.count("PROGRAM ENDED") != 1 or "ABORT" in text:
            raise RuntimeError(f"{row['case']}: incomplete run")
        oracle = [tuple(map(float, item)) for item in ORACLE_RE.findall(text)]
        maxima = tuple(max(item[index] for item in oracle) for index in range(3))
        expected_internal = tuple(
            float(row[key])
            for key in (
                "max_internal_dE_Ha",
                "max_internal_dVsh_Ha",
                "max_internal_dFfold_Ha",
            )
        )
        if not all(close(a, e, printed=True) for a, e in zip(maxima, expected_internal)):
            raise RuntimeError(f"{row['case']}: oracle residual drift")
        if max(maxima) > 1.0e-10:
            raise RuntimeError(f"{row['case']}: oracle gate failed")

        group = (int(row["k_groups"]), int(row["nred"]), int(row["nfull"]))
        if not GROUP_RE.findall(text) or any(tuple(map(int, item)) != group for item in GROUP_RE.findall(text)):
            raise RuntimeError(f"{row['case']}: group metadata drift")

        result = parse(output)
        deltas = (
            abs(result[0] - reference[0]),
            maximum_delta(result[1], reference[1]),
            maximum_delta(result[2], reference[2]),
        )
        expected_deltas = tuple(
            float(row[key])
            for key in (
                "dEtot_vs_P1_Ha",
                "dForce_vs_P1_Ha_per_bohr",
                "dStress_vs_P1_bar",
            )
        )
        if not all(close(a, e) for a, e in zip(deltas, expected_deltas)):
            raise RuntimeError(f"{row['case']}: final observable drift")

        timing = metrics(run / "time.txt")
        for key in (
            "wall_seconds",
            "child_user_seconds",
            "child_system_seconds",
            "child_maxrss_kb",
        ):
            if not close(timing[key], float(row[key])):
                raise RuntimeError(f"{row['case']}: timing metadata drift")
        if timing["returncode"] != 0:
            raise RuntimeError(f"{row['case']}: nonzero return code")
        memory = MEMORY_RE.findall(text)
        if not memory or int(memory[-1]) != int(row["cp2k_peak_process_MiB"]):
            raise RuntimeError(f"{row['case']}: CP2K memory metadata drift")
        if row["status"] != "PASS":
            raise RuntimeError(f"{row['case']}: non-PASS status")
        print(f"PASS\t{row['case']}")


if __name__ == "__main__":
    main()
