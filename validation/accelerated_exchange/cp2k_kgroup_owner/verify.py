#!/usr/bin/env python3
"""Re-evaluate the frozen KGROUP_OWNER qualification table."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

from compare import maximum_delta, parse


ROOT = Path(__file__).resolve().parent
CASES = {
    "CH4_full_2x2x2_P2": (
        "runs/ch4_full_p2_oracle_final/cp2k.out",
        "references/ch4_full_reference.out",
    ),
    "CH4_full_2x2x2_Release_P2": (
        "runs/ch4_full_release_p2_oracle_final/cp2k.out",
        "references/ch4_full_reference.out",
    ),
    "CH4_full_2x2x2_P4": (
        "runs/ch4_full_p4_oracle_final/cp2k.out",
        "references/ch4_full_reference.out",
    ),
    "CH4_Gamma_P2": (
        "runs/ch4_k1_p2_oracle_final/cp2k.out",
        "references/ch4_k1_reference.out",
    ),
    "H2_3x1x1_TR_P2": (
        "runs/h2_k311_tr_p2_oracle_final/cp2k.out",
        "references/h2_k311_tr_reference.out",
    ),
    "Ar2_1D_2x1x1_P2": (
        "runs/ar2_1d_p2_oracle_final/cp2k.out",
        "references/ar2_1d_reference.out",
    ),
    "Ar4_2D_2x2x1_P2": (
        "runs/ar4_2d_p2_oracle_final/cp2k.out",
        "references/ar4_2d_reference.out",
    ),
    "CH4_SPGLIB_2x2x2_P2": (
        "runs/ch4_spglib_p2_oracle_final/cp2k.out",
        "references/ch4_spglib_reference.out",
    ),
    "CH4_K290_2x2x2_P2": (
        "runs/ch4_k290_p2_oracle_final/cp2k.out",
        "references/ch4_k290_reference.out",
    ),
    "Si_shifted_2x2x2_P2": (
        "runs/si_shifted_p2_oracle_final/cp2k.out",
        "references/si_shifted_reference.out",
    ),
    "O2_UKS_3x1x1_TR_P2": (
        "runs/o2_uks_p2_oracle_final/cp2k.out",
        "references/o2_uks_reference.out",
    ),
}

ORACLE_RE = re.compile(
    r"KGROUP-ORACLE iter=1 dE=\s*([-+0-9.Ee]+) "
    r"dVsh=\s*([-+0-9.Ee]+) dFfold=\s*([-+0-9.Ee]+)"
)
GROUP_RE = re.compile(r"KGROUP-OWNER groups=(\d+), nred=(\d+), nfull=(\d+)")


def close(actual: float, expected: float, *, printed: bool = False) -> bool:
    if printed:
        return math.isclose(actual, expected, rel_tol=5.0e-7, abs_tol=5.0e-17)
    return math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-18)


def main() -> None:
    with (ROOT / "qualification_summary.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if {row["case"] for row in rows} != set(CASES):
        raise RuntimeError("table cases and verifier cases differ")

    for row in rows:
        candidate_rel, reference_rel = CASES[row["case"]]
        candidate_path = ROOT / candidate_rel
        text = candidate_path.read_text(errors="replace")
        if text.count("PROGRAM ENDED") != 1:
            raise RuntimeError(f"{row['case']}: expected one complete run")
        oracle = [tuple(map(float, match)) for match in ORACLE_RE.findall(text)]
        if not oracle:
            raise RuntimeError(f"{row['case']}: missing in-process oracle")
        maxima = tuple(max(values[index] for values in oracle) for index in range(3))
        expected_internal = tuple(
            float(row[name])
            for name in (
                "max_internal_dE_Ha",
                "max_internal_dVsh_Ha",
                "max_internal_dFfold_Ha",
            )
        )
        if not all(close(a, e, printed=True) for a, e in zip(maxima, expected_internal)):
            raise RuntimeError(f"{row['case']}: in-process residual differs from table")
        if max(maxima) > 1.0e-10:
            raise RuntimeError(f"{row['case']}: in-process gate failed")

        group_data = [tuple(map(int, match)) for match in GROUP_RE.findall(text)]
        expected_group = (int(row["k_groups"]), int(row["nred"]), int(row["nfull"]))
        if not group_data or any(value != expected_group for value in group_data):
            raise RuntimeError(f"{row['case']}: communicator/mesh metadata differs")

        candidate = parse(candidate_path)
        reference = parse(ROOT / reference_rel)
        external = (
            abs(candidate[0] - reference[0]),
            maximum_delta(candidate[1], reference[1]),
            maximum_delta(candidate[2], reference[2]),
        )
        expected_external = tuple(
            float(row[name])
            for name in (
                "external_dEtot_Ha",
                "external_dForce_Ha_per_bohr",
                "external_dStress_bar",
            )
        )
        if not all(close(a, e) for a, e in zip(external, expected_external)):
            raise RuntimeError(f"{row['case']}: final observable delta differs from table")
        if row["status"] != "PASS":
            raise RuntimeError(f"{row['case']}: non-PASS table status")
        print(f"PASS\t{row['case']}")


if __name__ == "__main__":
    main()
