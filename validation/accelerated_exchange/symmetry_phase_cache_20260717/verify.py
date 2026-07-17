#!/usr/bin/env python3
"""Recompute phase-cache/reference equivalence from archived CP2K logs."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FLOAT = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"
CASES = (
    "CH4_gxtb_kp_k290",
    "CH4_gxtb_kp_k290_force_stress",
    "CH4_gxtb_kp_spglib",
    "CH4_gxtb_kp_spglib_111_force_stress",
    "H2_gxtb_kp_311_tr_force_stress",
    "O2_gxtb_uks_kp_311_tr",
    "Si_prim_gxtb_kp_shifted_spglib",
)


def parse(
    path: Path,
) -> tuple[
    list[float],
    list[tuple[float, ...]],
    list[tuple[float, ...]],
    str,
    list[str],
]:
    text = path.read_text(errors="replace")
    if "PROGRAM ENDED AT" not in text:
        raise RuntimeError(f"incomplete CP2K output: {path}")
    energies = [
        float(value)
        for value in re.findall(
            rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})",
            text,
        )
    ]
    if not energies:
        raise RuntimeError(f"missing final energy: {path}")
    forces = [
        tuple(float(value) for value in match)
        for match in re.findall(
            rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
            text,
            re.MULTILINE,
        )
    ]
    blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    stress: list[tuple[float, ...]] = []
    if blocks:
        stress = [
            tuple(float(value) for value in match)
            for match in re.findall(
                rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
                blocks[-1],
                re.MULTILINE,
            )
        ]
        if len(stress) != 3:
            raise RuntimeError(f"incomplete stress block: {path}")
    revision_match = re.search(r"source code revision number:\s+(\S+)", text)
    if revision_match is None:
        raise RuntimeError(f"missing CP2K source revision: {path}")
    scientific_lines = [
        line
        for line in text.splitlines()
        if line.startswith(" ENERGY|") or line.startswith(" DEBUG|")
    ]
    return energies, forces, stress, revision_match.group(1), scientific_lines


def maximum_delta(
    left: list[tuple[float, ...]], right: list[tuple[float, ...]]
) -> float | None:
    if len(left) != len(right):
        raise RuntimeError("observable blocks have different lengths")
    if not left:
        return None
    return max(abs(x - y) for a, b in zip(left, right) for x, y in zip(a, b))


def format_optional(value: float | None) -> str:
    return "NA" if value is None else f"{value:.17e}"


def main() -> None:
    print(
        "case\tabs_energy_delta_hartree\tmax_force_delta_hartree_per_bohr\t"
        "max_stress_delta_bar\tscientific_line_count\tscientific_lines_identical\t"
        "reference_revision\tcandidate_revision"
    )
    for case in CASES:
        reference = parse(ROOT / "raw" / "reference" / case / "run.out")
        candidate = parse(ROOT / "raw" / "candidate" / case / "run.out")
        if len(candidate[0]) != len(reference[0]):
            raise RuntimeError(f"energy traces have different lengths: {case}")
        energy_delta = max(
            abs(value - oracle)
            for value, oracle in zip(candidate[0], reference[0], strict=True)
        )
        force_delta = maximum_delta(candidate[1], reference[1])
        stress_delta = maximum_delta(candidate[2], reference[2])
        if energy_delta > 1.0e-10:
            raise RuntimeError(f"energy equivalence failed: {case}")
        if force_delta is not None and force_delta > 1.0e-8:
            raise RuntimeError(f"force equivalence failed: {case}")
        if stress_delta is not None and stress_delta > 1.0e-3:
            raise RuntimeError(f"stress equivalence failed: {case}")
        scientific_lines_identical = candidate[4] == reference[4]
        if not scientific_lines_identical:
            raise RuntimeError(f"ENERGY/DEBUG scientific trace changed: {case}")
        print(
            f"{case}\t{energy_delta:.17e}\t{format_optional(force_delta)}\t"
            f"{format_optional(stress_delta)}\t{len(candidate[4])}\t"
            f"{str(scientific_lines_identical).lower()}\t"
            f"{reference[3]}\t{candidate[3]}"
        )


if __name__ == "__main__":
    main()
