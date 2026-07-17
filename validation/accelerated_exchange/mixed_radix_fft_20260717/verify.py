#!/usr/bin/env python3
"""Recompute dense/separable/mixed-radix equivalence from archived CP2K logs."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FLOAT = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"
CASES = {
    "ch4": "ch4",
    "ch4_efs": "ch4_efs",
    "h2_1d": "h2_1d",
    "h2_1d_efs": "h2_1d_efs",
    "h2_2d": "h2_2d",
    "h2_2d_efs": "h2_2d_efs",
    "o2": "o2",
    "si": "si",
}
MODES = ("separable", "fft")


def parse(
    path: Path,
) -> tuple[list[float], list[tuple[float, ...]], list[tuple[float, ...]], str]:
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
    return energies, forces, stress, revision_match.group(1)


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
        "case\tbackend\tabs_energy_delta_hartree\t"
        "max_printed_energy_delta_hartree\t"
        "max_force_delta_hartree_per_bohr\tmax_stress_delta_bar\t"
        "dense_revision\tcandidate_revision"
    )
    for case, stem in CASES.items():
        dense_path = ROOT / "raw" / f"{stem}_dense" / f"{stem}_dense.out"
        dense_energy, dense_forces, dense_stress, dense_revision = parse(dense_path)
        for mode in MODES:
            candidate_path = ROOT / "raw" / f"{stem}_{mode}" / f"{stem}_{mode}.out"
            energy, forces, stress, revision = parse(candidate_path)
            if len(energy) != len(dense_energy):
                raise RuntimeError(f"energy traces have different lengths: {case}/{mode}")
            energy_delta = abs(energy[-1] - dense_energy[-1])
            energy_trace_delta = max(
                abs(candidate - oracle)
                for candidate, oracle in zip(energy, dense_energy, strict=True)
            )
            force_delta = maximum_delta(forces, dense_forces)
            stress_delta = maximum_delta(stress, dense_stress)
            if energy_delta > 1.0e-10:
                raise RuntimeError(f"energy equivalence failed: {case}/{mode}")
            if energy_trace_delta > 1.0e-10:
                raise RuntimeError(f"energy-trace equivalence failed: {case}/{mode}")
            if force_delta is not None and force_delta > 1.0e-8:
                raise RuntimeError(f"force equivalence failed: {case}/{mode}")
            if stress_delta is not None and stress_delta > 1.0e-3:
                raise RuntimeError(f"stress equivalence failed: {case}/{mode}")
            print(
                f"{case}\t{mode}\t{energy_delta:.17e}\t{energy_trace_delta:.17e}\t"
                f"{format_optional(force_delta)}\t{format_optional(stress_delta)}\t"
                f"{dense_revision}\t{revision}"
            )


if __name__ == "__main__":
    main()
