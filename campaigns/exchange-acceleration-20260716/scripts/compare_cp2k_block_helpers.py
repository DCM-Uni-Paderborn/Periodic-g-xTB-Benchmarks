#!/usr/bin/env python3
"""Deterministically compare CP2K block-helper runs with their full-array baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parent.parent
RAW = CAMPAIGN / "raw" / "cp2k_block_helpers"
FLOAT_RE = re.compile(
    r"(?<![A-Za-z_])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
)
PAIRS = (
    ("observables", "k290_222_force_stress"),
    ("observables", "spglib_shifted_222_force_stress"),
    ("observables", "time_reversal_311"),
    ("derivatives", "k290_222_force_stress_debug"),
    ("derivatives", "time_reversal_311_force_stress"),
)


def floats(text: str) -> list[float]:
    return [float(value.replace("D", "E").replace("d", "e")) for value in FLOAT_RE.findall(text)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def energies(lines: list[str]) -> list[float]:
    return [floats(line)[-1] for line in lines if "ENERGY| Total FORCE_EVAL" in line]


def kpoint_metadata(lines: list[str]) -> tuple[list[int] | None, int | None]:
    grid = None
    irreducible = None
    for line in lines:
        if grid is None and "BRILLOUIN| K-Point grid" in line:
            values = floats(line)
            if len(values) >= 3:
                grid = [int(value) for value in values[-3:]]
        if irreducible is None and "BRILLOUIN| List of Kpoints" in line:
            values = floats(line)
            if values:
                irreducible = int(values[-1])
    return grid, irreducible


def force_components(lines: list[str]) -> list[float]:
    result = []
    row = re.compile(r"^\s*FORCES\|\s+\d+\s+")
    for line in lines:
        if row.match(line):
            values = floats(line.split("|", 1)[-1])
            if len(values) >= 5:
                result.extend(values[1:5])
    return result


def analytical_stress_components(lines: list[str]) -> list[float]:
    result = []
    active = False
    row = re.compile(r"\s*[xyz]\s+")
    for line in lines:
        if "STRESS| Analytical stress tensor" in line:
            active = True
            continue
        if active:
            if not line.strip():
                active = False
                continue
            part = line.split("|", 1)[-1]
            if row.match(part):
                values = floats(part)
                if len(values) == 3:
                    result.extend(values)
    return result


def debug_virial_components(lines: list[str], header: str) -> list[float]:
    result = []
    active = False
    for line in lines:
        if header in line:
            active = True
            continue
        if active:
            if not line.strip():
                active = False
                continue
            if "DEBUG|" in line:
                values = floats(line.split("|", 1)[-1])
                if len(values) == 3:
                    result.extend(values)
    return result


def debug_force_summary(lines: list[str]) -> dict[str, list[float]]:
    result = {"numerical": [], "analytical": [], "difference": []}
    active = False
    row = re.compile(r"DEBUG\|\s+\d+\s+[xyz]\s+")
    for line in lines:
        if "BEGIN OF SUMMARY" in line:
            active = True
            continue
        if "END OF SUMMARY" in line:
            active = False
        if active and row.search(line):
            values = floats(line.split("|", 1)[-1])
            if len(values) >= 4:
                result["numerical"].append(values[1])
                result["analytical"].append(values[2])
                result["difference"].append(values[3])
    return result


def debug_difference_sums(lines: list[str]) -> dict[str, list[float]]:
    result = {"stress": [], "periodic_stress": [], "force": []}
    for line in lines:
        if "DEBUG| Periodic-subspace sum of differences" in line:
            values = floats(line)
            if values:
                result["periodic_stress"].append(values[-1])
        elif "DEBUG| Sum of differences:" in line:
            values = floats(line)
            if values:
                result["force"].append(values[0])
        elif "DEBUG| Sum of differences" in line:
            values = floats(line)
            if values:
                result["stress"].append(values[-1])
    return result


def max_abs_difference(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right):
        return None
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0)


def finite_difference_max(sums: dict[str, list[float]]) -> float:
    values = [abs(value) for entries in sums.values() for value in entries]
    return max(values, default=0.0)


def compare_pair(kind: str, case: str) -> dict[str, object]:
    current = RAW / f"current_{kind}" / case
    baseline = RAW / f"baseline_{kind}" / case
    current_lines = load(current / "cp2k.out")
    baseline_lines = load(baseline / "cp2k.out")

    current_grid, current_nirreducible = kpoint_metadata(current_lines)
    baseline_grid, baseline_nirreducible = kpoint_metadata(baseline_lines)
    current_forces = force_components(current_lines)
    baseline_forces = force_components(baseline_lines)
    current_stress = analytical_stress_components(current_lines)
    baseline_stress = analytical_stress_components(baseline_lines)
    current_force_debug = debug_force_summary(current_lines)
    baseline_force_debug = debug_force_summary(baseline_lines)
    current_sums = debug_difference_sums(current_lines)
    baseline_sums = debug_difference_sums(baseline_lines)

    differences = {
        "energy_sequence_hartree_max_abs": max_abs_difference(
            energies(current_lines), energies(baseline_lines)
        ),
        "atomic_force_hartree_per_bohr_max_abs": max_abs_difference(
            current_forces, baseline_forces
        ),
        "analytical_stress_bar_max_abs": max_abs_difference(current_stress, baseline_stress),
        "numerical_virial_hartree_max_abs": max_abs_difference(
            debug_virial_components(current_lines, "DEBUG| Numerical pv_virial"),
            debug_virial_components(baseline_lines, "DEBUG| Numerical pv_virial"),
        ),
        "analytical_virial_hartree_max_abs": max_abs_difference(
            debug_virial_components(current_lines, "DEBUG| Analytical pv_virial"),
            debug_virial_components(baseline_lines, "DEBUG| Analytical pv_virial"),
        ),
        "debug_force_numerical_hartree_per_bohr_max_abs": max_abs_difference(
            current_force_debug["numerical"], baseline_force_debug["numerical"]
        ),
        "debug_force_analytical_hartree_per_bohr_max_abs": max_abs_difference(
            current_force_debug["analytical"], baseline_force_debug["analytical"]
        ),
    }
    comparison_gate = all(value is not None for value in differences.values()) and (
        differences["energy_sequence_hartree_max_abs"] <= 1.0e-10
        and differences["atomic_force_hartree_per_bohr_max_abs"] <= 1.0e-10
        and differences["analytical_stress_bar_max_abs"] <= 1.0e-6
        and differences["numerical_virial_hartree_max_abs"] <= 1.0e-9
        and differences["analytical_virial_hartree_max_abs"] <= 1.0e-10
        and differences["debug_force_numerical_hartree_per_bohr_max_abs"] <= 1.0e-10
        and differences["debug_force_analytical_hartree_per_bohr_max_abs"] <= 1.0e-10
    )
    current_fd = finite_difference_max(current_sums)
    baseline_fd = finite_difference_max(baseline_sums)
    run_gate = (
        (current / "returncode.txt").read_text().strip() == "0"
        and (baseline / "returncode.txt").read_text().strip() == "0"
        and any("PROGRAM ENDED AT" in line for line in current_lines)
        and any("PROGRAM ENDED AT" in line for line in baseline_lines)
        and current_grid == baseline_grid
        and current_nirreducible == baseline_nirreducible
        and sha256(current / "input.inp") == sha256(baseline / "input.inp")
    )
    derivative_gate = kind != "derivatives" or (
        current_fd <= 1.0e-6 and baseline_fd <= 1.0e-6
    )

    return {
        "case": case,
        "kind": kind,
        "grid": current_grid,
        "irreducible_kpoints": current_nirreducible,
        "input_sha256": sha256(current / "input.inp"),
        "current_output_sha256": sha256(current / "cp2k.out"),
        "baseline_output_sha256": sha256(baseline / "cp2k.out"),
        "counts": {
            "energies": [len(energies(current_lines)), len(energies(baseline_lines))],
            "force_components": [len(current_forces), len(baseline_forces)],
            "stress_components": [len(current_stress), len(baseline_stress)],
        },
        "differences": differences,
        "finite_difference_sum_hartree_max_abs": {
            "current": current_fd,
            "baseline": baseline_fd,
        },
        "run_gate": run_gate,
        "comparison_gate": comparison_gate,
        "derivative_gate": derivative_gate,
        "passed": run_gate and comparison_gate and derivative_gate,
    }


def build_result() -> dict[str, object]:
    comparisons = [compare_pair(kind, case) for kind, case in PAIRS]
    difference_names = comparisons[0]["differences"].keys()
    maxima = {}
    for name in difference_names:
        values = [entry["differences"][name] for entry in comparisons]
        maxima[name] = max(value for value in values if value is not None)
    return {
        "schema_version": 1,
        "component": "cp2k-block-expansion-and-weighted-adjoint-foldback-helpers",
        "oracle": "unchanged-full-array-expansion-and-foldback",
        "internal_runtime_gates": {
            "physical_overlap_expansion_relative": "max(1e-6, 100*eps_geo)",
            "weighted_real_adjoint_relative": 1.0e-10,
            "blockwise_fold_vs_full_array_relative": 1.0e-12,
            "probe": "deterministic-complex-Hermitian",
            "antiunitary_time_reversal": "covered-by-negative-rotp-case",
        },
        "baseline_comparison_gates": {
            "energy_hartree_abs": 1.0e-10,
            "force_hartree_per_bohr_abs": 1.0e-10,
            "analytical_stress_bar_abs": 1.0e-6,
            "numerical_virial_hartree_abs": 1.0e-9,
            "analytical_virial_hartree_abs": 1.0e-10,
            "finite_difference_sum_hartree_abs": 1.0e-6,
        },
        "comparisons": comparisons,
        "maxima": maxima,
        "all_passed": all(entry["passed"] for entry in comparisons),
        "scope_note": (
            "This qualifies the CP2K one-block expansion and weighted-adjoint foldback "
            "helpers only; the production streamed exchange module is not qualified here."
        ),
    }


def render_text(result: dict[str, object]) -> str:
    lines = [
        "CP2K_BLOCK_HELPER_QUALIFICATION",
        "case\tgrid\tnirr\tdE_Ha\tdF_HaBohr\tdStress_bar\tdVirial_Ha\tFDmax_Ha\tpass",
    ]
    for entry in result["comparisons"]:
        diff = entry["differences"]
        dvirial = max(
            diff["numerical_virial_hartree_max_abs"],
            diff["analytical_virial_hartree_max_abs"],
        )
        grid = "x".join(str(value) for value in entry["grid"])
        fdmax = entry["finite_difference_sum_hartree_max_abs"]["current"]
        lines.append(
            f"{entry['case']}\t{grid}\t{entry['irreducible_kpoints']}\t"
            f"{diff['energy_sequence_hartree_max_abs']:.3e}\t"
            f"{diff['atomic_force_hartree_per_bohr_max_abs']:.3e}\t"
            f"{diff['analytical_stress_bar_max_abs']:.3e}\t{dvirial:.3e}\t"
            f"{fdmax:.3e}\t{str(entry['passed']).lower()}"
        )
    lines.extend(
        [
            "",
            "MAXIMA",
            *(
                f"{name}\t{value:.12e}"
                for name, value in sorted(result["maxima"].items())
            ),
            "",
            f"ALL_PASSED\t{str(result['all_passed']).lower()}",
            f"SCOPE\t{result['scope_note']}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "text"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build_result()
    if args.format == "json":
        rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    else:
        rendered = render_text(result)
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
