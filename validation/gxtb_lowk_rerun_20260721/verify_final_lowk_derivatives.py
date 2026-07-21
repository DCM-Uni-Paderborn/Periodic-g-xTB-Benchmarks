#!/usr/bin/env python3
"""Verify the final-build low-k periodic g-xTB derivative qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


EXPECTED_BINARY = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
HARTREE_TO_KJ_MOL = 2625.4996394799
VIRIAL_HARTREE_TO_GPA_ANGSTROM3 = 4359.7447222071

DERIVATIVE_CASES = (
    "CH4_gxtb_gamma_force_stress",
    "CH4_gxtb_kp_111_force_stress",
    "CH4_gxtb_kp_full_222_force_stress",
    "CH4_gxtb_kp_k290_force_stress",
    "CH4_gxtb_kp_spglib_111_force_stress",
    "CH4_gxtb_kp_spglib_222_force_stress",
    "H2O_gxtb_gamma_force_stress",
    "H2O_gxtb_kp_gamma_force_stress",
    "H2_gxtb_gamma_supercell_311_force_stress",
    "H2_gxtb_kp_311_force_stress",
    "H2_gxtb_kp_311_tr_force_stress",
)

PARTIAL_CASES = (
    "gxtb_1d_native_gamma_centered_k211",
    "gxtb_1d_native_k211",
    "gxtb_1d_supercell_k111",
    "gxtb_1d_x_k211_force_stress",
    "gxtb_1d_x_k211_force_stress_full",
    "gxtb_1d_x_k211_force_stress_spglib",
    "gxtb_2d_native_gamma_centered_k212",
    "gxtb_2d_native_k212",
    "gxtb_2d_supercell_k111",
    "gxtb_2d_xz_k212_force_stress",
    "gxtb_2d_xz_k212_force_stress_full",
    "gxtb_2d_xz_k212_force_stress_spglib",
)

PARTIAL_DEBUG_CASES = tuple(
    case for case in PARTIAL_CASES if case.endswith(("force_stress", "force_stress_full", "force_stress_spglib"))
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_hash(path: Path) -> str:
    fields = path.read_text().split()
    if not fields:
        raise ValueError(f"empty hash record: {path}")
    return fields[0]


def first_float(pattern: str, text: str, label: str) -> float:
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        raise ValueError(f"missing {label}")
    return float(match.group(1))


def parse_output(path: Path, debug_required: bool) -> dict[str, float | bool | None]:
    text = path.read_text(errors="replace")
    energy = first_float(
        r"^\s*ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$",
        text,
        "total energy",
    )
    volume = first_float(
        r"^\s*CELL\| Volume \[angstrom\^3\]:\s*([-+0-9.Ee]+)\s*$",
        text,
        "cell volume",
    )
    force_matches = re.findall(
        r"^\s*DEBUG\| Sum of differences:\s*([-+0-9.Ee]+)", text, re.MULTILINE
    )
    virial_matches = re.findall(
        r"^\s*DEBUG\| Sum of differences\s+([-+0-9.Ee]+)\s*$", text, re.MULTILINE
    )
    if debug_required and (not force_matches or not virial_matches):
        raise ValueError(f"missing force/virial finite-difference summary in {path}")
    force = float(force_matches[-1]) if force_matches else None
    virial = float(virial_matches[-1]) if virial_matches else None
    stress = (
        virial * VIRIAL_HARTREE_TO_GPA_ANGSTROM3 / volume
        if virial is not None
        else None
    )
    return {
        "program_ended": "PROGRAM ENDED AT" in text,
        "energy_hartree": energy,
        "volume_angstrom3": volume,
        "force_sum_abs_hartree_per_bohr": force,
        "virial_sum_abs_hartree": virial,
        "stress_sum_abs_gpa": stress,
    }


def verify_case(
    root: Path, suite: str, case: str, debug_required: bool, allowed_cpus: set[int]
) -> dict[str, object]:
    result = root / "results" / suite / case
    input_path = root / "inputs" / suite / f"{case}.inp"
    output_path = result / "cp2k.out"
    if (result / "exit_status").read_text().strip() != "0":
        raise ValueError(f"nonzero exit status: {suite}/{case}")
    if first_hash(result / "binary.sha256") != EXPECTED_BINARY:
        raise ValueError(f"wrong binary identity: {suite}/{case}")
    if first_hash(result / "input.sha256") != sha256(input_path):
        raise ValueError(f"wrong input identity: {suite}/{case}")
    proof = (result / "affinity_preexec.txt").read_text()
    match = re.search(r"expected_cpu=(\d+) allowed=([0-9,-]+)", proof)
    if match is None or match.group(1) != match.group(2):
        raise ValueError(f"invalid singleton affinity proof: {suite}/{case}")
    cpu = int(match.group(1))
    if cpu not in allowed_cpus:
        raise ValueError(f"unexpected CPU {cpu}: {suite}/{case}")
    parsed = parse_output(output_path, debug_required)
    if not parsed["program_ended"]:
        raise ValueError(f"missing normal termination: {suite}/{case}")
    return {
        **parsed,
        "cpu": cpu,
        "input_sha256": sha256(input_path),
        "output_sha256": sha256(output_path),
    }


def delta(values: dict[str, dict[str, object]], left: str, right: str, divisor: float = 1.0) -> float:
    left_energy = float(values[left]["energy_hartree"])
    right_energy = float(values[right]["energy_hartree"])
    return left_energy / divisor - right_energy


def require_at_most(value: float, limit: float, label: str) -> None:
    if not math.isfinite(value) or abs(value) > limit:
        raise ValueError(f"{label}: |{value:.12e}| > {limit:.12e}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--legacy-partial-manifest", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    if (root / "controller.exit_status").read_text().strip() != "0":
        raise ValueError("low-k controller did not complete successfully")
    if args.binary is not None and sha256(args.binary) != EXPECTED_BINARY:
        raise ValueError("live CP2K executable does not match the qualified identity")

    identity = {}
    for line in (root / "run_identity.txt").read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            identity[key] = value
    try:
        allowed_cpus = {int(identity["cpu"])}
    except (KeyError, ValueError) as exc:
        raise ValueError("missing or invalid run_identity CPU") from exc
    derivative = {
        case: verify_case(root, "derivative", case, True, allowed_cpus)
        for case in DERIVATIVE_CASES
    }
    partial = {
        case: verify_case(
            root, "partial", case, case in PARTIAL_DEBUG_CASES, allowed_cpus
        )
        for case in PARTIAL_CASES
    }

    energy_checks = {
        "h2o_implicit_minus_explicit_gamma_hartree": delta(
            derivative,
            "H2O_gxtb_gamma_force_stress",
            "H2O_gxtb_kp_gamma_force_stress",
        ),
        "h2_supercell_per_primitive_minus_full_k311_hartree": delta(
            derivative,
            "H2_gxtb_gamma_supercell_311_force_stress",
            "H2_gxtb_kp_311_force_stress",
            3.0,
        ),
        "h2_full_minus_time_reversal_k311_hartree": delta(
            derivative,
            "H2_gxtb_kp_311_force_stress",
            "H2_gxtb_kp_311_tr_force_stress",
        ),
        "ch4_full_minus_k290_k222_hartree": delta(
            derivative,
            "CH4_gxtb_kp_full_222_force_stress",
            "CH4_gxtb_kp_k290_force_stress",
        ),
        "ch4_full_minus_spglib_k222_hartree": delta(
            derivative,
            "CH4_gxtb_kp_full_222_force_stress",
            "CH4_gxtb_kp_spglib_222_force_stress",
        ),
        "one_dimensional_supercell_per_primitive_minus_gamma_centered_hartree": delta(
            partial,
            "gxtb_1d_supercell_k111",
            "gxtb_1d_native_gamma_centered_k211",
            2.0,
        ),
        "two_dimensional_supercell_per_primitive_minus_gamma_centered_hartree": delta(
            partial,
            "gxtb_2d_supercell_k111",
            "gxtb_2d_native_gamma_centered_k212",
            4.0,
        ),
        "one_dimensional_full_minus_k290_hartree": delta(
            partial,
            "gxtb_1d_x_k211_force_stress_full",
            "gxtb_1d_x_k211_force_stress",
        ),
        "one_dimensional_full_minus_spglib_hartree": delta(
            partial,
            "gxtb_1d_x_k211_force_stress_full",
            "gxtb_1d_x_k211_force_stress_spglib",
        ),
        "two_dimensional_full_minus_k290_hartree": delta(
            partial,
            "gxtb_2d_xz_k212_force_stress_full",
            "gxtb_2d_xz_k212_force_stress",
        ),
        "two_dimensional_full_minus_spglib_hartree": delta(
            partial,
            "gxtb_2d_xz_k212_force_stress_full",
            "gxtb_2d_xz_k212_force_stress_spglib",
        ),
    }

    for label in (
        "h2o_implicit_minus_explicit_gamma_hartree",
        "h2_supercell_per_primitive_minus_full_k311_hartree",
    ):
        require_at_most(energy_checks[label], 1.0e-9, label)
    for label in (
        "h2_full_minus_time_reversal_k311_hartree",
        "ch4_full_minus_k290_k222_hartree",
        "ch4_full_minus_spglib_k222_hartree",
        "one_dimensional_full_minus_k290_hartree",
        "one_dimensional_full_minus_spglib_hartree",
        "two_dimensional_full_minus_k290_hartree",
        "two_dimensional_full_minus_spglib_hartree",
    ):
        require_at_most(energy_checks[label], 5.0e-12, label)
    for label in (
        "one_dimensional_supercell_per_primitive_minus_gamma_centered_hartree",
        "two_dimensional_supercell_per_primitive_minus_gamma_centered_hartree",
    ):
        require_at_most(energy_checks[label], 1.0e-8, label)

    derivative_force_max = max(
        float(value["force_sum_abs_hartree_per_bohr"]) for value in derivative.values()
    )
    derivative_virial_max = max(
        float(value["virial_sum_abs_hartree"]) for value in derivative.values()
    )
    derivative_stress_max = max(
        float(value["stress_sum_abs_gpa"]) for value in derivative.values()
    )
    partial_debug = [partial[case] for case in PARTIAL_DEBUG_CASES]
    partial_force_max = max(
        float(value["force_sum_abs_hartree_per_bohr"]) for value in partial_debug
    )
    partial_virial_max = max(
        float(value["virial_sum_abs_hartree"]) for value in partial_debug
    )
    partial_stress_max = max(float(value["stress_sum_abs_gpa"]) for value in partial_debug)
    require_at_most(derivative_force_max, 2.0e-8, "derivative force maximum")
    require_at_most(derivative_virial_max, 1.0e-7, "derivative virial maximum")
    require_at_most(partial_force_max, 3.0e-8, "partial-PBC force maximum")
    require_at_most(partial_virial_max, 2.0e-7, "partial-PBC virial maximum")

    bvk = {
        label: {
            "hartree": energy_checks[label],
            "kj_mol_per_primitive": energy_checks[label] * HARTREE_TO_KJ_MOL,
        }
        for label in (
            "one_dimensional_supercell_per_primitive_minus_gamma_centered_hartree",
            "two_dimensional_supercell_per_primitive_minus_gamma_centered_hartree",
        )
    }
    legacy_bvk_comparison = None
    if args.legacy_partial_manifest is not None:
        legacy = json.loads(args.legacy_partial_manifest.read_text())
        legacy_bvk = legacy["bvk_supercell_diagnostic"]
        legacy_values = {
            "one_dimensional": float(
                legacy_bvk["one_dimensional_supercell_per_primitive_minus_native_hartree"]
            ),
            "two_dimensional": float(
                legacy_bvk["two_dimensional_supercell_per_primitive_minus_native_hartree"]
            ),
        }
        current_values = {
            "one_dimensional": energy_checks[
                "one_dimensional_supercell_per_primitive_minus_gamma_centered_hartree"
            ],
            "two_dimensional": energy_checks[
                "two_dimensional_supercell_per_primitive_minus_gamma_centered_hartree"
            ],
        }
        legacy_bvk_comparison = {}
        for dimension in ("one_dimensional", "two_dimensional"):
            old = legacy_values[dimension]
            new = current_values[dimension]
            legacy_bvk_comparison[dimension] = {
                "legacy_hartree": old,
                "current_hartree": new,
                "absolute_residual_reduction_hartree": abs(old) - abs(new),
                "current_over_legacy_absolute_residual": (
                    abs(new) / abs(old) if old != 0.0 else None
                ),
                "legacy_over_current_absolute_residual": (
                    abs(old) / abs(new) if new != 0.0 else None
                ),
            }
    summary = {
        "schema_version": 1,
        "status": "passed",
        "qualified_binary_sha256": EXPECTED_BINARY,
        "counts": {
            "total_cases": len(derivative) + len(partial),
            "derivative_cases": len(derivative),
            "partial_pbc_cases": len(partial),
            "normal_terminations": sum(
                bool(value["program_ended"])
                for value in tuple(derivative.values()) + tuple(partial.values())
            ),
        },
        "maxima": {
            "derivative_force_sum_abs_hartree_per_bohr": derivative_force_max,
            "derivative_virial_sum_abs_hartree": derivative_virial_max,
            "derivative_stress_sum_abs_gpa": derivative_stress_max,
            "partial_pbc_force_sum_abs_hartree_per_bohr": partial_force_max,
            "partial_pbc_virial_sum_abs_hartree": partial_virial_max,
            "partial_pbc_stress_sum_abs_gpa": partial_stress_max,
        },
        "energy_route_differences": energy_checks,
        "bvk_supercell_equivalence": bvk,
        "legacy_bvk_comparison": legacy_bvk_comparison,
        "legacy_partial_manifest_sha256": (
            sha256(args.legacy_partial_manifest)
            if args.legacy_partial_manifest is not None
            else None
        ),
        "derivative_cases": derivative,
        "partial_pbc_cases": partial,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
