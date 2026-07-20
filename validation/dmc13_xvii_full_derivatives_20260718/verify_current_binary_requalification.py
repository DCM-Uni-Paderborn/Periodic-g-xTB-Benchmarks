#!/usr/bin/env python3
"""Verify the frozen-binary ice-XVII force and stress requalification.

The campaign is deliberately independent of the archived qualification: all
six calculations are rerun with one SHA-256-qualified CP2K executable.  This
verifier checks the run provenance, full-grid/SPGLIB equivalence, and central
finite differences reconstructed directly from the retained input geometry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

from verify_full_derivatives import (
    ANGSTROM_TO_BOHR,
    HARTREE_PER_BOHR3_TO_BAR,
    cell_volume_bohr3,
    matrix_differences,
    parse_cp2k,
)


HERE = Path(__file__).resolve().parent
ARCHIVED_SUMMARY = HERE / "summary.json"
DEFAULT_BINARY_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
CASES = (
    "full",
    "reduced",
    "force-plus",
    "force-minus",
    "strain-plus",
    "strain-minus",
)

THRESHOLDS = {
    "full_minus_reduced_energy_Ha": 1.0e-12,
    "full_minus_reduced_maximum_force_Ha_per_bohr": 1.0e-9,
    "full_minus_reduced_maximum_stress_bar": 1.0e-5,
    "force_finite_difference_absolute_difference_Ha_per_bohr": 5.0e-7,
    "stress_finite_difference_absolute_difference_Ha": 5.0e-5,
    "archived_energy_absolute_difference_Ha": 5.0e-10,
    "archived_force_absolute_difference_Ha_per_bohr": 5.0e-9,
    "archived_virial_absolute_difference_Ha": 5.0e-7,
    "input_geometry_absolute_tolerance_A": 5.0e-10,
    "force_displacement_span_absolute_tolerance_A": 5.0e-10,
    "strain_span_absolute_tolerance": 5.0e-10,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_sha256(path: Path) -> str:
    match = re.match(r"\s*([0-9a-f]{64})(?:\s|$)", path.read_text(encoding="utf-8"))
    if match is None:
        raise RuntimeError(f"cannot parse SHA-256 from {path}")
    return match.group(1)


def parse_manifest(path: Path, root: Path) -> int:
    count = 0
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})\s+\*?(.+?)\s*", raw)
        if match is None:
            raise RuntimeError(f"invalid manifest line {line_number} in {path}")
        expected, relative = match.groups()
        candidate = (root / relative).resolve()
        if root.resolve() not in candidate.parents or not candidate.is_file():
            raise RuntimeError(f"manifest artifact is unavailable: {relative}")
        actual = sha256_file(candidate)
        if actual != expected:
            raise RuntimeError(
                f"manifest mismatch for {relative}: {actual} != {expected}"
            )
        count += 1
    return count


def parse_cell_and_coordinates(path: Path) -> tuple[list[list[float]], list[list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    cell: dict[str, list[float]] = {}
    coordinates: list[list[float]] = []
    in_cell = False
    in_coordinates = False
    scaled = False
    for raw in lines:
        stripped = raw.strip()
        upper = stripped.upper()
        if upper == "&CELL":
            in_cell = True
            continue
        if upper == "&COORD":
            in_coordinates = True
            continue
        if upper in ("&END CELL", "&ENDCELL"):
            in_cell = False
            continue
        if upper in ("&END COORD", "&ENDCOORD"):
            in_coordinates = False
            continue
        if in_cell:
            fields = stripped.split()
            if len(fields) >= 4 and fields[0].upper() in ("A", "B", "C"):
                cell[fields[0].upper()] = [float(value) for value in fields[1:4]]
        elif in_coordinates:
            if upper == "SCALED":
                scaled = True
                continue
            if not stripped or stripped.startswith(("#", "!")):
                continue
            fields = stripped.split()
            if len(fields) >= 4 and not fields[0].startswith("&"):
                coordinates.append([float(value) for value in fields[1:4]])
    if set(cell) != {"A", "B", "C"}:
        raise RuntimeError(f"cannot parse all cell vectors from {path}")
    if not scaled or not coordinates:
        raise RuntimeError(f"expected non-empty SCALED coordinates in {path}")
    return [cell[label] for label in ("A", "B", "C")], coordinates


def fractional_to_cartesian(
    cell: list[list[float]], coordinates: list[list[float]]
) -> list[list[float]]:
    return [
        [
            sum(fractional[axis] * cell[axis][component] for axis in range(3))
            for component in range(3)
        ]
        for fractional in coordinates
    ]


def maximum_matrix_difference(left: list[list[float]], right: list[list[float]]) -> float:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise RuntimeError("matrix dimensions differ")
    return max(
        abs(a - b)
        for left_row, right_row in zip(left, right)
        for a, b in zip(left_row, right_row)
    )


def infer_force_displacement(
    base_path: Path, plus_path: Path, minus_path: Path
) -> dict[str, float | int]:
    base_cell, base_fractional = parse_cell_and_coordinates(base_path)
    plus_cell, plus_fractional = parse_cell_and_coordinates(plus_path)
    minus_cell, minus_fractional = parse_cell_and_coordinates(minus_path)
    tolerance = THRESHOLDS["input_geometry_absolute_tolerance_A"]
    if maximum_matrix_difference(base_cell, plus_cell) > tolerance:
        raise RuntimeError("positive-displacement cell differs from the base cell")
    if maximum_matrix_difference(base_cell, minus_cell) > tolerance:
        raise RuntimeError("negative-displacement cell differs from the base cell")
    base = fractional_to_cartesian(base_cell, base_fractional)
    plus = fractional_to_cartesian(plus_cell, plus_fractional)
    minus = fractional_to_cartesian(minus_cell, minus_fractional)
    if len(base) != len(plus) or len(base) != len(minus):
        raise RuntimeError("force finite-difference atom counts differ")

    changed: list[tuple[int, int, float]] = []
    midpoint_maximum = 0.0
    for atom in range(len(base)):
        for component in range(3):
            span = plus[atom][component] - minus[atom][component]
            midpoint = 0.5 * (plus[atom][component] + minus[atom][component])
            midpoint_maximum = max(midpoint_maximum, abs(midpoint - base[atom][component]))
            if abs(span) > tolerance:
                changed.append((atom, component, span))
    if len(changed) != 1:
        raise RuntimeError(
            "force finite difference must change exactly one Cartesian component; "
            f"found {changed}"
        )
    if midpoint_maximum > tolerance:
        raise RuntimeError(
            f"force finite-difference midpoint differs from base by {midpoint_maximum:.3e} A"
        )
    atom, component, span = changed[0]
    if not math.isclose(
        abs(span),
        0.001,
        rel_tol=0.0,
        abs_tol=THRESHOLDS["force_displacement_span_absolute_tolerance_A"],
    ):
        raise RuntimeError(f"unexpected force finite-difference span: {span:.17g} A")
    return {
        "atom_zero_based": atom,
        "component_zero_based": component,
        "signed_span_A": span,
        "midpoint_maximum_difference_A": midpoint_maximum,
    }


def infer_xx_strain(
    base_path: Path, plus_path: Path, minus_path: Path
) -> dict[str, float]:
    base_cell, base_coordinates = parse_cell_and_coordinates(base_path)
    plus_cell, plus_coordinates = parse_cell_and_coordinates(plus_path)
    minus_cell, minus_coordinates = parse_cell_and_coordinates(minus_path)
    tolerance = THRESHOLDS["input_geometry_absolute_tolerance_A"]
    if maximum_matrix_difference(base_coordinates, plus_coordinates) > tolerance:
        raise RuntimeError("positive-strain scaled coordinates differ from base")
    if maximum_matrix_difference(base_coordinates, minus_coordinates) > tolerance:
        raise RuntimeError("negative-strain scaled coordinates differ from base")

    plus_strains: list[float] = []
    minus_strains: list[float] = []
    unchanged_maximum = 0.0
    for vector in range(3):
        if abs(base_cell[vector][0]) < 1.0e-12:
            raise RuntimeError("cannot infer xx strain from a zero cell-vector x component")
        plus_strains.append(plus_cell[vector][0] / base_cell[vector][0] - 1.0)
        minus_strains.append(minus_cell[vector][0] / base_cell[vector][0] - 1.0)
        for component in (1, 2):
            unchanged_maximum = max(
                unchanged_maximum,
                abs(plus_cell[vector][component] - base_cell[vector][component]),
                abs(minus_cell[vector][component] - base_cell[vector][component]),
            )
    if unchanged_maximum > tolerance:
        raise RuntimeError("strain inputs modify a non-x Cartesian cell component")
    plus_spread = max(plus_strains) - min(plus_strains)
    minus_spread = max(minus_strains) - min(minus_strains)
    if max(plus_spread, minus_spread) > 5.0e-10:
        raise RuntimeError("cell vectors do not share one homogeneous xx strain")
    plus_strain = sum(plus_strains) / 3.0
    minus_strain = sum(minus_strains) / 3.0
    span = plus_strain - minus_strain
    if not math.isclose(
        span,
        1.0e-4,
        rel_tol=0.0,
        abs_tol=THRESHOLDS["strain_span_absolute_tolerance"],
    ):
        raise RuntimeError(f"unexpected strain span: {span:.17g}")
    return {
        "positive_strain": plus_strain,
        "negative_strain": minus_strain,
        "signed_span": span,
        "maximum_non_x_cell_change_A": unchanged_maximum,
    }


def require_at_most(actual: float, threshold_key: str) -> None:
    threshold = THRESHOLDS[threshold_key]
    if actual > threshold:
        raise RuntimeError(f"{threshold_key}: {actual:.17g} > {threshold:.17g}")


def validate_inputs(campaign: Path) -> dict:
    manifest = campaign / "INPUT_SHA256SUMS"
    if not manifest.is_file():
        raise RuntimeError(f"missing input manifest: {manifest}")
    manifest_count = parse_manifest(manifest, campaign)
    inputs = {case: campaign / "inputs" / case / "input.inp" for case in CASES}
    for case, path in inputs.items():
        if not path.is_file():
            raise RuntimeError(f"missing {case} input: {path}")
    force_geometry = infer_force_displacement(
        inputs["full"], inputs["force-plus"], inputs["force-minus"]
    )
    strain_geometry = infer_xx_strain(
        inputs["full"], inputs["strain-plus"], inputs["strain-minus"]
    )
    return {
        "manifest_entries_verified": manifest_count,
        "input_sha256": {case: sha256_file(path) for case, path in inputs.items()},
        "force_geometry": force_geometry,
        "strain_geometry": strain_geometry,
    }


def validate_run_provenance(campaign: Path, expected_binary_sha256: str) -> dict:
    provenance = {}
    for case in CASES:
        run = campaign / "runs" / case
        output = run / "cp2k.out"
        if (run / "exit_status").read_text(encoding="utf-8").strip() != "0":
            raise RuntimeError(f"{case} did not exit successfully")
        parsed_binary_sha256 = first_sha256(run / "binary.sha256")
        if parsed_binary_sha256 != expected_binary_sha256:
            raise RuntimeError(
                f"{case} binary SHA-256 {parsed_binary_sha256} != {expected_binary_sha256}"
            )
        input_path = campaign / "inputs" / case / "input.inp"
        parsed_input_sha256 = first_sha256(run / "input.sha256")
        actual_input_sha256 = sha256_file(input_path)
        if parsed_input_sha256 != actual_input_sha256:
            raise RuntimeError(f"{case} run/input SHA-256 mismatch")
        parsed = parse_cp2k(output.read_bytes(), case)
        provenance[case] = {
            "binary_sha256": parsed_binary_sha256,
            "input_sha256": actual_input_sha256,
            "output_sha256": sha256_file(output),
            "energy_Ha": parsed["energy"],
        }
    return provenance


def evaluate(campaign: Path, expected_binary_sha256: str) -> dict:
    input_validation = validate_inputs(campaign)
    provenance = validate_run_provenance(campaign, expected_binary_sha256)
    parsed = {
        case: parse_cp2k((campaign / "runs" / case / "cp2k.out").read_bytes(), case)
        for case in CASES
    }
    full = parsed["full"]
    reduced = parsed["reduced"]
    if len(full["forces"]) != 18 or len(reduced["forces"]) != 18:
        raise RuntimeError("full and reduced ice-XVII runs must each contain 18 atoms")
    if len(full["stress"]) != 3 or len(reduced["stress"]) != 3:
        raise RuntimeError("full and reduced runs must each contain a 3x3 stress tensor")

    energy_difference = full["energy"] - reduced["energy"]
    force_maximum, force_rms = matrix_differences(full["forces"], reduced["forces"])
    stress_maximum, stress_rms = matrix_differences(full["stress"], reduced["stress"])
    require_at_most(abs(energy_difference), "full_minus_reduced_energy_Ha")
    require_at_most(
        force_maximum, "full_minus_reduced_maximum_force_Ha_per_bohr"
    )
    require_at_most(stress_maximum, "full_minus_reduced_maximum_stress_bar")

    force_geometry = input_validation["force_geometry"]
    force_atom = int(force_geometry["atom_zero_based"])
    force_component = int(force_geometry["component_zero_based"])
    force_span_bohr = float(force_geometry["signed_span_A"]) * ANGSTROM_TO_BOHR
    finite_difference_force = -(
        parsed["force-plus"]["energy"] - parsed["force-minus"]["energy"]
    ) / force_span_bohr
    analytic_force = full["forces"][force_atom][force_component]
    force_residual = abs(analytic_force - finite_difference_force)
    require_at_most(
        force_residual, "force_finite_difference_absolute_difference_Ha_per_bohr"
    )

    strain_geometry = input_validation["strain_geometry"]
    strain_span = float(strain_geometry["signed_span"])
    finite_difference_virial = (
        parsed["strain-plus"]["energy"] - parsed["strain-minus"]["energy"]
    ) / strain_span
    volume_bohr3 = cell_volume_bohr3(campaign / "inputs" / "full" / "input.inp")
    analytic_virial = (
        -full["stress"][0][0] * volume_bohr3 / HARTREE_PER_BOHR3_TO_BAR
    )
    virial_residual = abs(analytic_virial - finite_difference_virial)
    require_at_most(
        virial_residual, "stress_finite_difference_absolute_difference_Ha"
    )

    archived = json.loads(ARCHIVED_SUMMARY.read_text(encoding="utf-8"))["dmc_xvii_k222"]
    archived_differences = {
        "energy_Ha": abs(full["energy"] - archived["energy_Ha_per_primitive"]),
        "analytic_force_Ha_per_bohr": abs(
            analytic_force - archived["atom_1_x_analytic_force_Ha_per_bohr"]
        ),
        "analytic_virial_Ha": abs(analytic_virial - archived["xx_analytic_virial_Ha"]),
    }
    require_at_most(
        archived_differences["energy_Ha"], "archived_energy_absolute_difference_Ha"
    )
    require_at_most(
        archived_differences["analytic_force_Ha_per_bohr"],
        "archived_force_absolute_difference_Ha_per_bohr",
    )
    require_at_most(
        archived_differences["analytic_virial_Ha"],
        "archived_virial_absolute_difference_Ha",
    )

    return {
        "schema": "periodic-gxtb-part-i-current-binary-derivatives-v1",
        "status": "PASS",
        "expected_binary_sha256": expected_binary_sha256,
        "thresholds": THRESHOLDS,
        "integrity": {
            "input_manifest_entries_verified": input_validation[
                "manifest_entries_verified"
            ],
            "runs_verified": len(provenance),
            "input_sha256": input_validation["input_sha256"],
            "run_provenance": provenance,
        },
        "input_geometry": {
            "force": force_geometry,
            "strain": strain_geometry,
        },
        "full_minus_reduced": {
            "energy_Ha": energy_difference,
            "maximum_force_Ha_per_bohr": force_maximum,
            "rms_force_Ha_per_bohr": force_rms,
            "maximum_stress_bar": stress_maximum,
            "rms_stress_bar": stress_rms,
            "force_component_count": 54,
            "stress_component_count": 9,
        },
        "finite_difference": {
            "atom_one_based": force_atom + 1,
            "component_zero_based": force_component,
            "analytic_force_Ha_per_bohr": analytic_force,
            "finite_difference_force_Ha_per_bohr": finite_difference_force,
            "force_absolute_difference_Ha_per_bohr": force_residual,
            "xx_analytic_virial_Ha": analytic_virial,
            "xx_finite_difference_virial_Ha": finite_difference_virial,
            "xx_absolute_difference_Ha": virial_residual,
        },
        "archived_qualification_absolute_differences": archived_differences,
        "interpretation": (
            "All six calculations use one SHA-256-qualified CP2K executable. "
            "The explicit full mesh and SPGLIB-reduced mesh agree for energy, "
            "forces, and stress, while independently reconstructed central "
            "differences validate the force and homogeneous-strain response."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path)
    parser.add_argument(
        "--expected-binary-sha256", default=DEFAULT_BINARY_SHA256
    )
    parser.add_argument("--inputs-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    if args.inputs_only:
        result = {
            "schema": "periodic-gxtb-part-i-current-binary-inputs-v1",
            "status": "PASS",
            "input_validation": validate_inputs(campaign),
        }
    else:
        result = evaluate(campaign, args.expected_binary_sha256)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
