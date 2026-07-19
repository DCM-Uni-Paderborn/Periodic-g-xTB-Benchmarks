#!/usr/bin/env python3
"""Qualify direct save_tblite and native CP2K energies on one BvK mesh."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from verify_k222_cli_native_requalification import (
    DIRECT_DISPERSION_RE,
    HARTREE_TO_KJ_MOL,
    NATIVE_DISPERSION_RE,
    PHASES,
    STRUCTURE_TOLERANCE_ANGSTROM,
    digest,
    inverse_3x3,
    parse_cp2k_structure,
    parse_poscar_structure,
    qualify_affinity,
    read_cli_energy,
    read_component,
    read_native_energy,
    recorded_digest,
    require_status_zero,
    row_times_matrix,
    vector_scale,
    verify_source_identity,
)


SCHEME_RE = re.compile(
    r"^\s*SCHEME\s+MACDONALD\s+(\d+)\s+(\d+)\s+(\d+)\s+"
    r"([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s*$",
    re.IGNORECASE | re.MULTILINE,
)


def parse_phase_list(value: str) -> tuple[str, ...]:
    phases = tuple(item.strip() for item in value.split(",") if item.strip())
    if not phases or len(set(phases)) != len(phases):
        raise argparse.ArgumentTypeError("phase list must be nonempty and unique")
    unknown = set(phases) - set(PHASES)
    if unknown:
        raise argparse.ArgumentTypeError(f"unknown phases: {sorted(unknown)}")
    if "Ih" not in phases:
        raise argparse.ArgumentTypeError("phase list must include Ih")
    return phases


def parse_cpu_set(value: str) -> frozenset[int]:
    try:
        cpus = frozenset(int(item.strip()) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("CPU list must contain integers") from error
    if not cpus or min(cpus) < 0:
        raise argparse.ArgumentTypeError("CPU list must be nonempty and nonnegative")
    return cpus


def macdonald_shift(mesh_size: int) -> float:
    return 0.0 if mesh_size % 2 else (mesh_size - 1.0) / (2.0 * mesh_size)


def verify_native_input(path: Path, mesh_size: int) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = list(SCHEME_RE.finditer(text))
    if len(matches) != 1:
        raise AssertionError(f"expected one MacDonald mesh: {path}")
    match = matches[0]
    dimensions = tuple(int(match.group(index)) for index in range(1, 4))
    shifts = tuple(
        float(match.group(index).replace("D", "E").replace("d", "e"))
        for index in range(4, 7)
    )
    expected_shift = macdonald_shift(mesh_size)
    if dimensions != (mesh_size,) * 3 or any(
        abs(value - expected_shift) > 5.0e-13 for value in shifts
    ):
        raise AssertionError(
            f"noncanonical {mesh_size}x{mesh_size}x{mesh_size} MacDonald mesh: {path}"
        )
    for keyword, expected in (("SYMMETRY", "T"), ("FULL_GRID", "F")):
        if (
            re.search(
                rf"^\s*{keyword}\s+{expected}\s*$",
                text,
                re.IGNORECASE | re.MULTILINE,
            )
            is None
        ):
            raise AssertionError(f"wrong {keyword} setting: {path}")
    return digest(path)


def verify_structure_mapping(
    cp2k_input: Path, poscar: Path, mesh_size: int
) -> tuple[float, float]:
    primitive_cell, primitive_coordinates = parse_cp2k_structure(cp2k_input)
    supercell, explicit_coordinates = parse_poscar_structure(poscar)
    expected_supercell = [
        vector_scale(vector, float(mesh_size)) for vector in primitive_cell
    ]
    cell_residual = max(
        abs(supercell[i][j] - expected_supercell[i][j])
        for i in range(3)
        for j in range(3)
    )
    if cell_residual > STRUCTURE_TOLERANCE_ANGSTROM:
        raise AssertionError(
            f"native/direct {mesh_size}^3 cell mismatch: "
            f"{cell_residual:.6e} Angstrom"
        )

    expected: dict[str, list[tuple[float, float, float]]] = {}
    for element, fractional in primitive_coordinates:
        wrapped = tuple(value % 1.0 for value in fractional)
        for iz in range(mesh_size):
            for iy in range(mesh_size):
                for ix in range(mesh_size):
                    expected.setdefault(element, []).append(
                        tuple(
                            (wrapped[axis] + shift) / mesh_size
                            for axis, shift in enumerate((ix, iy, iz))
                        )
                    )
    inverse_supercell = inverse_3x3(supercell)
    actual: dict[str, list[tuple[float, float, float]]] = {}
    for element, cartesian in explicit_coordinates:
        fractional = row_times_matrix(cartesian, inverse_supercell)
        actual.setdefault(element, []).append(
            tuple(value % 1.0 for value in fractional)
        )
    if {key: len(value) for key, value in expected.items()} != {
        key: len(value) for key, value in actual.items()
    }:
        raise AssertionError("native/direct species multiplicity mismatch")

    coordinate_residual = 0.0
    for element, expected_positions in expected.items():
        remaining = list(actual[element])
        for expected_position in expected_positions:
            best_index = -1
            best_distance = math.inf
            for index, actual_position in enumerate(remaining):
                delta_fractional = tuple(
                    (actual_position[j] - expected_position[j])
                    - round(actual_position[j] - expected_position[j])
                    for j in range(3)
                )
                delta_cartesian = row_times_matrix(delta_fractional, supercell)
                distance = math.sqrt(sum(value * value for value in delta_cartesian))
                if distance < best_distance:
                    best_distance = distance
                    best_index = index
            if best_index < 0 or best_distance > STRUCTURE_TOLERANCE_ANGSTROM:
                raise AssertionError(
                    f"native/direct coordinate mismatch for {element}: "
                    f"{best_distance:.6e} Angstrom"
                )
            coordinate_residual = max(coordinate_residual, best_distance)
            remaining.pop(best_index)
    return cell_residual, coordinate_residual


def water_count(structure: Path, replicas: int) -> int:
    lines = structure.read_text(encoding="utf-8").splitlines()
    if len(lines) < 7:
        raise AssertionError(f"incomplete POSCAR: {structure}")
    supercell_atoms = sum(int(value) for value in lines[6].split())
    denominator = 3 * replicas
    if supercell_atoms % denominator:
        raise AssertionError(f"atom count is incompatible with the mesh: {structure}")
    return supercell_atoms // denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive_root", type=Path)
    parser.add_argument("direct_root", type=Path)
    parser.add_argument("native_run_root", type=Path)
    parser.add_argument("native_input_root", type=Path)
    parser.add_argument("--mesh-size", type=int, required=True)
    parser.add_argument("--phases", type=parse_phase_list, default=PHASES)
    parser.add_argument("--direct-controller-status", type=Path, required=True)
    parser.add_argument("--native-controller-status", type=Path, required=True)
    parser.add_argument("--source-identity", type=Path, required=True)
    parser.add_argument("--expected-source-revision", required=True)
    parser.add_argument("--expected-direct-binary", required=True)
    parser.add_argument("--expected-native-binary", required=True)
    parser.add_argument("--expected-native-provider-archive", required=True)
    parser.add_argument("--expected-native-cp2k-revision", required=True)
    parser.add_argument("--expected-native-cmake-cache", required=True)
    parser.add_argument("--expected-native-build-ninja", required=True)
    parser.add_argument("--expected-direct-cpus", type=parse_cpu_set, required=True)
    parser.add_argument("--tolerance-ha", type=float, default=2.0e-7)
    parser.add_argument("--component-tolerance-ha", type=float, default=2.0e-7)
    parser.add_argument("--relative-tolerance-kj-mol", type=float, default=5.0e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.mesh_size < 1:
        parser.error("mesh size must be positive")
    for label, value, length in (
        ("direct binary", args.expected_direct_binary, 64),
        ("native binary", args.expected_native_binary, 64),
        ("native provider archive", args.expected_native_provider_archive, 64),
        ("native CMake cache", args.expected_native_cmake_cache, 64),
        ("native build.ninja", args.expected_native_build_ninja, 64),
        ("source revision", args.expected_source_revision, 40),
        ("native CP2K revision", args.expected_native_cp2k_revision, 40),
    ):
        if re.fullmatch(rf"[0-9a-f]{{{length}}}", value) is None:
            parser.error(f"invalid {label} digest")
    if min(
        args.tolerance_ha,
        args.component_tolerance_ha,
        args.relative_tolerance_kj_mol,
    ) <= 0.0:
        parser.error("energy tolerances must be positive")

    require_status_zero(args.direct_controller_status, "direct controller status")
    require_status_zero(args.native_controller_status, "native controller status")
    source_hash = verify_source_identity(
        args.source_identity,
        args.expected_source_revision,
        args.expected_direct_binary,
        args.expected_native_binary,
        args.expected_native_provider_archive,
        args.expected_native_cp2k_revision,
        args.expected_native_cmake_cache,
        args.expected_native_build_ninja,
    )

    mesh_id = f"k{args.mesh_size}{args.mesh_size}{args.mesh_size}"
    replicas = args.mesh_size**3
    rows: list[dict[str, object]] = []
    for phase in args.phases:
        structure = args.archive_root / "structures" / mesh_id / phase / "POSCAR"
        direct_dir = args.direct_root / mesh_id / phase
        native_dir = args.native_run_root / phase
        native_input = args.native_input_root / phase / "input.inp"
        if not structure.is_file() or not native_input.is_file():
            raise AssertionError(f"missing frozen input for phase {phase}")
        require_status_zero(direct_dir / "exit_status", f"direct phase {phase}")
        require_status_zero(native_dir / "exit_status", f"native phase {phase}")
        if recorded_digest(direct_dir / "binary.sha256") != args.expected_direct_binary:
            raise AssertionError(f"direct binary mismatch: {phase}")
        if recorded_digest(direct_dir / "input.sha256") != digest(structure):
            raise AssertionError(f"direct input mismatch: {phase}")
        if recorded_digest(native_dir / "binary.sha256") != args.expected_native_binary:
            raise AssertionError(f"native binary mismatch: {phase}")
        native_input_hash = verify_native_input(native_input, args.mesh_size)
        cell_residual, coordinate_residual = verify_structure_mapping(
            native_input, structure, args.mesh_size
        )
        if recorded_digest(native_dir / "input.sha256") != native_input_hash:
            raise AssertionError(f"native input mismatch: {phase}")
        direct_cpu = qualify_affinity(direct_dir / "affinity_preexec.txt")
        if direct_cpu not in args.expected_direct_cpus:
            raise AssertionError(f"unexpected direct CPU for {phase}: {direct_cpu}")
        native_cpu = qualify_affinity(native_dir / "affinity_preexec.txt")

        direct_text = (direct_dir / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        if (
            "total energy" not in direct_text
            or "JSON dump of results written" not in direct_text
        ):
            raise AssertionError(f"incomplete direct output: {phase}")
        direct_json = direct_dir / "result.json"
        direct_total = read_cli_energy(direct_json)
        direct_primitive = direct_total / replicas
        native_output = native_dir / "cp2k.out"
        native = read_native_energy(native_output)
        direct_dispersion = read_component(
            direct_text, DIRECT_DISPERSION_RE, f"direct dispersion energy: {phase}"
        ) / replicas
        native_text = native_output.read_text(encoding="utf-8", errors="replace")
        native_dispersion = read_component(
            native_text,
            NATIVE_DISPERSION_RE,
            f"native non-self-consistent dispersion energy: {phase}",
        )
        dispersion_delta = native_dispersion - direct_dispersion
        if abs(dispersion_delta) > args.component_tolerance_ha:
            raise AssertionError(
                f"native/direct dispersion mismatch {phase}: {dispersion_delta:+.6e} Ha"
            )
        delta = native - direct_primitive
        if abs(delta) > args.tolerance_ha:
            raise AssertionError(f"native/direct mismatch {phase}: {delta:+.6e} Ha")
        rows.append(
            {
                "phase": phase,
                "water_count_primitive": water_count(structure, replicas),
                "direct_supercell_energy_Ha": direct_total,
                "direct_primitive_energy_Ha": direct_primitive,
                "native_primitive_energy_Ha": native,
                "native_minus_direct_Ha": delta,
                "direct_primitive_dispersion_energy_Ha": direct_dispersion,
                "native_dispersion_energy_Ha": native_dispersion,
                "native_minus_direct_dispersion_Ha": dispersion_delta,
                "structure_sha256": digest(structure),
                "direct_json_sha256": digest(direct_json),
                "native_output_sha256": digest(native_output),
                "native_input_sha256": native_input_hash,
                "supercell_cell_residual_Angstrom": cell_residual,
                "supercell_coordinate_residual_Angstrom": coordinate_residual,
                "direct_cpu": direct_cpu,
                "native_cpu": native_cpu,
            }
        )

    by_phase = {str(row["phase"]): row for row in rows}
    ih = by_phase["Ih"]
    ih_delta_per_water = float(ih["native_minus_direct_Ha"]) / int(
        ih["water_count_primitive"]
    )
    relative_deltas = {
        phase: HARTREE_TO_KJ_MOL
        * (
            float(row["native_minus_direct_Ha"]) / int(row["water_count_primitive"])
            - ih_delta_per_water
        )
        for phase, row in by_phase.items()
        if phase != "Ih"
    }
    maximum_relative = max(
        (abs(value) for value in relative_deltas.values()), default=0.0
    )
    if maximum_relative > args.relative_tolerance_kj_mol:
        raise AssertionError(
            f"native/direct relative mismatch: {maximum_relative:.6e} kJ mol-1 per water"
        )
    signed = [float(row["native_minus_direct_Ha"]) for row in rows]
    dispersion_signed = [
        float(row["native_minus_direct_dispersion_Ha"]) for row in rows
    ]
    payload = {
        "status": "PASS",
        "mesh": f"{args.mesh_size}x{args.mesh_size}x{args.mesh_size}",
        "replicas": replicas,
        "phase_count": len(rows),
        "rows": rows,
        "relative_native_minus_direct_kJ_mol_per_water": relative_deltas,
        "statistics": {
            "max_abs_native_minus_direct_Ha": max(abs(value) for value in signed),
            "rms_native_minus_direct_Ha": math.sqrt(
                sum(value * value for value in signed) / len(signed)
            ),
            "max_abs_native_minus_direct_dispersion_Ha": max(
                abs(value) for value in dispersion_signed
            ),
            "rms_native_minus_direct_dispersion_Ha": math.sqrt(
                sum(value * value for value in dispersion_signed)
                / len(dispersion_signed)
            ),
            "max_supercell_cell_residual_Angstrom": max(
                float(row["supercell_cell_residual_Angstrom"]) for row in rows
            ),
            "max_supercell_coordinate_residual_Angstrom": max(
                float(row["supercell_coordinate_residual_Angstrom"]) for row in rows
            ),
            "max_abs_relative_native_minus_direct_kJ_mol_per_water": maximum_relative,
        },
        "provenance": {
            "direct_binary_sha256": args.expected_direct_binary,
            "native_binary_sha256": args.expected_native_binary,
            "native_provider_archive_sha256": args.expected_native_provider_archive,
            "native_cp2k_revision": args.expected_native_cp2k_revision,
            "native_cmake_cache_sha256": args.expected_native_cmake_cache,
            "native_build_ninja_sha256": args.expected_native_build_ninja,
            "direct_source_revision": args.expected_source_revision,
            "source_identity_sha256": source_hash,
            "allowed_direct_cpus": sorted(args.expected_direct_cpus),
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
