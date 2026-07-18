#!/usr/bin/env python3
"""Compare CP2K-native and explicit-supercell save_tblite derivatives."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


HARTREE_PER_BOHR3_TO_BAR = 294_210_156.96522176
ANGSTROM_TO_BOHR = 1.889_726_125_457_828


def numbers(text: str) -> list[float]:
    return [float(value.replace("D", "E")) for value in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[ED][-+]?\d+)?", text)]


def determinant(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def poscar_volume_bohr3(path: Path) -> float:
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    lattice = [[scale * float(value) for value in lines[index].split()[:3]] for index in range(2, 5)]
    return abs(determinant(lattice)) * ANGSTROM_TO_BOHR**3


def parse_cli_gradient(path: Path) -> tuple[float, list[list[float]], list[list[float]]]:
    text = path.read_text()
    energy_match = re.search(r"^energy\s+:real:0:\s*\n\s*([-+0-9.ED]+)", text, flags=re.MULTILINE)
    gradient_match = re.search(
        r"^gradient\s+:real:2:3,(\d+)\s*\n(.*?)^virial\s+:real:2:3,3\s*\n(.*)$",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not energy_match or not gradient_match:
        raise ValueError(f"cannot parse save_tblite gradient file {path}")
    natoms = int(gradient_match.group(1))
    gradient_values = numbers(gradient_match.group(2))
    virial_values = numbers(gradient_match.group(3))
    if len(gradient_values) != 3 * natoms or len(virial_values) != 9:
        raise ValueError(f"unexpected derivative dimensions in {path}")
    gradient = [gradient_values[3 * index : 3 * index + 3] for index in range(natoms)]
    virial = [virial_values[3 * index : 3 * index + 3] for index in range(3)]
    return float(energy_match.group(1).replace("D", "E")), gradient, virial


def parse_cp2k(path: Path) -> tuple[float, list[list[float]], list[list[float]]]:
    text = path.read_text()
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL.*?([-+]\d+\.\d+(?:[ED][-+]?\d+)?)", text)
    force_blocks = re.findall(
        r"FORCES\| Atomic forces \[hartree/bohr\](.*?)(?=\n\s*STRESS\||\n\s*ENERGY\||\Z)",
        text,
        flags=re.DOTALL,
    )
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\].*?\n\s*"
        r"STRESS\|\s+x\s+y\s+z\s*\n\s*"
        r"STRESS\|\s+x\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s*\n\s*"
        r"STRESS\|\s+y\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s*\n\s*"
        r"STRESS\|\s+z\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)",
        text,
        flags=re.DOTALL,
    )
    if not energies or not force_blocks or not stress_blocks:
        raise ValueError(f"cannot parse CP2K energy, forces, and stress from {path}")
    force_rows = re.findall(
        r"FORCES\|\s+\d+\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+[-+0-9.ED]+",
        force_blocks[-1],
    )
    forces = [[float(value.replace("D", "E")) for value in row] for row in force_rows]
    stress_values = [float(value.replace("D", "E")) for value in stress_blocks[-1]]
    stress = [stress_values[3 * index : 3 * index + 3] for index in range(3)]
    return float(energies[-1].replace("D", "E")), forces, stress


def fold_forces(gradient: list[list[float]], natoms: int, replicas: int) -> list[list[float]]:
    if len(gradient) != natoms * replicas:
        raise ValueError("CLI atom count does not equal primitive atoms times replicas")
    return [
        [
            -sum(gradient[atom * replicas + image][axis] for image in range(replicas)) / replicas
            for axis in range(3)
        ]
        for atom in range(natoms)
    ]


def differences(left: list[list[float]], right: list[list[float]]) -> tuple[float, float]:
    values = [left[i][j] - right[i][j] for i in range(len(left)) for j in range(len(left[i]))]
    return max(abs(value) for value in values), math.sqrt(sum(value * value for value in values) / len(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cli-gradient", type=Path, required=True)
    parser.add_argument("--poscar", type=Path, required=True)
    parser.add_argument("--replicas", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    cli_energy, cli_gradient, cli_virial = parse_cli_gradient(args.cli_gradient)
    cp2k_energy, cp2k_forces, cp2k_stress = parse_cp2k(args.cp2k)
    cli_forces = fold_forces(cli_gradient, len(cp2k_forces), args.replicas)
    volume = poscar_volume_bohr3(args.poscar)
    cli_stress = [
        [-cli_virial[i][j] / volume * HARTREE_PER_BOHR3_TO_BAR for j in range(3)]
        for i in range(3)
    ]
    force_max, force_rms = differences(cp2k_forces, cli_forces)
    stress_max, stress_rms = differences(cp2k_stress, cli_stress)
    result = {
        "replicas": args.replicas,
        "cp2k_native_energy_Ha_per_primitive": cp2k_energy,
        "save_tblite_cli_energy_Ha_supercell": cli_energy,
        "save_tblite_cli_energy_Ha_per_primitive": cli_energy / args.replicas,
        "native_minus_cli_energy_Ha_per_primitive": cp2k_energy - cli_energy / args.replicas,
        "force_max_abs_difference_Ha_per_bohr": force_max,
        "force_rms_difference_Ha_per_bohr": force_rms,
        "stress_max_abs_difference_bar": stress_max,
        "stress_rms_difference_bar": stress_rms,
        "supercell_volume_bohr3": volume,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
