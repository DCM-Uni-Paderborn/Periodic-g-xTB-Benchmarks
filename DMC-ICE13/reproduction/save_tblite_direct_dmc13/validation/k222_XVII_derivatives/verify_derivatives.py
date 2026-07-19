#!/usr/bin/env python3
"""Independently verify the archived ice-XVII energy/force/stress gates."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
sys.path.insert(0, str(PACKAGE / "tools"))

import compare_derivatives as comparison  # noqa: E402


ANGSTROM_TO_BOHR = comparison.ANGSTROM_TO_BOHR
REPLICAS = 8


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(
            f"{label}: actual={actual:.16g} expected={expected:.16g} "
            f"difference={actual - expected:+.3e} tolerance={tolerance:.3e}"
        )


def require_normal_cli(path: Path) -> None:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "JSON dump of results written" not in text:
        raise AssertionError(f"incomplete standalone calculation: {path}")


def energy(path: Path) -> float:
    value = float(json.loads(path.read_text())["energy"])
    if not math.isfinite(value):
        raise AssertionError(f"non-finite energy: {path}")
    return value


def parse_poscar(path: Path) -> tuple[list[list[float]], list[list[float]]]:
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    lattice = [
        [scale * float(value) for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]
    counts = [int(value) for value in lines[6].split()]
    mode_index = 7
    if lines[mode_index].strip().lower().startswith("s"):
        mode_index += 1
    if not lines[mode_index].strip().lower().startswith(("c", "k")):
        raise AssertionError(f"expected Cartesian POSCAR: {path}")
    start = mode_index + 1
    coordinates = [
        [scale * float(value) for value in lines[start + index].split()[:3]]
        for index in range(sum(counts))
    ]
    return lattice, coordinates


def verify_force_geometry(minus: Path, plus: Path, step: float) -> None:
    lattice_minus, xyz_minus = parse_poscar(minus)
    lattice_plus, xyz_plus = parse_poscar(plus)
    for i in range(3):
        for j in range(3):
            close(lattice_plus[i][j], lattice_minus[i][j], 2.0e-12, "force lattice")
    changed: list[tuple[int, int]] = []
    for atom in range(len(xyz_minus)):
        for axis in range(3):
            delta = xyz_plus[atom][axis] - xyz_minus[atom][axis]
            if abs(delta) > 2.0e-12:
                changed.append((atom, axis))
                close(delta, 2.0 * step, 2.0e-12, "force displacement")
    if changed != [(atom, 0) for atom in range(REPLICAS)]:
        raise AssertionError(f"wrong jointly displaced BvK images: {changed}")


def verify_strain_geometry(minus: Path, plus: Path, step: float) -> None:
    lattice_minus, xyz_minus = parse_poscar(minus)
    lattice_plus, xyz_plus = parse_poscar(plus)
    for i in range(3):
        for j in range(3):
            midpoint = 0.5 * (lattice_plus[i][j] + lattice_minus[i][j])
            delta = lattice_plus[i][j] - lattice_minus[i][j]
            expected = 2.0 * step * midpoint if j == 0 else 0.0
            close(delta, expected, 3.0e-12, "homogeneous lattice strain")
    for atom in range(len(xyz_minus)):
        for axis in range(3):
            midpoint = 0.5 * (xyz_plus[atom][axis] + xyz_minus[atom][axis])
            delta = xyz_plus[atom][axis] - xyz_minus[atom][axis]
            expected = 2.0 * step * midpoint if axis == 0 else 0.0
            close(delta, expected, 3.0e-12, "homogeneous coordinate strain")


def verify_native_cli() -> dict[str, float | int]:
    cp2k_path = ROOT / "native_reduced/cp2k_production.out"
    cli_gradient = ROOT / "cli_supercell/gradient.txt"
    poscar = ROOT / "cli_supercell/POSCAR"
    cp2k_text = cp2k_path.read_text(encoding="utf-8", errors="replace")
    if "PROGRAM ENDED AT" not in cp2k_text:
        raise AssertionError("incomplete CP2K derivative calculation")
    if (ROOT / "native_reduced/exit_status").read_text().strip() != "0":
        raise AssertionError("nonzero CP2K derivative exit status")
    if (ROOT / "cli_supercell/exit_status").read_text().strip() != "0":
        raise AssertionError("nonzero standalone derivative exit status")
    require_normal_cli(ROOT / "cli_supercell/process.out")

    cli_energy, cli_gradient_values, cli_virial = comparison.parse_cli_gradient(cli_gradient)
    cp2k_energy, cp2k_forces, cp2k_stress = comparison.parse_cp2k(cp2k_path)
    cli_forces = comparison.fold_forces(cli_gradient_values, len(cp2k_forces), REPLICAS)
    volume = comparison.poscar_volume_bohr3(poscar)
    cli_stress = [
        [
            -cli_virial[i][j] / volume * comparison.HARTREE_PER_BOHR3_TO_BAR
            for j in range(3)
        ]
        for i in range(3)
    ]
    force_max, force_rms = comparison.differences(cp2k_forces, cli_forces)
    stress_max, stress_rms = comparison.differences(cp2k_stress, cli_stress)
    result: dict[str, float | int] = {
        "replicas": REPLICAS,
        "cp2k_native_energy_Ha_per_primitive": cp2k_energy,
        "save_tblite_cli_energy_Ha_supercell": cli_energy,
        "save_tblite_cli_energy_Ha_per_primitive": cli_energy / REPLICAS,
        "native_minus_cli_energy_Ha_per_primitive": cp2k_energy - cli_energy / REPLICAS,
        "force_max_abs_difference_Ha_per_bohr": force_max,
        "force_rms_difference_Ha_per_bohr": force_rms,
        "stress_max_abs_difference_bar": stress_max,
        "stress_rms_difference_bar": stress_rms,
        "supercell_volume_bohr3": volume,
    }
    archived = json.loads((ROOT / "native_vs_cli_summary.json").read_text())
    if set(result) != set(archived):
        raise AssertionError("native/CLI summary field set changed")
    for key, value in result.items():
        close(float(value), float(archived[key]), 2.0e-12, f"native/CLI summary {key}")
    if abs(float(result["native_minus_cli_energy_Ha_per_primitive"])) > 2.0e-7:
        raise AssertionError("native/CLI energy parity exceeds 2e-7 Ha")
    if float(result["force_max_abs_difference_Ha_per_bohr"]) > 1.0e-6:
        raise AssertionError("native/CLI force parity exceeds 1e-6 Ha/bohr")
    if float(result["stress_max_abs_difference_bar"]) > 0.5:
        raise AssertionError("native/CLI stress parity exceeds 0.5 bar")
    return result


def verify_finite_differences() -> dict[str, object]:
    archived = json.loads((ROOT / "finite_difference_summary.json").read_text())
    force_rows = []
    for step, label in ((0.001, "h0001"), (0.0005, "h00005")):
        minus = ROOT / f"finite_difference_force/{label}_minus"
        plus = ROOT / f"finite_difference_force/{label}_plus"
        verify_force_geometry(minus / "POSCAR", plus / "POSCAR", step)
        require_normal_cli(minus / "process.out")
        require_normal_cli(plus / "process.out")
        value = -(
            energy(plus / "result.json") - energy(minus / "result.json")
        ) / (2.0 * step * ANGSTROM_TO_BOHR * REPLICAS)
        force_rows.append({"step_Angstrom": step, "force_Ha_per_bohr": value})

    stress_rows = []
    for step, label in (
        (0.00005, "h000005"),
        (0.0001, "h00001"),
        (0.0002, "h00002"),
        (0.0005, "h00005"),
    ):
        minus = ROOT / f"finite_difference_stress/{label}_minus"
        plus = ROOT / f"finite_difference_stress/{label}_plus"
        verify_strain_geometry(minus / "POSCAR", plus / "POSCAR", step)
        require_normal_cli(minus / "process.out")
        require_normal_cli(plus / "process.out")
        value = (
            energy(plus / "result.json") - energy(minus / "result.json")
        ) / (2.0 * step)
        stress_rows.append({"strain_step": step, "virial_Ha": value})

    for actual, expected in zip(force_rows, archived["force"]["central_difference"]):
        close(actual["step_Angstrom"], expected["step_Angstrom"], 1.0e-15, "force step")
        close(actual["force_Ha_per_bohr"], expected["force_Ha_per_bohr"], 2.0e-12, "force FD")
    for actual, expected in zip(stress_rows, archived["stress"]["central_difference"]):
        close(actual["strain_step"], expected["strain_step"], 1.0e-15, "strain step")
        close(actual["virial_Ha"], expected["virial_Ha"], 2.0e-12, "stress FD")

    analytic_force = float(archived["force"]["analytic_native_Ha_per_bohr"])
    analytic_virial = float(archived["stress"]["analytic_cli_supercell_virial_Ha"])
    if abs(force_rows[-1]["force_Ha_per_bohr"] - analytic_force) > 1.0e-6:
        raise AssertionError("force finite difference does not validate the analytical force")
    if abs(stress_rows[0]["virial_Ha"] - analytic_virial) > 2.0e-4:
        raise AssertionError("strain finite difference does not validate the analytical virial")
    return {"force": force_rows, "stress": stress_rows}


native_cli = verify_native_cli()
finite_difference = verify_finite_differences()
print(
    json.dumps(
        {
            "energy_delta_Ha": native_cli["native_minus_cli_energy_Ha_per_primitive"],
            "force_max_delta_Ha_per_bohr": native_cli[
                "force_max_abs_difference_Ha_per_bohr"
            ],
            "stress_max_delta_bar": native_cli["stress_max_abs_difference_bar"],
            "finest_force_fd_Ha_per_bohr": finite_difference["force"][-1][
                "force_Ha_per_bohr"
            ],
            "finest_stress_fd_virial_Ha": finite_difference["stress"][0]["virial_Ha"],
            "status": "PASS",
        },
        indent=2,
        sort_keys=True,
    )
)
