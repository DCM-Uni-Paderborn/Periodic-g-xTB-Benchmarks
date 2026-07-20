#!/usr/bin/env python3
"""Reconstruct the complete ice-XVII k=2 derivative qualification."""

from __future__ import annotations

import hashlib
import json
import math
import re
import tarfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "raw_diagnostics.tar.gz"
SUMMARY = HERE / "summary.json"
CLI_ROOT = HERE / "cli_supercell_validation"
OUTPUT = HERE / "verification.json"

ANGSTROM_TO_BOHR = 1.889_726_125_457_828
HARTREE_PER_BOHR3_TO_BAR = 294_210_156.965_221_76
EXPECTED_CP2K_REVISION = "8520b2e592cd04d35081ab4ad46d92c606071e23"
EXPECTED_SAVE_TBLITE_REVISION = "15915c9435644eb257178ca8f8bf7220c38b1a84"
EXPECTED_ARCHIVE_SHA256 = "f7cd8cbbedd8653150aa3f947b5978acec3954f5718e6c4d0bb9c0ebab6b4b0e"

ARCHIVE_MAP = {
    "ch4_k222_full/cp2k.out": "gxtb_acp_image_regression_20260718/CH4_full/cp2k.out",
    "ch4_k222_reduced/cp2k.out": "gxtb_acp_image_regression_20260718/CH4_reduced/cp2k.out",
    "dmc_full/cp2k.out": "dmc_xvii_fullmodel_acp_image_20260718/cp2k.out",
    "dmc_reduced/cp2k.out": "dmc_xvii_reduced_acp_image_20260718/cp2k.out",
    "h2_complex_k311/cp2k.out": "gxtb_acp_image_regression_20260718/H2_shifted/cp2k.out",
    "stress_fd_5e-5/minus/cp2k.out": (
        "dmc_xvii_stress_fd5e5_acp_image_20260718/minus/cp2k.out"
    ),
    "stress_fd_5e-5/plus/cp2k.out": (
        "dmc_xvii_stress_fd5e5_acp_image_20260718/plus/cp2k.out"
    ),
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise RuntimeError(f"{label}: {actual:.17g} != {expected:.17g}")


def matrix_differences(
    left: list[list[float]], right: list[list[float]]
) -> tuple[float, float]:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise RuntimeError("matrix dimensions differ")
    values = [
        left[row][column] - right[row][column]
        for row in range(len(left))
        for column in range(len(left[row]))
    ]
    return max(abs(value) for value in values), math.sqrt(
        sum(value * value for value in values) / len(values)
    )


def determinant(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def cell_volume_bohr3(path: Path) -> float:
    text = path.read_text(encoding="utf-8")
    cell = []
    for label in ("A", "B", "C"):
        match = re.search(
            rf"^\s*{label}\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            text,
            flags=re.MULTILINE,
        )
        if match is None:
            raise RuntimeError(f"cannot parse cell vector {label} from {path}")
        cell.append([float(value) for value in match.groups()])
    return abs(determinant(cell)) * ANGSTROM_TO_BOHR**3


def poscar_volume_bohr3(path: Path) -> float:
    lines = path.read_text(encoding="utf-8").splitlines()
    scale = float(lines[1])
    cell = [
        [scale * float(value) for value in lines[index].split()[:3]]
        for index in range(2, 5)
    ]
    return abs(determinant(cell)) * ANGSTROM_TO_BOHR**3


def parse_cp2k(data: bytes, label: str) -> dict:
    text = data.decode("utf-8")
    if "PROGRAM ENDED AT" not in text:
        raise RuntimeError(f"{label} is not a completed CP2K output")
    energies = re.findall(
        r"ENERGY\| Total FORCE_EVAL.*?([-+]\d+\.\d+(?:[ED][-+]?\d+)?)",
        text,
    )
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
    forces: list[list[float]] = []
    if force_blocks:
        force_rows = re.findall(
            r"FORCES\|\s+\d+\s+([-+0-9.ED]+)\s+([-+0-9.ED]+)\s+"
            r"([-+0-9.ED]+)\s+[-+0-9.ED]+",
            force_blocks[-1],
        )
        forces = [
            [float(value.replace("D", "E")) for value in row]
            for row in force_rows
        ]
    stress: list[list[float]] = []
    if stress_blocks:
        values = [float(value.replace("D", "E")) for value in stress_blocks[-1]]
        stress = [values[3 * row : 3 * row + 3] for row in range(3)]
    if not energies:
        raise RuntimeError(f"cannot parse an energy from {label}")
    return {
        "energy": float(energies[-1].replace("D", "E")),
        "forces": forces,
        "stress": stress,
        "text": text,
    }


def parse_cli_gradient(path: Path) -> tuple[float, list[list[float]], list[list[float]]]:
    text = path.read_text(encoding="utf-8")
    energy_match = re.search(
        r"^energy\s+:real:0:\s*\n\s*([-+0-9.ED]+)",
        text,
        flags=re.MULTILINE,
    )
    derivative_match = re.search(
        r"^gradient\s+:real:2:3,(\d+)\s*\n(.*?)^virial\s+:real:2:3,3\s*\n(.*)$",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if energy_match is None or derivative_match is None:
        raise RuntimeError("cannot parse the save_tblite gradient archive")
    natoms = int(derivative_match.group(1))

    def numbers(raw: str) -> list[float]:
        return [
            float(value.replace("D", "E"))
            for value in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[ED][-+]?\d+)?", raw)
        ]

    gradient_values = numbers(derivative_match.group(2))
    virial_values = numbers(derivative_match.group(3))
    if len(gradient_values) != 3 * natoms or len(virial_values) != 9:
        raise RuntimeError("unexpected save_tblite derivative dimensions")
    gradient = [gradient_values[3 * atom : 3 * atom + 3] for atom in range(natoms)]
    virial = [virial_values[3 * row : 3 * row + 3] for row in range(3)]
    return float(energy_match.group(1).replace("D", "E")), gradient, virial


def fold_cli_forces(
    gradient: list[list[float]], primitive_atoms: int, replicas: int
) -> list[list[float]]:
    if len(gradient) != primitive_atoms * replicas:
        raise RuntimeError("CLI atom count is inconsistent with the BvK replication")
    return [
        [
            -sum(
                gradient[primitive_atom * replicas + image][axis]
                for image in range(replicas)
            )
            / replicas
            for axis in range(3)
        ]
        for primitive_atom in range(primitive_atoms)
    ]


def parse_debug_residuals(text: str, label: str) -> tuple[float, float]:
    virials = re.findall(r"DEBUG\| Sum of differences\s+([-+0-9.Ee]+)", text)
    forces = re.findall(
        r"DEBUG\| Sum of differences:\s+([-+0-9.Ee]+)", text
    )
    if not virials or not forces:
        raise RuntimeError(f"cannot parse finite-difference residuals from {label}")
    return float(forces[-1]), float(virials[-1])


def main() -> None:
    if sha256_file(ARCHIVE) != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("the raw diagnostic archive SHA-256 does not match")
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if summary["cp2k_revision"] != EXPECTED_CP2K_REVISION:
        raise RuntimeError("unexpected CP2K source revision")
    if summary["save_tblite_revision"] != EXPECTED_SAVE_TBLITE_REVISION:
        raise RuntimeError("unexpected save_tblite source revision")

    cli_manifest_entries = 0
    for line_number, raw in enumerate(
        (HERE / "CLI_SHA256SUMS").read_text(encoding="utf-8").splitlines(), 1
    ):
        if not raw.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})\s+\./(.+)", raw)
        if match is None:
            raise RuntimeError(f"invalid CLI_SHA256SUMS line {line_number}")
        expected, relative = match.groups()
        candidate = (HERE / relative).resolve()
        if HERE.resolve() not in candidate.parents or not candidate.is_file():
            raise RuntimeError(f"CLI manifest artifact is unavailable: {relative}")
        if sha256_file(candidate) != expected:
            raise RuntimeError(f"CLI manifest mismatch: {relative}")
        cli_manifest_entries += 1

    with tarfile.open(ARCHIVE, "r:gz") as archive:
        members = {member.name: member for member in archive.getmembers() if member.isfile()}

        def archived(name: str) -> bytes:
            if name not in members:
                raise RuntimeError(f"archive member is missing: {name}")
            handle = archive.extractfile(members[name])
            if handle is None:
                raise RuntimeError(f"archive member is unreadable: {name}")
            return handle.read()

        manifest_entries = 0
        for line_number, raw in enumerate(
            (HERE / "SHA256SUMS").read_text(encoding="utf-8").splitlines(), 1
        ):
            if not raw.strip():
                continue
            match = re.fullmatch(r"([0-9a-f]{64})\s+\./(.+)", raw)
            if match is None:
                raise RuntimeError(f"invalid SHA256SUMS line {line_number}")
            expected, relative = match.groups()
            local = HERE / relative
            if local.is_file():
                actual = sha256_file(local)
            elif relative in ARCHIVE_MAP:
                actual = sha256_bytes(archived(ARCHIVE_MAP[relative]))
            else:
                raise RuntimeError(f"manifest artifact is unavailable: {relative}")
            if actual != expected:
                raise RuntimeError(f"manifest mismatch: {relative}")
            manifest_entries += 1

        full = parse_cp2k(archived(ARCHIVE_MAP["dmc_full/cp2k.out"]), "full mesh")
        reduced = parse_cp2k(
            archived(ARCHIVE_MAP["dmc_reduced/cp2k.out"]), "SPGLIB mesh"
        )
        if len(full["forces"]) != 18 or len(reduced["forces"]) != 18:
            raise RuntimeError("ice XVII must contain 18 primitive-cell atoms")
        energy_difference = full["energy"] - reduced["energy"]
        force_max, force_rms = matrix_differences(full["forces"], reduced["forces"])
        stress_max, stress_rms = matrix_differences(full["stress"], reduced["stress"])
        dmc = summary["dmc_xvii_k222"]
        close(full["energy"], dmc["energy_Ha_per_primitive"], 5.0e-14, "energy")
        close(
            energy_difference,
            dmc["full_minus_reduced_energy_Ha"],
            1.0e-15,
            "full/reduced energy",
        )
        close(
            force_max,
            dmc["full_minus_reduced_max_force_Ha_per_bohr"],
            1.0e-15,
            "full/reduced force",
        )
        close(
            stress_max,
            dmc["full_minus_reduced_max_stress_bar"],
            1.0e-15,
            "full/reduced stress",
        )

        analytic_force = full["forces"][0][0]
        close(
            analytic_force,
            dmc["atom_1_x_analytic_force_Ha_per_bohr"],
            5.0e-12,
            "analytical atom-1 x force",
        )
        force_plus = json.loads(
            (
                CLI_ROOT
                / "finite_difference_force/h00005_plus/result.json"
            ).read_text(encoding="utf-8")
        )["energy"]
        force_minus = json.loads(
            (
                CLI_ROOT
                / "finite_difference_force/h00005_minus/result.json"
            ).read_text(encoding="utf-8")
        )["energy"]
        displacement_bohr = 0.0005 * ANGSTROM_TO_BOHR
        finite_difference_force = -(
            force_plus - force_minus
        ) / (2.0 * displacement_bohr * 8.0)
        force_absolute_difference = abs(analytic_force - finite_difference_force)
        close(
            finite_difference_force,
            dmc["atom_1_x_finite_difference_force_Ha_per_bohr"],
            1.0e-14,
            "finite-difference atom-1 x force",
        )
        close(
            force_absolute_difference,
            dmc["atom_1_x_abs_difference_Ha_per_bohr"],
            1.0e-14,
            "force absolute difference",
        )

        stress_plus = parse_cp2k(
            archived(ARCHIVE_MAP["stress_fd_5e-5/plus/cp2k.out"]),
            "positive strain",
        )["energy"]
        stress_minus = parse_cp2k(
            archived(ARCHIVE_MAP["stress_fd_5e-5/minus/cp2k.out"]),
            "negative strain",
        )["energy"]
        finite_difference_virial = (stress_plus - stress_minus) / (2.0 * 5.0e-5)
        volume_bohr3 = cell_volume_bohr3(HERE / "dmc_full/input.inp")
        analytic_virial = (
            -full["stress"][0][0]
            * volume_bohr3
            / HARTREE_PER_BOHR3_TO_BAR
        )
        close(
            analytic_virial,
            dmc["xx_analytic_virial_Ha"],
            2.0e-10,
            "analytical xx virial",
        )
        close(
            finite_difference_virial,
            dmc["xx_finite_difference_virial_Ha"],
            1.0e-12,
            "finite-difference xx virial",
        )

        cli_energy, cli_gradient, cli_virial = parse_cli_gradient(
            CLI_ROOT / "cli_supercell/gradient.txt"
        )
        cli_forces = fold_cli_forces(cli_gradient, len(reduced["forces"]), 8)
        cli_volume = poscar_volume_bohr3(CLI_ROOT / "cli_supercell/POSCAR")
        cli_stress = [
            [
                -cli_virial[row][column]
                / cli_volume
                * HARTREE_PER_BOHR3_TO_BAR
                for column in range(3)
            ]
            for row in range(3)
        ]
        cli_energy_difference = reduced["energy"] - cli_energy / 8.0
        cli_force_max, cli_force_rms = matrix_differences(
            reduced["forces"], cli_forces
        )
        cli_stress_max, cli_stress_rms = matrix_differences(
            reduced["stress"], cli_stress
        )
        cli_summary = json.loads(
            (CLI_ROOT / "native_vs_cli_summary.json").read_text(encoding="utf-8")
        )
        for actual, key, tolerance in (
            (cli_energy_difference, "native_minus_cli_energy_Ha_per_primitive", 1.0e-13),
            (cli_force_max, "force_max_abs_difference_Ha_per_bohr", 1.0e-14),
            (cli_force_rms, "force_rms_difference_Ha_per_bohr", 1.0e-14),
            (cli_stress_max, "stress_max_abs_difference_bar", 1.0e-8),
            (cli_stress_rms, "stress_rms_difference_bar", 1.0e-8),
        ):
            close(actual, cli_summary[key], tolerance, f"CLI/native {key}")

        regression = {}
        for key, relative in (
            ("h2_complex_k311", "h2_complex_k311/cp2k.out"),
            ("ch4_full_k222", "ch4_k222_full/cp2k.out"),
            ("ch4_reduced_k222", "ch4_k222_reduced/cp2k.out"),
        ):
            parsed = parse_cp2k(archived(ARCHIVE_MAP[relative]), key)
            force_sum, virial_sum = parse_debug_residuals(parsed["text"], key)
            regression[key] = {
                "force_sum_difference_Ha_per_bohr": force_sum,
                "virial_sum_difference_Ha": virial_sum,
            }
            close(
                force_sum,
                summary["regression_gates"][
                    f"{key}_force_sum_difference_Ha_per_bohr"
                ],
                5.0e-13,
                f"{key} force residual",
            )
            close(
                virial_sum,
                summary["regression_gates"][
                    f"{key}_virial_sum_difference_Ha"
                ],
                5.0e-13,
                f"{key} virial residual",
            )

    result = {
        "schema": "periodic-gxtb-part-i-ice-xvii-full-derivatives-v1",
        "status": "PASS",
        "source_revisions": {
            "cp2k": EXPECTED_CP2K_REVISION,
            "save_tblite": EXPECTED_SAVE_TBLITE_REVISION,
        },
        "integrity": {
            "cli_manifest_entries_verified": cli_manifest_entries,
            "legacy_manifest_entries_verified": manifest_entries,
            "raw_archive_sha256": EXPECTED_ARCHIVE_SHA256,
        },
        "full_minus_reduced": {
            "energy_Ha": energy_difference,
            "maximum_force_Ha_per_bohr": force_max,
            "rms_force_Ha_per_bohr": force_rms,
            "maximum_stress_bar": stress_max,
            "rms_stress_bar": stress_rms,
            "force_component_count": sum(len(row) for row in full["forces"]),
            "stress_component_count": sum(len(row) for row in full["stress"]),
        },
        "finite_difference": {
            "atom_1_x_analytic_force_Ha_per_bohr": analytic_force,
            "atom_1_x_cli_supercell_force_Ha_per_bohr": finite_difference_force,
            "atom_1_x_absolute_difference_Ha_per_bohr": force_absolute_difference,
            "xx_analytic_virial_Ha": analytic_virial,
            "xx_native_finite_difference_virial_Ha": finite_difference_virial,
            "xx_absolute_difference_Ha": abs(
                analytic_virial - finite_difference_virial
            ),
        },
        "native_vs_direct_cli_supercell": {
            "energy_difference_Ha_per_primitive": cli_energy_difference,
            "maximum_force_difference_Ha_per_bohr": cli_force_max,
            "rms_force_difference_Ha_per_bohr": cli_force_rms,
            "maximum_stress_difference_bar": cli_stress_max,
            "rms_stress_difference_bar": cli_stress_rms,
        },
        "regression_gates": regression,
        "interpretation": (
            "The complete 18-atom ice-XVII 2x2x2 full-grid and SPGLIB "
            "derivatives are identical at printed precision. Independent collective "
            "supercell displacement and native homogeneous-strain differences reproduce "
            "the manuscript force and virial values, and the complex/self-inverse "
            "regression probes pass."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
