#!/usr/bin/env python3
"""Assemble the author-facing DMC-ICE13 branch and CP2K comparison tables."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
TABLES = PACKAGE / "tables"
RAW = PACKAGE / "raw"
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
HARTREE_TO_KJMOL = 2625.4996394798254
QUALIFIED_CP2K_SHA256 = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
QUALIFIED_PBC_CLI_SHA256 = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
QUALIFIED_MSTORE_CLI_SHA256 = "8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa"
PARITY_TOLERANCE_HA_PER_PRIMITIVE = 2.0e-7
PARITY_RELATIVE_TOLERANCE_KJMOL_PER_WATER = 5.0e-5
PARITY_CLI_ACCURACY = 0.1
PARITY_CLI_ENERGY_CONVERGENCE_HA = 1.0e-6 * PARITY_CLI_ACCURACY
PARITY_CLI_DENSITY_CONVERGENCE_E = 2.0e-5 * PARITY_CLI_ACCURACY
ENERGY_PATTERN = re.compile(r"ENERGY\| Total FORCE_EVAL .*?([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s*$")
CLI_ENERGY_CONVERGENCE_PATTERN = re.compile(
    r"energy convergence\s+([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+Eh"
)
CLI_DENSITY_CONVERGENCE_PATTERN = re.compile(
    r"density convergence\s+([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+e"
)
CP2K_REQUIRED_ACCURACY = 0.1
CP2K_REQUIRED_EPS_SCF = 1.0e-9
CP2K_REQUIRED_EPS_DEFAULT = 1.0e-12
CP2K_REQUIRED_MAX_SCF = 300
CP2K_REQUIRED_MIXER_ITERATIONS = 300
CP2K_REQUIRED_MIXING_ALPHA = 0.2


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_sidecar_hash(path: Path) -> str:
    if not path.is_file():
        return ""
    fields = path.read_text(encoding="utf-8").split()
    return fields[0].lower() if fields else ""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cp2k_energy(path: Path) -> float | None:
    if not path.is_file():
        return None
    value = None
    last_start = -1
    last_energy = -1
    last_end = -1
    for index, line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines()
    ):
        if "PROGRAM STARTED AT" in line:
            last_start = index
        match = ENERGY_PATTERN.search(line)
        if match:
            value = float(match.group(1))
            last_energy = index
        if "PROGRAM ENDED AT" in line:
            last_end = index
    normally_ended_latest_segment = (
        last_start >= 0
        and last_energy > last_start
        and last_end > last_energy
    )
    return value if normally_ended_latest_segment else None


def cp2k_exit_status(path: Path) -> int | None:
    """Return an archived integer process status, or None if it is unusable."""
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not re.fullmatch(r"[-+]?\d+", text):
        return None
    return int(text)


def cli_convergence(path: Path) -> tuple[float, float]:
    if not path.is_file():
        raise AssertionError(f"missing CLI process output: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    energy_match = CLI_ENERGY_CONVERGENCE_PATTERN.search(text)
    density_match = CLI_DENSITY_CONVERGENCE_PATTERN.search(text)
    if energy_match is None or density_match is None:
        raise AssertionError(f"missing CLI convergence thresholds: {path}")
    return float(energy_match.group(1)), float(density_match.group(1))


def cp2k_input_qualified(path: Path, mesh: int) -> bool:
    """Check every setting that defines the native parity calculation."""
    if not path.is_file():
        return False
    text = path.read_text(encoding="utf-8", errors="replace")

    def scalar(keyword: str) -> float | None:
        match = re.search(
            rf"^\s*{keyword}\s+([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$",
            text,
            re.MULTILINE | re.IGNORECASE,
        )
        return None if match is None else float(match.group(1))

    accuracy = scalar("ACCURACY")
    eps_scf = scalar("EPS_SCF")
    eps_default = scalar("EPS_DEFAULT")
    max_scf = scalar("MAX_SCF")
    mixer_iterations = scalar("ITERATIONS")
    mixing_alpha = scalar("ALPHA")
    expected_shift = 0.0 if mesh % 2 else (mesh - 1) / (2 * mesh)
    scheme = re.search(
        r"^\s*SCHEME\s+MACDONALD\s+(\d+)\s+(\d+)\s+(\d+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if scheme is None:
        return False
    mesh_edges = tuple(int(scheme.group(index)) for index in range(1, 4))
    shifts = tuple(float(scheme.group(index)) for index in range(4, 7))
    required_lines = (
        r"^\s*RUN_TYPE\s+ENERGY\s*$",
        r"^\s*METHOD\s+Quickstep\s*$",
        r"^\s*METHOD\s+xTB\s*$",
        r"^\s*GFN_TYPE\s+TBLITE\s*$",
        r"^\s*METHOD\s+GXTB\s*$",
        r"^\s*SCC_MIXER\s+TBLITE\s*$",
        r"^\s*SCF_GUESS\s+MOPAC\s*$",
        r"^\s*METHOD\s+DIRECT_P_MIXING\s*$",
        r"^\s*SYMMETRY\s+T\s*$",
        r"^\s*FULL_GRID\s+F\s*$",
        r"^\s*SYMMETRY_BACKEND\s+SPGLIB\s*$",
        r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$",
        r"^\s*CANONICALIZE\s+TRUE\s*$",
        r"^\s*PERIODIC\s+XYZ\s*$",
    )
    return (
        accuracy is not None
        and math.isclose(accuracy, CP2K_REQUIRED_ACCURACY, rel_tol=0.0, abs_tol=1.0e-15)
        and eps_scf is not None
        and math.isclose(eps_scf, CP2K_REQUIRED_EPS_SCF, rel_tol=0.0, abs_tol=1.0e-18)
        and eps_default is not None
        and math.isclose(eps_default, CP2K_REQUIRED_EPS_DEFAULT, rel_tol=0.0, abs_tol=1.0e-21)
        and max_scf is not None
        and math.isclose(max_scf, CP2K_REQUIRED_MAX_SCF, rel_tol=0.0, abs_tol=0.0)
        and mixer_iterations is not None
        and math.isclose(
            mixer_iterations,
            CP2K_REQUIRED_MIXER_ITERATIONS,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and mixing_alpha is not None
        and math.isclose(
            mixing_alpha,
            CP2K_REQUIRED_MIXING_ALPHA,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        and mesh_edges == (mesh, mesh, mesh)
        and all(math.isclose(value, expected_shift, rel_tol=0.0, abs_tol=1.0e-15) for value in shifts)
        and all(re.search(pattern, text, re.MULTILINE | re.IGNORECASE) for pattern in required_lines)
    )


def statistics(rows: list[dict[str, object]]) -> dict[str, float]:
    errors = [float(row["error_kj_mol_per_H2O"]) for row in rows]
    return {
        "me_kj_mol_per_H2O": sum(errors) / len(errors),
        "mae_kj_mol_per_H2O": sum(abs(value) for value in errors) / len(errors),
        "rmse_kj_mol_per_H2O": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "maxae_kj_mol_per_H2O": max(abs(value) for value in errors),
    }


def native_mesh_from_directory(name: str) -> int | None:
    """Return N for an exact kNNN-reduced directory name, including N >= 10."""
    if not name.startswith("k") or not name.endswith("-reduced"):
        return None
    digits = name[1:-len("-reduced")]
    if not digits.isdigit() or len(digits) % 3 != 0:
        return None
    width = len(digits) // 3
    components = (digits[:width], digits[width:2 * width], digits[2 * width:])
    if any(
        not component or (len(component) > 1 and component.startswith("0"))
        for component in components
    ):
        return None
    if len(set(components)) != 1:
        return None
    mesh = int(components[0])
    return mesh if mesh > 0 else None


def discover_native_meshes() -> tuple[int, ...]:
    native_root = RAW / "cp2k_native"
    if not native_root.is_dir():
        raise AssertionError(f"missing CP2K-native raw-data directory: {native_root}")
    meshes = {
        mesh
        for child in native_root.iterdir()
        if child.is_dir() and (mesh := native_mesh_from_directory(child.name)) is not None
    }
    if not meshes:
        raise AssertionError(f"no regular CP2K-native mesh directories found in {native_root}")
    return tuple(sorted(meshes))


def main() -> None:
    references = {
        row["phase"]: float(row["reference_relative_energy_kJmol_per_H2O"])
        for row in read_csv(TABLES / "dmc_reference_relative_energies.csv")
    }
    waters = {}
    for phase in PHASES:
        lines = (PACKAGE / "structures/primitive" / phase / "POSCAR").read_text(encoding="utf-8").splitlines()
        waters[phase] = sum(int(value) for value in lines[6].split()) // 3

    absolute_rows: list[dict[str, object]] = []
    relative_rows: list[dict[str, object]] = []
    for mesh in discover_native_meshes():
        mesh_id = f"k{mesh}{mesh}{mesh}-reduced"
        mesh_energies: dict[str, float] = {}
        for phase in PHASES:
            run = RAW / "cp2k_native" / mesh_id / phase
            output = run / "cp2k.out"
            energy = cp2k_energy(output)
            exit_status = cp2k_exit_status(run / "exit_status")
            binary_hash = read_sidecar_hash(run / "binary.sha256")
            input_hash = read_sidecar_hash(run / "input.sha256")
            input_file = run / "input.inp"
            input_hash_matches = (
                input_file.is_file() and input_hash == sha256(input_file)
            )
            input_settings_qualified = cp2k_input_qualified(input_file, mesh)
            qualified = (
                energy is not None
                and exit_status == 0
                and binary_hash == QUALIFIED_CP2K_SHA256
                and input_hash_matches
                and input_settings_qualified
            )
            if not qualified:
                continue
            mesh_energies[phase] = energy
            absolute_rows.append({
                "mesh_n": mesh,
                "mesh_id": mesh_id,
                "phase": phase,
                "water_molecules_primitive": waters[phase],
                "cp2k_native_energy_Ha_per_primitive": f"{energy:.15f}",
                "cp2k_binary_sha256": binary_hash,
                "input_sha256": input_hash,
                "input_file": input_file.relative_to(PACKAGE),
                "input_settings_qualification": "PASS",
                "exit_status": exit_status,
                "normal_termination_qualification": "PASS",
                "output_sha256": sha256(output),
                "raw_output": output.relative_to(PACKAGE),
                "qualification": "PASS",
            })
        if "Ih" not in mesh_energies:
            continue
        ih_per_water = mesh_energies["Ih"] / waters["Ih"]
        for phase in PHASES[1:]:
            if phase not in mesh_energies:
                continue
            relative = (mesh_energies[phase] / waters[phase] - ih_per_water) * HARTREE_TO_KJMOL
            error = relative - references[phase]
            relative_rows.append({
                "method": "CP2K-native pbc provider",
                "mesh_n": mesh,
                "phase": phase,
                "relative_energy_kj_mol_per_H2O": f"{relative:.12f}",
                "dmc_reference_kj_mol_per_H2O": f"{references[phase]:.6f}",
                "error_kj_mol_per_H2O": f"{error:.12f}",
                "absolute_error_kj_mol_per_H2O": f"{abs(error):.12f}",
                "qualification": "PASS",
            })

    write_csv(
        TABLES / "cp2k_native_absolute_energies_by_mesh.csv",
        absolute_rows,
        (
            "mesh_n", "mesh_id", "phase", "water_molecules_primitive",
            "cp2k_native_energy_Ha_per_primitive", "cp2k_binary_sha256",
            "input_sha256", "input_file", "input_settings_qualification",
            "exit_status", "normal_termination_qualification",
            "output_sha256", "raw_output", "qualification",
        ),
    )
    write_csv(
        TABLES / "cp2k_native_relative_energies_by_mesh.csv",
        relative_rows,
        (
            "method", "mesh_n", "phase", "relative_energy_kj_mol_per_H2O",
            "dmc_reference_kj_mol_per_H2O", "error_kj_mol_per_H2O",
            "absolute_error_kj_mol_per_H2O", "qualification",
        ),
    )

    current_rows = read_csv(TABLES / "current_absolute_energies_by_mesh.csv")
    current_fields = tuple(current_rows[0])
    expected_current_keys = {
        (mesh, phase) for mesh in range(1, 5) for phase in PHASES
    }
    current_keys = [(int(row["mesh_n"]), row["phase"]) for row in current_rows]
    if len(current_keys) != len(set(current_keys)):
        raise AssertionError("duplicate mesh/phase rows in current absolute-energy table")
    if set(current_keys) != expected_current_keys:
        missing = sorted(expected_current_keys - set(current_keys))
        extra = sorted(set(current_keys) - expected_current_keys)
        raise AssertionError(
            f"invalid current absolute-energy matrix: missing={missing} extra={extra}"
        )

    qualified_native = {
        (int(row["mesh_n"]), str(row["phase"])): row for row in absolute_rows
    }
    cli_convergence_rows: list[dict[str, object]] = []
    cli_convergence_by_key: dict[tuple[int, str], dict[str, object]] = {}
    for row in current_rows:
        mesh = int(row["mesh_n"])
        phase = row["phase"]
        expected_atoms = waters[phase] * 3
        if row["mesh_id"] != f"k{mesh}{mesh}{mesh}":
            raise AssertionError(f"invalid mesh identifier: mesh={mesh} phase={phase}")
        if int(row["natom_primitive"]) != expected_atoms:
            raise AssertionError(f"invalid primitive atom count: mesh={mesh} phase={phase}")
        if int(row["natom_cli_supercell"]) != expected_atoms * mesh**3:
            raise AssertionError(f"invalid CLI supercell atom count: mesh={mesh} phase={phase}")
        native = qualified_native.get((mesh, phase))
        if native is None:
            raise AssertionError(
                f"missing qualified CP2K-native raw result: mesh={mesh} phase={phase}"
            )
        # The qualified raw output is authoritative.  Rebuild both the energy
        # and its provenance hash instead of carrying either value forward
        # from a previous table-generation run.
        row["cp2k_native_energy_Ha_per_primitive"] = native[
            "cp2k_native_energy_Ha_per_primitive"
        ]
        row["cp2k_output_sha256"] = native["output_sha256"]
        run = RAW / "current_pbc_cli" / f"cli-k{mesh}{mesh}{mesh}" / phase
        result = run / "tblite.json"
        status = run / "exit_status"
        binary_hash = read_sidecar_hash(run / "binary.sha256")
        if (
            result.is_file()
            and status.is_file()
            and status.read_text(encoding="utf-8").strip() == "0"
            and binary_hash == QUALIFIED_PBC_CLI_SHA256
        ):
            energy_convergence, density_convergence = cli_convergence(run / "process.out")
            convergence_matches = (
                math.isclose(
                    energy_convergence,
                    PARITY_CLI_ENERGY_CONVERGENCE_HA,
                    rel_tol=0.0,
                    abs_tol=1.0e-18,
                )
                and math.isclose(
                    density_convergence,
                    PARITY_CLI_DENSITY_CONVERGENCE_E,
                    rel_tol=0.0,
                    abs_tol=1.0e-18,
                )
            )
            inferred_accuracy_energy = energy_convergence / 1.0e-6
            inferred_accuracy_density = density_convergence / 2.0e-5
            convergence_row = {
                "mesh_n": mesh,
                "phase": phase,
                "required_accuracy": f"{PARITY_CLI_ACCURACY:.16g}",
                "inferred_accuracy_from_energy": f"{inferred_accuracy_energy:.16g}",
                "inferred_accuracy_from_density": f"{inferred_accuracy_density:.16g}",
                "energy_convergence_Ha": f"{energy_convergence:.16e}",
                "density_convergence_e": f"{density_convergence:.16e}",
                "process_output_sha256": sha256(run / "process.out"),
                "qualification": "PASS" if convergence_matches else "WRONG_ACCURACY",
            }
            cli_convergence_rows.append(convergence_row)
            cli_convergence_by_key[(mesh, phase)] = convergence_row
            total = float(json.loads(result.read_text(encoding="utf-8"))["energy"])
            if not math.isfinite(total):
                raise AssertionError(
                    f"non-finite CLI energy: mesh={mesh} phase={phase}"
                )
            per_primitive = total / mesh**3
            native = float(row["cp2k_native_energy_Ha_per_primitive"])
            row["save_tblite_cli_energy_Ha_supercell"] = f"{total:.15f}"
            row["save_tblite_cli_energy_Ha_per_primitive"] = f"{per_primitive:.15f}"
            row["native_minus_cli_per_primitive_Ha"] = f"{native - per_primitive:+.12e}"
            row["poscar_sha256"] = read_sidecar_hash(run / "input.sha256")
            row["save_tblite_json_sha256"] = sha256(result)
    write_csv(
        TABLES / "current_cli_convergence_provenance.csv",
        cli_convergence_rows,
        (
            "mesh_n", "phase", "required_accuracy",
            "inferred_accuracy_from_energy", "inferred_accuracy_from_density",
            "energy_convergence_Ha", "density_convergence_e",
            "process_output_sha256", "qualification",
        ),
    )
    write_csv(TABLES / "current_absolute_energies_by_mesh.csv", current_rows, current_fields)
    parity_rows = []
    comparison_rows: list[dict[str, object]] = list(relative_rows)
    for row in current_rows:
        mesh = int(row["mesh_n"])
        phase = row["phase"]
        convergence = cli_convergence_by_key.get((mesh, phase), {})
        parity_rows.append({
            "mesh_n": mesh,
            "phase": phase,
            "water_molecules_primitive": waters[phase],
            "cp2k_native_energy_Ha_per_primitive": row["cp2k_native_energy_Ha_per_primitive"],
            "pbc_cli_energy_Ha_per_primitive": row["save_tblite_cli_energy_Ha_per_primitive"],
            "native_minus_cli_Ha_per_primitive": row["native_minus_cli_per_primitive_Ha"],
            "cp2k_output_sha256": row["cp2k_output_sha256"],
            "pbc_cli_json_sha256": row["save_tblite_json_sha256"],
            "poscar_sha256": row["poscar_sha256"],
            "cli_accuracy": convergence.get("inferred_accuracy_from_energy", ""),
            "cli_accuracy_qualification": convergence.get("qualification", "MISSING"),
        })
    write_csv(
        TABLES / "pbc_cli_vs_cp2k_native_absolute_parity.csv",
        parity_rows,
        (
            "mesh_n", "phase", "water_molecules_primitive",
            "cp2k_native_energy_Ha_per_primitive", "pbc_cli_energy_Ha_per_primitive",
            "native_minus_cli_Ha_per_primitive", "cp2k_output_sha256",
            "pbc_cli_json_sha256", "poscar_sha256", "cli_accuracy",
            "cli_accuracy_qualification",
        ),
    )

    parity_by_mesh = []
    for mesh in sorted({int(row["mesh_n"]) for row in parity_rows}):
        compared = [
            row for row in parity_rows
            if int(row["mesh_n"]) == mesh and row["native_minus_cli_Ha_per_primitive"]
        ]
        absolute_differences = [
            (abs(float(row["native_minus_cli_Ha_per_primitive"])), str(row["phase"]))
            for row in compared
        ]
        maximum_entry = max(absolute_differences, default=None)
        maximum = None if maximum_entry is None else maximum_entry[0]
        maximum_phase = None if maximum_entry is None else maximum_entry[1]
        by_phase = {str(row["phase"]): row for row in compared}
        relative_differences = []
        if "Ih" in by_phase:
            ih = by_phase["Ih"]
            ih_per_water = (
                float(ih["native_minus_cli_Ha_per_primitive"])
                / float(ih["water_molecules_primitive"])
            )
            for phase, row in by_phase.items():
                if phase == "Ih":
                    continue
                phase_per_water = (
                    float(row["native_minus_cli_Ha_per_primitive"])
                    / float(row["water_molecules_primitive"])
                )
                relative_differences.append(
                    (abs((phase_per_water - ih_per_water) * HARTREE_TO_KJMOL), phase)
                )
        maximum_relative_entry = max(relative_differences, default=None)
        maximum_relative = (
            None if maximum_relative_entry is None else maximum_relative_entry[0]
        )
        maximum_relative_phase = (
            None if maximum_relative_entry is None else maximum_relative_entry[1]
        )
        complete = len(compared) == len(PHASES)
        accuracy_complete = complete and all(
            str(row["cli_accuracy_qualification"]) == "PASS" for row in compared
        )
        if not complete:
            status = "INCOMPLETE"
        elif not accuracy_complete:
            status = "WRONG_ACCURACY"
        elif (
            maximum is not None
            and maximum <= PARITY_TOLERANCE_HA_PER_PRIMITIVE
            and maximum_relative is not None
            and maximum_relative <= PARITY_RELATIVE_TOLERANCE_KJMOL_PER_WATER
        ):
            status = "PASS"
        else:
            status = "FAIL"
        parity_by_mesh.append({
            "mesh_n": mesh,
            "compared_phase_count": len(compared),
            "expected_phase_count": len(PHASES),
            "required_cli_accuracy": PARITY_CLI_ACCURACY,
            "accuracy_qualified_phase_count": sum(
                str(row["cli_accuracy_qualification"]) == "PASS" for row in compared
            ),
            "maximum_absolute_difference_Ha_per_primitive": maximum,
            "maximum_absolute_difference_phase": maximum_phase,
            "tolerance_Ha_per_primitive": PARITY_TOLERANCE_HA_PER_PRIMITIVE,
            "maximum_relative_difference_kJ_mol_per_H2O": maximum_relative,
            "maximum_relative_difference_phase": maximum_relative_phase,
            "relative_tolerance_kJ_mol_per_H2O": PARITY_RELATIVE_TOLERANCE_KJMOL_PER_WATER,
            "status": status,
        })
    required_parity_meshes = {1, 2, 3, 4}
    parity_status = {int(row["mesh_n"]): str(row["status"]) for row in parity_by_mesh}
    failed_parity_meshes = sorted(
        mesh for mesh in required_parity_meshes if parity_status.get(mesh) != "PASS"
    )
    for method, energy_column in (
        ("current pbc CLI", "save_tblite_cli_energy_Ha_per_primitive"),
    ):
        by_mesh = {}
        for row in current_rows:
            if row[energy_column]:
                by_mesh.setdefault(int(row["mesh_n"]), {})[row["phase"]] = float(row[energy_column])
        for mesh, values in sorted(by_mesh.items()):
            if "Ih" not in values:
                continue
            ih_per_water = values["Ih"] / waters["Ih"]
            for phase in PHASES[1:]:
                if phase not in values:
                    continue
                relative = (values[phase] / waters[phase] - ih_per_water) * HARTREE_TO_KJMOL
                error = relative - references[phase]
                comparison_rows.append({
                    "method": method,
                    "mesh_n": mesh,
                    "phase": phase,
                    "relative_energy_kj_mol_per_H2O": f"{relative:.12f}",
                    "dmc_reference_kj_mol_per_H2O": f"{references[phase]:.6f}",
                    "error_kj_mol_per_H2O": f"{error:.12f}",
                    "absolute_error_kj_mol_per_H2O": f"{abs(error):.12f}",
                    "qualification": "PASS",
                })

    author_rows = read_csv(TABLES / "author_pbc_relative_energies.csv")
    for method, column in (
        ("author pbc CLI", "author_pbc_kj_mol_per_water"),
    ):
        for row in author_rows:
            phase = row["phase"]
            relative = float(row[column])
            error = relative - references[phase]
            comparison_rows.append({
                "method": method,
                "mesh_n": int(row["mesh_n"]),
                "phase": phase,
                "relative_energy_kj_mol_per_H2O": f"{relative:.12f}",
                "dmc_reference_kj_mol_per_H2O": f"{references[phase]:.6f}",
                "error_kj_mol_per_H2O": f"{error:.12f}",
                "absolute_error_kj_mol_per_H2O": f"{abs(error):.12f}",
                "qualification": "PASS",
            })

    mstore_absolute_rows = []
    mstore_rows = []
    mstore_accuracy_by_key: dict[tuple[int, str], float] = {}
    mstore_root = RAW / "mstore_inorganic_cli"
    for mesh in (1, 2, 3):
        energies = {}
        binary_hashes = set()
        qualified_phases = set()
        for phase in PHASES:
            run = mstore_root / f"k{mesh}{mesh}{mesh}" / phase
            result = run / "result.json"
            status = run / "exit_status"
            if not result.is_file() or not status.is_file() or status.read_text(encoding="utf-8").strip() != "0":
                continue
            energy = float(json.loads(result.read_text(encoding="utf-8"))["energy"])
            binary_hash = read_sidecar_hash(run / "binary.sha256")
            input_hash = read_sidecar_hash(run / "input.sha256")
            archived_input = run / "POSCAR"
            input_hash_consistent = (
                archived_input.is_file() and input_hash == sha256(archived_input)
            )
            energy_convergence, density_convergence = cli_convergence(run / "process.out")
            accuracy_from_energy = energy_convergence / 1.0e-6
            accuracy_from_density = density_convergence / 2.0e-5
            accuracy_consistent = math.isclose(
                accuracy_from_energy,
                accuracy_from_density,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            energies[phase] = energy
            binary_hashes.add(binary_hash)
            exact_binary = binary_hash == QUALIFIED_MSTORE_CLI_SHA256
            if exact_binary and input_hash_consistent and accuracy_consistent:
                qualified_phases.add(phase)
            mstore_accuracy_by_key[(mesh, phase)] = accuracy_from_energy
            mstore_absolute_rows.append({
                "mesh_n": mesh,
                "phase": phase,
                "water_molecules_primitive": waters[phase],
                "mstore_energy_Ha_supercell": f"{energy:.15f}",
                "mstore_energy_Ha_per_primitive": f"{energy / mesh**3:.15f}",
                "binary_sha256": binary_hash,
                "input_sha256": input_hash,
                "cli_accuracy": f"{accuracy_from_energy:.16g}",
                "energy_convergence_Ha": f"{energy_convergence:.16e}",
                "density_convergence_e": f"{density_convergence:.16e}",
                "result_sha256": sha256(result),
                "process_output_sha256": sha256(run / "process.out") if (run / "process.out").is_file() else "",
                "raw_result": result.relative_to(PACKAGE),
                "qualification": (
                    "PASS"
                    if exact_binary and input_hash_consistent and accuracy_consistent
                    else "WRONG_BINARY_HASH"
                    if not exact_binary
                    else "INPUT_HASH_MISMATCH"
                    if not input_hash_consistent
                    else "INCONSISTENT_ACCURACY"
                ),
            })
        if "Ih" not in energies:
            continue
        ih_per_water = energies["Ih"] / (waters["Ih"] * mesh**3)
        for phase in PHASES[1:]:
            if phase not in energies:
                continue
            relative = (
                energies[phase] / (waters[phase] * mesh**3) - ih_per_water
            ) * HARTREE_TO_KJMOL
            error = relative - references[phase]
            row = {
                "method": "historical mstore-inorganic CLI",
                "mesh_n": mesh,
                "phase": phase,
                "relative_energy_kj_mol_per_H2O": f"{relative:.12f}",
                "dmc_reference_kj_mol_per_H2O": f"{references[phase]:.6f}",
                "error_kj_mol_per_H2O": f"{error:.12f}",
                "absolute_error_kj_mol_per_H2O": f"{abs(error):.12f}",
                "qualification": (
                    "PASS"
                    if binary_hashes == {QUALIFIED_MSTORE_CLI_SHA256}
                    and qualified_phases == set(energies)
                    else "UNQUALIFIED_INPUTS"
                ),
            }
            mstore_rows.append(row)
            comparison_rows.append(row)
    write_csv(
        TABLES / "mstore_inorganic_absolute_energies.csv",
        mstore_absolute_rows,
        (
            "mesh_n", "phase", "water_molecules_primitive",
            "mstore_energy_Ha_supercell", "mstore_energy_Ha_per_primitive",
            "binary_sha256", "input_sha256", "cli_accuracy",
            "energy_convergence_Ha", "density_convergence_e", "result_sha256",
            "process_output_sha256", "raw_result", "qualification",
        ),
    )
    write_csv(
        TABLES / "mstore_inorganic_relative_energies_by_mesh.csv",
        mstore_rows,
        (
            "method", "mesh_n", "phase", "relative_energy_kj_mol_per_H2O",
            "dmc_reference_kj_mol_per_H2O", "error_kj_mol_per_H2O",
            "absolute_error_kj_mol_per_H2O", "qualification",
        ),
    )

    comparison_rows.sort(key=lambda row: (str(row["method"]), int(row["mesh_n"]), PHASES.index(str(row["phase"]))))
    write_csv(
        TABLES / "all_branch_relative_energy_comparison.csv",
        comparison_rows,
        (
            "method", "mesh_n", "phase", "relative_energy_kj_mol_per_H2O",
            "dmc_reference_kj_mol_per_H2O", "error_kj_mol_per_H2O",
            "absolute_error_kj_mol_per_H2O", "qualification",
        ),
    )

    overlapping_rows = []
    comparison_index = {
        (str(row["method"]), int(row["mesh_n"]), str(row["phase"])): row
        for row in comparison_rows
    }
    for mesh in (1, 2, 3):
        for phase in PHASES[1:]:
            pbc = comparison_index.get(("current pbc CLI", mesh, phase))
            mstore = comparison_index.get(("historical mstore-inorganic CLI", mesh, phase))
            if pbc is None or mstore is None:
                continue
            pbc_relative = float(pbc["relative_energy_kj_mol_per_H2O"])
            mstore_relative = float(mstore["relative_energy_kj_mol_per_H2O"])
            pbc_accuracy = float(
                cli_convergence_by_key[(mesh, phase)]["inferred_accuracy_from_energy"]
            )
            mstore_accuracy = mstore_accuracy_by_key[(mesh, phase)]
            overlapping_rows.append({
                "mesh_n": mesh,
                "phase": phase,
                "dmc_reference_kj_mol_per_H2O": f"{references[phase]:.6f}",
                "pbc_relative_kj_mol_per_H2O": f"{pbc_relative:.12f}",
                "mstore_inorganic_relative_kj_mol_per_H2O": f"{mstore_relative:.12f}",
                "mstore_minus_pbc_kj_mol_per_H2O": f"{mstore_relative - pbc_relative:.12f}",
                "pbc_absolute_error_kj_mol_per_H2O": pbc["absolute_error_kj_mol_per_H2O"],
                "mstore_absolute_error_kj_mol_per_H2O": mstore["absolute_error_kj_mol_per_H2O"],
                "pbc_cli_accuracy": f"{pbc_accuracy:.16g}",
                "mstore_cli_accuracy": f"{mstore_accuracy:.16g}",
                "same_cli_accuracy": str(
                    math.isclose(
                        pbc_accuracy,
                        mstore_accuracy,
                        rel_tol=0.0,
                        abs_tol=1.0e-15,
                    )
                ).lower(),
            })
    write_csv(
        TABLES / "mstore_vs_pbc_relative_differences.csv",
        overlapping_rows,
        (
            "mesh_n", "phase", "dmc_reference_kj_mol_per_H2O",
            "pbc_relative_kj_mol_per_H2O",
            "mstore_inorganic_relative_kj_mol_per_H2O",
            "mstore_minus_pbc_kj_mol_per_H2O",
            "pbc_absolute_error_kj_mol_per_H2O",
            "mstore_absolute_error_kj_mol_per_H2O",
            "pbc_cli_accuracy", "mstore_cli_accuracy", "same_cli_accuracy",
        ),
    )

    summary_rows = []
    grouped: dict[tuple[str, int], list[dict[str, object]]] = {}
    for row in comparison_rows:
        grouped.setdefault((str(row["method"]), int(row["mesh_n"])), []).append(row)
    for (method, mesh), rows in sorted(grouped.items()):
        if (
            len(rows) != 12
            or {str(row["phase"]) for row in rows} != set(PHASES[1:])
            or any(str(row["qualification"]) != "PASS" for row in rows)
        ):
            continue
        stats = statistics(rows)
        summary_rows.append({
            "method": method,
            "mesh_n": mesh,
            "phase_count": len(rows),
            **{key: f"{value:.12f}" for key, value in stats.items()},
        })
    write_csv(
        TABLES / "branch_comparison_statistics.csv",
        summary_rows,
        (
            "method", "mesh_n", "phase_count", "me_kj_mol_per_H2O",
            "mae_kj_mol_per_H2O", "rmse_kj_mol_per_H2O", "maxae_kj_mol_per_H2O",
        ),
    )

    closure_path = PACKAGE / "evidence" / "three_route_k333_closure" / "summary.json"
    oracle_path = PACKAGE / "evidence" / "cp2k_gamma_supercell_oracle" / "verification.json"
    mic_path = PACKAGE / "evidence" / "second_order_mic_attribution" / "verification.json"
    closure = json.loads(closure_path.read_text(encoding="utf-8")) if closure_path.is_file() else {}
    oracle = json.loads(oracle_path.read_text(encoding="utf-8")) if oracle_path.is_file() else {}
    mic = json.loads(mic_path.read_text(encoding="utf-8")) if mic_path.is_file() else {}
    source_states = json.loads((PACKAGE / "sources.json").read_text(encoding="utf-8"))
    if mic:
        source_states.update({
            "second_order_mic_evidence": "evidence/second_order_mic_attribution/verification.json",
            "second_order_mic_reverted_commit": mic["source"]["reverted_commit"],
            "second_order_mic_test_tree": mic["source"]["hybrid_tree"],
            "second_order_mic_test_cli_sha256": read_sidecar_hash(
                PACKAGE
                / "evidence"
                / "second_order_mic_attribution"
                / "build"
                / "pbc-without-mic-binary.sha256"
            ),
        })
    package_summary = {
        "schema_version": 1,
        "classification": (
            "The exact current pbc-derived CLI and CP2K-native implementation "
            "are numerically equivalent for all 52 required parity points on "
            "the complete 1x1x1 through 4x4x4 matrices. Reciprocal source-patch "
            "tests attribute the historical mstore-inorganic/pbc sparse-mesh "
            "separation to the Wigner--Seitz exchange self-image correction plus "
            "the later minimum-image second-order Coulomb form, not to the CP2K "
            "interface."
            if not failed_parity_meshes
            else "The current pbc-derived CLI and CP2K-native parity matrix is "
            "not yet complete on every required mesh; mstore-inorganic and pbc "
            "remain distinct model source states."
        ),
        "source_states": source_states,
        "current_pbc_cli_vs_cp2k_native_k333": {
            "maximum_absolute_difference_Ha_per_primitive": closure.get(
                "maximum_absolute_native_minus_current_Ha"
            ),
            "maximum_relative_difference_kJ_mol_per_H2O": closure.get(
                "maximum_absolute_native_minus_current_relative_kJ_mol_per_H2O"
            ),
            "status": closure.get("status"),
        },
        "current_pbc_cli_vs_cp2k_native_absolute_parity_by_mesh": parity_by_mesh,
        "cp2k_native_vs_explicit_gamma_bvk_oracle": {
            "difference_Ha_per_primitive": oracle.get("deltas_hartree_per_primitive", {}).get(
                "native_minus_gamma"
            ),
            "status": oracle.get("status"),
        },
        "historical_branch_causality_k222_ice_VII_minus_Ih": ({
            "mstore_wsc_corrected_kj_mol_per_H2O": mic["relative_energies_kj_mol_per_H2O"]["mstore-wsc-corrected"],
            "pbc_correct_kj_mol_per_H2O": mic["relative_energies_kj_mol_per_H2O"]["pbc-correct"],
            "pbc_correct_minus_mstore_wsc_corrected_kj_mol_per_H2O": mic["pbc_correct_minus_mstore_wsc_corrected_kj_mol_per_H2O"],
            "pbc_without_minimum_image_second_order_kj_mol_per_H2O": mic["relative_energies_kj_mol_per_H2O"]["pbc-without-mic"],
            "pbc_without_minimum_image_second_order_minus_mstore_wsc_corrected_kj_mol_per_H2O": mic["pbc_without_mic_minus_mstore_wsc_corrected_kj_mol_per_H2O"],
            "post_wsc_residual_explained_percent": mic["residual_explained_percent"],
            "independent_ice_XVII_crosscheck": {
                "pbc_correct_minus_mstore_wsc_corrected_kj_mol_per_H2O": mic["phase_resolved_attribution"]["XVII"]["pbc_correct_minus_mstore_wsc_corrected_kj_mol_per_H2O"],
                "pbc_without_minimum_image_second_order_minus_mstore_wsc_corrected_kj_mol_per_H2O": mic["phase_resolved_attribution"]["XVII"]["pbc_without_mic_minus_mstore_wsc_corrected_kj_mol_per_H2O"],
                "post_wsc_residual_explained_percent": mic["phase_resolved_attribution"]["XVII"]["residual_explained_percent"],
                "status": mic["status"],
            },
            "second_order_mic_evidence": "evidence/second_order_mic_attribution/verification.json",
            "status": mic["status"],
            "wigner_seitz_evidence": "evidence/wigner_seitz_self_image_attribution/verification.json",
        } if mic else {}),
        "complete_branch_statistics": summary_rows,
    }
    (PACKAGE / "comparison_summary.json").write_text(
        json.dumps(package_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        f"cp2k_absolute={len(absolute_rows)} cp2k_relative={len(relative_rows)} "
        f"mstore_relative={len(mstore_rows)} comparison={len(comparison_rows)} "
        f"statistics={len(summary_rows)}"
    )
    if failed_parity_meshes:
        raise AssertionError(
            "incomplete or failed current-pbc CLI/CP2K-native parity meshes: "
            + ", ".join(str(mesh) for mesh in failed_parity_meshes)
        )


if __name__ == "__main__":
    main()
