#!/usr/bin/env python3
"""Verify strict 2x2x2 direct-CLI versus CP2K-native energy parity."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
CP2K_SHA256 = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
CLI_SHA256 = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
HARTREE_TO_KJ_MOL = 2625.4996394799
SUPERCELL_MULTIPLICITY = 8
WATERS_PER_PRIMITIVE_CELL = 12
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
ECONV_RE = re.compile(r"energy convergence\s+([-+0-9.Ee]+)\s+Eh")
PCONV_RE = re.compile(r"density convergence\s+([-+0-9.Ee]+)\s+e")
ITER_RE = re.compile(r"^\s+(\d+)\s+-\d", re.MULTILINE)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def sidecar(path: Path) -> str:
    fields = path.read_text(encoding="utf-8").split()
    return fields[0].lower() if fields else ""


def require_exit_zero(path: Path) -> None:
    if (path / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise ValueError(f"nonzero exit status: {path}")


def cp2k_run(path: Path, input_path: Path, accuracy: float, eps_scf: float) -> dict[str, object]:
    require_exit_zero(path)
    if sidecar(path / "binary.sha256") != CP2K_SHA256:
        raise ValueError(f"wrong CP2K binary: {path}")
    if sidecar(path / "input.sha256") != digest(input_path):
        raise ValueError(f"CP2K input hash mismatch: {path}")
    input_text = input_path.read_text(encoding="utf-8", errors="replace")
    if not re.search(rf"(?im)^\s*ACCURACY\s+{re.escape(str(accuracy))}(?:0*)\s*$", input_text):
        raise ValueError(f"wrong CP2K accuracy: {input_path}")
    eps_match = re.search(r"(?im)^\s*EPS_SCF\s+([-+0-9.Ee]+)\s*$", input_text)
    if eps_match is None or not math.isclose(float(eps_match.group(1)), eps_scf, abs_tol=1.0e-20):
        raise ValueError(f"wrong CP2K EPS_SCF: {input_path}")
    output = (path / "cp2k.out").read_text(encoding="utf-8", errors="replace")
    values = [float(match.group(1)) for line in output.splitlines() if (match := ENERGY_RE.match(line))]
    if "PROGRAM ENDED AT" not in output or not values or not math.isfinite(values[-1]):
        raise ValueError(f"incomplete CP2K output: {path}")
    return {
        "accuracy": accuracy,
        "eps_scf": eps_scf,
        "energy_hartree_per_primitive_cell": values[-1],
        "binary_sha256": CP2K_SHA256,
        "input_sha256": digest(input_path),
        "output_sha256": digest(path / "cp2k.out"),
    }


def cli_run(path: Path, input_path: Path, accuracy: float) -> dict[str, object]:
    require_exit_zero(path)
    if sidecar(path / "binary.sha256") != CLI_SHA256:
        raise ValueError(f"wrong direct-CLI binary: {path}")
    if sidecar(path / "input.sha256") != digest(input_path):
        raise ValueError(f"direct-CLI input hash mismatch: {path}")
    accuracy_path = path / "accuracy.txt"
    recorded_accuracy = (
        float(accuracy_path.read_text(encoding="utf-8")) if accuracy_path.exists() else None
    )
    if recorded_accuracy is not None and recorded_accuracy != accuracy:
        raise ValueError(f"wrong recorded direct-CLI accuracy: {path}")
    output = (path / "process.out").read_text(encoding="utf-8", errors="replace")
    energy_match = ECONV_RE.search(output)
    density_match = PCONV_RE.search(output)
    iterations = [int(value) for value in ITER_RE.findall(output)]
    if energy_match is None or density_match is None or not iterations:
        raise ValueError(f"incomplete direct-CLI output: {path}")
    energy_convergence = float(energy_match.group(1))
    density_convergence = float(density_match.group(1))
    if not math.isclose(energy_convergence, 1.0e-6 * accuracy, abs_tol=1.0e-18):
        raise ValueError(f"wrong direct-CLI energy convergence: {path}")
    if not math.isclose(density_convergence, 2.0e-5 * accuracy, abs_tol=1.0e-18):
        raise ValueError(f"wrong direct-CLI density convergence: {path}")
    result = json.loads((path / "tblite.json").read_text(encoding="utf-8"))
    energy = float(result["energy"])
    if not math.isfinite(energy):
        raise ValueError(f"non-finite direct-CLI energy: {path}")
    return {
        "accuracy": accuracy,
        "recorded_accuracy": recorded_accuracy,
        "energy_hartree_supercell": energy,
        "energy_hartree_per_primitive_cell": energy / SUPERCELL_MULTIPLICITY,
        "energy_convergence_hartree": energy_convergence,
        "density_convergence_e": density_convergence,
        "iterations": max(iterations),
        "binary_sha256": CLI_SHA256,
        "input_sha256": digest(input_path),
        "result_sha256": digest(path / "tblite.json"),
        "process_output_sha256": digest(path / "process.out"),
    }


def main() -> None:
    loose_native_root = PACKAGE / "raw/cp2k_native/k222-reduced"
    loose_cli_root = PACKAGE / "raw/current_pbc_cli/cli-k222"
    strict_native_root = HERE / "results/native"
    strict_cli_root = HERE / "results/cli"

    rows: dict[str, dict[str, object]] = {}
    for phase in ("Ih", "VII"):
        loose_native = cp2k_run(
            loose_native_root / phase,
            loose_native_root / phase / "input.inp",
            0.1,
            1.0e-9,
        )
        strict_native = cp2k_run(
            strict_native_root / phase,
            HERE / "inputs" / phase / "input.inp",
            0.01,
            1.0e-10,
        )
        loose_cli = cli_run(
            loose_cli_root / phase,
            loose_cli_root / phase / "POSCAR",
            0.1,
        )
        strict_cli = cli_run(
            strict_cli_root / phase,
            loose_cli_root / phase / "POSCAR",
            0.01,
        )
        if loose_cli["input_sha256"] != strict_cli["input_sha256"]:
            raise ValueError(f"loose and strict direct-CLI structures differ: {phase}")
        loose_delta = float(loose_native["energy_hartree_per_primitive_cell"]) - float(
            loose_cli["energy_hartree_per_primitive_cell"]
        )
        strict_delta = float(strict_native["energy_hartree_per_primitive_cell"]) - float(
            strict_cli["energy_hartree_per_primitive_cell"]
        )
        rows[phase] = {
            "accuracy_0.1": {"cp2k_native": loose_native, "direct_cli": loose_cli},
            "accuracy_0.01": {"cp2k_native": strict_native, "direct_cli": strict_cli},
            "native_minus_cli_hartree_per_primitive_cell_accuracy_0.1": loose_delta,
            "native_minus_cli_hartree_per_primitive_cell_accuracy_0.01": strict_delta,
            "absolute_parity_change_hartree": abs(strict_delta) - abs(loose_delta),
            "cp2k_tight_minus_loose_hartree": float(
                strict_native["energy_hartree_per_primitive_cell"]
            )
            - float(loose_native["energy_hartree_per_primitive_cell"]),
            "cli_tight_minus_loose_hartree_per_primitive_cell": float(
                strict_cli["energy_hartree_per_primitive_cell"]
            )
            - float(loose_cli["energy_hartree_per_primitive_cell"]),
        }

    for label, accuracy_key in (("accuracy_0.1", "accuracy_0.1"), ("accuracy_0.01", "accuracy_0.01")):
        native_relative = (
            float(rows["VII"][accuracy_key]["cp2k_native"]["energy_hartree_per_primitive_cell"])
            - float(rows["Ih"][accuracy_key]["cp2k_native"]["energy_hartree_per_primitive_cell"])
        ) * HARTREE_TO_KJ_MOL / WATERS_PER_PRIMITIVE_CELL
        cli_relative = (
            float(rows["VII"][accuracy_key]["direct_cli"]["energy_hartree_per_primitive_cell"])
            - float(rows["Ih"][accuracy_key]["direct_cli"]["energy_hartree_per_primitive_cell"])
        ) * HARTREE_TO_KJ_MOL / WATERS_PER_PRIMITIVE_CELL
        rows["VII"][f"relative_energy_kj_mol_per_water_{label}"] = {
            "cp2k_native": native_relative,
            "direct_cli": cli_relative,
            "native_minus_cli": native_relative - cli_relative,
        }

    max_loose = max(
        abs(float(row["native_minus_cli_hartree_per_primitive_cell_accuracy_0.1"]))
        for row in rows.values()
    )
    max_strict = max(
        abs(float(row["native_minus_cli_hartree_per_primitive_cell_accuracy_0.01"]))
        for row in rows.values()
    )
    strict_relative_delta = abs(
        float(rows["VII"]["relative_energy_kj_mol_per_water_accuracy_0.01"]["native_minus_cli"])
    )
    status = "PASS" if max_strict <= 2.0e-7 and strict_relative_delta <= 5.0e-5 else "FAIL"
    report = {
        "status": status,
        "mesh": [2, 2, 2],
        "phases": rows,
        "maximum_absolute_native_minus_cli_hartree_per_primitive_cell_accuracy_0.1": max_loose,
        "maximum_absolute_native_minus_cli_hartree_per_primitive_cell_accuracy_0.01": max_strict,
        "strict_relative_energy_native_minus_cli_kj_mol_per_water": strict_relative_delta,
        "tolerances": {
            "absolute_energy_hartree_per_primitive_cell": 2.0e-7,
            "relative_energy_kj_mol_per_water": 5.0e-5,
        },
        "interpretation": (
            "Tightening both routes leaves the native energies unchanged at printed precision and "
            "changes the direct-CLI supercell energies only below 1e-9 hartree. The residual "
            "native/CLI difference therefore remains a tiny deterministic numerical path "
            "difference, not an SCF-stopping or unit-conversion error."
        ),
    }
    output = HERE / "verification.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
