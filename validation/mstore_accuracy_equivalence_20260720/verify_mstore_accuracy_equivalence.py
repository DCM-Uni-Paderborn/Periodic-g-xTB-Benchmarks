#!/usr/bin/env python3
"""Verify the controlled mstore-inorganic ACCURACY 0.1/0.01 comparison."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
HARTREE_TO_KJMOL = Decimal("2625.4996394798254")
EXPECTED_BINARY_SHA256 = "8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa"
TOLERANCE_HARTREE = 1.0e-10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_hash(path: Path) -> str:
    return path.read_text(encoding="utf-8").split()[0]


def read_case(accuracy_dir: str, phase: str) -> tuple[dict[str, object], Decimal]:
    case = ROOT / accuracy_dir / "runs" / "k333" / phase
    poscar = case / "POSCAR"
    if not poscar.is_file():
        poscar = ROOT / accuracy_dir / "structures" / "k333" / phase / "POSCAR"
    required = ("binary.sha256", "exit_status", "input.sha256", "process.out", "result.json")
    missing = [name for name in required if not (case / name).is_file()]
    if not poscar.is_file():
        missing.append("POSCAR")
    if missing:
        raise RuntimeError(f"{accuracy_dir}/{phase}: missing {', '.join(missing)}")

    if (case / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise RuntimeError(f"{accuracy_dir}/{phase}: nonzero exit status")
    if parse_hash(case / "binary.sha256") != EXPECTED_BINARY_SHA256:
        raise RuntimeError(f"{accuracy_dir}/{phase}: unexpected executable hash")
    if parse_hash(case / "input.sha256") != sha256(poscar):
        raise RuntimeError(f"{accuracy_dir}/{phase}: input hash mismatch")

    payload = json.loads(
        (case / "result.json").read_text(encoding="utf-8"), parse_float=Decimal
    )
    energy_hartree = Decimal(payload["energy"])
    atom_counts = [int(value) for value in poscar.read_text(encoding="utf-8").splitlines()[6].split()]
    water_molecules = sum(atom_counts) // 3
    return {
        "energy_hartree_supercell": float(energy_hartree),
        "water_molecules_supercell": water_molecules,
        "input_sha256": sha256(poscar),
        "result_sha256": sha256(case / "result.json"),
        "process_output_sha256": sha256(case / "process.out"),
    }, energy_hartree


def main() -> None:
    raw_cases = {}
    maximum = Decimal(0)
    maximum_phase = ""

    for phase in PHASES:
        loose, loose_energy = read_case("acc010", phase)
        tight, tight_energy = read_case("acc001", phase)
        if loose["input_sha256"] != tight["input_sha256"]:
            raise RuntimeError(f"{phase}: the two accuracy cases do not use the same POSCAR")
        delta = tight_energy - loose_energy
        if abs(delta) > maximum:
            maximum = abs(delta)
            maximum_phase = phase
        raw_cases[phase] = {
            "accuracy_0.1": loose,
            "accuracy_0.01": tight,
            "tight_minus_loose_hartree_supercell": float(delta),
            "loose_energy": loose_energy,
            "tight_energy": tight_energy,
        }

    ih = raw_cases["Ih"]
    ih_delta_per_water = (
        ih["tight_energy"] - ih["loose_energy"]
    ) / Decimal(ih["accuracy_0.1"]["water_molecules_supercell"])
    phases: dict[str, object] = {}
    maximum_relative = Decimal(0)
    maximum_relative_phase = ""
    for phase in PHASES:
        case = raw_cases[phase]
        relative_delta = (
            (case["tight_energy"] - case["loose_energy"])
            / Decimal(case["accuracy_0.1"]["water_molecules_supercell"])
            - ih_delta_per_water
        ) * HARTREE_TO_KJMOL
        if abs(relative_delta) > maximum_relative:
            maximum_relative = abs(relative_delta)
            maximum_relative_phase = phase
        phases[phase] = {
            "accuracy_0.1": case["accuracy_0.1"],
            "accuracy_0.01": case["accuracy_0.01"],
            "tight_minus_loose_hartree_supercell": case[
                "tight_minus_loose_hartree_supercell"
            ],
            "tight_minus_loose_relative_kj_mol_per_H2O": float(relative_delta),
        }

    status = "PASS" if maximum <= TOLERANCE_HARTREE else "FAIL"
    report = {
        "status": status,
        "mesh": [3, 3, 3],
        "phase_count_including_Ih": len(PHASES),
        "executable_sha256": EXPECTED_BINARY_SHA256,
        "maximum_accuracy_sensitivity_hartree_supercell": float(maximum),
        "maximum_accuracy_sensitivity_phase": maximum_phase,
        "maximum_relative_accuracy_sensitivity_kj_mol_per_H2O": float(
            maximum_relative
        ),
        "maximum_relative_accuracy_sensitivity_phase": maximum_relative_phase,
        "tolerance_hartree_supercell": TOLERANCE_HARTREE,
        "interpretation": (
            "The independently rebuilt historical mstore-inorganic executable is energetically "
            "unchanged on the DMC energy scale when its direct CLI accuracy is tightened "
            "from 0.1 to 0.01. The result.json energy is interpreted directly in hartree; "
            "the maximum same-mesh Ih-referenced response per water is also reported. The "
            "mstore-inorganic/pbc DMC difference is therefore a model-source difference, "
            "not an SCC stopping-threshold artifact."
        ),
        "phases": phases,
    }
    (ROOT / "verification.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
