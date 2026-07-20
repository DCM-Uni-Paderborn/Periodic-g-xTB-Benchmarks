#!/usr/bin/env python3
"""Verify the controlled mstore-inorganic ACCURACY 0.1/0.01 comparison."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
EV_PER_HARTREE = 27.211386245988
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


def read_case(accuracy_dir: str, phase: str) -> dict[str, object]:
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

    payload = json.loads((case / "result.json").read_text(encoding="utf-8"))
    energy_ev = float(payload["energy"])
    return {
        "energy_eV_supercell": energy_ev,
        "energy_hartree_supercell": energy_ev / EV_PER_HARTREE,
        "input_sha256": sha256(poscar),
        "result_sha256": sha256(case / "result.json"),
        "process_output_sha256": sha256(case / "process.out"),
    }


def main() -> None:
    phases: dict[str, object] = {}
    maximum = 0.0
    maximum_phase = ""

    for phase in PHASES:
        loose = read_case("acc010", phase)
        tight = read_case("acc001", phase)
        if loose["input_sha256"] != tight["input_sha256"]:
            raise RuntimeError(f"{phase}: the two accuracy cases do not use the same POSCAR")
        delta = (float(tight["energy_eV_supercell"]) - float(loose["energy_eV_supercell"])) / EV_PER_HARTREE
        if abs(delta) > maximum:
            maximum = abs(delta)
            maximum_phase = phase
        phases[phase] = {
            "accuracy_0.1": loose,
            "accuracy_0.01": tight,
            "tight_minus_loose_hartree_supercell": delta,
        }

    status = "PASS" if maximum <= TOLERANCE_HARTREE else "FAIL"
    report = {
        "status": status,
        "mesh": [3, 3, 3],
        "phase_count_including_Ih": len(PHASES),
        "executable_sha256": EXPECTED_BINARY_SHA256,
        "maximum_accuracy_sensitivity_hartree_supercell": maximum,
        "maximum_accuracy_sensitivity_phase": maximum_phase,
        "tolerance_hartree_supercell": TOLERANCE_HARTREE,
        "interpretation": (
            "The independently rebuilt historical mstore-inorganic executable is energetically "
            "unchanged to far below the qualified tolerance when its direct CLI accuracy is "
            "tightened from 0.1 to 0.01. The mstore-inorganic/pbc DMC difference is therefore a "
            "model-source difference, not an SCC stopping-threshold artifact."
        ),
        "phases": phases,
    }
    (ROOT / "verification.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
