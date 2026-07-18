#!/usr/bin/env python3
"""Verify the DMC-ICE13 g-xTB ACCURACY 0.1/0.01 comparison."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
          "XIII", "XIV", "XV", "XVII")
SELECTED_K2 = ("Ih", "VII", "XI", "XIV")
N_WATER = {
    "Ih": 12, "II": 12, "III": 12, "IV": 16, "VI": 10, "VII": 12,
    "VIII": 8, "IX": 12, "XI": 8, "XIII": 28, "XIV": 12, "XV": 10,
    "XVII": 6,
}
REFERENCE = {
    "II": 0.31, "III": 1.25, "IV": 3.83, "VI": 1.78, "VII": 4.99,
    "VIII": 4.23, "IX": 0.60, "XI": 0.16, "XIII": 2.12, "XIV": 1.70,
    "XV": 1.74, "XVII": 1.75,
}
HARTREE_TO_KJMOL = 2625.499638


def cli_energy(path: Path, replication: int = 1) -> float:
    if not path.is_file():
        raise SystemExit(f"missing direct-CLI result: {path}")
    data = json.loads(path.read_text())
    energy = float(data["energy"]) / replication
    if not math.isfinite(energy):
        raise SystemExit(f"non-finite direct-CLI energy: {path}")
    return energy


def native_energy(path: Path) -> float:
    if not path.is_file():
        raise SystemExit(f"missing CP2K output: {path}")
    text = path.read_text()
    if "PROGRAM ENDED" not in text:
        raise SystemExit(f"CP2K did not terminate normally: {path}")
    values = re.findall(
        r"^ ENERGY\|.*?([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s*$", text, re.M
    )
    if not values:
        raise SystemExit(f"missing CP2K energy: {path}")
    return float(values[-1])


def relative(energies: dict[str, float], phase: str) -> float:
    return (
        energies[phase] / N_WATER[phase]
        - energies["Ih"] / N_WATER["Ih"]
    ) * HARTREE_TO_KJMOL


cli_01_gamma = {
    phase: cli_energy(
        PACKAGE / "results/current_save_tblite_cli/k111" / phase / "result.json"
    )
    for phase in PHASES
}
cli_001_gamma = {
    phase: cli_energy(ROOT / "current_save_tblite_cli_acc001/k111" / phase / "result.json")
    for phase in PHASES
}
native_01_gamma = {
    phase: native_energy(
        PACKAGE / "results/current_cp2k_native/k111" / phase / "cp2k.out"
    )
    for phase in PHASES
}
native_001_gamma = {
    phase: native_energy(
        ROOT / "current_cp2k_native_acc001/k111" / phase / "cp2k.out"
    )
    for phase in PHASES
}

relative_shifts_gamma = [
    relative(cli_001_gamma, phase) - relative(cli_01_gamma, phase)
    for phase in PHASES[1:]
]
mae_01 = sum(
    abs(relative(cli_01_gamma, phase) - REFERENCE[phase]) for phase in PHASES[1:]
) / 12
mae_001 = sum(
    abs(relative(cli_001_gamma, phase) - REFERENCE[phase]) for phase in PHASES[1:]
) / 12

cli_01_k2 = {
    phase: cli_energy(
        PACKAGE / "results/current_save_tblite_cli/k222" / phase / "result.json", 8
    )
    for phase in SELECTED_K2
}
cli_001_k2 = {
    phase: cli_energy(
        ROOT / "current_save_tblite_cli_acc001/k222" / phase / "result.json", 8
    )
    for phase in SELECTED_K2
}
relative_shifts_k2 = [
    relative(cli_001_k2, phase) - relative(cli_01_k2, phase)
    for phase in SELECTED_K2[1:]
]

native_cli_delta = [
    native_001_gamma[phase] - cli_001_gamma[phase] for phase in PHASES
]
native_accuracy_delta = [
    native_001_gamma[phase] - native_01_gamma[phase] for phase in PHASES
]

max_gamma_shift = max(map(abs, relative_shifts_gamma))
max_k2_shift = max(map(abs, relative_shifts_k2))
max_native_cli = max(map(abs, native_cli_delta))
rms_native_cli = math.sqrt(sum(value * value for value in native_cli_delta) / len(PHASES))
max_native_accuracy = max(map(abs, native_accuracy_delta))

if max_gamma_shift > 1.0e-6 or max_k2_shift > 1.0e-6:
    raise SystemExit("ACCURACY sensitivity exceeds the documented tolerance")
if max_native_cli > 2.0e-7:
    raise SystemExit("native/direct-CLI Gamma parity exceeds tolerance")
if max_native_accuracy > 5.0e-13:
    raise SystemExit("native ACCURACY change exceeds printed precision")

print(f"Gamma max relative shift (kJ/mol): {max_gamma_shift:.12e}")
print(f"Gamma MAE at 0.1 (kJ/mol): {mae_01:.12f}")
print(f"Gamma MAE at 0.01 (kJ/mol): {mae_001:.12f}")
print(f"Gamma MAE shift (kJ/mol): {mae_001 - mae_01:+.12e}")
print(f"k222 max relative shift (kJ/mol): {max_k2_shift:.12e}")
print(f"Gamma native/CLI max delta (Ha): {max_native_cli:.12e}")
print(f"Gamma native/CLI RMS delta (Ha): {rms_native_cli:.12e}")
print(f"Gamma native 0.01/0.1 max delta (Ha): {max_native_accuracy:.12e}")
print("accuracy-sensitivity validation: pass")
