#!/usr/bin/env python3
"""Verify the phase-VII H0-anisotropy source attribution."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "results" / "historical_h0" / "result.json"
OUTPUT = ROOT / "results" / "historical_h0" / "process.out"
CURRENT_ENERGY = -7352.465349096680
FINAL_PBC_ENERGY = -7352.468520192124
EXPECTED_RESULT_SHA256 = "54cb4bce295fb06f08cdc5191dfe09238116df6411438d4ad7b09b4b8f12e7ca"
EXPECTED_OUTPUT_SHA256 = "2332c76b86dfb5e39c97fe6f7696bb4053463a4963e8811611303de3c90f77b9"
INVARIANCE_RESULT_SHA256 = {
    "current_original": "75810b402b24921b278c3b8b081697972a9adcb1d3de7a760ca0a7c1dfaced6e",
    "current_shifted": "b146f8b729da569fd7b886dcb34dc027e0558d244da282c8e7f15cc2083a457f",
    "historical_original": "54cb4bce295fb06f08cdc5191dfe09238116df6411438d4ad7b09b4b8f12e7ca",
    "historical_shifted": "e6195867bff164cc3d912bb45108a2cfac60564a7c7321bee0600855146ae274",
}


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


assert digest(RESULT) == EXPECTED_RESULT_SHA256
assert digest(OUTPUT) == EXPECTED_OUTPUT_SHA256
assert "SCC did not converge" not in OUTPUT.read_text(errors="replace")
energy = float(json.loads(RESULT.read_text())["energy"])
assert math.isfinite(energy)
assert abs(energy - FINAL_PBC_ENERGY) < 5.0e-11
full_gap = FINAL_PBC_ENERGY - CURRENT_ENERGY
unexplained = FINAL_PBC_ENERGY - energy
invariance_energies = {}
for name, expected_sha256 in INVARIANCE_RESULT_SHA256.items():
    result_path = ROOT / "invariance" / name / "result.json"
    output_path = ROOT / "invariance" / name / "process.out"
    assert digest(result_path) == expected_sha256
    assert "SCC did not converge" not in output_path.read_text(errors="replace")
    invariance_energies[name] = float(json.loads(result_path.read_text())["energy"])
current_invariance = invariance_energies["current_shifted"] - invariance_energies["current_original"]
historical_invariance = (
    invariance_energies["historical_shifted"] - invariance_energies["historical_original"]
)
assert abs(current_invariance) < 1.0e-10
assert abs(historical_invariance) > 1.0e-9
summary = {
    "status": "PASS",
    "phase": "VII",
    "mesh": "2x2x2 explicit BvK supercell",
    "energy_hartree_supercell": {
        "current_provider": CURRENT_ENERGY,
        "historical_h0_diagnostic": energy,
        "final_pbc_provider": FINAL_PBC_ENERGY,
    },
    "historical_h0_minus_current_hartree_supercell": energy - CURRENT_ENERGY,
    "historical_h0_minus_final_pbc_hartree_supercell": energy - FINAL_PBC_ENERGY,
    "fraction_of_provider_gap_accounted_for": 1.0 - abs(unexplained / full_gap),
    "equivalent_image_invariance_hartree_supercell": {
        "periodic_neighbour_h0": current_invariance,
        "historical_central_cell_h0": historical_invariance,
    },
}
(ROOT / "verification.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
