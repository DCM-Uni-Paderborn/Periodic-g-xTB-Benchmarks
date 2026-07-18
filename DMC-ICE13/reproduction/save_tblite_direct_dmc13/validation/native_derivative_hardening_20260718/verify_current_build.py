#!/usr/bin/env python3
"""Verify the final-build energy, force, and stress retention gate."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
TOOLS = ROOT.parent.parent / "tools" / "compare_derivatives.py"
SPEC = importlib.util.spec_from_file_location("derivative_comparison", TOOLS)
if SPEC is None or SPEC.loader is None:
    raise SystemExit(f"cannot load {TOOLS}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if abs(actual - expected) > tolerance:
        raise SystemExit(
            f"{label}: actual={actual:.15g} expected={expected:.15g} "
            f"difference={actual - expected:+.3e}"
        )


reference_path = ROOT / "dmc_reduced" / "cp2k.out"
current_path = ROOT / "current_build_gate" / "cp2k.out"
summary = json.loads((ROOT / "current_build_gate" / "comparison.json").read_text())
reference_energy, reference_forces, reference_stress = MODULE.parse_cp2k(reference_path)
current_energy, current_forces, current_stress = MODULE.parse_cp2k(current_path)
force_max, force_rms = MODULE.differences(current_forces, reference_forces)
stress_max, stress_rms = MODULE.differences(current_stress, reference_stress)

close(current_energy, summary["energy_Ha_per_primitive"], 5.0e-13, "energy")
close(current_energy - reference_energy, summary["energy_difference_Ha"],
      5.0e-14, "energy difference")
close(force_max, summary["maximum_force_difference_Ha_per_bohr"],
      5.0e-12, "maximum force difference")
close(force_rms, summary["rms_force_difference_Ha_per_bohr"],
      5.0e-12, "RMS force difference")
close(stress_max, summary["maximum_stress_difference_bar"],
      5.0e-8, "maximum stress difference")
close(stress_rms, summary["rms_stress_difference_bar"],
      5.0e-8, "RMS stress difference")
if ("PROGRAM ENDED AT" in current_path.read_text()) != summary["program_ended"]:
    raise SystemExit("program completion state does not match the summary")

print("current CP2K build retention gate: all archived values verified")
