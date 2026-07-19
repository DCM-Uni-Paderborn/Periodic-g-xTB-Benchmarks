#!/usr/bin/env python3
"""Integration test for the archived Part-I implementation audit."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
VERIFY = TOOLS / "verify_part_i_implementation.py"
EXPECTED_GATES = {
    "absolute_energy_parity",
    "accuracy_sensitivity",
    "response_fix",
    "energy_force_stress",
    "kpoint_grid_bvk",
    "model_revision",
    "native_derivative_hardening",
    "provider_revision",
    "qualified_energy_sentinels",
    "wigner_seitz_branch_diagnosis",
    "final_lowk_derivatives",
    "phase_viii_component_ablation",
    "phase_xvii_derivative_component_ablation",
    "provider_component_attribution",
    "pbc_h0_anisotropy_attribution",
    "save_tblite_periodic_source_tests",
    "three_route_k333_closure",
    "seidler_recalculation_package",
    "archive_sha256",
}


class PartIImplementationTests(unittest.TestCase):
    def test_complete_archived_audit_passes(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(VERIFY)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["gate_count"], len(EXPECTED_GATES))
        self.assertEqual(payload["passed_gate_count"], len(EXPECTED_GATES))
        self.assertEqual(
            {row["name"] for row in payload["gates"]}, EXPECTED_GATES
        )
        self.assertEqual(payload["failed_gates"], [])


if __name__ == "__main__":
    unittest.main()
