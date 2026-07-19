#!/usr/bin/env python3
"""Independently verify the archived provider-component comparison."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
COMPARISON = ROOT / "validation" / "comparison.json"
PROVIDERS = ("current", "pbc")
MODES = ("full", "no_exchange", "no_acp", "no_exchange_no_acp")
EXPECTED_EXECUTABLES = {
    "current": "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a",
    "pbc": "692b7da28bdef43cf8795fb33a2b60ccca49a116115715ffda2bda85fefe565d",
}


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


payload = json.loads(COMPARISON.read_text())
assert payload["phase"] == "VII"
assert payload["primitive_repetitions"] == 8
rows = {row["mode"]: row for row in payload["rows"]}
assert tuple(rows) == MODES

cpus: set[int] = set()
recomputed: dict[str, float] = {}
for mode in MODES:
    row = rows[mode]
    energies: dict[str, float] = {}
    for provider in PROVIDERS:
        directory = ROOT / "results" / provider / mode
        result_path = directory / "result.json"
        output_path = directory / "process.out"
        affinity_path = directory / "affinity_preexec.txt"
        data = json.loads(result_path.read_text())
        record = row[provider]
        energy = float(data["energy"])
        assert math.isfinite(energy)
        assert int((directory / "exit_status").read_text()) == 0
        assert "SCC did not converge" not in output_path.read_text(errors="replace")
        assert math.isclose(energy, record["energy_hartree_supercell"], rel_tol=0.0, abs_tol=1.0e-12)
        assert digest(result_path) == record["result_sha256"]
        assert digest(output_path) == record["output_sha256"]
        assert record["executable_sha256"] == EXPECTED_EXECUTABLES[provider]
        match = re.search(r"allowed=(\d+)", affinity_path.read_text())
        assert match is not None
        cpu = int(match.group(1))
        assert cpu == record["singleton_cpu"]
        assert cpu not in cpus
        cpus.add(cpu)
        energies[provider] = energy
    delta = energies["pbc"] - energies["current"]
    assert math.isclose(delta, row["pbc_minus_current_hartree_supercell"], rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(delta / 8.0, row["pbc_minus_current_hartree_primitive"], rel_tol=0.0, abs_tol=1.0e-12)
    recomputed[mode] = delta

full_gap = abs(recomputed["full"])
no_exchange_gap = abs(recomputed["no_exchange"])
summary = {
    "status": "PASS",
    "phase": payload["phase"],
    "mesh": payload["mesh"],
    "unique_singleton_cpus": sorted(cpus),
    "pbc_minus_current_hartree_supercell": recomputed,
    "gap_reduction_when_exchange_disabled_percent": 100.0 * (1.0 - no_exchange_gap / full_gap),
    "interpretation_limit": (
        "Each ablation was reconverged self-consistently. Component differences are therefore "
        "non-additive and do not constitute a fixed-density energy decomposition."
    ),
}
(ROOT / "verification.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
