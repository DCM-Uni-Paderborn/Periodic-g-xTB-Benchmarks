#!/usr/bin/env python3
"""Recompute the Wigner--Seitz branch-diagnosis relative energies."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HARTREE_TO_KJMOL = 2625.4996394799
WATER_COUNT = {"Ih": 96, "VII": 96, "XVII": 48}
ROOT = Path(__file__).resolve().parent
VALIDATION = ROOT.parent


def energy(path: Path) -> float:
    return float(json.loads(path.read_text(encoding="utf-8"))["energy"])


def relative(root: Path, phase: str) -> float:
    ih = energy(root / "k222" / "Ih" / "result.json") / WATER_COUNT["Ih"]
    value = energy(root / "k222" / phase / "result.json") / WATER_COUNT[phase]
    return (value - ih) * HARTREE_TO_KJMOL


full_roots = {
    "current_integration": (
        VALIDATION
        / "provider_revision_bvk_ab_20260718"
        / "current_save_tblite_cli"
    ),
    "seidler_pbc": (
        VALIDATION
        / "provider_revision_bvk_ab_20260718"
        / "seidler_pbc_cli_linux"
    ),
    "seidler_mstore-inorganic": (
        VALIDATION
        / "model_revision_coarse_grid_ab_20260718"
        / "raw"
        / "authors_exchange"
    ),
    "controlled_compact_index_only": (
        ROOT / "raw" / "controlled_current_source" / "compact_index_only" / "full"
    ),
    "controlled_threshold_only": (
        ROOT / "raw" / "controlled_current_source" / "threshold_only" / "full"
    ),
}

no_coulomb_roots = {
    "current_integration": ROOT / "raw" / "no_coulomb_current",
    "seidler_mstore-inorganic": ROOT / "raw" / "no_coulomb_mstore-inorganic",
    "controlled_compact_index_only": (
        ROOT
        / "raw"
        / "controlled_current_source"
        / "compact_index_only"
        / "no_coulomb"
    ),
    "controlled_threshold_only": (
        ROOT
        / "raw"
        / "controlled_current_source"
        / "threshold_only"
        / "no_coulomb"
    ),
}


def verify_table(path: Path, roots: dict[str, Path]) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        name = row["implementation"]
        for phase in ("VII", "XVII"):
            column = f"{phase}_k222_kj_mol_per_water"
            observed = relative(roots[name], phase)
            reported = float(row[column])
            residual = abs(observed - reported)
            if residual > 5.0e-7:
                raise AssertionError(
                    f"{path.name}: {name} {phase} residual {residual:.3e}"
                )
            print(
                f"PASS {path.name} {name} {phase}: "
                f"{observed:.12f} kJ mol-1 per water"
            )


verify_table(ROOT / "relative_energy_comparison.csv", full_roots)
verify_table(ROOT / "no_coulomb_relative_energy_comparison.csv", no_coulomb_roots)

legacy_exit = int(
    (ROOT / "raw" / "unit_tests" / "legacy_index_wignerseitz_ctest.exitcode")
    .read_text(encoding="utf-8")
    .strip()
)
fixed_exit = int(
    (ROOT / "raw" / "unit_tests" / "fixed_index_wignerseitz_ctest.exitcode")
    .read_text(encoding="utf-8")
    .strip()
)
if legacy_exit == 0 or fixed_exit != 0:
    raise AssertionError(
        f"unexpected unit-test status: legacy={legacy_exit}, fixed={fixed_exit}"
    )
print(f"PASS unit-test statuses: legacy={legacy_exit}, fixed={fixed_exit}")
