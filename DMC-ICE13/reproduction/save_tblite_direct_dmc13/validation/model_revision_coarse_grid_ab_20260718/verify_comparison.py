#!/usr/bin/env python3
"""Recompute and validate the g-XTB model-revision DMC-ICE13 diagnostic."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parents[1]
PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
    "XIII", "XIV", "XV", "XVII",
)
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
MESH_PROVIDERS = {
    1: ("current", "legacy_mstore_inorganic", "gxtb_v201", "dcm_main"),
    2: ("current", "legacy_mstore_inorganic", "gxtb_v201", "dcm_main"),
    3: ("current", "final_pbc", "legacy_mstore_inorganic"),
}
EXPECTED_MAE = {
    (1, "current"): 163.834465930395,
    (1, "legacy_mstore_inorganic"): 130.899820086448,
    (1, "gxtb_v201"): 41.390727637315,
    (1, "dcm_main"): 162.626025116661,
    (2, "current"): 88.681375524804,
    (2, "legacy_mstore_inorganic"): 48.710763687881,
    (2, "gxtb_v201"): 37.201256501330,
    (2, "dcm_main"): 89.396043594685,
    (3, "current"): 34.04849186179439,
    (3, "final_pbc"): 34.070588090850464,
    (3, "legacy_mstore_inorganic"): 17.83062347616431,
}
EXPECTED_EXECUTABLE_HASH = {
    "legacy_mstore_inorganic_macos": "324c2c1e4968eab579fae1bd8571a467d62a8eaf372f2b88906bb0d9f7ba7549",
    "legacy_mstore_inorganic_linux_k333": "4fa6fd99e1b0de2d0aa76b80cc9089a0ceeefdaf1bc787042221c7fb63479ffd",
    "gxtb_v201": "c87471101170b506dae7f54700d5724aad9ce3dc5923e48d5317a4fd8f6cac60",
    "dcm_main": "2af03fdc70875df823038e49319f69751ae4a94dada58ce2960d09d358884bf0",
}


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(
            f"{label}: actual={actual:.15g}, expected={expected:.15g}, "
            f"tolerance={tolerance:.3g}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def result_path(provider: str, mesh: int, phase: str) -> Path:
    mesh_id = f"k{mesh}{mesh}{mesh}"
    if provider == "current":
        return PACKAGE / "results/current_save_tblite_cli" / mesh_id / phase / "result.json"
    if provider == "legacy_mstore_inorganic":
        if mesh == 3:
            return ROOT / "raw/authors_exchange_linux_k333" / phase / "result.json"
        return ROOT / "raw/authors_exchange" / mesh_id / phase / "result.json"
    if provider == "final_pbc":
        return (
            ROOT.parent / "provider_revision_bvk_ab_20260718"
            / "seidler_pbc_cli_linux" / mesh_id / phase / "result.json"
        )
    return ROOT / "raw" / provider / mesh_id / phase / "result.json"


def energy(provider: str, mesh: int, phase: str) -> float:
    path = result_path(provider, mesh, phase)
    if not path.is_file():
        raise AssertionError(f"missing result: {path}")
    value = float(json.loads(path.read_text())["energy"]) / mesh**3
    if not math.isfinite(value):
        raise AssertionError(f"non-finite energy: {path}")
    return value


def relative(provider: str, mesh: int, phase: str) -> float:
    return (
        energy(provider, mesh, phase) / N_WATER[phase]
        - energy(provider, mesh, "Ih") / N_WATER["Ih"]
    ) * HARTREE_TO_KJMOL


computed: dict[tuple[int, str], float] = {}
for mesh, providers in MESH_PROVIDERS.items():
    for provider in providers:
        errors = [
            abs(relative(provider, mesh, phase) - REFERENCE[phase])
            for phase in PHASES[1:]
        ]
        computed[(mesh, provider)] = sum(errors) / len(errors)
        close(
            computed[(mesh, provider)],
            EXPECTED_MAE[(mesh, provider)],
            5.0e-10,
            f"MAE mesh={mesh} provider={provider}",
        )

with (ROOT / "coarse_mae_summary.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
if len(rows) != len(EXPECTED_MAE):
    raise AssertionError("MAE summary has the wrong row count")
for row in rows:
    key = (int(row["mesh"]), row["provider"])
    if key not in EXPECTED_MAE:
        raise AssertionError(f"unexpected MAE row: {key}")
    close(float(row["mae_kj_mol_per_water"]), computed[key], 5.0e-10, f"summary {key}")

for mesh in (2, 3):
    path = ROOT / f"relative_energies_k{mesh}{mesh}{mesh}.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if [row["phase"] for row in rows] != list(PHASES[1:]):
        raise AssertionError(f"relative-energy phase order is incomplete: {path}")
    for row in rows:
        phase = row["phase"]
        close(
            float(row["dmc_reference_kj_mol_per_water"]),
            REFERENCE[phase],
            5.0e-13,
            f"DMC reference mesh={mesh} phase={phase}",
        )
        for provider in MESH_PROVIDERS[mesh]:
            close(
                float(row[f"{provider}_kj_mol_per_water"]),
                relative(provider, mesh, phase),
                7.0e-10,
                f"relative energy mesh={mesh} provider={provider} phase={phase}",
            )

manifest_pairs: set[tuple[str, str]] = set()
with (ROOT / "k333_input_manifest.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
for row in rows:
    if int(row["mesh"]) != 3:
        raise AssertionError("k333 manifest contains the wrong mesh")
    relative_path = Path(row["path_from_package_root"])
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise AssertionError(f"unsafe k333 manifest path: {relative_path}")
    path = PACKAGE / relative_path
    if sha256(path) != row["sha256"]:
        raise AssertionError(f"k333 manifest hash mismatch: {relative_path}")
    manifest_pairs.add((row["provider"], row["phase"]))
expected_manifest_pairs = {
    (provider, phase) for provider in MESH_PROVIDERS[3] for phase in PHASES
}
if manifest_pairs != expected_manifest_pairs or len(rows) != len(expected_manifest_pairs):
    raise AssertionError("k333 manifest has incomplete or duplicate provider/phase coverage")

linux_root = ROOT / "raw/authors_exchange_linux_k333"
identity = (linux_root / "build_identity.txt").read_text()
if "source_commit=be87ef681acd880705d83b8b1f7c19b58ca5ea85" not in identity:
    raise AssertionError("wrong legacy mstore-inorganic source commit")
if "source_tree=3acace864498c38e25449ec6a93350b4ab35aa4c" not in identity:
    raise AssertionError("wrong legacy mstore-inorganic source tree")
if EXPECTED_EXECUTABLE_HASH["legacy_mstore_inorganic_linux_k333"] not in identity:
    raise AssertionError("wrong Linux k333 legacy executable hash")
if {path.name for path in linux_root.iterdir() if path.is_dir()} != set(PHASES):
    raise AssertionError("legacy Linux k333 phase directory set is incomplete")
for phase in PHASES:
    run = linux_root / phase
    if (run / "exit_status").read_text().strip() != "0":
        raise AssertionError(f"nonzero legacy Linux k333 exit status: {phase}")
    if "JSON dump of results written" not in (run / "process.out").read_text():
        raise AssertionError(f"missing normal legacy Linux k333 completion marker: {phase}")
    if "expected_cpu=90 allowed=90" not in (run / "affinity_preexec.txt").read_text():
        raise AssertionError(f"wrong legacy Linux k333 affinity: {phase}")
    checksum_lines = (run / "SHA256SUMS").read_text().splitlines()
    if len(checksum_lines) != 2:
        raise AssertionError(f"wrong raw checksum count: {phase}")
    structure_hash = checksum_lines[0].split()[0]
    result_hash = checksum_lines[1].split()[0]
    if structure_hash != sha256(PACKAGE / f"structures/k333/{phase}/POSCAR"):
        raise AssertionError(f"legacy Linux k333 structure hash mismatch: {phase}")
    if result_hash != sha256(run / "result.json"):
        raise AssertionError(f"legacy Linux k333 result hash mismatch: {phase}")

for provider, expected_hash in (
    ("authors_exchange", EXPECTED_EXECUTABLE_HASH["legacy_mstore_inorganic_macos"]),
    ("gxtb_v201", EXPECTED_EXECUTABLE_HASH["gxtb_v201"]),
    ("dcm_main", EXPECTED_EXECUTABLE_HASH["dcm_main"]),
):
    actual = (ROOT / "raw" / provider / "executable.sha256").read_text().split()[0]
    if actual != expected_hash:
        raise AssertionError(f"executable hash {provider}: actual={actual} expected={expected_hash}")

platform_deltas = []
for phase in ("Ih", "VII", "XVII"):
    macos = energy("legacy_mstore_inorganic", 2, phase) * 8
    linux = float(json.loads(
        (ROOT / "raw/authors_exchange_linux_k222" / phase / "result.json").read_text()
    )["energy"])
    platform_deltas.append(linux - macos)
max_platform_delta = max(map(abs, platform_deltas))
if max_platform_delta > 2.0e-11:
    raise AssertionError(f"legacy Linux/macOS mismatch: {max_platform_delta:.3e} Ha")

summary = json.loads((ROOT / "k333_revision_summary.json").read_text())
for key, expected in (
    ("current_mae_kj_mol_per_water", computed[(3, "current")]),
    ("final_pbc_mae_kj_mol_per_water", computed[(3, "final_pbc")]),
    ("legacy_mstore_inorganic_mae_kj_mol_per_water", computed[(3, "legacy_mstore_inorganic")]),
    (
        "final_pbc_minus_current_mae_kj_mol_per_water",
        computed[(3, "final_pbc")] - computed[(3, "current")],
    ),
    (
        "legacy_minus_final_pbc_mae_kj_mol_per_water",
        computed[(3, "legacy_mstore_inorganic")] - computed[(3, "final_pbc")],
    ),
):
    close(float(summary[key]), expected, 5.0e-10, f"k333 summary {key}")
shifts = [
    (
        abs(relative("legacy_mstore_inorganic", 3, phase) - relative("final_pbc", 3, phase)),
        phase,
    )
    for phase in PHASES[1:]
]
maximum_shift, maximum_phase = max(shifts)
close(
    float(summary["maximum_legacy_final_relative_energy_shift_kj_mol_per_water"]),
    maximum_shift,
    5.0e-10,
    "k333 maximum legacy/final relative-energy shift",
)
if summary["maximum_shift_phase"] != maximum_phase:
    raise AssertionError("wrong maximum legacy/final shift phase")

print("mesh provider MAE_kJ_mol_per_water")
for key in sorted(computed):
    print(key[0], key[1], f"{computed[key]:.12f}")
print(f"legacy Linux/macOS max delta (Ha): {max_platform_delta:.12e}")
print("model-revision coarse-grid validation: pass")
