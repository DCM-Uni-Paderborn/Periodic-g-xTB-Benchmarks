#!/usr/bin/env python3
"""Qualify the mstore-inorganic versus pbc component ablation matrix."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from decimal import Decimal, getcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
PACKAGE = next(
    (
        candidate
        for candidate in (HERE, *HERE.parents)
        if (candidate / "structures/structure_manifest.csv").is_file()
    ),
    HERE.parents[1] / "DMC-ICE13/reproduction/seidler_dmc13_recalculation",
)
PROVIDERS = ("mstore", "pbc")
MODES = ("full", "no_exchange", "no_acp", "no_exchange_no_acp")
PHASES = ("Ih", "VII")
EXPECTED_BINARY = {
    "mstore": "8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa",
    "pbc": "81f1d9690ff040836c2f40cfe0eaf6aa33822681ec029479c5633785537d1aee",
}
WATERS = Decimal(12 * 8)
HARTREE_TO_KJ_MOL = Decimal("2625.4996394798254")
getcontext().prec = 40


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recorded_hash(path: Path) -> str:
    value = path.read_text(encoding="utf-8").split()[0]
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise RuntimeError(f"invalid SHA-256 record: {path}")
    return value


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def printed_threshold(text: str, label: str) -> Decimal:
    match = re.search(rf"^{re.escape(label)}\s+([0-9.Ee+-]+)\s+", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {label}")
    return Decimal(match.group(1))


def main() -> None:
    energies: dict[tuple[str, str, str], Decimal] = {}
    raw_hashes: dict[str, dict[str, str]] = {}
    checks: dict[str, bool] = {}
    for provider in PROVIDERS:
        for mode in MODES:
            for phase in PHASES:
                label = f"{provider}/{mode}/{phase}"
                run = HERE / "results" / provider / mode / phase
                required = (
                    run / "result.json",
                    run / "process.out",
                    run / "process.err",
                    run / "exit_status",
                    run / "binary.sha256",
                    run / "input.sha256",
                    run / "command.txt",
                )
                if not all(path.is_file() for path in required):
                    raise RuntimeError(f"incomplete component result: {label}")
                result = json.loads(
                    (run / "result.json").read_text(encoding="utf-8"), parse_float=Decimal
                )
                output = (run / "process.out").read_text(encoding="utf-8")
                energy = Decimal(result["energy"])
                energies[(provider, mode, phase)] = energy
                mode_checks = {
                    "exit_zero": (run / "exit_status").read_text(encoding="utf-8").strip() == "0",
                    "binary_hash": recorded_hash(run / "binary.sha256") == EXPECTED_BINARY[provider],
                    "input_hash": recorded_hash(run / "input.sha256")
                    == sha256(PACKAGE / f"raw/mstore_inorganic_cli/k222/{phase}/POSCAR"),
                    "energy_threshold": printed_threshold(output, "energy convergence")
                    == Decimal("1e-7"),
                    "density_threshold": printed_threshold(output, "density convergence")
                    == Decimal("2e-6"),
                    "finite_energy": math.isfinite(float(energy)),
                    "normal_scc": "SCC did not converge" not in output,
                    "parameter_mode": ("--param" in (run / "command.txt").read_text())
                    == (mode != "full"),
                }
                checks.update({f"{label}:{key}": value for key, value in mode_checks.items()})
                raw_hashes[label] = {path.name: sha256(path) for path in required}

    archived_mstore = {
        row["phase"]: Decimal(row["relative_energy_kj_mol_per_H2O"])
        for row in load_rows(PACKAGE / "tables/mstore_inorganic_relative_energies_by_mesh.csv")
        if row["mesh_n"] == "2"
    }
    archived_pbc = {
        row["phase"]: Decimal(row["author_pbc_kj_mol_per_water"])
        for row in load_rows(PACKAGE / "tables/author_pbc_relative_energies.csv")
        if row["mesh_n"] == "2"
    }
    rows: list[dict[str, str]] = []
    relative: dict[tuple[str, str], Decimal] = {}
    for mode in MODES:
        for provider in PROVIDERS:
            value = (
                energies[(provider, mode, "VII")] - energies[(provider, mode, "Ih")]
            ) / WATERS * HARTREE_TO_KJ_MOL
            relative[(provider, mode)] = value
            rows.append(
                {
                    "provider": provider,
                    "mode": mode,
                    "relative_energy_kj_mol_per_H2O": format(value, ".12f"),
                    "pbc_minus_mstore_kj_mol_per_H2O": "",
                }
            )
        gap = relative[("pbc", mode)] - relative[("mstore", mode)]
        rows[-2]["pbc_minus_mstore_kj_mol_per_H2O"] = format(gap, ".12f")
        rows[-1]["pbc_minus_mstore_kj_mol_per_H2O"] = format(gap, ".12f")

    full_gap = relative[("pbc", "full")] - relative[("mstore", "full")]
    reductions = {
        mode: (Decimal(1) - abs(relative[("pbc", mode)] - relative[("mstore", mode)]) / abs(full_gap))
        * Decimal(100)
        for mode in MODES[1:]
    }
    checks["mstore_full_reproduces_archive"] = abs(
        relative[("mstore", "full")] - archived_mstore["VII"]
    ) < Decimal("1e-9")
    checks["pbc_full_reproduces_archive"] = abs(
        relative[("pbc", "full")] - archived_pbc["VII"]
    ) < Decimal("1e-6")
    passed = all(checks.values())

    with (HERE / "component_relative_energies.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "schema": "periodic-gxtb-mstore-pbc-component-ablation-v1",
        "status": "PASS" if passed else "FAIL",
        "phase": "VII relative to same-mode Ih",
        "mesh": "2x2x2 explicit BvK supercell",
        "relative_energies_kj_mol_per_H2O": {
            f"{provider}/{mode}": float(relative[(provider, mode)])
            for provider in PROVIDERS
            for mode in MODES
        },
        "pbc_minus_mstore_full_kj_mol_per_H2O": float(full_gap),
        "gap_reduction_percent": {mode: float(value) for mode, value in reductions.items()},
        "checks": checks,
        "raw_sha256": raw_hashes,
        "interpretation_limit": (
            "Every mode was reconverged self-consistently. The ablations identify which "
            "coupled model terms control the branch gap, but they are not an additive "
            "fixed-density energy decomposition."
        ),
    }
    (HERE / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
