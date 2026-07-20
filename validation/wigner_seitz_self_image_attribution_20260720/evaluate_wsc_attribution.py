#!/usr/bin/env python3
"""Verify the causal Wigner--Seitz self-image-index attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from decimal import Decimal, getcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
PHASES = ("Ih", "VII")
HARTREE_TO_KJ_MOL = Decimal("2625.4996394798254")
WATERS_IN_K222_SUPERCELL = Decimal(96)
getcontext().prec = 40

VARIANTS = {
    "pbc_old_wsc_same_build": {
        "directory": "old-wsc",
        "binary": "old-wsc-tblite",
        "binary_sha256": "12113f242f3ee3985dbe12b839c733e99cce7bbdf924ba8c326c9eb6aa7f37c0",
        "version": "0.6.0",
    },
    "pbc_correct_wsc_same_build": {
        "directory": "correct-wsc-same-build",
        "binary": "correct-wsc-same-build-tblite",
        "binary_sha256": "b1fc71ad2302323218324c59440099fd0609ce3760a39558f09ed8d77913f15b",
        "version": "0.6.0",
    },
    "mstore_correct_wsc": {
        "directory": "mstore-corrected-wsc",
        "binary": "mstore-corrected-wsc-tblite",
        "binary_sha256": "59ba4d5a646c6f6f1d1b5d6be008edea1007c49382b8a389247a8debc065f38a",
        "version": "0.5.0",
    },
}

HISTORICAL_MSTORE_BINARY_SHA256 = (
    "8df9fcc990f15600f0b99316602d1d6adfad43f85a2b0203ae14aad44ad4b1aa"
)
ARCHIVED_AUTHOR_PBC_BINARY_SHA256 = (
    "81f1d9690ff040836c2f40cfe0eaf6aa33822681ec029479c5633785537d1aee"
)


def find_package() -> Path:
    candidates = [HERE, *HERE.parents]
    candidates.append(
        Path("/Users/tkuehne/.cache/gxtb-part-i-working")
        / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
    )
    for candidate in candidates:
        if (candidate / "tables/mstore_inorganic_absolute_energies.csv").is_file():
            return candidate
        nested = (
            candidate
            / "DMC-ICE13/reproduction/seidler_dmc13_recalculation"
        )
        if (nested / "tables/mstore_inorganic_absolute_energies.csv").is_file():
            return nested
    raise RuntimeError("cannot locate the Seidler DMC-ICE13 reproduction package")


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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def printed_threshold(text: str, label: str) -> Decimal:
    match = re.search(rf"^{re.escape(label)}\s+([0-9.Ee+-]+)\s+", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {label}")
    return Decimal(match.group(1))


def relative_energy(energies: dict[str, Decimal]) -> Decimal:
    return (
        (energies["VII"] - energies["Ih"])
        / WATERS_IN_K222_SUPERCELL
        * HARTREE_TO_KJ_MOL
    )


def load_variant(
    name: str,
    metadata: dict[str, str],
    package: Path,
    checks: dict[str, bool],
    raw_hashes: dict[str, dict[str, str]],
) -> dict[str, Decimal]:
    staged_binary = HERE / metadata["binary"]
    checks[f"{name}:staged_binary_hash_if_available"] = (
        not staged_binary.exists()
        or (staged_binary.is_file() and sha256(staged_binary) == metadata["binary_sha256"])
    )

    energies: dict[str, Decimal] = {}
    for phase in PHASES:
        label = f"{name}/{phase}"
        run = HERE / "results" / metadata["directory"] / phase
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
            raise RuntimeError(f"incomplete WSC result: {label}")

        result = json.loads(
            (run / "result.json").read_text(encoding="utf-8"), parse_float=Decimal
        )
        output = (run / "process.out").read_text(encoding="utf-8")
        energy = Decimal(result["energy"])
        energies[phase] = energy
        input_path = package / f"raw/mstore_inorganic_cli/k222/{phase}/POSCAR"
        run_checks = {
            "exit_zero": (run / "exit_status").read_text(encoding="utf-8").strip()
            == "0",
            "binary_hash_record": recorded_hash(run / "binary.sha256")
            == metadata["binary_sha256"],
            "input_hash_record": recorded_hash(run / "input.sha256")
            == sha256(input_path),
            "energy_threshold": printed_threshold(output, "energy convergence")
            == Decimal("1e-7"),
            "density_threshold": printed_threshold(output, "density convergence")
            == Decimal("2e-6"),
            "finite_energy": math.isfinite(float(energy)),
            "normal_scc": "SCC did not converge" not in output,
            "version": str(result["version"]).strip() == metadata["version"],
        }
        checks.update({f"{label}:{key}": value for key, value in run_checks.items()})
        raw_hashes[label] = {path.name: sha256(path) for path in required}
    return energies


def load_historical_mstore(
    package: Path,
    checks: dict[str, bool],
    raw_hashes: dict[str, dict[str, str]],
) -> dict[str, Decimal]:
    energies: dict[str, Decimal] = {}
    rows = {
        row["phase"]: row
        for row in read_rows(package / "tables/mstore_inorganic_absolute_energies.csv")
        if row["mesh_n"] == "2" and row["phase"] in PHASES
    }
    for phase in PHASES:
        label = f"historical_mstore/{phase}"
        run = package / f"raw/mstore_inorganic_cli/k222/{phase}"
        result_path = run / "result.json"
        output_path = run / "process.out"
        binary_record = run / "binary.sha256"
        input_path = run / "POSCAR"
        result = json.loads(result_path.read_text(encoding="utf-8"), parse_float=Decimal)
        output = output_path.read_text(encoding="utf-8")
        energy = Decimal(result["energy"])
        energies[phase] = energy
        checks.update(
            {
                f"{label}:binary_hash_record": recorded_hash(binary_record)
                == HISTORICAL_MSTORE_BINARY_SHA256,
                f"{label}:input_hash_table": sha256(input_path)
                == rows[phase]["input_sha256"],
                f"{label}:result_hash_table": sha256(result_path)
                == rows[phase]["result_sha256"],
                f"{label}:energy_table": abs(
                    energy - Decimal(rows[phase]["mstore_energy_Ha_supercell"])
                )
                <= Decimal("2e-12"),
                f"{label}:energy_threshold": printed_threshold(
                    output, "energy convergence"
                )
                == Decimal("1e-7"),
                f"{label}:density_threshold": printed_threshold(
                    output, "density convergence"
                )
                == Decimal("2e-6"),
                f"{label}:finite_energy": math.isfinite(float(energy)),
                f"{label}:normal_scc": "SCC did not converge" not in output,
            }
        )
        raw_hashes[label] = {
            "result.json": sha256(result_path),
            "process.out": sha256(output_path),
            "binary.sha256": sha256(binary_record),
            "POSCAR": sha256(input_path),
        }
    return energies


def load_archived_pbc(package: Path, checks: dict[str, bool]) -> dict[str, Decimal]:
    table = package / "tables/author_pbc_absolute_energies.csv"
    rows = {
        row["phase"]: row
        for row in read_rows(table)
        if row["mesh_n"] == "2" and row["phase"] in PHASES
    }
    binary = HERE / "pbc-correct-tblite"
    checks["archived_author_pbc:staged_binary_hash_if_available"] = (
        not binary.exists()
        or (binary.is_file() and sha256(binary) == ARCHIVED_AUTHOR_PBC_BINARY_SHA256)
    )
    checks["archived_author_pbc:table_has_both_phases"] = set(rows) == set(PHASES)
    return {phase: Decimal(rows[phase]["author_pbc_total_Ha"]) for phase in PHASES}


def main() -> None:
    package = find_package()
    checks: dict[str, bool] = {}
    raw_hashes: dict[str, dict[str, str]] = {}

    energies = {
        name: load_variant(name, metadata, package, checks, raw_hashes)
        for name, metadata in VARIANTS.items()
    }
    energies["historical_mstore"] = load_historical_mstore(
        package, checks, raw_hashes
    )
    energies["archived_author_pbc"] = load_archived_pbc(package, checks)

    relative = {name: relative_energy(values) for name, values in energies.items()}
    full_gap = relative["archived_author_pbc"] - relative["historical_mstore"]
    pbc_wsc_shift = (
        relative["pbc_correct_wsc_same_build"]
        - relative["pbc_old_wsc_same_build"]
    )
    mstore_wsc_shift = (
        relative["mstore_correct_wsc"] - relative["historical_mstore"]
    )
    pbc_side_residual = (
        relative["pbc_old_wsc_same_build"] - relative["historical_mstore"]
    )
    corrected_branch_residual = (
        relative["mstore_correct_wsc"] - relative["pbc_correct_wsc_same_build"]
    )
    pbc_side_explained = (
        Decimal(1) - abs(pbc_side_residual) / abs(full_gap)
    ) * Decimal(100)
    inverse_side_explained = (
        Decimal(1) - abs(corrected_branch_residual) / abs(full_gap)
    ) * Decimal(100)
    reciprocal_shift_difference = abs(pbc_wsc_shift - mstore_wsc_shift)

    for phase in PHASES:
        checks[f"correct_same_build_matches_archived_pbc:{phase}"] = abs(
            energies["pbc_correct_wsc_same_build"][phase]
            - energies["archived_author_pbc"][phase]
        ) <= Decimal("2e-12")
    checks["pbc_side_wsc_explains_more_than_95_percent"] = (
        pbc_side_explained > Decimal(95)
    )
    checks["inverse_side_wsc_explains_more_than_95_percent"] = (
        inverse_side_explained > Decimal(95)
    )
    checks["reciprocal_wsc_shifts_agree_within_0p5_kj_mol"] = (
        reciprocal_shift_difference <= Decimal("0.5")
    )
    checks["corrected_branch_residual_below_10_kj_mol"] = (
        abs(corrected_branch_residual) < Decimal(10)
    )
    checks["old_wsc_pbc_is_closer_to_historical_mstore_than_correct_pbc"] = abs(
        pbc_side_residual
    ) < abs(full_gap)

    rows = [
        {
            "source_state": name,
            "Ih_total_Ha": format(energies[name]["Ih"], ".15f"),
            "VII_total_Ha": format(energies[name]["VII"], ".15f"),
            "VII_minus_Ih_kJ_mol_per_H2O": format(relative[name], ".12f"),
        }
        for name in (
            "historical_mstore",
            "pbc_old_wsc_same_build",
            "pbc_correct_wsc_same_build",
            "archived_author_pbc",
            "mstore_correct_wsc",
        )
    ]
    with (HERE / "wsc_relative_energies.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    passed = all(checks.values())
    payload = {
        "schema": "periodic-gxtb-wsc-self-image-attribution-v1",
        "status": "PASS" if passed else "FAIL",
        "system": "DMC-ICE13 Ih and VII at 2x2x2 explicit BvK sampling",
        "water_molecules_per_supercell": int(WATERS_IN_K222_SUPERCELL),
        "relative_energies_kj_mol_per_H2O": {
            name: float(value) for name, value in relative.items()
        },
        "attribution": {
            "correct_pbc_minus_historical_mstore_kj_mol_per_H2O": float(full_gap),
            "pbc_old_to_correct_wsc_shift_kj_mol_per_H2O": float(pbc_wsc_shift),
            "mstore_old_to_correct_wsc_shift_kj_mol_per_H2O": float(
                mstore_wsc_shift
            ),
            "reciprocal_shift_difference_kj_mol_per_H2O": float(
                reciprocal_shift_difference
            ),
            "pbc_old_wsc_minus_historical_mstore_kj_mol_per_H2O": float(
                pbc_side_residual
            ),
            "mstore_correct_wsc_minus_pbc_correct_wsc_kj_mol_per_H2O": float(
                corrected_branch_residual
            ),
            "pbc_side_gap_explained_percent": float(pbc_side_explained),
            "inverse_side_gap_explained_percent": float(inverse_side_explained),
        },
        "source_history": {
            "merge_base": "84274bd621502365ab62fad486162300d5534469",
            "historical_mstore_head": "be87ef681acd880705d83b8b1f7c19b58ca5ea85",
            "author_pbc_head": "c932120d2580811901de6a1fe3f89b943c251766",
            "wsc_fix_commit": "30b04691e0afd1e89d7d74977e679a54fc32f288",
            "wsc_fix_subject": "fix wsc self-image indexing",
        },
        "dependency_history": {
            "pbc_same_build": {
                "dftd4": "99d64ee8383",
                "s_dftd3": "6f0b06f",
                "toml_f": "d5e9270",
                "multicharge": "8abf156",
                "mctc_lib": "e9de066",
                "mstore": "a9070de",
            },
            "historical_mstore_build": {
                "dftd4": "e56749a",
                "s_dftd3": "3425201",
                "toml_f": "d37d83f",
                "multicharge": "8abf156",
                "mctc_lib": "e9de066",
                "mstore": "35e76c9",
                "jonquil": "4d43ffe",
                "test_drive": "9c3401",
            },
        },
        "checks": checks,
        "raw_sha256": raw_hashes,
        "interpretation_limit": (
            "The reciprocal one-patch tests causally attribute more than 95% of the "
            "sparse-mesh branch gap to Wigner--Seitz self-image indexing. The residual "
            "contains all other source and dependency differences and is not assigned "
            "to a single term by this test. Every compared endpoint was reconverged "
            "self-consistently; the shifts are not a fixed-density decomposition."
        ),
        "build_note": (
            "The inverse mstore build used the historical source and pinned dependency "
            "objects with only the WSC orig-index patch. Its final macOS link command "
            "was issued manually because the generated command contained a malformed "
            "Accelerate-framework token; no compiled physics object was replaced."
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
