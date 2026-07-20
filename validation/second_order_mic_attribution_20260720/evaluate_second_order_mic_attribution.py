#!/usr/bin/env python3
"""Verify the reciprocal one-patch attribution of the post-WSC branch residual."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from decimal import Decimal, getcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
HARTREE_TO_KJ_MOL = Decimal("2625.4996394798254")
WATERS = Decimal(96)
PHASES = ("Ih", "VII")
PROVIDERS = ("pbc-correct", "mstore-wsc-corrected", "pbc-without-mic")
EXPECTED_BINARY = {
    "pbc-correct": "81f1d9690ff040836c2f40cfe0eaf6aa33822681ec029479c5633785537d1aee",
    "mstore-wsc-corrected": "59ba4d5a646c6f6f1d1b5d6be008edea1007c49382b8a389247a8debc065f38a",
    "pbc-without-mic": "1355ba96a0b5ae491a4dbe5015eda072b355f70e675ccb95c3411b16498a6827",
}
EXPECTED_INPUT = {
    "Ih": "cc6ec119078ecd4769f50893aa7cd6128390cab0eb6f069d58d9972b0a1904b7",
    "VII": "4de281e3ab3632f443b22f99162e80a3a327c0b8124489c015d4b68fa62d1d91",
}
EXPECTED_SOURCE = {
    "pbc_head": "c932120d2580811901de6a1fe3f89b943c251766",
    "reverted_commit": "083f22030f0be7abff3e3d27b35c141c04c2aa6d",
    "hybrid_tree": "d56a0195324b25c99a8971987d13a10d9f9aa678",
    "patch_sha256": "548858893ab921230852d2c005ec7196ea2ad326d9877b92d4d10785da395ac6",
}
EXPECTED_CHANGED_FILES = (
    "M\tsrc/tblite/coulomb/charge/effective.f90\n"
    "M\tsrc/tblite/coulomb/charge/type.f90\n"
    "M\tsrc/tblite/coulomb/thirdorder/twobody.f90\n"
    "M\ttest/cli/06a-gfn1xtb.json\n"
    "M\ttest/unit/test_coulomb_charge.f90\n"
)
CACHE_KEYS = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_C_COMPILER",
    "CMAKE_Fortran_COMPILER",
    "BUILD_SHARED_LIBS",
    "WITH_API",
    "WITH_DDX",
    "WITH_HDF5",
    "WITH_OpenMP",
    "WITH_TESTS",
    "WITH_TREXIO",
    "mctc-lib_DIR",
    "dftd4_DIR",
    "multicharge_DIR",
    "s-dftd3_DIR",
    "toml-f_DIR",
)
getcontext().prec = 50


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


def poscar_composition(path: Path) -> dict[str, int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    symbols = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    return dict(zip(symbols, counts, strict=True))


def cmake_cache(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line or ":" not in line:
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        values[key] = value
    return {key: values.get(key, "<missing>") for key in CACHE_KEYS}


def main() -> None:
    checks: dict[str, bool] = {}
    raw_sha256: dict[str, dict[str, str]] = {}
    energies: dict[tuple[str, str], Decimal] = {}

    for phase in PHASES:
        structure = HERE / "inputs" / phase / "POSCAR"
        checks[f"input/{phase}:hash"] = sha256(structure) == EXPECTED_INPUT[phase]
        checks[f"input/{phase}:composition"] = poscar_composition(structure) == {
            "H": 192,
            "O": 96,
        }

    for provider in PROVIDERS:
        for phase in PHASES:
            label = f"{provider}/{phase}"
            run = HERE / "results" / provider / phase
            required = tuple(
                run / name
                for name in (
                    "result.json",
                    "process.out",
                    "process.err",
                    "command.txt",
                    "exit_status",
                    "binary.sha256",
                    "input.sha256",
                )
            )
            if not all(path.is_file() for path in required):
                raise RuntimeError(f"incomplete raw result: {label}")
            result = json.loads(
                (run / "result.json").read_text(encoding="utf-8"),
                parse_float=Decimal,
            )
            output = (run / "process.out").read_text(encoding="utf-8")
            command = (run / "command.txt").read_text(encoding="utf-8")
            energy_match = re.search(
                r"^energy convergence\s+([0-9.Ee+-]+)\s+", output, re.MULTILINE
            )
            density_match = re.search(
                r"^density convergence\s+([0-9.Ee+-]+)\s+", output, re.MULTILINE
            )
            energy = Decimal(result["energy"])
            energies[(provider, phase)] = energy
            checks.update(
                {
                    f"{label}:exit_zero": (run / "exit_status").read_text().strip()
                    == "0",
                    f"{label}:binary_hash": recorded_hash(run / "binary.sha256")
                    == EXPECTED_BINARY[provider],
                    f"{label}:input_hash": recorded_hash(run / "input.sha256")
                    == EXPECTED_INPUT[phase],
                    f"{label}:command_settings": all(
                        token in command
                        for token in (
                            "run --method gxtb",
                            "--acc 0.1",
                            "--iterations 300",
                            "--no-restart",
                            "--json result.json",
                        )
                    ),
                    f"{label}:energy_threshold": energy_match is not None
                    and Decimal(energy_match.group(1)) == Decimal("1e-7"),
                    f"{label}:density_threshold": density_match is not None
                    and Decimal(density_match.group(1)) == Decimal("2e-6"),
                    f"{label}:scc_converged": "SCC did not converge" not in output,
                    f"{label}:finite_energy": math.isfinite(float(energy)),
                    f"{label}:empty_stderr": (run / "process.err").stat().st_size == 0,
                }
            )
            raw_sha256[label] = {path.name: sha256(path) for path in required}

    pbc_cache = cmake_cache(HERE / "build/pbc-correct-CMakeCache.txt")
    no_mic_cache = cmake_cache(HERE / "build/pbc-without-mic-CMakeCache.txt")
    checks["build:relevant_cmake_options_identical"] = pbc_cache == no_mic_cache
    checks["source:pbc_head"] = (
        (HERE / "source/pbc-head.txt").read_text().strip()
        == EXPECTED_SOURCE["pbc_head"]
    )
    checks["source:reverted_commit"] = (
        (HERE / "source/reverted-commit.txt").read_text().strip()
        == EXPECTED_SOURCE["reverted_commit"]
    )
    checks["source:hybrid_tree"] = (
        (HERE / "source/hybrid-tree.txt").read_text().strip()
        == EXPECTED_SOURCE["hybrid_tree"]
    )
    checks["source:changed_files"] = (
        (HERE / "source/changed-files.txt").read_text(encoding="utf-8")
        == EXPECTED_CHANGED_FILES
    )
    checks["source:inverse_patch_hash"] = (
        sha256(HERE / "source/revert-083f220.patch")
        == EXPECTED_SOURCE["patch_sha256"]
    )
    checks["source:inverse_patch_mentions_mic_routine"] = (
        "get_amat_wsc_3d"
        in (HERE / "source/revert-083f220.patch").read_text(encoding="utf-8")
    )

    relative: dict[str, Decimal] = {}
    for provider in PROVIDERS:
        relative[provider] = (
            energies[(provider, "VII")] - energies[(provider, "Ih")]
        ) / WATERS * HARTREE_TO_KJ_MOL

    original_residual = relative["pbc-correct"] - relative["mstore-wsc-corrected"]
    no_mic_residual = (
        relative["pbc-without-mic"] - relative["mstore-wsc-corrected"]
    )
    mic_shift = relative["pbc-correct"] - relative["pbc-without-mic"]
    explained_percent = (
        Decimal(1) - abs(no_mic_residual) / abs(original_residual)
    ) * Decimal(100)
    checks["numerics:no_mic_residual_below_5e-5"] = (
        abs(no_mic_residual) < Decimal("5e-5")
    )
    checks["numerics:explained_above_99_999_percent"] = (
        explained_percent > Decimal("99.999")
    )

    with (HERE / "relative_energies.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("provider", "relative_energy_kj_mol_per_H2O"),
        )
        writer.writeheader()
        for provider in PROVIDERS:
            writer.writerow(
                {
                    "provider": provider,
                    "relative_energy_kj_mol_per_H2O": format(relative[provider], ".12f"),
                }
            )

    payload = {
        "schema": "periodic-gxtb-second-order-mic-attribution-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "test": "ice VII minus same-source ice Ih on an explicit 2x2x2 BvK supercell",
        "waters_per_cell": 96,
        "relative_energies_kj_mol_per_H2O": {
            key: float(value) for key, value in relative.items()
        },
        "pbc_correct_minus_mstore_wsc_corrected_kj_mol_per_H2O": float(
            original_residual
        ),
        "pbc_without_mic_minus_mstore_wsc_corrected_kj_mol_per_H2O": float(
            no_mic_residual
        ),
        "mic_variant_shift_kj_mol_per_H2O": float(mic_shift),
        "residual_explained_percent": float(explained_percent),
        "source": EXPECTED_SOURCE,
        "relevant_cmake_options": pbc_cache,
        "checks": checks,
        "raw_sha256": raw_sha256,
        "interpretation": (
            "After independently correcting the historical Wigner--Seitz self-image "
            "index, the remaining pbc/mstore relative-energy difference is reproduced "
            "by the later minimum-image second-order Coulomb variant. Reverting only "
            "that source change on pbc removes the residual to numerical SCC noise."
        ),
    }
    (HERE / "verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
