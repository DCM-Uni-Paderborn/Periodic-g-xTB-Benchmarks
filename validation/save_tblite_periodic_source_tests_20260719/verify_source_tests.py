#!/usr/bin/env python3
"""Verify the archived periodic save_tblite source tests."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path


CURRENT_REVISION = "15915c9435644eb257178ca8f8bf7220c38b1a84"
PBC_REVISION_PREFIX = "c932120d25"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def main() -> int:
    root = Path(__file__).resolve().parent
    for line in (root / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        path = root / name
        require(path.is_file(), f"missing archived file: {name}")
        require(
            hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            f"SHA-256 mismatch: {name}",
        )
    data = json.loads((root / "results.json").read_text())
    tests = data["tests"]

    require(data["current"]["revision"] == CURRENT_REVISION, "current revision")
    require(
        data["pbc_baseline"]["revision"].startswith(PBC_REVISION_PREFIX),
        "pbc baseline revision",
    )
    require(data["current"]["git_status"] == "", "current source is not clean")
    require(
        data["pbc_baseline"]["build_only_repair_files"]
        == ["test/unit/CMakeLists.txt", "test/unit/test_gxtb.f90"],
        "baseline repairs must be restricted to its broken test registration",
    )

    for name in ("current_h0_diamond", "current_h0_supercell", "current_h0_gradient"):
        require(tests[name]["returncode"] == 0, f"{name} return code")
        require(tests[name]["passed_count"] == 1, f"{name} pass count")
        require(tests[name]["failed_count"] == 0, f"{name} failure count")
    require(tests["current_wignerseitz"]["returncode"] == 0, "Wigner-Seitz")
    require(tests["current_wignerseitz"]["passed_count"] == 6, "Wigner-Seitz count")
    require(tests["current_exchange"]["returncode"] == 0, "exchange")
    require(tests["current_exchange"]["passed_count"] == 36, "exchange count")
    component_counts = {
        "current_acp": 22,
        "current_coulomb_charge": 58,
        "current_coulomb_multipole": 51,
        "current_dispersion": 15,
        "current_repulsion": 19,
        "current_gxtb": 40,
    }
    for name, expected_count in component_counts.items():
        require(tests[name]["returncode"] == 0, f"{name} return code")
        require(tests[name]["passed_count"] == expected_count, f"{name} pass count")
        require(tests[name]["failed_count"] == 0, f"{name} failure count")
    require(
        data["exchange_residuals_hartree"]["whole_mesh_vs_explicit_bvk"] <= 1.0e-9,
        "whole-mesh exchange residual",
    )
    require(
        data["exchange_residuals_hartree"]["mixed_radix_fft_vs_dense"] <= 1.0e-11,
        "FFT exchange residual",
    )

    require(tests["current_hamiltonian"]["passed_count"] == 74, "Hamiltonian pass count")
    require(tests["current_hamiltonian"]["failed_count"] == 1, "Hamiltonian failure count")
    require(tests["current_hamiltonian"]["returncode"] == 1, "Hamiltonian return code")
    for name in ("current_cecl3", "pbc_cecl3"):
        require(tests[name]["returncode"] == 1, f"{name} expected threshold exit")
        require(tests[name]["failed_count"] == 1, f"{name} failure count")
    cecl3 = data["cecl3_nonperiodic_finite_difference"]
    require(
        cecl3["identical_difference_components"],
        "current and pbc CeCl3 gradient differences differ",
    )
    require(
        cecl3["current_max_abs_hartree_per_bohr"]
        == cecl3["pbc_max_abs_hartree_per_bohr"],
        "current and pbc CeCl3 maxima differ",
    )
    require(
        cecl3["current_max_abs_hartree_per_bohr"] <= 3.0e-10,
        "CeCl3 compiler-sensitive residual is unexpectedly large",
    )

    print("PASS: periodic save_tblite source qualification")
    print("  H0 anisotropy: 3/3")
    print("  Wigner-Seitz: 6/6")
    print("  exchange: 36/36")
    print("  ACP/Coulomb/dispersion/repulsion/full g-xTB: 205/205")
    print("  Hamiltonian: 74/75; sole pbc-identical nonperiodic CeCl3 threshold case")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
