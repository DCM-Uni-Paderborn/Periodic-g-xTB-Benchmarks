#!/usr/bin/env python3
"""Prove exact MacDonald/native-Bloch and Gamma-supercell BvK equivalence."""

from __future__ import annotations

import hashlib
import itertools
import json
import re
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = ROOT / "DMC-ICE13/reproduction/seidler_dmc13_recalculation/raw/cp2k_native"
SOURCE_REVISION = "8520b2e592cd04d35081ab4ad46d92c606071e23"
SCHEME_RE = re.compile(
    r"^\s*SCHEME\s+MACDONALD\s+(\d+)\s+(\d+)\s+(\d+)\s+"
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*$",
    re.IGNORECASE | re.MULTILINE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decimal_fraction(token: str) -> Fraction:
    return Fraction(token)


def expected_shift(n: int) -> Fraction:
    return Fraction(0) if n % 2 else Fraction(n - 1, 2 * n)


def native_coordinate(n: int, i: int, shift: Fraction) -> Fraction:
    # Production inputs do not enable GAMMA_CENTERED, so the ELSE branch in
    # full_grid_gen is the authoritative coordinate rule.
    return (Fraction(2 * i - n - 1, 2 * n) + shift) % 1


def bvk_coordinate(n: int, index: int) -> Fraction:
    return Fraction(index, n)


def coordinate_set(n: int, shift: Fraction) -> tuple[Fraction, ...]:
    return tuple(sorted(native_coordinate(n, i, shift) for i in range(1, n + 1)))


def bvk_set(n: int) -> tuple[Fraction, ...]:
    return tuple(bvk_coordinate(n, i) for i in range(n))


def parse_input(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    match = SCHEME_RE.search(text)
    if match is None:
        raise AssertionError(f"missing MACDONALD scheme: {path}")
    mesh = tuple(int(match.group(i)) for i in range(1, 4))
    shifts = tuple(decimal_fraction(match.group(i)) for i in range(4, 7))
    gamma_enabled = bool(
        re.search(r"^\s*GAMMA_CENTERED\s+(?:T|TRUE|ON|YES)\s*$", text, re.I | re.M)
    )
    spglib = bool(
        re.search(
            r"^\s*SYMMETRY_REDUCTION_METHOD\s+SPGLIB\s*$", text, re.I | re.M
        )
    )
    return {
        "mesh": mesh,
        "shifts": shifts,
        "gamma_centered_enabled": gamma_enabled,
        "spglib_reduction_requested": spglib,
    }


def main() -> None:
    excerpt = HERE / "cryssym_full_grid_gen_excerpt.F90"
    excerpt_text = excerpt.read_text(encoding="utf-8")
    source_formula_present = all(
        token in excerpt_text
        for token in (
            SOURCE_REVISION,
            "2*ik(idim) - nk(idim) - 1",
            "xkp(1:3, i) = xkp(1:3, i) + shift(1:3)",
        )
    )

    analytic_meshes = {}
    all_analytic_pass = True
    for n in range(1, 10):
        shift = expected_shift(n)
        native = coordinate_set(n, shift)
        bvk = bvk_set(n)
        one_dimensional_pass = native == bvk
        native_3d = set(itertools.product(native, repeat=3))
        bvk_3d = set(itertools.product(bvk, repeat=3))
        three_dimensional_pass = native_3d == bvk_3d and len(native_3d) == n**3
        passed = one_dimensional_pass and three_dimensional_pass
        all_analytic_pass = all_analytic_pass and passed
        analytic_meshes[str(n)] = {
            "expected_shift": str(shift),
            "point_count_3d": len(native_3d),
            "one_dimensional_native_modulo_one": [str(value) for value in native],
            "one_dimensional_bvk": [str(value) for value in bvk],
            "one_dimensional_pass": one_dimensional_pass,
            "three_dimensional_pass": three_dimensional_pass,
        }

    input_paths = sorted(RAW.glob("k???-reduced/*/input.inp"))
    if not input_paths:
        raise AssertionError(f"no archived production inputs under {RAW}")
    input_results = {}
    input_counts = {}
    all_inputs_pass = True
    for path in input_paths:
        relative = str(path.relative_to(ROOT))
        parsed = parse_input(path)
        nx, ny, nz = parsed["mesh"]
        directory_code = path.parents[1].name[1:4]
        directory_is_cubic = len(directory_code) == 3 and len(set(directory_code)) == 1
        mesh_from_directory = int(directory_code[0]) if directory_is_cubic else -1
        cubic = nx == ny == nz == mesh_from_directory
        shifts = parsed["shifts"]
        expected = expected_shift(nx) if cubic else None
        # The input is parsed into binary64 by CP2K.  For repeating fractions
        # (5/12 at N=6), the written decimal is not rationally identical to the
        # ideal shift, but it rounds to exactly the same binary64 value.
        shifts_round_to_expected = cubic and all(
            float(value) == float(expected) for value in shifts
        )
        coordinate_residuals = []
        for n, shift in zip(parsed["mesh"], shifts):
            native = coordinate_set(n, shift)
            target = bvk_set(n)
            coordinate_residuals.extend(abs(a - b) for a, b in zip(native, target))
        max_coordinate_residual = max(coordinate_residuals, default=Fraction(0))
        coordinates_match = cubic and float(max_coordinate_residual) <= 1.0e-15
        passed = (
            cubic
            and shifts_round_to_expected
            and coordinates_match
            and not parsed["gamma_centered_enabled"]
            and parsed["spglib_reduction_requested"]
        )
        all_inputs_pass = all_inputs_pass and passed
        input_counts[str(nx)] = input_counts.get(str(nx), 0) + 1
        input_results[relative] = {
            **parsed,
            "mesh": list(parsed["mesh"]),
            "shifts": [str(value) for value in shifts],
            "expected_shift": str(expected) if expected is not None else None,
            "input_shifts_round_to_expected_binary64": shifts_round_to_expected,
            "max_text_rational_coordinate_residual": str(max_coordinate_residual),
            "max_text_rational_coordinate_residual_float": float(max_coordinate_residual),
            "coordinates_equal_gamma_bvk": coordinates_match,
            "passed": passed,
            "sha256": sha256(path),
        }

    passed = source_formula_present and all_analytic_pass and all_inputs_pass
    output = {
        "schema": "periodic-gxtb-macdonald-bvk-equivalence-v1",
        "status": "PASS" if passed else "FAIL",
        "qualified_cp2k_source_revision": SOURCE_REVISION,
        "source_excerpt": {
            "file": str(excerpt.relative_to(ROOT)),
            "sha256": sha256(excerpt),
            "authoritative_formula_present": source_formula_present,
        },
        "proof_domain": "exact rational arithmetic modulo reciprocal lattice vectors",
        "analytic_meshes": analytic_meshes,
        "archived_input_count": len(input_paths),
        "archived_input_count_by_mesh": input_counts,
        "all_archived_inputs_pass": all_inputs_pass,
        "archived_inputs": input_results,
        "interpretation": (
            "Every archived native DMC-ICE13 MacDonald mesh is exactly the "
            "Gamma-supercell BvK folding grid before symmetry reduction."
        ),
    }
    output_path = HERE / "verification.json"
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
