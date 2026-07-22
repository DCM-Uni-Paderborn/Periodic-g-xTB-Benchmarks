#!/usr/bin/env python3
"""Verify the compact range-local/mixed-radix combination archive."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ENERGY_RE = re.compile(
    r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+([-+0-9.Ee]+)"
)
FORWARD_RE = re.compile(
    r"KGROUP-PARTIAL-ROOT iter=\d+ dE=\s*([-+0-9.Ee]+) "
    r"dVsh=\s*([-+0-9.Ee]+) dFfold=\s*([-+0-9.Ee]+)"
)
REVERSE_RE = re.compile(
    r"KGROUP-PARTIAL-ROOT-REVERSE dOverlap=\s*([-+0-9.Ee]+) "
    r"dForce=\s*([-+0-9.Ee]+) dStress=\s*([-+0-9.Ee]+)"
)
ACP_RE = re.compile(r"ACP-SPARSE-REVERSE residual=\s*([-+0-9.Ee]+)")
MIXER_RE = re.compile(
    r"denseFullComplex=(\d+) streamedPeakComplex=(\d+)"
)
MATRIX_CASES = (
    "h2_1d_tr_3x1x1",
    "h2_2d_tr_3x1x3",
    "o2_uks_tr_3x1x1",
    "si_full_shifted_2x2x2",
    "ch4_spglib_shifted_4x4x4",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def verify_pair(platform: str, production: Path, qualification: Path) -> list[str]:
    prod = production.read_text(errors="replace")
    qual = qualification.read_text(errors="replace")
    for label, text in (("production", prod), ("qualification", qual)):
        require("PROGRAM ENDED AT" in text, f"{platform} {label} did not terminate normally")
        require(
            "stage=FORWARD backend=MIXED_RADIX_FFT fallback=0" in text,
            f"{platform} {label} lacks the exact forward transform marker",
        )
        require(
            "stage=REVERSE backend=MIXED_RADIX_FFT fallback=0" in text,
            f"{platform} {label} lacks the exact reverse transform marker",
        )
        require(
            "GXTB-KGROUP-PARTIAL-DISTRIBUTED-IMAGES importers=2, empty=0" in text,
            f"{platform} {label} lacks the two-importer range marker",
        )

    prod_energy = float(ENERGY_RE.findall(prod)[-1])
    qual_energy = float(ENERGY_RE.findall(qual)[-1])
    require(abs(prod_energy - qual_energy) <= 1.0e-12, f"{platform} energy mismatch")

    forward = [tuple(float(value) for value in match) for match in FORWARD_RE.findall(qual)]
    require(forward, f"{platform} qualification lacks forward residuals")
    max_de = max(values[0] for values in forward)
    max_dvsh = max(values[1] for values in forward)
    max_dfold = max(values[2] for values in forward)
    require(max_de <= 1.0e-12, f"{platform} forward energy residual")
    require(max_dvsh <= 1.0e-12, f"{platform} forward shell residual")
    require(max_dfold <= 1.0e-12, f"{platform} forward folded-Fock residual")

    reverse_match = REVERSE_RE.search(qual)
    require(reverse_match is not None, f"{platform} qualification lacks reverse residuals")
    reverse = tuple(float(value) for value in reverse_match.groups())
    require(max(reverse) <= 1.0e-12, f"{platform} reverse residual")

    acp_match = ACP_RE.search(qual)
    require(acp_match is not None, f"{platform} qualification lacks ACP residual")
    acp = float(acp_match.group(1))
    require(acp <= 1.0e-7, f"{platform} ACP residual")

    mixer = [tuple(int(value) for value in match) for match in MIXER_RE.findall(qual)]
    require(mixer and set(mixer) == {(512, 192)}, f"{platform} mixer storage counters")

    return [
        platform,
        f"{prod_energy:.15f}",
        f"{qual_energy:.15f}",
        f"{max_de:.16e}",
        f"{max_dvsh:.16e}",
        f"{max_dfold:.16e}",
        *(f"{value:.16e}" for value in reverse),
        f"{acp:.16e}",
        "512",
        "192",
    ]


def verify_matrix_case(platform: str, case: str, output: Path) -> list[str]:
    text = output.read_text(errors="replace")
    require("PROGRAM ENDED AT" in text, f"{platform} {case} did not terminate normally")
    require(
        "stage=FORWARD backend=MIXED_RADIX_FFT fallback=0" in text,
        f"{platform} {case} lacks the exact forward transform marker",
    )
    require(
        "stage=REVERSE backend=MIXED_RADIX_FFT fallback=0" in text,
        f"{platform} {case} lacks the exact reverse transform marker",
    )
    require(
        "GXTB-KGROUP-PARTIAL-DISTRIBUTED-IMAGES importers=2" in text,
        f"{platform} {case} lacks the two-importer range marker",
    )
    energies = ENERGY_RE.findall(text)
    require(energies, f"{platform} {case} lacks a final energy")
    forward = [tuple(float(value) for value in match) for match in FORWARD_RE.findall(text)]
    require(forward, f"{platform} {case} lacks forward residuals")
    forward_maxima = tuple(max(values[i] for values in forward) for i in range(3))
    require(max(forward_maxima) <= 1.0e-12, f"{platform} {case} forward residual")
    reverse_match = REVERSE_RE.search(text)
    require(reverse_match is not None, f"{platform} {case} lacks reverse residuals")
    reverse = tuple(float(value) for value in reverse_match.groups())
    require(max(reverse) <= 1.0e-12, f"{platform} {case} reverse residual")
    return [
        platform,
        case,
        f"{float(energies[-1]):.15f}",
        *(f"{value:.16e}" for value in forward_maxima),
        *(f"{value:.16e}" for value in reverse),
    ]


def main() -> None:
    pre_guard = (ROOT / "local/pre_guard_fix/cp2k.out").read_text(errors="replace")
    require(
        "Compact g-XTB transform has an incompatible exchange/gradient" in pre_guard
        and "combination" in pre_guard,
        "missing guard history",
    )
    require("PROGRAM ENDED AT" not in pre_guard, "guard history unexpectedly terminated normally")

    missing_spglib = (
        ROOT / "matrix/local_missing_spglib/ch4_spglib_shifted_4x4x4/cp2k.out"
    ).read_text(errors="replace")
    require(
        "SPGLIB k-point symmetry was requested, but SPGLIB is not available" in missing_spglib,
        "missing local SPGLIB fail-closed record",
    )
    fullgrid_guard = (
        ROOT / "matrix/local_pre_fullgrid_guard/si_full_shifted_2x2x2/cp2k.out"
    ).read_text(errors="replace")
    require(
        "SYMMETRY_STAR_CONTRACTION requires an actually" in fullgrid_guard,
        "missing full-grid streamed-star guard record",
    )
    acp_gate = (
        ROOT / "matrix/local_pre_acp_gate/h2_1d_tr_3x1x1/cp2k.out"
    ).read_text(errors="replace")
    require(
        "sparse ACP reverse exceeded its 1.0E-7" in acp_gate,
        "missing independent H2 ACP qualification finding",
    )

    rows = [
        verify_pair(
            "macOS",
            ROOT / "local/production/cp2k.out",
            ROOT / "local/qualification/cp2k.out",
        ),
        verify_pair(
            "Terok-Linux",
            ROOT / "linux/results/production/cp2k.out",
            ROOT / "linux/results/qualification/cp2k.out",
        ),
    ]
    header = [
        "platform",
        "production_energy_hartree",
        "qualification_energy_hartree",
        "max_forward_dE",
        "max_forward_dVsh",
        "max_forward_dFfold",
        "reverse_dOverlap",
        "reverse_dForce",
        "reverse_dStress",
        "acp_sparse_reverse_residual",
        "mixer_dense_full_complex",
        "mixer_streamed_peak_complex",
    ]
    rendered = "\n".join("\t".join(row) for row in (header, *rows)) + "\n"
    expected = (ROOT / "summary.tsv").read_text()
    require(rendered == expected, "summary.tsv does not match the raw outputs")
    print(rendered, end="")

    matrix_rows = []
    for case in MATRIX_CASES[:-1]:
        matrix_rows.append(
            verify_matrix_case("macOS", case, ROOT / f"matrix/local/{case}/cp2k.out")
        )
    for case in MATRIX_CASES:
        matrix_rows.append(
            verify_matrix_case(
                "Terok-Linux", case, ROOT / f"linux/results/matrix-{case}/cp2k.out"
            )
        )
    local_energy = {row[1]: float(row[2]) for row in matrix_rows if row[0] == "macOS"}
    linux_energy = {row[1]: float(row[2]) for row in matrix_rows if row[0] == "Terok-Linux"}
    for case, energy in local_energy.items():
        require(abs(energy - linux_energy[case]) <= 1.0e-10, f"cross-platform {case} energy")
    matrix_header = [
        "platform",
        "case",
        "energy_hartree",
        "max_forward_dE",
        "max_forward_dVsh",
        "max_forward_dFfold",
        "reverse_dOverlap",
        "reverse_dForce",
        "reverse_dStress",
    ]
    matrix_rendered = "\n".join(
        "\t".join(row) for row in (matrix_header, *matrix_rows)
    ) + "\n"
    expected_matrix = (ROOT / "matrix_summary.tsv").read_text()
    require(matrix_rendered == expected_matrix, "matrix_summary.tsv does not match raw outputs")
    print(matrix_rendered, end="")
    print("PASS: numerical, marker, termination, and storage gates")


if __name__ == "__main__":
    main()
