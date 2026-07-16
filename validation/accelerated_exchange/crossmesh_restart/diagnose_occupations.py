#!/usr/bin/env python3
"""Diagnose metric density eigenvalues of a transferred v2 restart.

The restart P(R) is evaluated at every irreducible target k point exactly as
`apply_bvk_transfer` does.  S(k) and target weights are read from CP2K's
high-precision MO_KP/OVERLAP_MATRIX output.  We report eigenvalues of
S(k)^(1/2) P_sigma(k) S(k)^(1/2), before and after CP2K's scalar electron-count
normalization rule.
"""

from __future__ import annotations

import argparse
import pathlib
import re

import numpy as np

from parse_restart import parse


def parse_mokp(path: pathlib.Path) -> tuple[list[tuple[int, np.ndarray, float]], dict[int, np.ndarray]]:
    lines = path.read_text().splitlines()
    dims = next(line for line in lines if line.startswith("# DIMENSIONS:"))
    _, _, nao, nkp = map(int, re.findall(r"[-+]?\d+", dims.split("=")[-1]))
    kpoints: list[tuple[int, np.ndarray, float]] = []
    overlaps: dict[int, np.ndarray] = {}
    i = lines.index("# KPOINT_LIST: ikp  kx  ky  kz  weight") + 1
    while i < len(lines) and not lines[i].startswith("#"):
        fields = lines[i].split()
        if fields:
            kpoints.append((int(fields[0]), np.array(list(map(float, fields[1:4]))), float(fields[4])))
        i += 1
    if len(kpoints) != nkp:
        raise ValueError(f"expected {nkp} k points, got {len(kpoints)}")
    i = 0
    while i < len(lines):
        if not lines[i].startswith("# BEGIN_OVERLAP"):
            i += 1
            continue
        ikp = int(re.findall(r"\d+", lines[i])[0])
        matrix = np.zeros((nao, nao), dtype=complex)
        i += 2  # skip OVERLAP_RE heading
        while i < len(lines) and not lines[i].startswith("# OVERLAP_IM"):
            fields = lines[i].split()
            if fields:
                row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
                matrix[row - 1, col - 1] = value
            i += 1
        i += 1
        while i < len(lines) and not lines[i].startswith("# END_OVERLAP"):
            fields = lines[i].split()
            if fields:
                row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
                matrix[row - 1, col - 1] += 1j * value
            i += 1
        overlaps[ikp] = matrix
        i += 1
    if set(overlaps) != {item[0] for item in kpoints}:
        raise ValueError("missing overlap blocks")
    return kpoints, overlaps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("restart", type=pathlib.Path, help="source-mesh v2 density restart")
    parser.add_argument("mokp", type=pathlib.Path, help="target-mesh MO_KP overlap dump")
    args = parser.parse_args()
    restart = parse(args.restart)
    kpoints, overlaps = parse_mokp(args.mokp)
    nspin = restart["nspin"]
    expected = [sum(restart["nelectron"])] if nspin == 1 else list(restart["nelectron"][:nspin])
    density: dict[tuple[int, int], np.ndarray] = {}
    hermiticity = 0.0
    traces = np.zeros(nspin)
    for ikp, kvec, weight in kpoints:
        overlap = 0.5 * (overlaps[ikp] + overlaps[ikp].conj().T)
        for ispin in range(1, nspin + 1):
            pmat = np.zeros((restart["nao"], restart["nao"]), dtype=complex)
            for (payload_spin, cell), block in restart["payload"].items():
                if payload_spin == ispin:
                    pmat += np.exp(2j * np.pi * np.dot(kvec, cell)) * block
            scale = max(1.0, float(np.max(np.abs(pmat))))
            hermiticity = max(hermiticity, float(np.max(np.abs(pmat - pmat.conj().T))) / scale)
            pmat = 0.5 * (pmat + pmat.conj().T)
            density[ispin, ikp] = pmat
            traces[ispin - 1] += weight * np.trace(pmat @ overlap).real
    factors = np.ones(nspin)
    for ispin in range(nspin):
        error = abs(traces[ispin] - expected[ispin])
        if error > 1.0e-10:
            factors[ispin] = expected[ispin] / traces[ispin]
    raw_min = np.inf
    raw_max = -np.inf
    normalized_min = np.inf
    normalized_max = -np.inf
    overlap_min = np.inf
    for ikp, _, _ in kpoints:
        overlap = 0.5 * (overlaps[ikp] + overlaps[ikp].conj().T)
        seig, svec = np.linalg.eigh(overlap)
        overlap_min = min(overlap_min, float(seig.min()))
        if seig.min() <= 0.0:
            raise ValueError(f"non-positive target overlap at k point {ikp}: {seig.min()}")
        shalf = (svec * np.sqrt(seig)) @ svec.conj().T
        for ispin in range(1, nspin + 1):
            qmat = shalf @ density[ispin, ikp] @ shalf
            qeig = np.linalg.eigvalsh(0.5 * (qmat + qmat.conj().T))
            raw_min = min(raw_min, float(qeig.min()))
            raw_max = max(raw_max, float(qeig.max()))
            qeig *= factors[ispin - 1]
            normalized_min = min(normalized_min, float(qeig.min()))
            normalized_max = max(normalized_max, float(qeig.max()))
    upper = 2.0 if nspin == 1 else 1.0
    print(f"restart={args.restart}")
    print(f"mokp={args.mokp}")
    print(f"nspin={nspin} target_irreducible_kpoints={len(kpoints)} occupation_upper_bound={upper:.1f}")
    print(f"hermiticity_relative_max={hermiticity:.17e}")
    print(f"overlap_eigenvalue_min={overlap_min:.17e}")
    for ispin in range(nspin):
        print(
            f"spin={ispin + 1} trace_before={traces[ispin]:.17e} "
            f"expected={expected[ispin]:.17e} normalization_factor={factors[ispin]:.17e}"
        )
    print(f"metric_density_raw_min={raw_min:.17e}")
    print(f"metric_density_raw_max={raw_max:.17e}")
    print(f"metric_density_normalized_min={normalized_min:.17e}")
    print(f"metric_density_normalized_max={normalized_max:.17e}")
    print(f"lower_bound_violation={max(0.0, -normalized_min):.17e}")
    print(f"upper_bound_violation={max(0.0, normalized_max - upper):.17e}")


if __name__ == "__main__":
    main()
