#!/usr/bin/env python3
"""Fail-closed verifier for distributed-image partial importers."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MATRIX = json.loads((ROOT / "matrix.json").read_text())
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
FORWARD_RE = re.compile(
    rf"GXTB-QUALIFICATION_ONLY KGROUP-PARTIAL-ROOT iter=(\d+)"
    rf"\s+dE=\s*({FLOAT})\s+dVsh=\s*({FLOAT})\s+dFfold=\s*({FLOAT})"
)
REVERSE_RE = re.compile(
    rf"GXTB-QUALIFICATION_ONLY KGROUP-PARTIAL-ROOT-REVERSE"
    rf"\s+dOverlap=\s*({FLOAT})\s+dForce=\s*({FLOAT})\s+dStress=\s*({FLOAT})"
)
MODE_RE = re.compile(
    r"GXTB-KGROUP-PARTIAL-DISTRIBUTED-IMAGES importers=(\d+), empty=(\d+), "
    r"nred=(\d+), nfull=(\d+), batch=(\d+);"
)
REVERSE_MODE_RE = re.compile(
    r"GXTB KGROUP-PARTIAL-DISTRIBUTED-IMAGES-REVERSE nfull=(\d+), importers=(\d+), "
    r"empty=(\d+), batch=(\d+), "
    r"(?:provider_)?peak_complex=(\d+);"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expanded_cases() -> list[dict]:
    result = []
    for case in MATRIX["cases"]:
        for ranks in case["ranks"]:
            result.append({**case, "ranks": int(ranks)})
    return result


def parse_observables(path: Path) -> tuple[float, list[tuple[float, float, float]], list[tuple[float, float, float]]]:
    text = path.read_text(errors="replace")
    if text.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"incomplete CP2K output: {path}")
    energies = [float(value) for value in re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )]
    force_matches = re.findall(
        rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
        text,
        re.MULTILINE,
    )
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    if not energies or not force_matches or not stress_blocks:
        raise RuntimeError(f"missing energy, force, or stress block: {path}")
    stress_matches = re.findall(
        rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        stress_blocks[-1],
        re.MULTILINE,
    )
    if len(stress_matches) != 3:
        raise RuntimeError(f"malformed final stress tensor: {path}")
    force = [tuple(map(float, row)) for row in force_matches]
    stress = [tuple(map(float, row)) for row in stress_matches]
    values = [energies[-1], *(value for row in force for value in row), *(value for row in stress for value in row)]
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"non-finite final observable: {path}")
    return energies[-1], force, stress


def max_delta(left: list[tuple[float, float, float]], right: list[tuple[float, float, float]]) -> float:
    if len(left) != len(right):
        raise RuntimeError(f"block lengths differ: {len(left)} != {len(right)}")
    return max(abs(a - b) for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b))


def checked_metadata(run_dir: Path, variant: str, case: dict) -> dict:
    meta_path = run_dir / "run.json"
    rc_path = run_dir / "returncode.txt"
    output_path = run_dir / "cp2k.out"
    stderr_path = run_dir / "cp2k.err"
    for path in (meta_path, rc_path, output_path, stderr_path):
        if not path.is_file():
            raise RuntimeError(f"missing result file: {path}")
    if rc_path.read_text().strip() != "0":
        raise RuntimeError(f"nonzero return code: {run_dir}")
    meta = json.loads(meta_path.read_text())
    if meta.get("returncode") != 0 or meta.get("variant") != variant:
        raise RuntimeError(f"metadata result/variant mismatch: {run_dir}")
    if meta.get("ranks") != case["ranks"] or meta.get("case") != case["name"]:
        raise RuntimeError(f"metadata case/rank mismatch: {run_dir}")
    if meta.get("expected_nfull") != case["nfull"]:
        raise RuntimeError(f"metadata nfull mismatch: {run_dir}")
    if meta.get("input_sha256") != sha256(ROOT / "inputs" / case["input"]):
        raise RuntimeError(f"input hash mismatch: {run_dir}")
    if meta.get("output_sha256") != sha256(output_path) or meta.get("stderr_sha256") != sha256(stderr_path):
        raise RuntimeError(f"raw-output hash mismatch: {run_dir}")
    return meta


def internal_residuals(path: Path, case: dict, batch_size: int) -> tuple[float, ...]:
    text = path.read_text(errors="replace")
    forwards = [tuple(map(float, match[1:])) for match in FORWARD_RE.findall(text)]
    reverses = [tuple(map(float, match)) for match in REVERSE_RE.findall(text)]
    modes = [tuple(map(int, match)) for match in MODE_RE.findall(text)]
    reverse_modes = [tuple(map(int, match)) for match in REVERSE_MODE_RE.findall(text)]
    if not forwards or not reverses or not modes or not reverse_modes:
        raise RuntimeError(f"missing replicated-importer qualification marker: {path}")
    expected_batch = min(batch_size, case["nfull"])
    if not all(nfull == case["nfull"] and batch == expected_batch
               for _, _, _, nfull, batch in modes):
        raise RuntimeError(f"forward mode metadata mismatch: {path}")
    if not all(nfull == case["nfull"] and batch == expected_batch
               for nfull, _, _, batch, _ in reverse_modes):
        raise RuntimeError(f"reverse mode metadata mismatch: {path}")
    importer_values = {importers for importers, _, _, _, _ in modes} | {
        importers for _, importers, _, _, _ in reverse_modes
    }
    if importer_values != {case["ranks"]}:
        raise RuntimeError(f"inconsistent/invalid importer count: {path}")
    expected_empty = max(0, case["ranks"] - case["nfull"])
    empty_values = {empty for _, empty, _, _, _ in modes} | {
        empty for _, _, empty, _, _ in reverse_modes
    }
    if empty_values != {expected_empty}:
        raise RuntimeError(f"inconsistent/invalid empty-importer count: {path}")
    nred_values = {nred for _, _, nred, _, _ in modes}
    if len(nred_values) != 1 or not (1 <= next(iter(nred_values)) <= case["nfull"]):
        raise RuntimeError(f"inconsistent/invalid irreducible mesh size: {path}")
    values = [value for row in forwards + reverses for value in row]
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise RuntimeError(f"invalid internal residual: {path}")
    maxima = tuple(max(row[index] for row in forwards) for index in range(3)) + tuple(
        max(row[index] for row in reverses) for index in range(3)
    )
    forward_gate = float(MATRIX["gates"]["internal_forward"])
    reverse_gate = float(MATRIX["gates"]["internal_reverse"])
    if max(maxima[:3]) > forward_gate or max(maxima[3:]) > reverse_gate:
        raise RuntimeError(f"internal oracle gate failed: {path}: {maxima}")
    return maxima


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-partial", action="store_true", help="verify only pairs present on disk")
    args = parser.parse_args()
    rows = []
    for case in expanded_cases():
        stem = f"{case['name']}_p{case['ranks']}"
        dense_dir = ROOT / "runs" / f"{stem}_dense"
        partial_dir = ROOT / "runs" / f"{stem}_partial_distributed_images"
        if args.allow_partial and (not dense_dir.is_dir() or not partial_dir.is_dir()):
            continue
        dense_meta = checked_metadata(dense_dir, "DENSE", case)
        partial_meta = checked_metadata(partial_dir, "PARTIAL_DISTRIBUTED_IMAGES", case)
        if dense_meta["cp2k_sha256"] != partial_meta["cp2k_sha256"]:
            raise RuntimeError(f"DENSE/partial executable hashes differ: {stem}")

        dense = parse_observables(dense_dir / "cp2k.out")
        partial = parse_observables(partial_dir / "cp2k.out")
        d_energy = abs(dense[0] - partial[0])
        d_force = max_delta(dense[1], partial[1])
        d_stress = max_delta(dense[2], partial[2])
        internal = internal_residuals(partial_dir / "cp2k.out", case, int(partial_meta["batch_size"]))
        gates = MATRIX["gates"]
        if d_energy > gates["external_energy_Ha"]:
            raise RuntimeError(f"external energy gate failed: {stem}: {d_energy}")
        if d_force > gates["external_force_Ha_per_bohr"]:
            raise RuntimeError(f"external force gate failed: {stem}: {d_force}")
        if d_stress > gates["external_stress_bar"]:
            raise RuntimeError(f"external stress gate failed: {stem}: {d_stress}")
        rows.append({
            "case": case["name"],
            "ranks": case["ranks"],
            "features": ",".join(case["features"]),
            "nfull": case["nfull"],
            "dense_wall_s": f"{dense_meta['wall_seconds']:.9f}",
            "partial_wall_s": f"{partial_meta['wall_seconds']:.9f}",
            "external_dE_Ha": f"{d_energy:.16e}",
            "external_dForce_Ha_per_bohr": f"{d_force:.16e}",
            "external_dStress_bar": f"{d_stress:.16e}",
            "internal_dE_Ha": f"{internal[0]:.16e}",
            "internal_dVsh_Ha": f"{internal[1]:.16e}",
            "internal_dFfold_Ha": f"{internal[2]:.16e}",
            "internal_dOverlap": f"{internal[3]:.16e}",
            "internal_dForce_Ha_per_bohr": f"{internal[4]:.16e}",
            "internal_dStress_Ha": f"{internal[5]:.16e}",
            "status": "PASS",
        })

    expected = sum(len(case["ranks"]) for case in MATRIX["cases"])
    if not args.allow_partial and len(rows) != expected:
        raise RuntimeError(f"matrix incomplete: {len(rows)} != {expected}")
    if not rows:
        raise RuntimeError("no complete DENSE/partial pairs found")
    fields = list(rows[0])
    with (ROOT / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)}/{expected} DENSE/KGROUP_PARTIAL_DISTRIBUTED_IMAGES pairs")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
