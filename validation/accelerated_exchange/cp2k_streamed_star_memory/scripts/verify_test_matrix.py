#!/usr/bin/env python3
"""Fail-closed verifier for mixer symmetry-star storage qualification."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
import tarfile
import tempfile
from pathlib import Path


ARCHIVE_ROOT = Path(__file__).resolve().parent.parent
RAW_ARCHIVE = (
    ARCHIVE_ROOT
    / "raw_archive"
    / "gxtb_streamed_star_memory_evidence_20260717.tar.gz"
)
_EXTRACTED = tempfile.TemporaryDirectory(prefix="gxtb-streamed-star-")
_EXTRACTED_ROOT = Path(_EXTRACTED.name).resolve()
with tarfile.open(RAW_ARCHIVE, "r:gz") as archive:
    members = archive.getmembers()
    for member in members:
        if member.issym() or member.islnk():
            raise RuntimeError(f"archive link is not permitted: {member.name}")
        destination = (_EXTRACTED_ROOT / member.name).resolve()
        if destination != _EXTRACTED_ROOT and _EXTRACTED_ROOT not in destination.parents:
            raise RuntimeError(f"archive path escapes extraction root: {member.name}")
    archive.extractall(_EXTRACTED_ROOT, members=members)

ROOT = _EXTRACTED_ROOT / "gxtb_streamed_star_memory_qualification_v4_20260717"
MATRIX = json.loads((ROOT / "test_matrix.json").read_text())
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
QUALIFY_RE = re.compile(
    rf"GXTB-QUALIFICATION_ONLY MIXER-STAR iter=(\d+)"
    rf"\s+denseCov=\s*({FLOAT})\s+streamCov=\s*({FLOAT})"
    rf"\s+streamRoundtrip=\s*({FLOAT})\s+covDelta=\s*({FLOAT})"
    rf"\s+denseFullComplex=(\d+)\s+streamedPeakComplex=(\d+)"
)
STREAM_RE = re.compile(
    rf"GXTB-MIXER-STAR-STREAMED denseFullComplexAvoided=(\d+),"
    rf" peakComplex=(\d+), covariance=\s*({FLOAT}), roundtrip=\s*({FLOAT})"
)
MODE_RE = re.compile(
    r"GXTB-KGROUP-PARTIAL-ROOT groups=(\d+), nred=(\d+), nfull=(\d+), batch=(\d+);"
)


def dependency_hash(path_suffix: str) -> str:
    matches = [
        line.split()[0]
        for manifest in ("dependencies.sha256", "terok_provenance.sha256")
        for line in (ARCHIVE_ROOT / manifest).read_text().splitlines()
        if line.split(maxsplit=1)[1].endswith(path_suffix)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one dependency hash for {path_suffix}")
    return matches[0]


EXPECTED_CP2K_SHA256 = dependency_hash("/bin/cp2k.psmp")
EXPECTED_CP2K_LIB_SHA256 = dependency_hash("/src/libcp2k.so.2026.2")


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


def checked_run(case: dict, variant: str) -> tuple[dict, str]:
    stem = f"{case['name']}_p{case['ranks']}_{variant.lower()}"
    run_dir = ROOT / "runs" / stem
    required = [run_dir / name for name in ("run.json", "returncode.txt", "cp2k.out", "cp2k.err")]
    if not all(path.is_file() for path in required):
        raise RuntimeError(f"missing result file in {run_dir}")
    if (run_dir / "returncode.txt").read_text().strip() != "0":
        raise RuntimeError(f"nonzero return code: {run_dir}")
    metadata = json.loads((run_dir / "run.json").read_text())
    if metadata.get("returncode") != 0 or metadata.get("variant") != variant:
        raise RuntimeError(f"metadata variant/result mismatch: {run_dir}")
    if metadata.get("case") != case["name"] or metadata.get("ranks") != case["ranks"]:
        raise RuntimeError(f"metadata case/rank mismatch: {run_dir}")
    if metadata.get("input_sha256") != sha256(ROOT / "inputs" / case["input"]):
        raise RuntimeError(f"input hash mismatch: {run_dir}")
    if metadata.get("cp2k_sha256") != EXPECTED_CP2K_SHA256:
        raise RuntimeError(f"CP2K executable hash mismatch: {run_dir}")
    if metadata.get("cp2k_lib_sha256") != EXPECTED_CP2K_LIB_SHA256:
        raise RuntimeError(f"CP2K shared-library hash mismatch: {run_dir}")
    affinity = metadata.get("affinity_proof")
    if not isinstance(affinity, list) or len(affinity) != case["ranks"]:
        raise RuntimeError(f"missing live rank-affinity proof: {run_dir}")
    first_cpu, last_cpu = map(int, metadata["cpu_set"].split("-"))
    expected_cpus = set(range(first_cpu, last_cpu + 1))
    if any(set(item.get("cpus_allowed", [])) != expected_cpus or
           item.get("processor") not in expected_cpus for item in affinity):
        raise RuntimeError(f"invalid live rank-affinity proof: {run_dir}")
    if metadata.get("output_sha256") != sha256(run_dir / "cp2k.out"):
        raise RuntimeError(f"output hash mismatch: {run_dir}")
    if metadata.get("stderr_sha256") != sha256(run_dir / "cp2k.err"):
        raise RuntimeError(f"stderr hash mismatch: {run_dir}")
    text = (run_dir / "cp2k.out").read_text(errors="replace")
    if text.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"incomplete output: {run_dir}")
    modes = [tuple(map(int, match)) for match in MODE_RE.findall(text)]
    if not modes or any(nfull != case["nfull"] for _, _, nfull, _ in modes):
        raise RuntimeError(f"missing/mismatched partial-root mode marker: {run_dir}")
    return metadata, text


def observables(text: str) -> tuple[float, list[float], list[float]]:
    energies = [float(value) for value in re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )]
    forces = [tuple(map(float, row)) for row in re.findall(
        rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
        text,
        re.MULTILINE,
    )]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    if not energies or not forces or not stress_blocks:
        raise RuntimeError("missing energy, force, or stress observable")
    stress_rows = re.findall(
        rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        stress_blocks[-1],
        re.MULTILINE,
    )
    if len(stress_rows) != 3:
        raise RuntimeError("malformed stress block")
    force_values = [value for row in forces for value in row]
    stress_values = [float(value) for row in stress_rows for value in row]
    values = [energies[-1], *force_values, *stress_values]
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError("non-finite observable")
    return energies[-1], force_values, stress_values


def max_delta(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise RuntimeError("observable block lengths differ")
    return max(abs(a - b) for a, b in zip(left, right))


def star_residuals(case: dict, variant: str, text: str) -> tuple[float, float, float, int, int]:
    if variant == "STREAMED":
        matches = STREAM_RE.findall(text)
        if len(matches) != 1:
            raise RuntimeError("STREAMED run lacks exactly one selector marker")
        dense_full, streamed_peak, covariance, roundtrip = matches[0]
        dense_residual = float("nan")
        stream_residual = float(covariance)
        roundtrip_residual = float(roundtrip)
    elif variant == "QUALIFY":
        matches = QUALIFY_RE.findall(text)
        if not matches:
            raise RuntimeError("QUALIFY run lacks selector markers")
        dense_residual = max(float(match[1]) for match in matches)
        stream_residual = max(float(match[2]) for match in matches)
        roundtrip_residual = max(float(match[3]) for match in matches)
        covariance_delta = max(float(match[4]) for match in matches)
        if covariance_delta > MATRIX["gates"]["internal_covariance"]:
            raise RuntimeError(f"dense/stream covariance delta failed: {covariance_delta}")
        dense_values = {int(match[5]) for match in matches}
        stream_values = {int(match[6]) for match in matches}
        if len(dense_values) != 1 or len(stream_values) != 1:
            raise RuntimeError("memory counters changed during QUALIFY run")
        dense_full = str(next(iter(dense_values)))
        streamed_peak = str(next(iter(stream_values)))
    else:
        return float("nan"), float("nan"), float("nan"), 0, 0
    dense_full_i = int(dense_full)
    streamed_peak_i = int(streamed_peak)
    nspin = 2 if "UKS" in case["features"] else 1
    if dense_full_i * 3 != streamed_peak_i * nspin * case["nfull"]:
        raise RuntimeError("reported memory counters violate exact allocation formula")
    for value in (stream_residual, roundtrip_residual):
        if not math.isfinite(value) or value < 0.0:
            raise RuntimeError("invalid streamed residual")
    if stream_residual > MATRIX["gates"]["internal_covariance"]:
        raise RuntimeError(f"streamed covariance gate failed: {stream_residual}")
    if roundtrip_residual > MATRIX["gates"]["internal_roundtrip"]:
        raise RuntimeError(f"streamed roundtrip gate failed: {roundtrip_residual}")
    return dense_residual, stream_residual, roundtrip_residual, dense_full_i, streamed_peak_i


def main() -> int:
    rows = []
    executable_hashes = set()
    library_hashes = set()
    for case in expanded_cases():
        runs = {}
        for variant in ("DENSE", "STREAMED", "QUALIFY"):
            metadata, text = checked_run(case, variant)
            executable_hashes.add(metadata["cp2k_sha256"])
            library_hashes.add(metadata["cp2k_lib_sha256"])
            runs[variant] = (metadata, text, observables(text))
        if len(executable_hashes) != 1:
            raise RuntimeError("matrix used more than one CP2K executable")
        if len(library_hashes) != 1:
            raise RuntimeError("matrix used more than one CP2K shared library")
        dense_obs = runs["DENSE"][2]
        streamed_obs = runs["STREAMED"][2]
        qualify_obs = runs["QUALIFY"][2]
        d_energy = max(abs(dense_obs[0] - streamed_obs[0]), abs(dense_obs[0] - qualify_obs[0]))
        d_force = max(max_delta(dense_obs[1], streamed_obs[1]), max_delta(dense_obs[1], qualify_obs[1]))
        d_stress = max(max_delta(dense_obs[2], streamed_obs[2]), max_delta(dense_obs[2], qualify_obs[2]))
        gates = MATRIX["gates"]
        if d_energy > gates["external_energy_Ha"]:
            raise RuntimeError(f"energy gate failed for {case['name']}_p{case['ranks']}: {d_energy}")
        if d_force > gates["external_force_Ha_per_bohr"]:
            raise RuntimeError(f"force gate failed for {case['name']}_p{case['ranks']}: {d_force}")
        if d_stress > gates["external_stress_bar"]:
            raise RuntimeError(f"stress gate failed for {case['name']}_p{case['ranks']}: {d_stress}")
        _, stream_cov, stream_roundtrip, dense_full, streamed_peak = star_residuals(
            case, "STREAMED", runs["STREAMED"][1]
        )
        dense_cov, qualify_cov, qualify_roundtrip, qualify_dense, qualify_peak = star_residuals(
            case, "QUALIFY", runs["QUALIFY"][1]
        )
        if (dense_full, streamed_peak) != (qualify_dense, qualify_peak):
            raise RuntimeError("STREAMED/QUALIFY memory counters differ")
        rows.append({
            "case": case["name"],
            "ranks": case["ranks"],
            "features": ",".join(case["features"]),
            "nfull": case["nfull"],
            "external_dE_Ha": f"{d_energy:.16e}",
            "external_dForce_Ha_per_bohr": f"{d_force:.16e}",
            "external_dStress_bar": f"{d_stress:.16e}",
            "dense_covariance": f"{dense_cov:.16e}",
            "stream_covariance": f"{max(stream_cov, qualify_cov):.16e}",
            "stream_roundtrip": f"{max(stream_roundtrip, qualify_roundtrip):.16e}",
            "dense_full_complex": dense_full,
            "streamed_peak_complex": streamed_peak,
            "status": "PASS",
        })
    expected = sum(len(case["ranks"]) for case in MATRIX["cases"])
    if len(rows) != expected:
        raise RuntimeError(f"matrix incomplete: {len(rows)} != {expected}")
    with (ROOT / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)}/{expected} DENSE/STREAMED/QUALIFY triples")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        raise
