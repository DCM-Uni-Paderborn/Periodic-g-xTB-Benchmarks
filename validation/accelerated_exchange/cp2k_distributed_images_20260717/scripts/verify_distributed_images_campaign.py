#!/usr/bin/env python3
"""Verify every distributed-image run, fault gate, and affinity proof."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[1]
MATRIX = json.loads((CAMPAIGN / "harness" / "campaign_matrix.json").read_text())
RUNS = CAMPAIGN / "formal_runs"
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
    r"empty=(\d+), batch=(\d+), provider_peak_complex=(\d+);"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_observables(path: Path) -> tuple[float, list[tuple[float, float, float]], list[tuple[float, float, float]]]:
    text = path.read_text(errors="replace")
    if text.count("PROGRAM ENDED") != 1:
        raise RuntimeError(f"incomplete CP2K output: {path}")
    energies = [
        float(value)
        for value in re.findall(
            rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
        )
    ]
    forces = [
        tuple(map(float, row))
        for row in re.findall(
            rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
            text,
            re.MULTILINE,
        )
    ]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    if not energies or not forces or not stress_blocks:
        raise RuntimeError(f"missing energy, force, or analytical stress: {path}")
    stress = [
        tuple(map(float, row))
        for row in re.findall(
            rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
            stress_blocks[-1],
            re.MULTILINE,
        )
    ]
    if len(stress) != 3:
        raise RuntimeError(f"malformed stress tensor: {path}")
    values = [energies[-1], *(v for row in forces for v in row), *(v for row in stress for v in row)]
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"non-finite observable: {path}")
    return energies[-1], forces, stress


def max_delta(left: list[tuple[float, float, float]], right: list[tuple[float, float, float]]) -> float:
    if len(left) != len(right):
        raise RuntimeError(f"observable block lengths differ: {len(left)} != {len(right)}")
    return max(abs(a - b) for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b))


def checked_run(run_dir: Path, job: dict, variant: str) -> dict:
    paths = [
        run_dir / name
        for name in ("run.json", "cp2k.out", "cp2k.err", "proc_affinity.tsv", "preexec_affinity.tsv")
    ]
    if not all(path.is_file() for path in paths):
        raise RuntimeError(f"missing evidence in {run_dir}")
    meta = json.loads(paths[0].read_text())
    if meta.get("returncode") != 0 or meta.get("variant") != variant:
        raise RuntimeError(f"bad result/variant metadata: {run_dir}")
    if (
        meta.get("case") != job["case"]
        or meta.get("ranks") != job["ranks"]
        or meta.get("batch_size") != job["batch"]
        or meta.get("expected_nfull") != job["nfull"]
    ):
        raise RuntimeError(f"case metadata mismatch: {run_dir}")
    if meta.get("input_sha256") != sha256(CAMPAIGN / "harness" / "inputs" / job["input"]):
        raise RuntimeError(f"input hash mismatch: {run_dir}")
    if meta.get("output_sha256") != sha256(paths[1]) or meta.get("stderr_sha256") != sha256(paths[2]):
        raise RuntimeError(f"raw output hash mismatch: {run_dir}")
    if meta.get("affinity_sha256") != sha256(paths[3]) or meta.get("bindings_sha256") != sha256(paths[2]):
        raise RuntimeError(f"affinity/binding hash mismatch: {run_dir}")
    if meta.get("preexec_affinity_sha256") != sha256(paths[4]):
        raise RuntimeError(f"pre-exec affinity hash mismatch: {run_dir}")
    cores = meta.get("cores", [])
    masks = meta.get("observed_singleton_masks", [])
    expected_masks = [str(core) for core in cores]
    if (
        len(cores) != job["ranks"]
        or meta.get("preexec_singleton_masks") != expected_masks
        or (masks and masks != sorted(expected_masks, key=int))
    ):
        raise RuntimeError(f"singleton affinity proof mismatch: {run_dir}")
    return meta


def internal_residuals(path: Path, job: dict) -> tuple[float, ...]:
    text = path.read_text(errors="replace")
    forwards = [tuple(map(float, match[1:])) for match in FORWARD_RE.findall(text)]
    reverses = [tuple(map(float, match)) for match in REVERSE_RE.findall(text)]
    modes = [tuple(map(int, match)) for match in MODE_RE.findall(text)]
    reverse_modes = [tuple(map(int, match)) for match in REVERSE_MODE_RE.findall(text)]
    if not forwards or not reverses or not modes or not reverse_modes:
        raise RuntimeError(f"missing distributed-image qualification marker: {path}")
    actual_batch = min(job["batch"], job["nfull"])
    expected_empty = max(0, job["ranks"] - job["nfull"])
    if not all(
        importers == job["ranks"]
        and empty == expected_empty
        and 1 <= nred <= job["nfull"]
        and nfull == job["nfull"]
        and batch == actual_batch
        for importers, empty, nred, nfull, batch in modes
    ):
        raise RuntimeError(f"forward mode metadata mismatch: {path}")
    if not all(
        nfull == job["nfull"]
        and importers == job["ranks"]
        and empty == expected_empty
        and batch == actual_batch
        and peak > 0
        for nfull, importers, empty, batch, peak in reverse_modes
    ):
        raise RuntimeError(f"reverse mode metadata mismatch: {path}")
    values = [value for row in forwards + reverses for value in row]
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise RuntimeError(f"invalid internal residual: {path}")
    maxima = tuple(max(row[index] for row in forwards) for index in range(3)) + tuple(
        max(row[index] for row in reverses) for index in range(3)
    )
    if max(maxima[:3]) > MATRIX["gates"]["internal_forward"]:
        raise RuntimeError(f"forward oracle gate failed: {path}: {maxima[:3]}")
    if max(maxima[3:]) > MATRIX["gates"]["internal_reverse"]:
        raise RuntimeError(f"reverse oracle gate failed: {path}: {maxima[3:]}")
    return maxima


def verify_faults(binary_hash: str) -> None:
    for fault in MATRIX["faults"]:
        run_dir = RUNS / f"fault_{fault['name']}"
        meta_path = run_dir / "run.json"
        output_path = run_dir / "cp2k.out"
        stderr_path = run_dir / "cp2k.err"
        affinity_path = run_dir / "proc_affinity.tsv"
        preexec_path = run_dir / "preexec_affinity.tsv"
        if not all(
            path.is_file()
            for path in (meta_path, output_path, stderr_path, affinity_path, preexec_path)
        ):
            raise RuntimeError(f"missing fault evidence: {run_dir}")
        meta = json.loads(meta_path.read_text())
        combined = output_path.read_text(errors="replace") + "\n" + stderr_path.read_text(errors="replace")
        if (
            meta.get("returncode") == 0
            or meta.get("injection") != fault["injection"]
            or meta.get("expected_diagnostic") != fault["diagnostic"]
            or meta.get("diagnostic_count", 0) < 1
            or fault["diagnostic"] not in combined
            or "PROGRAM ENDED" in combined
            or meta.get("cp2k_sha256") != binary_hash
        ):
            raise RuntimeError(f"fault gate mismatch: {fault['name']}")
        if (
            meta.get("output_sha256") != sha256(output_path)
            or meta.get("stderr_sha256") != sha256(stderr_path)
            or meta.get("affinity_sha256") != sha256(affinity_path)
            or meta.get("preexec_affinity_sha256") != sha256(preexec_path)
        ):
            raise RuntimeError(f"fault raw hash mismatch: {fault['name']}")
        cores = meta.get("cores", [])
        expected_masks = [str(core) for core in cores]
        observed_masks = meta.get("observed_singleton_masks", [])
        if (
            meta.get("preexec_singleton_masks") != expected_masks
            or (observed_masks and observed_masks != sorted(expected_masks, key=int))
        ):
            raise RuntimeError(f"fault affinity mismatch: {fault['name']}")


def verify_concurrent_core_ownership() -> dict:
    """Reject any two wall-time-overlapping runs that claim the same CPU core."""
    intervals = []
    for meta_path in sorted(RUNS.glob("*/run.json")):
        meta = json.loads(meta_path.read_text())
        if "started_unix" not in meta or "finished_unix" not in meta:
            raise RuntimeError(f"missing run interval: {meta_path}")
        started = float(meta["started_unix"])
        finished = float(meta["finished_unix"])
        cores = tuple(int(core) for core in meta.get("cores", []))
        if (
            not math.isfinite(started)
            or not math.isfinite(finished)
            or finished < started
            or not cores
            or len(cores) != len(set(cores))
        ):
            raise RuntimeError(f"invalid run interval/core ownership: {meta_path}")
        intervals.append(
            {
                "run": meta_path.parent.name,
                "started_unix": started,
                "finished_unix": finished,
                "cores": list(cores),
            }
        )
    expected = 2 * len(MATRIX["jobs"]) + len(MATRIX["faults"])
    if len(intervals) != expected:
        raise RuntimeError(f"affinity interval count mismatch: {len(intervals)} != {expected}")
    conflicts = []
    for index, left in enumerate(intervals):
        for right in intervals[index + 1 :]:
            wall_overlap = max(left["started_unix"], right["started_unix"]) < min(
                left["finished_unix"], right["finished_unix"]
            )
            shared = sorted(set(left["cores"]) & set(right["cores"]))
            if wall_overlap and shared:
                conflicts.append(
                    {"left": left["run"], "right": right["run"], "shared_cores": shared}
                )
    audit = {
        "runs": len(intervals),
        "allowed_core_pool": sorted({core for row in intervals for core in row["cores"]}),
        "concurrent_core_conflicts": conflicts,
        "status": "PASS" if not conflicts else "FAIL",
    }
    (CAMPAIGN / "affinity_concurrency_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    if conflicts:
        raise RuntimeError(f"concurrent CPU-core ownership conflicts: {conflicts}")
    return audit


def main() -> None:
    rows = []
    binary_hashes = set()
    for job in MATRIX["jobs"]:
        stem = f"{job['case']}_p{job['ranks']}_b{job['batch']}"
        dense_dir = RUNS / f"{stem}_dense"
        partial_dir = RUNS / f"{stem}_distributed_images"
        dense_meta = checked_run(dense_dir, job, "DENSE")
        partial_meta = checked_run(partial_dir, job, "DISTRIBUTED_IMAGES")
        if dense_meta["cp2k_sha256"] != partial_meta["cp2k_sha256"]:
            raise RuntimeError(f"binary hash differs inside pair: {stem}")
        binary_hashes.add(dense_meta["cp2k_sha256"])
        dense = parse_observables(dense_dir / "cp2k.out")
        partial = parse_observables(partial_dir / "cp2k.out")
        d_energy = abs(dense[0] - partial[0])
        d_force = max_delta(dense[1], partial[1])
        d_stress = max_delta(dense[2], partial[2])
        internal = internal_residuals(partial_dir / "cp2k.out", job)
        gates = MATRIX["gates"]
        if d_energy > gates["external_energy_Ha"]:
            raise RuntimeError(f"external energy gate failed: {stem}: {d_energy}")
        if d_force > gates["external_force_Ha_per_bohr"]:
            raise RuntimeError(f"external force gate failed: {stem}: {d_force}")
        if d_stress > gates["external_stress_bar"]:
            raise RuntimeError(f"external stress gate failed: {stem}: {d_stress}")
        dense_wall = float(dense_meta["wall_seconds"])
        partial_wall = float(partial_meta["wall_seconds"])
        dense_rss = int(dense_meta["peak_sampled_cp2k_rank_rss_kb"])
        partial_rss = int(partial_meta["peak_sampled_cp2k_rank_rss_kb"])
        rows.append(
            {
                "case": job["case"],
                "ranks": job["ranks"],
                "batch": job["batch"],
                "features": ",".join(job["features"]),
                "nfull": job["nfull"],
                "dense_wall_s": f"{dense_wall:.9f}",
                "distributed_wall_s": f"{partial_wall:.9f}",
                "wall_ratio_distributed_over_dense": f"{partial_wall / dense_wall:.9f}",
                "dense_peak_rank_rss_kb": dense_rss,
                "distributed_peak_rank_rss_kb": partial_rss,
                "rss_ratio_distributed_over_dense": (
                    f"{partial_rss / dense_rss:.9f}" if dense_rss > 0 and partial_rss > 0 else "NA"
                ),
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
            }
        )
    if len(binary_hashes) != 1:
        raise RuntimeError(f"campaign used multiple CP2K binaries: {binary_hashes}")
    verify_faults(next(iter(binary_hashes)))
    affinity_audit = verify_concurrent_core_ownership()
    summary = CAMPAIGN / "summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    maxima = {
        "pairs": len(rows),
        "faults": len(MATRIX["faults"]),
        "max_external_dE_Ha": max(float(row["external_dE_Ha"]) for row in rows),
        "max_external_dForce_Ha_per_bohr": max(float(row["external_dForce_Ha_per_bohr"]) for row in rows),
        "max_external_dStress_bar": max(float(row["external_dStress_bar"]) for row in rows),
        "max_internal_forward": max(
            max(float(row[key]) for key in ("internal_dE_Ha", "internal_dVsh_Ha", "internal_dFfold_Ha"))
            for row in rows
        ),
        "max_internal_reverse": max(
            max(float(row[key]) for key in ("internal_dOverlap", "internal_dForce_Ha_per_bohr", "internal_dStress_Ha"))
            for row in rows
        ),
        "affinity_runs": affinity_audit["runs"],
        "concurrent_core_conflicts": len(affinity_audit["concurrent_core_conflicts"]),
        "rss_pairs_sampled": sum(
            row["dense_peak_rank_rss_kb"] > 0 and row["distributed_peak_rank_rss_kb"] > 0
            for row in rows
        ),
        "cp2k_sha256": next(iter(binary_hashes)),
    }
    (CAMPAIGN / "maxima.json").write_text(json.dumps(maxima, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {len(rows)} DENSE/distributed-image pairs and {len(MATRIX['faults'])} faults")
    print(json.dumps(maxima, sort_keys=True))


if __name__ == "__main__":
    main()
