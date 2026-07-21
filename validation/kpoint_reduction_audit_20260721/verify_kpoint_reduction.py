#!/usr/bin/env python3
"""Verify numerical identity and measured K290/SPGLIB reduction behavior."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPECTED_BINARY_SHA256 = (
    "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
)
EXPECTED_ROUNDS = (1, 2, 3)
ENERGY_TOLERANCE_HARTREE = 1.0e-10
CASES = {
    "ch4_k666_full_complete": {
        "kpoints": 216,
        "backend": "FULL",
        "fused": False,
    },
    "ch4_k666_k290_complete": {
        "kpoints": 10,
        "backend": "K290",
        "fused": False,
    },
    "ch4_k666_spglib_complete": {
        "kpoints": 10,
        "backend": "SPGLIB",
        "fused": False,
    },
    "ch4_k666_k290_fused": {
        "kpoints": 10,
        "backend": "K290",
        "fused": True,
    },
    "ch4_k666_spglib_fused": {
        "kpoints": 10,
        "backend": "SPGLIB",
        "fused": True,
    },
}
THREAD_VARIABLES = {
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_hash(path: Path) -> str:
    return path.read_text().split()[0]


def parse_profile(path: Path) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in path.read_text().splitlines():
        key, value = line.split("=", 1)
        if key == "elapsed_seconds":
            values[key] = float(value)
        else:
            values[key] = int(value)
    return values


def parse_affinity(path: Path) -> dict[str, object]:
    text = path.read_text()
    header = re.search(r"expected_cpu=(\d+) allowed=([^\s]+)", text)
    if not header:
        raise ValueError(f"Missing affinity header in {path}")
    threads = {}
    for name, value in re.findall(r"^([A-Z_]+)=(\d+)\s*$", text, re.MULTILINE):
        if name in THREAD_VARIABLES:
            threads[name] = int(value)
    return {
        "expected_cpu": int(header.group(1)),
        "allowed_cpu_list": header.group(2),
        "threads": threads,
        "singleton": header.group(2) == header.group(1),
        "all_threads_one": threads.keys() == THREAD_VARIABLES
        and all(value == 1 for value in threads.values()),
    }


def parse_output(path: Path, expected: dict[str, object]) -> dict[str, object]:
    text = path.read_text(errors="replace")
    energy_matches = re.findall(
        r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+([-+0-9.Ee]+)",
        text,
    )
    kpoint_matches = re.findall(
        r"BRILLOUIN\| List of Kpoints \[2 Pi/Bohr\]\s+(\d+)", text
    )
    scf_matches = re.findall(r"SCF run converged in\s+(\d+) steps", text)
    if not energy_matches or not kpoint_matches or not scf_matches:
        raise ValueError(f"Incomplete CP2K output: {path}")

    fused_matches = re.findall(
        r"GXTB-SYMMETRY-FUSED nred=(\d+), nfull=(\d+), "
        r"requestedBatch=(\d+), effectiveBatch=(\d+), "
        r"cachedRepresentativeOverlaps=(\d+), "
        r"persistentFullExchangeStarMatrices=(\d+)",
        text,
    )
    fused_records = [
        {
            "nred": int(nred),
            "nfull": int(nfull),
            "requested_batch": int(requested),
            "effective_batch": int(effective),
            "cached_representative_overlaps": int(cached),
            "persistent_full_exchange_star_matrices": int(persistent),
        }
        for nred, nfull, requested, effective, cached, persistent in fused_matches
    ]
    fused_marker_pass = all(
        record
        == {
            "nred": 10,
            "nfull": 216,
            "requested_batch": 8,
            "effective_batch": 8,
            "cached_representative_overlaps": 10,
            "persistent_full_exchange_star_matrices": 0,
        }
        for record in fused_records
    ) and bool(fused_records)

    backend = expected["backend"]
    if backend == "FULL":
        backend_pass = "K-Point point group symmetrization" in text and bool(
            re.search(r"K-Point point group symmetrization\s+OFF", text)
        )
    else:
        backend_pass = bool(
            re.search(rf"BRILLOUIN\| Symmetry backend\s+{backend}", text)
            and re.search(
                rf"BRILLOUIN\| Symmetry reduction method\s+{backend}", text
            )
        )

    return {
        "energy_hartree": float(energy_matches[-1]),
        "kpoints": int(kpoint_matches[-1]),
        "scf_steps": int(scf_matches[-1]),
        "normal_termination": "PROGRAM ENDED AT" in text,
        "backend_pass": backend_pass,
        "fused_records": fused_records,
        "fused_marker_pass": fused_marker_pass
        if expected["fused"]
        else not fused_records,
        "output_sha256": sha256(path),
    }


def parse_case(round_number: int, name: str, expected: dict[str, object]) -> dict[str, object]:
    result = ROOT / "results" / f"round{round_number}" / name
    input_path = ROOT / "inputs" / f"{name}.inp"
    output_path = result / "cp2k.out"
    profile = parse_profile(result / "time_verbose.txt")
    affinity = parse_affinity(result / "affinity_preexec.txt")
    parsed = parse_output(output_path, expected)
    exit_status = int((result / "exit_status").read_text().strip())
    recorded_binary = first_hash(result / "binary.sha256")
    recorded_input = first_hash(result / "input.sha256")
    checks = {
        "exit_status_zero": exit_status == 0,
        "normal_termination": parsed["normal_termination"],
        "qualified_binary": recorded_binary == EXPECTED_BINARY_SHA256,
        "input_hash_matches": recorded_input == sha256(input_path),
        "singleton_affinity": affinity["singleton"],
        "all_threads_one": affinity["all_threads_one"],
        "expected_backend": parsed["backend_pass"],
        "expected_kpoint_count": parsed["kpoints"] == expected["kpoints"],
        "expected_fused_storage_marker": parsed["fused_marker_pass"],
        "profile_exit_status_zero": profile["calculation_exit_status"] == 0,
    }
    return {
        "round": round_number,
        "case": name,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "binary_sha256": recorded_binary,
        "input_sha256": recorded_input,
        "affinity": affinity,
        "profile": profile,
        **parsed,
    }


def median(values: list[float | int]) -> float:
    return float(statistics.median(values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=ROOT / "verification.json"
    )
    args = parser.parse_args()

    records = []
    for round_number in EXPECTED_ROUNDS:
        for name, expected in CASES.items():
            records.append(parse_case(round_number, name, expected))

    reference_by_round = {
        record["round"]: record
        for record in records
        if record["case"] == "ch4_k666_full_complete"
    }
    for record in records:
        reference = reference_by_round[record["round"]]
        record["energy_minus_round_full_hartree"] = (
            record["energy_hartree"] - reference["energy_hartree"]
        )
        record["within_energy_tolerance"] = (
            abs(record["energy_minus_round_full_hartree"])
            <= ENERGY_TOLERANCE_HARTREE
        )

    summaries = {}
    for name in CASES:
        selected = [record for record in records if record["case"] == name]
        summaries[name] = {
            "median_elapsed_seconds": median(
                [record["profile"]["elapsed_seconds"] for record in selected]
            ),
            "median_peak_rss_kib": median(
                [
                    max(
                        record["profile"]["peak_sampled_rss_kib"],
                        record["profile"]["peak_sampled_hwm_kib"],
                    )
                    for record in selected
                ]
            ),
            "energies_hartree": [record["energy_hartree"] for record in selected],
            "scf_steps": [record["scf_steps"] for record in selected],
            "kpoints": selected[0]["kpoints"],
        }

    full_time = summaries["ch4_k666_full_complete"]["median_elapsed_seconds"]
    full_rss = summaries["ch4_k666_full_complete"]["median_peak_rss_kib"]
    for name, summary in summaries.items():
        summary["elapsed_change_vs_full_percent"] = 100.0 * (
            summary["median_elapsed_seconds"] - full_time
        ) / full_time
        summary["peak_rss_change_vs_full_percent"] = 100.0 * (
            summary["median_peak_rss_kib"] - full_rss
        ) / full_rss

    checks = {
        "all_provenance_and_execution_checks_pass": all(
            record["all_checks_pass"] for record in records
        ),
        "all_energies_match_full": all(
            record["within_energy_tolerance"] for record in records
        ),
        "all_scf_step_counts_match": len(
            {record["scf_steps"] for record in records}
        )
        == 1,
        "full_mesh_has_216_points": all(
            record["kpoints"] == 216
            for record in records
            if record["case"] == "ch4_k666_full_complete"
        ),
        "reduced_meshes_have_10_points": all(
            record["kpoints"] == 10
            for record in records
            if record["case"] != "ch4_k666_full_complete"
        ),
        "fused_paths_prove_zero_persistent_full_exchange_star_matrices": all(
            record["fused_marker_pass"]
            for record in records
            if CASES[record["case"]]["fused"]
        ),
    }
    kpoint_reduction_percent = 100.0 * (216 - 10) / 216
    status = "PASS" if all(checks.values()) else "FAIL"
    result = {
        "schema": "periodic-gxtb-kpoint-reduction-audit-v1",
        "status": status,
        "qualified_cp2k_binary_sha256": EXPECTED_BINARY_SHA256,
        "mesh": [6, 6, 6],
        "full_kpoints": 216,
        "irreducible_kpoints": 10,
        "kpoint_count_reduction_percent": kpoint_reduction_percent,
        "energy_tolerance_hartree": ENERGY_TOLERANCE_HARTREE,
        "checks": checks,
        "summaries": summaries,
        "records": records,
        "interpretation": (
            "K290 and SPGLIB reduce the SCF/diagonalization set from 216 to 10 "
            "points. With the default COMPLETE_ARRAY exchange backend, however, "
            "the coupled 216-point exchange mesh is still reconstructed and "
            "materialized, so the measured whole-process RSS and wall-time changes "
            "must be reported rather than inferred from the 95.37% point-count "
            "reduction. SYMMETRY_FUSED proves bounded batch storage and zero "
            "persistent full exchange-star matrices, but its end-to-end RSS and "
            "timing remain empirical properties of this small single-core case."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
