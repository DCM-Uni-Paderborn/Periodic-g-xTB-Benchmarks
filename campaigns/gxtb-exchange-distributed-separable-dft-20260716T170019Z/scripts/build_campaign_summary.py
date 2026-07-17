#!/usr/bin/env python3
"""Build and validate the immutable exchange-acceleration campaign index."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "raw" / "distributed_cp2k"
PROVIDER = ROOT / "raw" / "provider"
CONSUMER = ROOT / "raw" / "cp2k_separable_consumer"
DERIVED = ROOT / "derived"
MANIFEST = ROOT / "SHA256SUMS"

CONSUMER_PRODUCTION_PAIRS = {
    "ch4_full_debug": (
        "runs/ch4_full_debug_dense/cp2k.out",
        "runs/ch4_full_debug_separable/cp2k.out",
    ),
    "ch4_spglib": (
        "runs/ch4_spglib_dense_printed/cp2k.out",
        "runs/ch4_spglib_separable_printed/cp2k.out",
    ),
    "ch4_k290": (
        "runs/ch4_k290_dense_printed/cp2k.out",
        "runs/ch4_k290_separable_printed/cp2k.out",
    ),
    "ar2_1d": (
        "runs/ar2_1d_dense/cp2k.out",
        "runs/ar2_1d_separable/cp2k.out",
    ),
    "ar4_2d": (
        "runs/ar4_2d_dense/cp2k.out",
        "runs/ar4_2d_separable/cp2k.out",
    ),
    "si_shifted": (
        "runs/si_shifted_dense/cp2k.out",
        "runs/si_shifted_separable/cp2k.out",
    ),
    "o2_uks_time_reversal": (
        "runs/o2_uks_dense/cp2k.out",
        "runs/o2_uks_separable/cp2k.out",
    ),
    "ch4_distributed_smoke": (
        "runs/ch4_full_dense_nooracle/cp2k.out",
        "runs/ch4_full_distributed_p2/cp2k.out",
    ),
    "ch4_distributed_plus_separable": (
        "runs/ch4_full_distributed_p2/cp2k.out",
        "runs/ch4_full_distributed_separable_p2/cp2k.out",
    ),
}

CONSUMER_IDENTITY_DIAGNOSTICS = {
    "runs/ch4_full_dense/cp2k.out",
    "runs/ch4_full_separable/cp2k.out",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(value: str) -> float | None:
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_source_manifest(root: Path) -> dict[str, object]:
    failures: list[dict[str, str]] = []
    checked = 0
    with (root / "SHA256SUMS").open(encoding="utf-8") as handle:
        for line in handle:
            expected, relpath = line.rstrip("\n").split(maxsplit=1)
            path = root / relpath.removeprefix("./")
            checked += 1
            if not path.is_file():
                failures.append({"path": relpath, "error": "missing"})
            else:
                actual = sha256(path)
                if actual != expected:
                    failures.append(
                        {"path": relpath, "expected": expected, "actual": actual}
                    )
    return {"checked": checked, "passed": not failures, "failures": failures}


def output_inventory(root: Path) -> dict[str, object]:
    outputs = sorted(root.rglob("*.out"))
    ended = []
    diagnostic = []
    for path in outputs:
        record = {
            "path": path.relative_to(ROOT).as_posix(),
            "sha256": sha256(path),
        }
        if b"PROGRAM ENDED" in path.read_bytes():
            ended.append(record)
        else:
            diagnostic.append(record)
    return {
        "total_outputs": len(outputs),
        "normally_ended": len(ended),
        "preserved_diagnostic_or_incomplete": len(diagnostic),
        "diagnostic_or_incomplete_outputs": diagnostic,
    }


def aggregate_maxima(rows: list[dict[str, str]]) -> dict[str, float | None]:
    keys = [
        "max_abs_energy_Eh",
        "max_abs_force_Eh_per_bohr",
        "max_abs_analytical_stress_bar",
        "oracle_max_dE_Eh",
        "oracle_max_dVsh_Eh",
        "oracle_max_dFfold_Eh",
        "oracle_max_hermiticity",
        "oracle_max_cov_full",
        "oracle_max_cov_fold",
        "oracle_max_duality_fold",
    ]
    maxima: dict[str, float | None] = {}
    for key in keys:
        values = [value for row in rows if (value := numeric(row[key])) is not None]
        maxima[key] = max(values) if values else None
    return maxima


def aggregate_consumer_maxima(rows: list[dict[str, str]]) -> dict[str, float | None]:
    keys = [
        "max_abs_energy_delta_Ha",
        "max_abs_force_delta_Ha_per_bohr",
        "max_abs_stress_delta_bar",
    ]
    return {
        key: max(float(row[key]) for row in rows)
        for key in keys
    }


def build_summary() -> dict[str, object]:
    distributed_rows = read_tsv(DIST / "qualification_summary.tsv")
    transform_rows = read_tsv(PROVIDER / "separable_dft_benchmark_20260716.tsv")
    consumer_rows = read_tsv(CONSUMER / "qualification_summary.tsv")
    release_timing_rows = read_tsv(CONSUMER / "release_bench" / "timing_summary.tsv")
    distributed_manifest = verify_source_manifest(DIST)
    consumer_manifest = verify_source_manifest(CONSUMER)

    if not distributed_manifest["passed"]:
        raise SystemExit("copied distributed evidence failed its source manifest")
    if not consumer_manifest["passed"]:
        raise SystemExit("copied CP2K consumer evidence failed its source manifest")

    listed_outputs = []
    for row in distributed_rows:
        for field in ("legacy_output", "distributed_output"):
            path = DIST / row[field]
            if not path.is_file():
                raise SystemExit(f"missing output listed in qualification_summary.tsv: {path}")
            if b"PROGRAM ENDED" not in path.read_bytes():
                raise SystemExit(f"listed production output did not end normally: {path}")
            listed_outputs.append(path.relative_to(ROOT).as_posix())

    consumer_cases = {row["case"] for row in consumer_rows}
    if consumer_cases != set(CONSUMER_PRODUCTION_PAIRS):
        raise SystemExit(
            "consumer qualification cases do not match the explicit production-pair map"
        )

    consumer_listed_outputs = []
    for row in consumer_rows:
        if not row["status"].startswith("pass"):
            raise SystemExit(f"non-pass row in consumer qualification table: {row['case']}")
        for relpath in CONSUMER_PRODUCTION_PAIRS[row["case"]]:
            path = CONSUMER / relpath
            if not path.is_file():
                raise SystemExit(f"missing CP2K consumer production output: {path}")
            if b"PROGRAM ENDED" not in path.read_bytes():
                raise SystemExit(
                    f"tabulated CP2K consumer output did not end normally: {path}"
                )
            consumer_listed_outputs.append(path.relative_to(ROOT).as_posix())

    consumer_inventory = output_inventory(CONSUMER)
    observed_consumer_diagnostics = {
        Path(record["path"]).relative_to(
            Path("raw") / "cp2k_separable_consumer"
        ).as_posix()
        for record in consumer_inventory["diagnostic_or_incomplete_outputs"]
    }
    if observed_consumer_diagnostics != CONSUMER_IDENTITY_DIAGNOSTICS:
        raise SystemExit(
            "unexpected CP2K consumer output without PROGRAM ENDED: "
            f"{sorted(observed_consumer_diagnostics)}"
        )

    expected_release_runs = {
        (backend, str(repeat))
        for backend in ("dense", "separable_dft")
        for repeat in (1, 2, 3)
    }
    observed_release_runs = {
        (row["backend"], row["repeat"]) for row in release_timing_rows
    }
    if observed_release_runs != expected_release_runs:
        raise SystemExit("release timing table does not contain the expected 3+3 runs")

    release_outputs = []
    for row in release_timing_rows:
        if row["program_ended"].lower() != "true":
            raise SystemExit(f"release timing row is not a pass: {row}")
        backend_label = "dense" if row["backend"] == "dense" else "separable"
        path = (
            CONSUMER
            / "release_bench"
            / f"ch4_k666_{backend_label}_r{row['repeat']}"
            / "cp2k.out"
        )
        if not path.is_file() or b"PROGRAM ENDED" not in path.read_bytes():
            raise SystemExit(f"release benchmark output did not end normally: {path}")
        release_outputs.append(path.relative_to(ROOT).as_posix())

    def release_median(backend: str, field: str) -> float:
        return statistics.median(
            float(row[field])
            for row in release_timing_rows
            if row["backend"] == backend
        )

    release_medians = {
        backend: {
            field: release_median(backend, field)
            for field in ("real_s", "exclusive_s", "inclusive_s")
        }
        for backend in ("dense", "separable_dft")
    }
    release_speedups = {
        field: release_medians["dense"][field]
        / release_medians["separable_dft"][field]
        for field in ("real_s", "exclusive_s", "inclusive_s")
    }

    return {
        "schema_version": 1,
        "campaign_id": ROOT.name,
        "archive_created_utc": "2026-07-16T17:00:19Z",
        "archive_extended_with_cp2k_consumer_utc": "2026-07-16T18:27:48Z",
        "scope": "periodic g-xTB Brillouin-zone-coupled nonlocal exchange",
        "qualification_boundary": {
            "distributed_image_ranges": "qualified_end_to_end_forward",
            "separable_direct_dft_provider": "qualified_provider",
            "separable_direct_dft_cp2k_consumer": "qualified_end_to_end",
            "bounded_derivative_distribution": "not_implemented",
            "cross_mesh_restart": "not_covered_by_this_archive",
        },
        "methods": {
            "distributed_image_ranges": {
                "algorithm": "contiguous_disjoint_BvK_image_range_contraction",
                "parallelism": "caller-owned MPI rank partition with exact-once coverage",
                "reduction_targets": [
                    "exchange_energy",
                    "shell_potential",
                    "folded_Fock_response",
                ],
                "reference": "legacy explicit full-mesh CP2K path plus per-iteration full-mesh oracle",
                "derivative_boundary": (
                    "Analytical forces and stress were compared and agree, but their reverse/gradient "
                    "path still evaluates the complete mesh and is not distributed or batch-bounded."
                ),
            },
            "separable_direct_dft": {
                "algorithm": "factorized direct DFT over regular-grid pencils",
                "is_fft": False,
                "arithmetic_complexity": "O(Nrow * Nk * sum(nmesh))",
                "reference": "selectable dense phase-table transform",
                "dense_plan_complex_elements": "2 * Nk^2",
                "compact_plan_complex_elements": "Nk + 3 * max(nmesh)",
                "compact_plan_integer_elements": "2 * Nk",
                "per_call_scratch_complex_elements": "2 * Nrow * Nk",
                "applicability": [
                    "machine-regular complete uniform product mesh",
                    "common grid shift",
                    "arbitrary external k and image order",
                    "1D, 2D, and 3D meshes",
                ],
            },
        },
        "distributed_qualification": {
            "case_count": len(distributed_rows),
            "normally_ended_listed_outputs": len(listed_outputs),
            "cases": distributed_rows,
            "aggregate_maxima": aggregate_maxima(distributed_rows),
            "copied_source_manifest": distributed_manifest,
            "all_output_inventory": output_inventory(DIST),
        },
        "separable_direct_dft_qualification": {
            "debug_exchange": {"passed": 31, "total": 31},
            "debug_gxtb": {
                "ordinary_passes": 40,
                "expected_failure_diagnostics": 4,
                "total": 44,
            },
            "release_exchange": {"passed": 31, "total": 31},
            "focused_dense_oracle_absolute_tolerance": 1.0e-11,
            "benchmark_protocol": {
                "host": "local Apple-silicon host",
                "build": "release -O3",
                "rows": 64,
                "paired_calls_per_sample": "k-to-R and R-to-k",
                "repetitions": 20,
            },
            "transform_benchmarks": transform_rows,
        },
        "cp2k_separable_direct_dft_consumer_qualification": {
            "case_count": len(consumer_rows),
            "normally_ended_tabulated_outputs": len(consumer_listed_outputs),
            "distinct_normally_ended_tabulated_outputs": len(
                set(consumer_listed_outputs)
            ),
            "production_pairs": {
                case: list(paths)
                for case, paths in sorted(CONSUMER_PRODUCTION_PAIRS.items())
            },
            "cases": consumer_rows,
            "aggregate_maxima": aggregate_consumer_maxima(consumer_rows),
            "copied_source_manifest": consumer_manifest,
            "all_output_inventory": consumer_inventory,
            "preserved_nonpass_identity_diagnostics": [
                {
                    "path": "runs/ch4_full_dense/cp2k.out",
                    "classification": "qualification-only identity-gate diagnostic",
                    "termination": "intentional ABORT; not a production pass",
                },
                {
                    "path": "runs/ch4_full_separable/cp2k.out",
                    "classification": "qualification-only identity-duality diagnostic",
                    "termination": "intentional ABORT; not a production pass",
                },
            ],
            "derivative_boundary": (
                "Energy, force, and analytical-stress equivalence is qualified, including "
                "a 23-energy DEBUG trajectory. The derivative contraction itself still uses "
                "the complete mesh and has no distributed or batch-bounded scaling claim."
            ),
            "release_k666_timing": {
                "input": "inputs/ch4_full_k666_energy.inp",
                "repeat_count_per_backend": 3,
                "normally_ended_outputs": release_outputs,
                "raw_rows": release_timing_rows,
                "medians_s": release_medians,
                "dense_over_separable_speedups": release_speedups,
                "interpretation": (
                    "The exclusive build_tblite_ks_matrix kernel is 1.20x faster, "
                    "but the inclusive timer changes by only 1.03x and the approximately "
                    "7.18 s total wall time is indistinguishable. This does not establish "
                    "an end-to-end SCF speedup."
                ),
            },
        },
        "archive_notes": [
            "All copied source inputs and outputs are preserved byte-for-byte.",
            "Seven extra CP2K outputs without PROGRAM ENDED are retained as diagnostic/incomplete raw evidence and are not qualification rows.",
            "Two CP2K transform-consumer identity-diagnostic outputs terminate through qualification-only ABORT gates; they are preserved but are not production passes or tabulated comparison pairs.",
            "Provider transform timing is transform-only and must not be interpreted as whole-SCF speedup.",
            "The six-run CP2K release probe resolves a 1.20x exclusive kernel speedup but no measurable end-to-end speedup.",
            "Serial image batching is a bounded-memory/recomputation trade-off, not a serial acceleration claim.",
            "No paper or supporting-information file, PDF, or paper hash is contained in this campaign.",
        ],
    }


def write_manifest() -> None:
    files = [
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and path != MANIFEST
        and "__pycache__" not in path.parts
        and path.suffix != ".pyc"
    ]
    lines = [f"{sha256(path)}  ./{path.relative_to(ROOT).as_posix()}\n" for path in sorted(files)]
    MANIFEST.write_text("".join(lines), encoding="utf-8")


def validate_manifest() -> None:
    with MANIFEST.open(encoding="utf-8") as handle:
        for line in handle:
            expected, relpath = line.rstrip("\n").split(maxsplit=1)
            path = ROOT / relpath.removeprefix("./")
            if not path.is_file() or sha256(path) != expected:
                raise SystemExit(f"campaign manifest mismatch: {relpath}")


def main() -> None:
    DERIVED.mkdir(exist_ok=True)
    summary = build_summary()
    (DERIVED / "qualification_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_manifest()
    validate_manifest()
    print(
        "validated: "
        f"{summary['distributed_qualification']['case_count']} distributed cases, "
        f"{len(summary['separable_direct_dft_qualification']['transform_benchmarks'])} transform benchmarks, "
        f"{summary['cp2k_separable_direct_dft_consumer_qualification']['case_count']} CP2K consumer cases"
    )


if __name__ == "__main__":
    main()
