#!/usr/bin/env python3
"""Create deterministic summaries and checksums for this runtime archive."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
ORACLE_KEYS = (
    "dE",
    "dVsh",
    "dFfold",
    "hermFull",
    "covOracleFull",
    "covOracleFold",
    "dualOracleFold",
    "covStream",
    "dualStream",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_matrix_after(lines: list[str], marker: str) -> list[float] | None:
    for index, line in enumerate(lines):
        if marker not in line:
            continue
        values: list[float] = []
        for candidate in lines[index + 1 : index + 6]:
            if not candidate.startswith(" DEBUG|"):
                if values:
                    break
                continue
            numbers = re.findall(FLOAT, candidate.split("|", 1)[1])
            if len(numbers) != 3:
                if values:
                    break
                continue
            values.extend(float(number) for number in numbers)
            if len(values) == 9:
                return values
    return None


def force_summary(lines: list[str]) -> list[dict[str, object]]:
    inside = False
    rows: list[dict[str, object]] = []
    pattern = re.compile(
        rf"^ DEBUG\|\s+(\d+)\s+([xyz])\s+({FLOAT})\s+({FLOAT})\s+"
        rf"({FLOAT})\s+({FLOAT})"
    )
    for line in lines:
        if "BEGIN OF SUMMARY" in line:
            inside = True
            continue
        if "END OF SUMMARY" in line:
            break
        if not inside:
            continue
        match = pattern.match(line)
        if match:
            rows.append(
                {
                    "atom": int(match.group(1)),
                    "coordinate": match.group(2),
                    "numerical_au": float(match.group(3)),
                    "analytical_au": float(match.group(4)),
                    "difference_au": float(match.group(5)),
                    "error_percent": float(match.group(6)),
                }
            )
    return rows


def parse_output(path: Path) -> dict[str, object]:
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    version_match = re.search(r"CP2K\| version string:\s*(.+)", text)
    run_type_match = re.search(r"GLOBAL\| Run type\s+(\S+)", text)
    mpi_match = re.search(r"DBCSR\| MPI: Number of processes\s+(\d+)", text)
    thread_match = re.search(
        r"GLOBAL\| Number of threads for this process\s+(\d+)", text
    )
    energy_text = re.findall(
        rf"ENERGY\| Total FORCE_EVAL .*?energy \[hartree\]\s+({FLOAT})", text
    )
    energies = [float(value) for value in energy_text]
    timing = None
    for line in lines:
        fields = line.split()
        if len(fields) >= 7 and fields[0] == "CP2K" and fields[1] == "1":
            try:
                timing = float(fields[-1])
            except ValueError:
                pass
    oracle_max = {key: 0.0 for key in ORACLE_KEYS}
    oracle_lines = 0
    for line in lines:
        if "GXTB-QUALIFICATION_ONLY STREAM-ORACLE" not in line:
            continue
        oracle_lines += 1
        for key, value in re.findall(rf"(\w+)=\s*({FLOAT})", line):
            if key in oracle_max:
                oracle_max[key] = max(oracle_max[key], abs(float(value)))
    virial_sum_match = re.search(
        rf"DEBUG\| Sum of differences\s+({FLOAT})", text
    )
    periodic_sum_match = re.search(
        rf"DEBUG\| Periodic-subspace sum of differences\s+({FLOAT})", text
    )
    force_sum_match = re.search(
        rf"DEBUG\| Sum of differences:\s+({FLOAT})", text
    )
    return {
        "program_ended": "PROGRAM ENDED AT" in text,
        "cp2k_version": version_match.group(1).strip() if version_match else None,
        "run_type": run_type_match.group(1) if run_type_match else None,
        "mpi_processes": int(mpi_match.group(1)) if mpi_match else None,
        "omp_threads": int(thread_match.group(1)) if thread_match else None,
        "energy_evaluation_count": len(energies),
        "first_energy_hartree_text": energy_text[0] if energy_text else None,
        "final_energy_hartree_text": energy_text[-1] if energy_text else None,
        "first_energy_hartree": energies[0] if energies else None,
        "final_energy_hartree": energies[-1] if energies else None,
        "cp2k_total_time_seconds": timing,
        "oracle_line_count": oracle_lines,
        "oracle_max_abs": oracle_max if oracle_lines else None,
        "analytical_pv_virial_au": first_matrix_after(
            lines, "Analytical pv_virial [a.u.]"
        ),
        "debug_virial_sum_abs_au": (
            abs(float(virial_sum_match.group(1))) if virial_sum_match else None
        ),
        "debug_periodic_subspace_virial_sum_abs_au": (
            abs(float(periodic_sum_match.group(1)))
            if periodic_sum_match
            else None
        ),
        "debug_force_sum_abs_au": (
            abs(float(force_sum_match.group(1))) if force_sum_match else None
        ),
        "force_summary": force_summary(lines),
    }


def max_delta(left: list[float] | None, right: list[float] | None) -> float | None:
    if left is None or right is None or not left or len(left) != len(right):
        return None
    return max(abs(a - b) for a, b in zip(left, right, strict=True))


def analytical_forces(row: dict[str, object]) -> list[float]:
    return [entry["analytical_au"] for entry in row["force_summary"]]


def write_sha256sums() -> None:
    manifest = ROOT / "SHA256SUMS"
    files = sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and path != manifest
        and "__pycache__" not in path.parts
    )
    manifest.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(ROOT)}\n" for path in files)
    )


def main() -> None:
    with (ROOT / "run_manifest.tsv").open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

    runs: list[dict[str, object]] = []
    for manifest_row in manifest_rows:
        input_path = ROOT / manifest_row["input_file"]
        output_path = ROOT / manifest_row["output_file"]
        parsed = parse_output(output_path)
        row: dict[str, object] = dict(manifest_row)
        row["batch_size"] = (
            int(manifest_row["batch_size"]) if manifest_row["batch_size"] else None
        )
        row["input_sha256"] = sha256(input_path)
        row["output_sha256"] = sha256(output_path)
        row.update(parsed)
        runs.append(row)

    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in runs:
        groups[str(row["comparison_group"])].append(row)

    comparisons: list[dict[str, object]] = []
    for group, rows in groups.items():
        references = [row for row in rows if row["mode"] == "legacy"]
        reference = references[0] if references else rows[0]
        for row in rows:
            comparison: dict[str, object] = {
                "comparison_group": group,
                "reference_case_id": reference["case_id"],
                "case_id": row["case_id"],
            }
            for key in ("first_energy_hartree", "final_energy_hartree"):
                if reference[key] is not None and row[key] is not None:
                    comparison[f"delta_{key}"] = row[key] - reference[key]
                else:
                    comparison[f"delta_{key}"] = None
            comparison["max_abs_delta_analytical_pv_virial_au"] = max_delta(
                reference["analytical_pv_virial_au"],
                row["analytical_pv_virial_au"],
            )
            comparison["max_abs_delta_analytical_force_au"] = max_delta(
                analytical_forces(reference), analytical_forces(row)
            )
            if (
                reference["cp2k_total_time_seconds"] is not None
                and row["cp2k_total_time_seconds"] is not None
                and reference["cp2k_total_time_seconds"] != 0.0
            ):
                comparison["time_ratio_vs_reference"] = (
                    row["cp2k_total_time_seconds"]
                    / reference["cp2k_total_time_seconds"]
                )
            else:
                comparison["time_ratio_vs_reference"] = None
            comparisons.append(comparison)

    summary = {
        "schema_version": 1,
        "campaign": ROOT.name,
        "run_count": len(runs),
        "all_program_ended": all(bool(row["program_ended"]) for row in runs),
        "runs": runs,
        "comparisons": comparisons,
    }
    (ROOT / "derived" / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    tsv_fields = [
        "case_id",
        "phase",
        "comparison_group",
        "mode",
        "batch_size",
        "symmetry",
        "pbc",
        "spin",
        "build",
        "program_ended",
        "mpi_processes",
        "omp_threads",
        "energy_evaluation_count",
        "first_energy_hartree_text",
        "final_energy_hartree_text",
        "first_energy_hartree",
        "final_energy_hartree",
        "cp2k_total_time_seconds",
        "oracle_line_count",
        "debug_virial_sum_abs_au",
        "debug_periodic_subspace_virial_sum_abs_au",
        "debug_force_sum_abs_au",
        "input_sha256",
        "output_sha256",
    ]
    with (ROOT / "derived" / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=tsv_fields)
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in tsv_fields} for row in runs)

    write_sha256sums()


if __name__ == "__main__":
    main()
