#!/usr/bin/env python3
"""Compare two direct save_tblite BvK campaigns on one complete DMC-ICE13 mesh."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


HARTREE_TO_KJMOL = 2625.4996394799
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
N_WATER = {
    "Ih": 12,
    "II": 12,
    "III": 12,
    "IV": 16,
    "VI": 10,
    "VII": 12,
    "VIII": 8,
    "IX": 12,
    "XI": 8,
    "XIII": 28,
    "XIV": 12,
    "XV": 10,
    "XVII": 6,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_energy(root: Path, mesh: int, phase: str) -> tuple[float, Path]:
    path = root / f"k{mesh}{mesh}{mesh}" / phase / "result.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        value = float(payload["energy"])
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid or incomplete direct result: {path}: {exc}") from exc
    if not math.isfinite(value):
        raise RuntimeError(f"non-finite energy in {path}")
    return value, path


def write_csv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=int, required=True)
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--author-root", type=Path, required=True)
    parser.add_argument("--reference-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.mesh <= 0:
        raise SystemExit("mesh must be positive")

    with args.reference_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        reference_field = next(
            (
                candidate
                for candidate in (
                    "DMC_relative_kJmol",
                    "dmc_reference_kj_mol_per_water",
                    "dmc_reference_kJ_mol_per_water",
                )
                if candidate in (reader.fieldnames or ())
            ),
            None,
        )
        if reference_field is None:
            raise RuntimeError(
                f"no supported DMC reference column in {args.reference_csv}"
            )
        references = {row["phase"]: float(row[reference_field]) for row in reader}
    missing_references = sorted(set(PHASES[1:]) - references.keys())
    if missing_references:
        raise RuntimeError(f"missing DMC references: {', '.join(missing_references)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    replicas = args.mesh**3
    current: dict[str, float] = {}
    author: dict[str, float] = {}
    manifest_rows: list[dict[str, object]] = []
    absolute_rows: list[dict[str, object]] = []
    for phase in PHASES:
        current_total, current_path = load_energy(args.current_root, args.mesh, phase)
        author_total, author_path = load_energy(args.author_root, args.mesh, phase)
        current[phase] = current_total / replicas
        author[phase] = author_total / replicas
        absolute_rows.append(
            {
                "mesh_n": args.mesh,
                "phase": phase,
                "current_total_Ha": f"{current_total:.15f}",
                "author_total_Ha": f"{author_total:.15f}",
                "current_Ha_per_primitive": f"{current[phase]:.15f}",
                "author_Ha_per_primitive": f"{author[phase]:.15f}",
                "author_minus_current_Ha_per_primitive": f"{author[phase] - current[phase]:+.15e}",
            }
        )
        for implementation, path in (("current", current_path), ("author_pbc", author_path)):
            manifest_rows.append(
                {
                    "mesh_n": args.mesh,
                    "phase": phase,
                    "implementation": implementation,
                    "path": str(path),
                    "sha256": sha256(path),
                }
            )

    relative_rows: list[dict[str, object]] = []
    current_absolute_errors: list[float] = []
    author_absolute_errors: list[float] = []
    shifts: list[tuple[float, str]] = []
    for phase in PHASES[1:]:
        current_relative = (
            current[phase] / N_WATER[phase] - current["Ih"] / N_WATER["Ih"]
        ) * HARTREE_TO_KJMOL
        author_relative = (
            author[phase] / N_WATER[phase] - author["Ih"] / N_WATER["Ih"]
        ) * HARTREE_TO_KJMOL
        reference = references[phase]
        current_error = current_relative - reference
        author_error = author_relative - reference
        shift = author_relative - current_relative
        current_absolute_errors.append(abs(current_error))
        author_absolute_errors.append(abs(author_error))
        shifts.append((abs(shift), phase))
        relative_rows.append(
            {
                "mesh_n": args.mesh,
                "phase": phase,
                "dmc_reference_kj_mol_per_water": f"{reference:.12f}",
                "current_save_tblite_kj_mol_per_water": f"{current_relative:.12f}",
                "author_pbc_kj_mol_per_water": f"{author_relative:.12f}",
                "author_minus_current_kj_mol_per_water": f"{shift:+.12f}",
                "current_error_kj_mol_per_water": f"{current_error:+.12f}",
                "author_error_kj_mol_per_water": f"{author_error:+.12f}",
            }
        )

    maximum_shift, maximum_phase = max(shifts)
    current_mae = statistics.mean(current_absolute_errors)
    author_mae = statistics.mean(author_absolute_errors)
    summary = {
        "mesh_n": args.mesh,
        "phase_count_including_Ih": len(PHASES),
        "relative_phase_count": len(PHASES) - 1,
        "current_save_tblite_mae_kj_mol_per_water": current_mae,
        "author_pbc_mae_kj_mol_per_water": author_mae,
        "author_minus_current_mae_kj_mol_per_water": author_mae - current_mae,
        "maximum_absolute_relative_energy_shift_kj_mol_per_water": maximum_shift,
        "maximum_shift_phase": maximum_phase,
        "status": "PASS",
    }

    write_csv(
        args.output_dir / "absolute_energy_comparison.csv",
        (
            "mesh_n",
            "phase",
            "current_total_Ha",
            "author_total_Ha",
            "current_Ha_per_primitive",
            "author_Ha_per_primitive",
            "author_minus_current_Ha_per_primitive",
        ),
        absolute_rows,
    )
    write_csv(
        args.output_dir / "relative_energy_comparison.csv",
        (
            "mesh_n",
            "phase",
            "dmc_reference_kj_mol_per_water",
            "current_save_tblite_kj_mol_per_water",
            "author_pbc_kj_mol_per_water",
            "author_minus_current_kj_mol_per_water",
            "current_error_kj_mol_per_water",
            "author_error_kj_mol_per_water",
        ),
        relative_rows,
    )
    write_csv(
        args.output_dir / "input_manifest.csv",
        ("mesh_n", "phase", "implementation", "path", "sha256"),
        manifest_rows,
    )
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
